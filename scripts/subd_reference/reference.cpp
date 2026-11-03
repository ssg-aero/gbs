// Values of reference of the Catmull-Clark subdivision computed by OpenSubdiv, for the tests of
// gbs/subd (tests/tests_subd.cpp). Not part of the build of gbs: built by hand against the
// OpenSubdiv of the development environment (libopensubdiv on conda-forge), run when the cages
// (cages.py) or the rules change, its output committed:
//
//     cmake -S scripts/subd_reference -B build-subd-reference -G Ninja -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
//     cmake --build build-subd-reference
//     build-subd-reference/subd_reference tests/in/subd
//
// For each <name>.cage.json of the directory, writes <name>.osd.json:
// - level1, level2: the cage refined uniformly once and twice (vertices, faces);
// - limit2: OpenSubdiv's limit position and tangents at the vertices of level 2 (its masks,
//   exact except at darts and across irregular creases, see gbs/subd/limit.h);
// - patches: the number of ptex faces and the limit surface at random (ptex face, u, v), its
//   first and second derivatives, from the patch table (adaptive refinement to level 10,
//   infinitely sharp patches, double precision), kept only on regular patches, which are exact
//   (the Gregory end caps around irregular vertices are not).
//
// Rules: Catmark, boundary interpolation "edge and corner", infinitely sharp creases and corners.
#include <opensubdiv/far/patchMap.h>
#include <opensubdiv/far/patchTable.h>
#include <opensubdiv/far/patchTableFactory.h>
#include <opensubdiv/far/primvarRefiner.h>
#include <opensubdiv/far/ptexIndices.h>
#include <opensubdiv/far/topologyDescriptor.h>
#include <opensubdiv/version.h>

#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>
#include <rapidjson/ostreamwrapper.h>
#include <rapidjson/prettywriter.h>

#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

using namespace OpenSubdiv;
namespace fs = std::filesystem;
using Writer = rapidjson::PrettyWriter<rapidjson::OStreamWrapper>;

namespace
{
    struct Point
    {
        double p[3]{};
        void Clear(void * = nullptr) { p[0] = p[1] = p[2] = 0; }
        void AddWithWeight(const Point &s, double w)
        {
            for (int i = 0; i < 3; ++i)
                p[i] += w * s.p[i];
        }
    };

    struct Cage
    {
        std::vector<Point> vertices;
        std::vector<int> sizes, indices, creases, corners;
    };

    auto fnv1a(const std::string &s) -> std::uint64_t
    {
        std::uint64_t h = 14695981039346656037ull;
        for (const unsigned char c : s)
            h = (h ^ c) * 1099511628211ull;
        return h;
    }

    auto read_cage(const fs::path &path) -> Cage
    {
        std::ifstream in{path};
        rapidjson::IStreamWrapper isw{in};
        rapidjson::Document d;
        d.ParseStream(isw);
        Cage c;
        for (const auto &v : d["vertices"].GetArray())
            c.vertices.push_back({{v[0].GetDouble(), v[1].GetDouble(), v[2].GetDouble()}});
        for (const auto &f : d["faces"].GetArray())
        {
            c.sizes.push_back(int(f.Size()));
            for (const auto &i : f.GetArray())
                c.indices.push_back(i.GetInt());
        }
        for (const auto &e : d["creases"].GetArray())
            c.creases.insert(c.creases.end(), {e[0].GetInt(), e[1].GetInt()});
        for (const auto &v : d["corners"].GetArray())
            c.corners.push_back(v.GetInt());
        return c;
    }

    auto make_refiner(const Cage &c) -> std::unique_ptr<Far::TopologyRefiner>
    {
        Far::TopologyDescriptor desc;
        desc.numVertices = int(c.vertices.size());
        desc.numFaces = int(c.sizes.size());
        desc.numVertsPerFace = c.sizes.data();
        desc.vertIndicesPerFace = c.indices.data();
        std::vector<float> crease_weights(c.creases.size() / 2, Sdc::Crease::SHARPNESS_INFINITE);
        std::vector<float> corner_weights(c.corners.size(), Sdc::Crease::SHARPNESS_INFINITE);
        desc.numCreases = int(crease_weights.size());
        desc.creaseVertexIndexPairs = c.creases.data();
        desc.creaseWeights = crease_weights.data();
        desc.numCorners = int(corner_weights.size());
        desc.cornerVertexIndices = c.corners.data();
        desc.cornerWeights = corner_weights.data();
        Sdc::Options opts;
        opts.SetVtxBoundaryInterpolation(Sdc::Options::VTX_BOUNDARY_EDGE_AND_CORNER);
        using Factory = Far::TopologyRefinerFactory<Far::TopologyDescriptor>;
        return std::unique_ptr<Far::TopologyRefiner>{Factory::Create(desc, Factory::Options(Sdc::SCHEME_CATMARK, opts))};
    }

    void write_point(Writer &w, const Point &p)
    {
        w.StartArray();
        for (double x : p.p)
            w.Double(x);
        w.EndArray();
    }

    void write_level(Writer &w, const Far::TopologyLevel &level, const Point *vertices)
    {
        w.StartObject();
        w.Key("vertices");
        w.StartArray();
        for (int v = 0; v < level.GetNumVertices(); ++v)
            write_point(w, vertices[v]);
        w.EndArray();
        w.Key("faces");
        w.StartArray();
        for (int f = 0; f < level.GetNumFaces(); ++f)
        {
            w.StartArray();
            for (const int v : level.GetFaceVertices(f))
                w.Int(v);
            w.EndArray();
        }
        w.EndArray();
        w.EndObject();
    }

    void uniform(Writer &w, const Cage &c)
    {
        auto refiner = make_refiner(c);
        Far::TopologyRefiner::UniformOptions uo{2};
        uo.fullTopologyInLastLevel = true; // needed by Limit()
        refiner->RefineUniform(uo);
        std::vector<Point> buffer(refiner->GetNumVerticesTotal());
        std::copy(c.vertices.begin(), c.vertices.end(), buffer.begin());
        Far::PrimvarRefinerReal<double> primvar{*refiner};
        Point *src = buffer.data();
        std::array<Point *, 3> level{};
        level[0] = src;
        for (int l = 1; l <= 2; ++l)
        {
            Point *dst = level[l - 1] + refiner->GetLevel(l - 1).GetNumVertices();
            primvar.Interpolate(l, level[l - 1], dst);
            level[l] = dst;
        }
        w.Key("level1");
        write_level(w, refiner->GetLevel(1), level[1]);
        w.Key("level2");
        write_level(w, refiner->GetLevel(2), level[2]);

        const auto n = refiner->GetLevel(2).GetNumVertices();
        std::vector<Point> pos(n), tan1(n), tan2(n);
        Point *pp = pos.data(), *pt1 = tan1.data(), *pt2 = tan2.data();
        primvar.Limit(level[2], pp, pt1, pt2);
        w.Key("limit2");
        w.StartObject();
        for (const auto &[key, values] : {std::pair{"positions", &pos}, std::pair{"tangent1", &tan1}, std::pair{"tangent2", &tan2}})
        {
            w.Key(key);
            w.StartArray();
            for (const auto &p : *values)
                write_point(w, p);
            w.EndArray();
        }
        w.EndObject();
    }

    void patches(Writer &w, const Cage &c, std::uint64_t seed)
    {
        auto refiner = make_refiner(c);
        Far::PatchTableFactory::Options po{10};
        po.SetEndCapType(Far::PatchTableFactory::Options::ENDCAP_GREGORY_BASIS);
        po.useInfSharpPatch = true;
        po.SetPatchPrecision<double>();
        refiner->RefineAdaptive(po.GetRefineAdaptiveOptions());
        std::unique_ptr<Far::PatchTable> table{Far::PatchTableFactory::Create(*refiner, po)};

        const int total = refiner->GetNumVerticesTotal();
        std::vector<Point> buffer(total + table->GetNumLocalPoints());
        std::copy(c.vertices.begin(), c.vertices.end(), buffer.begin());
        Far::PrimvarRefinerReal<double> primvar{*refiner};
        Point *src = buffer.data();
        for (int l = 1; l <= refiner->GetMaxLevel(); ++l)
        {
            Point *dst = src + refiner->GetLevel(l - 1).GetNumVertices();
            primvar.Interpolate(l, src, dst);
            src = dst;
        }
        if (const auto *local = table->GetLocalPointStencilTable<double>())
            local->UpdateValues(buffer.data(), buffer.data() + total);

        const Far::PtexIndices ptex{*refiner};
        const Far::PatchMap map{*table};
        std::mt19937_64 rng{seed};
        std::uniform_real_distribution<double> unit{0.0, 1.0};
        std::uniform_int_distribution<int> face{0, ptex.GetNumFaces() - 1};

        w.Key("patches");
        w.StartObject();
        w.Key("count");
        w.Int(ptex.GetNumFaces());
        w.Key("samples");
        w.StartArray();
        int kept = 0;
        for (int s = 0; s < 4000 && kept < 300; ++s)
        {
            const int f = face(rng);
            const double u = unit(rng), v = unit(rng);
            const auto *h = map.FindPatch(f, u, v);
            if (!h || table->GetPatchDescriptor(*h).GetType() != Far::PatchDescriptor::REGULAR)
                continue;
            double wP[20], wDu[20], wDv[20], wDuu[20], wDuv[20], wDvv[20];
            table->EvaluateBasis(*h, u, v, wP, wDu, wDv, wDuu, wDuv, wDvv);
            const auto cvs = table->GetPatchVertices(*h);
            std::array<Point, 6> r{};
            for (int i = 0; i < cvs.size(); ++i)
            {
                const auto &q = buffer[cvs[i]];
                r[0].AddWithWeight(q, wP[i]);
                r[1].AddWithWeight(q, wDu[i]);
                r[2].AddWithWeight(q, wDv[i]);
                r[3].AddWithWeight(q, wDuu[i]);
                r[4].AddWithWeight(q, wDuv[i]);
                r[5].AddWithWeight(q, wDvv[i]);
            }
            w.StartObject();
            w.Key("patch");
            w.Int(f);
            w.Key("u");
            w.Double(u);
            w.Key("v");
            w.Double(v);
            const char *keys[] = {"p", "du", "dv", "duu", "duv", "dvv"};
            for (int i = 0; i < 6; ++i)
            {
                w.Key(keys[i]);
                write_point(w, r[i]);
            }
            w.EndObject();
            ++kept;
        }
        w.EndArray();
        w.EndObject();
    }
}

int main(int argc, char **argv)
{
    const fs::path dir = argc > 1 ? argv[1] : "tests/in/subd";
    for (const auto &entry : fs::directory_iterator{dir})
    {
        const auto name = entry.path().filename().string();
        const std::string suffix = ".cage.json";
        if (name.size() <= suffix.size() || name.compare(name.size() - suffix.size(), suffix.size(), suffix) != 0)
            continue;
        const auto stem = name.substr(0, name.size() - suffix.size());
        const auto cage = read_cage(entry.path());
        std::ofstream out{dir / (stem + ".osd.json")};
        rapidjson::OStreamWrapper osw{out};
        Writer w{osw};
        w.SetIndent(' ', 1);
        w.SetFormatOptions(rapidjson::kFormatSingleLineArray);
        w.StartObject();
        w.Key("opensubdiv");
        w.String((std::to_string(OPENSUBDIV_VERSION_MAJOR) + "." + std::to_string(OPENSUBDIV_VERSION_MINOR) + "." + std::to_string(OPENSUBDIV_VERSION_PATCH)).c_str());
        uniform(w, cage);
        patches(w, cage, fnv1a(stem)); // the samples of a cage depend on its name only
        w.EndObject();
        std::cout << "written " << (dir / (stem + ".osd.json")).string() << '\n';
    }
}
