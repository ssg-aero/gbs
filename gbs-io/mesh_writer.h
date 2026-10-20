#pragma once
/**
 * @file mesh_writer.h
 * @brief STL (ASCII and binary) and Wavefront OBJ files of a triangle mesh of the native BREP
 * (`gbs-brep/mesh.h`, stage 3, `docs/sources/design/brep_palier3_mesh.md`).
 *
 * - STL: one facet per triangle, its normal that of the triangle (the corners in the sense of
 *   the mesh, counter-clockwise seen from outside); single precision in the binary form, little
 *   endian whatever the host.
 * - OBJ: the vertices once (`v`), the normal of each corner (`vn`, those of the surfaces), one
 *   group per face of the model (`g face<index>`); the vertices of a shared edge are shared.
 */
#include <gbs-brep/mesh.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <string_view>

namespace gbs::io
{
    namespace mesh_writer_detail
    {
        /// Appends the bytes of x, little endian whatever the host.
        template <typename X>
        void put(std::string &out, X x)
        {
            using U = std::conditional_t<sizeof(X) == 4, std::uint32_t, std::uint16_t>;
            U u;
            std::memcpy(&u, &x, sizeof(X));
            for (std::size_t k = 0; k < sizeof(U); ++k)
                out.push_back(static_cast<char>((u >> (8 * k)) & 0xFF));
        }

        /// Appends printf-formatted text (a point as decimal separator, in the default C locale).
        template <typename... A>
        void print(std::string &out, const char *fmt, A... a)
        {
            char b[256];
            const int n = std::snprintf(b, sizeof b, fmt, a...);
            out.append(b, static_cast<std::size_t>(std::max(n, 0)));
        }

        template <std::floating_point T>
        point<T, 3> facet_normal(const brep::Mesh<T> &mesh, const std::array<std::uint32_t, 3> &t)
        {
            const auto n = cross(mesh.vertices[t[1]] - mesh.vertices[t[0]], mesh.vertices[t[2]] - mesh.vertices[t[0]]);
            const T l = norm(n);
            return l > T(0) ? n / l : n;
        }

        inline auto write_file(const std::string &text, const std::filesystem::path &path) -> brep::BuildResult<void>
        {
            std::ofstream f(path, std::ios::binary);
            if (!f || !f.write(text.data(), static_cast<std::streamsize>(text.size())))
                return std::unexpected(brep::BuildError{brep::BuildErrc::InvalidFile, "cannot write " + path.string(), {}});
            return {};
        }
    }

    /// The STL text of a mesh: binary (80-byte header, then 50 bytes per triangle) or ASCII.
    template <std::floating_point T>
    [[nodiscard]] std::string write_stl(const brep::Mesh<T> &mesh, bool binary = true, std::string_view name = "gbs")
    {
        using namespace mesh_writer_detail;
        std::string out;
        if (binary)
        {
            out.reserve(84 + 50 * mesh.triangles.size());
            std::string header(80, ' ');
            header.replace(0, std::min<std::size_t>(name.size(), 80), name.substr(0, 80));
            out += header;
            put(out, static_cast<std::uint32_t>(mesh.triangles.size()));
            for (const auto &t : mesh.triangles)
            {
                for (T x : facet_normal(mesh, t))
                    put(out, static_cast<float>(x));
                for (auto v : t)
                    for (T x : mesh.vertices[v])
                        put(out, static_cast<float>(x));
                put(out, std::uint16_t{0});
            }
            return out;
        }
        out.append("solid ").append(name).append("\n");
        for (const auto &t : mesh.triangles)
        {
            const auto n = facet_normal(mesh, t);
            print(out, "  facet normal %.9g %.9g %.9g\n    outer loop\n", double(n[0]), double(n[1]), double(n[2]));
            for (auto v : t)
            {
                const auto &p = mesh.vertices[v];
                print(out, "      vertex %.17g %.17g %.17g\n", double(p[0]), double(p[1]), double(p[2]));
            }
            out += "    endloop\n  endfacet\n";
        }
        out.append("endsolid ").append(name).append("\n");
        return out;
    }

    /// The Wavefront OBJ text of a mesh.
    template <std::floating_point T>
    [[nodiscard]] std::string write_obj(const brep::Mesh<T> &mesh)
    {
        using namespace mesh_writer_detail;
        std::string out = "# gbs mesh\n";
        for (const auto &p : mesh.vertices)
            print(out, "v %.17g %.17g %.17g\n", double(p[0]), double(p[1]), double(p[2]));
        for (const auto &n : mesh.normals)
            print(out, "vn %.9g %.9g %.9g\n", double(n[0]), double(n[1]), double(n[2]));
        for (const auto &fr : mesh.faces)
        {
            print(out, "g face%u\n", unsigned(fr.face.index));
            for (std::uint32_t k = fr.first; k < fr.first + fr.count; ++k)
            {
                const auto &t = mesh.triangles[k];
                // (indices from 1; the normals of the triangle k are 3 k .. 3 k + 2)
                print(out, "f %u//%u %u//%u %u//%u\n", t[0] + 1, 3 * k + 1, t[1] + 1, 3 * k + 2, t[2] + 1, 3 * k + 3);
            }
        }
        return out;
    }

    /// Writes the STL file of a mesh (`write_stl`); `InvalidFile` when it cannot be written.
    template <std::floating_point T>
    auto write_stl_file(const brep::Mesh<T> &mesh, const std::filesystem::path &path, bool binary = true) -> brep::BuildResult<void>
    {
        return mesh_writer_detail::write_file(write_stl(mesh, binary, path.stem().string()), path);
    }

    /// Writes the OBJ file of a mesh (`write_obj`); `InvalidFile` when it cannot be written.
    template <std::floating_point T>
    auto write_obj_file(const brep::Mesh<T> &mesh, const std::filesystem::path &path) -> brep::BuildResult<void>
    {
        return mesh_writer_detail::write_file(write_obj(mesh), path);
    }
}
