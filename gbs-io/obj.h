#pragma once
/**
 * @file obj.h
 * @brief Wavefront OBJ files of a polygon mesh: its vertices, its faces of any number of sides,
 * its line and point elements; read and written exactly, whatever the locale.
 *
 * A polygon mesh is what a modeller (Blender, Maya, Rhino) exports as OBJ: faces of three sides
 * or more given by their vertices, in the order of the file (counterclockwise seen from outside
 * in the usual convention). `l` (polylines) and `p` (points) are the other elements of the
 * format that refer to vertices only. The triangle meshes of the BREP are written by
 * `mesh_writer.h`, with their normals.
 *
 * Reading (`read_obj`):
 * - `v x y z`: a vertex; numbers after z (a weight, a colour) are allowed and not kept;
 * - `f`, `l`, `p`: elements of vertex references `i`, `i/t`, `i//n` or `i/t/n`, of which only i is
 *   kept: from 1, or negative (relative to the vertices read so far: -1 is the last one); a
 *   reference beyond the vertices read is allowed if the file defines it later;
 * - `#` starts a comment, a line ending with `\` goes on with the next one, a line may end
 *   with `\r\n`;
 * - every other statement (`vt`, `vn`, `g`, `o`, `s`, `usemtl`, `mtllib`, the free-form curves
 *   and surfaces...) is skipped, and counted by keyword (`ObjRead::skipped`): nothing is dropped
 *   without a word;
 * - a face of fewer than 3 vertices, a line of fewer than 2, a reference of 0, out of range or
 *   malformed, a coordinate missing, malformed or not finite throw `ObjError`, with the line.
 *
 * Writing (`write_obj`): comment lines, then `v` with the shortest decimal form of each
 * coordinate that reads back exactly (`append_real`), then `f`, `l`, `p` from 1. A file is written
 * whole or not at all (`write_text_atomic`).
 */
#include <gbs/gbslib.h>
#include <gbs-io/files.h>
#include <gbs-io/real_format.h>

#include <algorithm>
#include <cerrno>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <format>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

namespace gbs::io
{
    /// A polygon mesh: vertices, faces of any number of sides, polylines and points, by the index of their vertices (from 0).
    template <std::floating_point T>
    struct PolygonMesh
    {
        points_vector<T, 3> vertices;
        std::vector<std::vector<std::size_t>> faces; ///< 3 vertices or more each
        std::vector<std::vector<std::size_t>> lines; ///< `l` elements: polylines of 2 vertices or more
        std::vector<std::size_t> points;             ///< `p` elements

        auto operator==(const PolygonMesh &) const -> bool = default;
    };

    /// What `read_obj` read, and the statements it skipped, counted by keyword in the order met.
    template <std::floating_point T>
    struct ObjRead
    {
        PolygonMesh<T> mesh;
        std::vector<std::pair<std::string, std::size_t>> skipped;
    };

    /// An OBJ text that cannot be read: the line (from 1) and why.
    class ObjError : public std::runtime_error
    {
        std::size_t reason_at_{};

    public:
        std::size_t line{};

        ObjError(std::string_view source, std::size_t line_, std::string_view reason_)
            : std::runtime_error{std::format("{}: line {}: {}", source, line_, reason_)},
              reason_at_{std::string_view{what()}.size() - reason_.size()}, line{line_}
        {
        }

        /// why the line is not valid
        [[nodiscard]] auto reason() const noexcept -> std::string_view { return std::string_view{what()}.substr(reason_at_); }
    };

    namespace obj_detail
    {
        constexpr bool is_space(char c) noexcept { return c == ' ' || c == '\t' || c == '\r' || c == '\f' || c == '\v'; }

        /// The next word of a statement, from `pos`; empty at its end.
        inline auto word(std::string_view s, std::size_t &pos) -> std::string_view
        {
            while (pos < s.size() && is_space(s[pos]))
                ++pos;
            const auto start = pos;
            while (pos < s.size() && !is_space(s[pos]))
                ++pos;
            return s.substr(start, pos - start);
        }

        /// A reference to a vertex (its first number), resolved against the n vertices read so
        /// far; `none` and a reason when it is malformed or 0.
        inline auto reference(std::string_view w, std::size_t n, std::string &why) -> std::ptrdiff_t
        {
            const auto first = w.substr(0, w.find('/'));
            long long i{};
            const auto [p, ec] = std::from_chars(first.data(), first.data() + first.size(), i);
            if (first.empty() || ec != std::errc{} || p != first.data() + first.size())
            {
                why = std::format("malformed vertex reference \"{}\"", w);
                return -1;
            }
            if (i == 0)
            {
                why = std::format("vertex reference 0 in \"{}\" (they count from 1)", w);
                return -1;
            }
            if (i < 0)
            {
                if (static_cast<unsigned long long>(-i) > n)
                {
                    why = std::format("vertex reference {} before the first vertex ({} read so far)", i, n);
                    return -1;
                }
                return static_cast<std::ptrdiff_t>(n) + static_cast<std::ptrdiff_t>(i);
            }
            return static_cast<std::ptrdiff_t>(i - 1);
        }
    }

    /**
     * @brief The polygon mesh of an OBJ text (see the file comment).
     *
     * @param source what the text is, for the message of an error (a file name)
     * @throws ObjError naming the line of the first fault
     */
    template <std::floating_point T = double>
    [[nodiscard]] auto read_obj(std::string_view text, std::string_view source = "OBJ") -> ObjRead<T>
    {
        using namespace obj_detail;
        ObjRead<T> out;
        auto &m = out.mesh;
        // the line each forward reference was met on, to name it if it stays undefined
        std::vector<std::pair<std::size_t, std::size_t>> forward; // {vertex index, line}
        std::string statement;
        std::size_t line = 0, first_line = 1;
        std::size_t pos = 0;
        while (pos <= text.size())
        {
            const auto end = std::min(text.find('\n', pos), text.size());
            auto l = text.substr(pos, end - pos);
            pos = end + 1;
            ++line;
            if (l.ends_with('\r'))
                l.remove_suffix(1);
            if (statement.empty())
                first_line = line;
            // a comment ends the line; a backslash at its end joins the next one
            if (const auto hash = l.find('#'); hash != std::string_view::npos)
                l = l.substr(0, hash);
            const bool goes_on = l.ends_with('\\');
            if (goes_on)
                l.remove_suffix(1);
            statement.append(l);
            if (goes_on && pos <= text.size())
            {
                statement += ' ';
                continue;
            }
            const std::string_view s{statement};
            std::size_t at = 0;
            const auto key = word(s, at);
            const auto fail = [&](std::string_view why) { return ObjError{source, first_line, why}; };
            if (key.empty())
            {
                // a blank line, a comment
            }
            else if (key == "v")
            {
                point<T, 3> p{};
                for (std::size_t c = 0; c < 3; ++c)
                {
                    const auto w = word(s, at);
                    double x{};
                    if (w.empty())
                        throw fail("a vertex needs three coordinates");
                    if (parse_real(w, x) != std::errc{} || !std::isfinite(x))
                        throw fail(std::format("coordinate \"{}\" is not a finite number", w));
                    p[c] = static_cast<T>(x);
                }
                m.vertices.push_back(p);
            }
            else if (key == "f" || key == "fo" || key == "l" || key == "p")
            {
                std::vector<std::size_t> refs;
                std::string why;
                for (auto w = word(s, at); !w.empty(); w = word(s, at))
                {
                    const auto i = reference(w, m.vertices.size(), why);
                    if (i < 0)
                        throw fail(why);
                    if (static_cast<std::size_t>(i) >= m.vertices.size())
                        forward.emplace_back(static_cast<std::size_t>(i), first_line);
                    refs.push_back(static_cast<std::size_t>(i));
                }
                if (key == "p")
                {
                    if (refs.empty())
                        throw fail("a point element needs a vertex");
                    m.points.insert(m.points.end(), refs.begin(), refs.end());
                }
                else if (key == "l")
                {
                    if (refs.size() < 2)
                        throw fail("a line needs 2 vertices or more");
                    m.lines.push_back(std::move(refs));
                }
                else
                {
                    if (refs.size() < 3)
                        throw fail("a face needs 3 vertices or more");
                    m.faces.push_back(std::move(refs));
                }
            }
            else
            {
                const auto k = std::string{key};
                auto it = std::ranges::find(out.skipped, k, &std::pair<std::string, std::size_t>::first);
                if (it == out.skipped.end())
                    out.skipped.emplace_back(k, 1);
                else
                    ++it->second;
            }
            statement.clear();
        }
        for (const auto &[i, at] : forward)
            if (i >= m.vertices.size())
                throw ObjError{source, at, std::format("vertex {} referred to, the file has {}", i + 1, m.vertices.size())};
        return out;
    }

    /**
     * @brief The polygon mesh of an OBJ file (`read_obj`).
     *
     * @throws std::filesystem::filesystem_error if the file cannot be read, ObjError on a fault
     */
    template <std::floating_point T = double>
    [[nodiscard]] auto read_obj_file(const std::filesystem::path &path) -> ObjRead<T>
    {
        errno = 0;
        const auto text = read_text(path);
        if (!text)
        {
            const int reason = errno;
            std::error_code e;
            throw std::filesystem::filesystem_error{
                "cannot read the OBJ file", path,
                std::error_code{std::filesystem::is_directory(path, e) ? EISDIR : reason != 0 ? reason : EIO, std::generic_category()}};
        }
        return read_obj<T>(*text, utf8_string(path));
    }

    /**
     * @brief The OBJ text of a polygon mesh (see the file comment), each line of `comment` as a
     * comment line first.
     *
     * @throws std::invalid_argument for a vertex that is not finite, an element of too few
     * vertices or a reference out of range: what read_obj would refuse
     */
    template <std::floating_point T>
    [[nodiscard]] auto write_obj(const PolygonMesh<T> &mesh, std::string_view comment = {}) -> std::string
    {
        const auto n = mesh.vertices.size();
        const auto check = [&](const std::vector<std::size_t> &e, std::size_t least, std::string_view what, std::size_t k) {
            if (e.size() < least)
                throw std::invalid_argument{std::format("OBJ: {} {} has {} vertices, fewer than {}", what, k, e.size(), least)};
            for (const auto v : e)
                if (v >= n)
                    throw std::invalid_argument{std::format("OBJ: {} {} refers to vertex {}, the mesh has {}", what, k, v, n)};
        };
        for (std::size_t k = 0; k < mesh.faces.size(); ++k)
            check(mesh.faces[k], 3, "face", k);
        for (std::size_t k = 0; k < mesh.lines.size(); ++k)
            check(mesh.lines[k], 2, "line", k);
        check(mesh.points, 0, "points", 0);
        std::string out;
        for (std::size_t pos = 0; pos < comment.size();)
        {
            const auto end = std::min(comment.find('\n', pos), comment.size());
            out.append("# ").append(comment.substr(pos, end - pos)).append("\n");
            pos = end + 1;
        }
        for (std::size_t i = 0; i < n; ++i)
        {
            out += 'v';
            for (const T x : mesh.vertices[i])
            {
                if (!std::isfinite(x))
                    throw std::invalid_argument{std::format("OBJ: vertex {} is not finite", i)};
                out += ' ';
                append_real(out, double(x));
            }
            out += '\n';
        }
        const auto element = [&](char tag, const std::vector<std::size_t> &e) {
            out += tag;
            for (const auto v : e)
                out.append(" ").append(std::to_string(v + 1));
            out += '\n';
        };
        for (const auto &f : mesh.faces)
            element('f', f);
        for (const auto &l : mesh.lines)
            element('l', l);
        if (!mesh.points.empty())
            element('p', mesh.points);
        return out;
    }

    /**
     * @brief Writes the OBJ file of a polygon mesh (`write_obj`), whole or not at all.
     *
     * @throws std::invalid_argument as write_obj, std::filesystem::filesystem_error if the file
     * cannot be written
     */
    template <std::floating_point T>
    void write_obj_file(const PolygonMesh<T> &mesh, const std::filesystem::path &path, std::string_view comment = {})
    {
        write_text_atomic(path, write_obj(mesh, comment));
    }
}
