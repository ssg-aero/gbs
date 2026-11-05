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
 * - `v x y z`: a vertex; up to four numbers after z (a weight, a colour) are allowed and not
 *   kept, any other word is a fault;
 * - `f`, `l`, `p`: elements of vertex references `i`, `i/t`, `i//n` or `i/t/n`, of which only i is
 *   kept: from 1, or negative (relative to the vertices read so far: -1 is the last one); a
 *   reference beyond the vertices read is allowed if the file defines it later;
 * - `#` starts a comment, a line ending with `\` goes on with the next one; a line ends with
 *   `\n`, `\r\n` or a lone `\r` (the files of the old Mac OS); a UTF-8 byte order mark at the
 *   start is skipped;
 * - every other statement (`vt`, `vn`, `g`, `o`, `s`, `usemtl`, `mtllib`, the free-form curves
 *   and surfaces...) is skipped, and counted by keyword (`ObjRead::skipped`): nothing is dropped
 *   without a word;
 * - a face of fewer than 3 vertices, a line of fewer than 2, a reference of 0, out of range or
 *   malformed, a coordinate missing, malformed or not finite in T throw `ObjError`, with the line.
 *
 * Writing (`write_obj`): comment lines, then `v` with the shortest decimal form of each
 * coordinate that reads back exactly in T (`append_real`), then `f`, `l`, `p` from 1. A file is
 * written whole or not at all (`write_text_atomic`).
 */
#include <gbs/gbslib.h>
#include <gbs-io/files.h>
#include <gbs-io/real_format.h>

#include <algorithm>
#include <charconv>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <filesystem>
#include <format>
#include <initializer_list>
#include <iterator>
#include <limits>
#include <ranges>
#include <span>
#include <type_traits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

namespace gbs::io
{
    /// The coordinates an OBJ file is read into and written from: those that `parse_real` and
    /// `append_real` give exactly.
    template <typename T>
    concept obj_real = std::same_as<T, float> || std::same_as<T, double>;

    /**
     * Polygons (or polylines) by the index of their vertices, flat: the vertices of the k-th are
     * `indices[offsets[k], offsets[k + 1])`, one allocation for all of them rather than one each.
     * `p[k]` is a span of them; a range of spans.
     */
    struct Polygons
    {
        std::vector<std::uint32_t> indices;
        std::vector<std::size_t> offsets{0}; ///< size() + 1, from 0

        Polygons() = default;
        Polygons(std::initializer_list<std::initializer_list<std::uint32_t>> list)
        {
            for (const auto &p : list)
                push_back(p);
        }
        /// Each range of indices a polygon.
        template <std::ranges::input_range R>
            requires(!std::same_as<std::remove_cvref_t<R>, Polygons>) && std::ranges::input_range<std::ranges::range_reference_t<R>>
        explicit Polygons(R &&list)
        {
            for (auto &&p : list)
                push_back(p);
        }

        [[nodiscard]] auto size() const noexcept -> std::size_t { return offsets.size() - 1; }
        [[nodiscard]] auto empty() const noexcept -> bool { return offsets.size() == 1; }
        [[nodiscard]] auto operator[](std::size_t k) const -> std::span<const std::uint32_t>
        {
            return {indices.data() + offsets[k], offsets[k + 1] - offsets[k]};
        }
        /// Appends a polygon.
        void push_back(std::initializer_list<std::uint32_t> p) { push_back(std::span<const std::uint32_t>{p.begin(), p.size()}); }
        template <std::ranges::input_range R>
        void push_back(R &&p)
        {
            for (auto &&i : p)
                indices.push_back(static_cast<std::uint32_t>(i));
            offsets.push_back(indices.size());
        }
        void clear() noexcept
        {
            indices.clear();
            offsets.assign(1, 0);
        }

        /// The polygons in order, as spans.
        class const_iterator
        {
            const Polygons *p_{};
            std::size_t k_{};

        public:
            using iterator_concept = std::forward_iterator_tag;
            using value_type = std::span<const std::uint32_t>;
            using difference_type = std::ptrdiff_t;
            const_iterator() = default;
            const_iterator(const Polygons *p, std::size_t k) : p_{p}, k_{k} {}
            auto operator*() const -> value_type { return (*p_)[k_]; }
            auto operator++() -> const_iterator &
            {
                ++k_;
                return *this;
            }
            auto operator++(int) -> const_iterator
            {
                auto t = *this;
                ++k_;
                return t;
            }
            auto operator==(const const_iterator &) const -> bool = default;
        };
        [[nodiscard]] auto begin() const -> const_iterator { return {this, 0}; }
        [[nodiscard]] auto end() const -> const_iterator { return {this, size()}; }

        auto operator==(const Polygons &) const -> bool = default;
    };

    /// A polygon mesh: vertices, faces of any number of sides, polylines and points, by the index of their vertices (from 0).
    template <obj_real T>
    struct PolygonMesh
    {
        points_vector<T, 3> vertices;
        Polygons faces;                    ///< 3 vertices or more each
        Polygons lines;                    ///< `l` elements: polylines of 2 vertices or more
        std::vector<std::uint32_t> points; ///< `p` elements

        auto operator==(const PolygonMesh &) const -> bool = default;
    };

    /// What `read_obj` read, and the statements it skipped, counted by keyword in the order met.
    template <obj_real T>
    struct ObjRead
    {
        PolygonMesh<T> mesh;
        std::vector<std::pair<std::string, std::size_t>> skipped;
    };

    /// An OBJ text that cannot be read: the line (from 1) and why (`reason()`, the end of what()).
    class ObjError : public TextError
    {
    public:
        std::size_t line{};

        ObjError(std::string_view source, std::size_t line_, std::string_view reason_)
            : TextError{std::format("{}: line {}: {}", source, line_, reason_), reason_}, line{line_}
        {
        }
    };

    namespace obj_detail
    {
        /// The most vertices a mesh holds: its indices are 32 bits, as those of the BREP meshes.
        inline constexpr std::size_t max_vertices = std::numeric_limits<std::uint32_t>::max();

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
        /// far (n <= max_vertices); why it is not one: malformed, 0, before the first vertex or
        /// beyond the vertices a mesh holds. Beyond n, it refers to a vertex defined later.
        inline auto reference(std::string_view w, std::size_t n) -> std::expected<std::uint32_t, std::string>
        {
            const auto first = w.substr(0, w.find('/'));
            long long i{};
            const auto [p, ec] = std::from_chars(first.data(), first.data() + first.size(), i);
            if (first.empty() || ec != std::errc{} || p != first.data() + first.size())
                return std::unexpected{std::format("malformed vertex reference \"{}\"", w)};
            if (i == 0)
                return std::unexpected{std::format("vertex reference 0 in \"{}\" (they count from 1)", w)};
            if (i < 0)
            {
                if (i < -static_cast<long long>(n)) // (no negation of i: -LLONG_MIN overflows)
                    return std::unexpected{std::format("vertex reference {} before the first vertex ({} read so far)", i, n)};
                return static_cast<std::uint32_t>(static_cast<long long>(n) + i);
            }
            if (static_cast<unsigned long long>(i) > max_vertices)
                return std::unexpected{std::format("vertex reference {} beyond the {} vertices a mesh holds", i, max_vertices)};
            return static_cast<std::uint32_t>(i - 1);
        }

        /// The lines of a text, with their number from 1: ends at `\n`, `\r\n` or a lone `\r`.
        inline auto next_line(std::string_view text, std::size_t &pos) -> std::string_view
        {
            const auto end = std::min(text.find_first_of("\r\n", pos), text.size());
            const auto l = text.substr(pos, end - pos);
            pos = end + (end + 1 < text.size() && text[end] == '\r' && text[end + 1] == '\n' ? 2 : 1);
            return l;
        }
    }

    /**
     * @brief The polygon mesh of an OBJ text (see the file comment).
     *
     * @param source what the text is, for the message of an error (a file name)
     * @throws ObjError naming the line of the first fault
     */
    template <obj_real T = double>
    [[nodiscard]] auto read_obj(std::string_view text, std::string_view source = "OBJ") -> ObjRead<T>
    {
        using namespace obj_detail;
        ObjRead<T> out;
        auto &m = out.mesh;
        if (text.starts_with("\xEF\xBB\xBF")) // a UTF-8 byte order mark (Windows tools)
            text.remove_prefix(3);
        // the line each forward reference was met on, to name it if it stays undefined
        std::vector<std::pair<std::uint32_t, std::size_t>> forward; // {vertex index, line}
        std::vector<std::uint32_t> refs;                            // (reused from one element to the next)
        std::string statement;
        std::size_t line = 0, first_line = 1;
        std::size_t pos = 0;
        while (pos <= text.size())
        {
            auto l = next_line(text, pos);
            ++line;
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
            const auto real = [&](std::string_view w) {
                T x{};
                if (parse_real(w, x) != std::errc{} || !std::isfinite(x))
                    throw fail(std::format("coordinate \"{}\" is not a finite number{}", w, std::same_as<T, float> ? " (as a float)" : ""));
                return x;
            };
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
                    if (w.empty())
                        throw fail("a vertex needs three coordinates");
                    p[c] = real(w);
                }
                // a weight, a colour: numbers, not kept
                std::size_t more = 0;
                for (auto w = word(s, at); !w.empty(); w = word(s, at), ++more)
                {
                    if (more == 4)
                        throw fail(std::format("\"{}\": a vertex has at most 4 numbers after its coordinates", w));
                    (void)real(w);
                }
                if (m.vertices.size() == max_vertices)
                    throw fail(std::format("more than {} vertices", max_vertices));
                m.vertices.push_back(p);
            }
            else if (key == "f" || key == "fo" || key == "l" || key == "p")
            {
                refs.clear();
                for (auto w = word(s, at); !w.empty(); w = word(s, at))
                {
                    const auto i = reference(w, m.vertices.size());
                    if (!i)
                        throw fail(i.error());
                    if (*i >= m.vertices.size())
                        forward.emplace_back(*i, first_line);
                    refs.push_back(*i);
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
                    m.lines.push_back(refs);
                }
                else
                {
                    if (refs.size() < 3)
                        throw fail("a face needs 3 vertices or more");
                    m.faces.push_back(refs);
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
                throw ObjError{source, at, std::format("vertex {} referred to, the file has {}", std::size_t{i} + 1, m.vertices.size())};
        return out;
    }

    /**
     * @brief The polygon mesh of an OBJ file (`read_obj`).
     *
     * @throws std::filesystem::filesystem_error if the file cannot be read, ObjError on a fault
     */
    template <obj_real T = double>
    [[nodiscard]] auto read_obj_file(const std::filesystem::path &path) -> ObjRead<T>
    {
        return read_obj<T>(read_text_or_throw(path, "cannot read the OBJ file"), utf8_string(path));
    }

    /**
     * @brief The OBJ text of a polygon mesh (see the file comment), each line of `comment` as a
     * comment line first.
     *
     * @throws std::invalid_argument for a vertex that is not finite, an element of too few
     * vertices, a reference out of range or more vertices than an index holds: what read_obj
     * would refuse
     */
    template <obj_real T>
    [[nodiscard]] auto write_obj(const PolygonMesh<T> &mesh, std::string_view comment = {}) -> std::string
    {
        const auto n = mesh.vertices.size();
        if (n > obj_detail::max_vertices)
            throw std::invalid_argument{std::format("OBJ: {} vertices, more than the {} a mesh holds", n, obj_detail::max_vertices)};
        const auto check = [&](std::span<const std::uint32_t> e, std::size_t least, std::string_view what, std::size_t k) {
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
        for (std::size_t i = 0; i < n; ++i)
            if (!std::ranges::all_of(mesh.vertices[i], [](T x) { return std::isfinite(x); }))
                throw std::invalid_argument{std::format("OBJ: vertex {} is not finite", i)};
        std::string out;
        if (comment.ends_with('\n'))
            comment.remove_suffix(1); // (no empty comment line for it)
        if (!comment.empty())
            for (const auto line : std::views::split(comment, '\n'))
                out.append("# ").append(std::string_view(line.begin(), line.end())).append("\n");
        for (const auto &p : mesh.vertices)
        {
            out += 'v';
            for (const T x : p)
            {
                out += ' ';
                append_real(out, x);
            }
            out += '\n';
        }
        const auto element = [&](char tag, std::span<const std::uint32_t> e) {
            out += tag;
            char buf[16];
            for (const auto v : e)
            {
                out += ' ';
                out.append(buf, std::to_chars(buf, buf + sizeof buf, std::uint64_t{v} + 1).ptr);
            }
            out += '\n';
        };
        for (const auto f : mesh.faces)
            element('f', f);
        for (const auto l : mesh.lines)
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
    template <obj_real T>
    void write_obj_file(const PolygonMesh<T> &mesh, const std::filesystem::path &path, std::string_view comment = {})
    {
        write_text_atomic(path, write_obj(mesh, comment));
    }
}
