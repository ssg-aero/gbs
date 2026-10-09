#pragma once

/**
 * @file p21.h
 * @brief Reader of the ISO 10303-21 ("Part 21") exchange structure, the text
 * format of STEP files.
 *
 * Design: docs/sources/design/step_reader.md, section 2; architecture note
 * docs/sources/design/step_pr01_p21.md.
 *
 * This layer only reads the syntax: it turns the text into a table of entity
 * instances and their parameters. It knows nothing of the STEP schemas: a
 * CYLINDRICAL_SURFACE is just an instance with that type name. The mapping to
 * gbs geometry and topology is a separate layer.
 *
 * Values live in contiguous arenas (std::vector) referenced by index, like the
 * gbs::brep::Model; views (InstanceView, ValueView, ListView) are cheap handles
 * into a P21File and stay valid as long as it lives.
 */

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cstdint>
#include <cstdlib>
#include <expected>
#include <filesystem>
#include <fstream>
#include <limits>
#include <optional>
#include <ranges>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#if defined(_LIBCPP_VERSION)
#include <clocale>
#include <locale.h>
#if defined(__APPLE__)
#include <xlocale.h>
#endif
#endif

namespace gbs::step
{
    // =========================================================================
    // Values
    // =========================================================================

    enum class ValueKind : std::uint8_t
    {
        Reference, ///< #123
        Integer,   ///< 42
        Real,      ///< 1.E-07
        String,    ///< 'text' (decoded to UTF-8)
        Enum,      ///< .T., .UNSPECIFIED. (name without dots)
        Binary,    ///< "0123" (raw hexadecimal text)
        Omitted,   ///< $
        Derived,   ///< *
        List,      ///< ( ... )
        Typed,     ///< LENGTH_MEASURE(1.E-07)
    };

    /// Syntax error: position in the text and what went wrong.
    struct P21Error
    {
        std::size_t line{}, column{};
        std::string message;
    };

    /// Raised by the views when a value is not of the requested kind (content error, not syntax).
    class P21AccessError : public std::runtime_error
    {
    public:
        using std::runtime_error::runtime_error;
    };

    class P21File;

    namespace detail
    {
        struct Node
        {
            ValueKind kind{ValueKind::Omitted};
            std::uint32_t a{}; ///< String/Enum/Binary: offset in chars ; List: first item ; Typed: type id
            std::uint32_t n{}; ///< String/Enum/Binary: length ; List: count ; Typed: child node
            union
            {
                std::int64_t i;
                double r;
                std::uint64_t ref;
            };
            Node() : i{0} {}
        };

        struct Part
        {
            std::uint32_t type{}; ///< interned type name
            std::uint32_t args{}; ///< List node
        };

        struct Instance
        {
            std::uint64_t id{}; ///< 0 for header entities
            std::uint32_t first_part{}, n_parts{};
            std::uint32_t line{};
            bool complex{false};
        };
    } // namespace detail

    class ListView;

    /// Read-only view of one parameter value.
    class ValueView
    {
        const P21File *f_{};
        std::uint32_t node_{};

    public:
        ValueView() = default;
        ValueView(const P21File *f, std::uint32_t node) : f_{f}, node_{node} {}

        [[nodiscard]] ValueKind kind() const;
        [[nodiscard]] bool is(ValueKind k) const { return kind() == k; }
        [[nodiscard]] bool is_null() const { return kind() == ValueKind::Omitted || kind() == ValueKind::Derived; }

        [[nodiscard]] std::uint64_t as_ref() const;
        [[nodiscard]] std::int64_t as_int() const;
        [[nodiscard]] double as_real() const; ///< accepts an Integer too
        [[nodiscard]] std::string_view as_string() const;
        [[nodiscard]] std::string_view as_enum() const;
        [[nodiscard]] bool as_bool() const;  ///< .T. / .F.
        [[nodiscard]] std::optional<bool> as_logical() const; ///< .T. / .F. / .U. (nullopt)
        [[nodiscard]] std::string_view as_binary() const;
        [[nodiscard]] ListView as_list() const;
        [[nodiscard]] std::string_view typed_name() const; ///< LENGTH_MEASURE in LENGTH_MEASURE(1.E-07)
        [[nodiscard]] ValueView typed_value() const;
        /// Unwraps a typed value (if any) then reads a real.
        [[nodiscard]] double as_measure() const { return is(ValueKind::Typed) ? typed_value().as_real() : as_real(); }
    };

    /// Read-only view of a list of values.
    class ListView
    {
        const P21File *f_{};
        std::uint32_t first_{}, n_{};

    public:
        ListView() = default;
        ListView(const P21File *f, std::uint32_t first, std::uint32_t n) : f_{f}, first_{first}, n_{n} {}
        [[nodiscard]] std::size_t size() const noexcept { return n_; }
        [[nodiscard]] bool empty() const noexcept { return n_ == 0; }
        [[nodiscard]] ValueView operator[](std::size_t i) const;
        [[nodiscard]] ValueView at(std::size_t i) const
        {
            if (i >= n_)
                throw P21AccessError("parameter index " + std::to_string(i) + " out of " + std::to_string(n_));
            return (*this)[i];
        }

        class iterator
        {
            const ListView *l_{};
            std::size_t i_{};

        public:
            using value_type = ValueView;
            using difference_type = std::ptrdiff_t;
            iterator() = default;
            iterator(const ListView *l, std::size_t i) : l_{l}, i_{i} {}
            ValueView operator*() const { return (*l_)[i_]; }
            iterator &operator++() { ++i_; return *this; }
            iterator operator++(int) { auto t = *this; ++i_; return t; }
            bool operator==(const iterator &o) const { return i_ == o.i_; }
        };
        [[nodiscard]] iterator begin() const { return {this, 0}; }
        [[nodiscard]] iterator end() const { return {this, n_}; }
    };

    /// One partial entity of an instance (the only one for a simple instance).
    class PartView
    {
        const P21File *f_{};
        std::uint32_t part_{};

    public:
        PartView(const P21File *f, std::uint32_t part) : f_{f}, part_{part} {}
        [[nodiscard]] std::string_view type() const;
        [[nodiscard]] ListView args() const;
    };

    /// Read-only view of one entity instance.
    class InstanceView
    {
        const P21File *f_{};
        std::uint32_t inst_{};

    public:
        InstanceView(const P21File *f, std::uint32_t inst) : f_{f}, inst_{inst} {}
        [[nodiscard]] std::uint64_t id() const;
        [[nodiscard]] bool is_complex() const;
        [[nodiscard]] std::size_t line() const;
        /// Type of a simple instance; first partial type of a complex one.
        [[nodiscard]] std::string_view type() const;
        [[nodiscard]] std::size_t part_count() const;
        [[nodiscard]] PartView part(std::size_t i) const;
        /// True if the instance is of this type, or has it as a partial type.
        [[nodiscard]] bool has_type(std::string_view type) const;
        /// Partial entity of this type, if any.
        [[nodiscard]] std::optional<PartView> part(std::string_view type) const;
        /// Arguments of a simple instance (of the first part of a complex one).
        [[nodiscard]] ListView args() const;
        /// Arguments of the partial entity `type` (throws if absent).
        [[nodiscard]] ListView args(std::string_view type) const;
    };

    // =========================================================================
    // File
    // =========================================================================

    /**
     * @brief Content of a Part 21 exchange structure: header entities and data
     * instances, indexed by their #id and by type name.
     */
    class P21File
    {
        friend class ValueView;
        friend class ListView;
        friend class PartView;
        friend class InstanceView;
        friend class P21Parser;

        std::vector<detail::Node> nodes_;
        std::vector<std::uint32_t> items_; ///< list children, contiguous per list
        std::string chars_;                ///< decoded strings, enums, binaries
        std::vector<detail::Part> parts_;
        std::vector<detail::Instance> instances_; ///< data section instances
        std::vector<detail::Instance> header_;    ///< header entities (id 0)
        std::vector<std::string> types_;          ///< interned type names
        std::unordered_map<std::string, std::uint32_t> type_ids_;
        std::vector<std::vector<std::uint32_t>> by_type_; ///< type id -> instances (every partial type of a complex one)
        std::vector<std::uint32_t> dense_;                ///< #id -> instance index, when ids are dense enough
        std::unordered_map<std::uint64_t, std::uint32_t> sparse_; ///< otherwise
        std::vector<std::string> schemas_;

        static constexpr std::uint32_t npos = std::numeric_limits<std::uint32_t>::max();

        std::string_view chars(std::uint32_t off, std::uint32_t n) const { return std::string_view{chars_}.substr(off, n); }

        std::uint32_t find_index(std::uint64_t id) const
        {
            if (!dense_.empty())
                return id < dense_.size() ? dense_[id] : npos;
            auto it = sparse_.find(id);
            return it == sparse_.end() ? npos : it->second;
        }

        std::optional<std::uint32_t> type_id(std::string_view type) const
        {
            std::string up{type};
            std::ranges::transform(up, up.begin(), [](unsigned char c) { return char(std::toupper(c)); });
            auto it = type_ids_.find(up);
            if (it == type_ids_.end())
                return std::nullopt;
            return it->second;
        }

    public:
        /// Number of data instances.
        [[nodiscard]] std::size_t size() const noexcept { return instances_.size(); }
        [[nodiscard]] bool contains(std::uint64_t id) const { return find_index(id) != npos; }

        /// Instance #id, nullopt if absent.
        [[nodiscard]] std::optional<InstanceView> find(std::uint64_t id) const
        {
            const auto i = find_index(id);
            if (i == npos)
                return std::nullopt;
            return InstanceView{this, i};
        }

        /// Instance #id; throws P21AccessError if absent (dangling reference).
        [[nodiscard]] InstanceView instance(std::uint64_t id) const
        {
            if (auto v = find(id))
                return *v;
            throw P21AccessError("reference to missing instance #" + std::to_string(id));
        }

        /// Every data instance, in file order.
        [[nodiscard]] auto instances() const
        {
            return std::views::iota(std::uint32_t{0}, std::uint32_t(instances_.size())) |
                   std::views::transform([this](std::uint32_t i) { return InstanceView{this, i}; });
        }

        /// Instances of a type (a complex instance is listed under each of its partial types), in file order.
        [[nodiscard]] std::vector<InstanceView> instances_of(std::string_view type) const
        {
            std::vector<InstanceView> r;
            if (auto t = type_id(type))
                for (auto i : by_type_[*t])
                    r.emplace_back(this, i);
            return r;
        }

        /// Header entities (FILE_DESCRIPTION, FILE_NAME, FILE_SCHEMA, ...), in file order.
        [[nodiscard]] std::vector<InstanceView> header() const
        {
            std::vector<InstanceView> r;
            for (std::uint32_t i{}; i < header_.size(); ++i)
                r.emplace_back(this, npos - 1 - i); // header instances are addressed from the top of the index range
            return r;
        }

        /// Schema names of FILE_SCHEMA, e.g. "AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }".
        [[nodiscard]] std::span<const std::string> schemas() const noexcept { return schemas_; }

    private:
        const detail::Instance &inst(std::uint32_t i) const
        {
            return i >= instances_.size() ? header_[npos - 1 - i] : instances_[i];
        }
    };

    // ---- view implementations ---------------------------------------------------

    inline ValueKind ValueView::kind() const { return f_->nodes_[node_].kind; }

    namespace detail
    {
        inline const char *kind_name(ValueKind k)
        {
            switch (k)
            {
            case ValueKind::Reference: return "reference";
            case ValueKind::Integer: return "integer";
            case ValueKind::Real: return "real";
            case ValueKind::String: return "string";
            case ValueKind::Enum: return "enumeration";
            case ValueKind::Binary: return "binary";
            case ValueKind::Omitted: return "$";
            case ValueKind::Derived: return "*";
            case ValueKind::List: return "list";
            case ValueKind::Typed: return "typed value";
            }
            std::unreachable();
        }

        inline void expect(ValueKind got, ValueKind want)
        {
            if (got != want)
                throw P21AccessError(std::string("expected ") + kind_name(want) + ", got " + kind_name(got));
        }
    } // namespace detail

    inline std::uint64_t ValueView::as_ref() const { detail::expect(kind(), ValueKind::Reference); return f_->nodes_[node_].ref; }
    inline std::int64_t ValueView::as_int() const { detail::expect(kind(), ValueKind::Integer); return f_->nodes_[node_].i; }
    inline double ValueView::as_real() const
    {
        if (kind() == ValueKind::Integer)
            return double(f_->nodes_[node_].i);
        detail::expect(kind(), ValueKind::Real);
        return f_->nodes_[node_].r;
    }
    inline std::string_view ValueView::as_string() const
    {
        detail::expect(kind(), ValueKind::String);
        const auto &n = f_->nodes_[node_];
        return f_->chars(n.a, n.n);
    }
    inline std::string_view ValueView::as_enum() const
    {
        detail::expect(kind(), ValueKind::Enum);
        const auto &n = f_->nodes_[node_];
        return f_->chars(n.a, n.n);
    }
    inline bool ValueView::as_bool() const
    {
        const auto e = as_enum();
        if (e == "T")
            return true;
        if (e == "F")
            return false;
        throw P21AccessError("expected .T. or .F., got ." + std::string(e) + ".");
    }
    inline std::optional<bool> ValueView::as_logical() const
    {
        if (as_enum() == "U")
            return std::nullopt;
        return as_bool();
    }
    inline std::string_view ValueView::as_binary() const
    {
        detail::expect(kind(), ValueKind::Binary);
        const auto &n = f_->nodes_[node_];
        return f_->chars(n.a, n.n);
    }
    inline ListView ValueView::as_list() const
    {
        detail::expect(kind(), ValueKind::List);
        const auto &n = f_->nodes_[node_];
        return {f_, n.a, n.n};
    }
    inline std::string_view ValueView::typed_name() const
    {
        detail::expect(kind(), ValueKind::Typed);
        return f_->types_[f_->nodes_[node_].a];
    }
    inline ValueView ValueView::typed_value() const
    {
        detail::expect(kind(), ValueKind::Typed);
        return {f_, f_->nodes_[node_].n};
    }

    inline ValueView ListView::operator[](std::size_t i) const { return {f_, f_->items_[first_ + i]}; }

    inline std::string_view PartView::type() const { return f_->types_[f_->parts_[part_].type]; }
    inline ListView PartView::args() const
    {
        const auto &n = f_->nodes_[f_->parts_[part_].args];
        return {f_, n.a, n.n};
    }

    inline std::uint64_t InstanceView::id() const { return f_->inst(inst_).id; }
    inline bool InstanceView::is_complex() const { return f_->inst(inst_).complex; }
    inline std::size_t InstanceView::line() const { return f_->inst(inst_).line; }
    inline std::size_t InstanceView::part_count() const { return f_->inst(inst_).n_parts; }
    inline PartView InstanceView::part(std::size_t i) const
    {
        if (i >= part_count())
            throw P21AccessError("part index out of range");
        return {f_, f_->inst(inst_).first_part + std::uint32_t(i)};
    }
    inline std::string_view InstanceView::type() const { return part(0).type(); }
    inline std::optional<PartView> InstanceView::part(std::string_view type) const
    {
        const auto t = f_->type_id(type);
        if (!t)
            return std::nullopt;
        const auto &in = f_->inst(inst_);
        for (std::uint32_t k{}; k < in.n_parts; ++k)
            if (f_->parts_[in.first_part + k].type == *t)
                return PartView{f_, in.first_part + k};
        return std::nullopt;
    }
    inline bool InstanceView::has_type(std::string_view type) const { return part(type).has_value(); }
    inline ListView InstanceView::args() const { return part(0).args(); }
    inline ListView InstanceView::args(std::string_view type) const
    {
        if (auto p = part(type))
            return p->args();
        throw P21AccessError("instance #" + std::to_string(id()) + " has no partial type " + std::string(type));
    }

    // =========================================================================
    // Lexer
    // =========================================================================

    namespace detail
    {
        enum class Tok : std::uint8_t
        {
            Keyword, Ref, Integer, Real, String, Enum, Binary,
            LParen, RParen, Comma, Semicolon, Equal, Dollar, Star, End
        };

        /// Locale-independent conversion of a decimal real.
        inline bool parse_real(std::string_view s, double &out)
        {
#if defined(_LIBCPP_VERSION)
            // libc++ (conda-forge on macOS) marks floating-point from_chars unavailable:
            // use strtod_l with the "C" locale, as exact and locale independent.
            static const locale_t c_locale = newlocale(LC_ALL_MASK, "C", locale_t(0));
            char buf[128];
            if (s.size() >= sizeof(buf))
                return false;
            std::ranges::copy(s, buf);
            buf[s.size()] = '\0';
            char *end = nullptr;
            out = strtod_l(buf, &end, c_locale);
            return end == buf + s.size();
#else
            auto [p, ec] = std::from_chars(s.data(), s.data() + s.size(), out);
            return ec == std::errc{} && p == s.data() + s.size();
#endif
        }

        inline void append_utf8(std::string &out, std::uint32_t cp)
        {
            if (cp < 0x80)
                out += char(cp);
            else if (cp < 0x800)
            {
                out += char(0xC0 | (cp >> 6));
                out += char(0x80 | (cp & 0x3F));
            }
            else if (cp < 0x10000)
            {
                out += char(0xE0 | (cp >> 12));
                out += char(0x80 | ((cp >> 6) & 0x3F));
                out += char(0x80 | (cp & 0x3F));
            }
            else
            {
                out += char(0xF0 | (cp >> 18));
                out += char(0x80 | ((cp >> 12) & 0x3F));
                out += char(0x80 | ((cp >> 6) & 0x3F));
                out += char(0x80 | (cp & 0x3F));
            }
        }

        inline int hex(char c)
        {
            if (c >= '0' && c <= '9') return c - '0';
            if (c >= 'A' && c <= 'F') return c - 'A' + 10;
            if (c >= 'a' && c <= 'f') return c - 'a' + 10;
            return -1;
        }

        /**
         * Decodes the body of a Part 21 string (between the quotes) to UTF-8:
         * '' -> ', \\ -> \, \S\c -> ISO 8859 upper half, \X\hh -> 8-bit code,
         * \X2\hhhh...\X0\ -> UTF-16 (surrogates combined), \X4\hhhhhhhh...\X0\ -> UTF-32,
         * \Pc\ code page switches ignored, end-of-line characters ignored.
         * Returns false on a malformed escape.
         */
        inline bool decode_string(std::string_view raw, std::string &out)
        {
            for (std::size_t i = 0; i < raw.size();)
            {
                const char c = raw[i];
                if (c == '\'' && i + 1 < raw.size() && raw[i + 1] == '\'')
                {
                    out += '\'';
                    i += 2;
                }
                else if (c == '\r' || c == '\n')
                    ++i;
                else if (c != '\\')
                {
                    out += c;
                    ++i;
                }
                else
                {
                    const auto rest = raw.substr(i);
                    if (rest.starts_with("\\\\"))
                    {
                        out += '\\';
                        i += 2;
                    }
                    else if (rest.starts_with("\\S\\") && rest.size() >= 4)
                    {
                        append_utf8(out, std::uint32_t(static_cast<unsigned char>(rest[3])) + 128u);
                        i += 4;
                    }
                    else if (rest.starts_with("\\X\\") && rest.size() >= 5)
                    {
                        const int h = hex(rest[3]), l = hex(rest[4]);
                        if (h < 0 || l < 0)
                            return false;
                        append_utf8(out, std::uint32_t(h * 16 + l));
                        i += 5;
                    }
                    else if (rest.starts_with("\\X2\\") || rest.starts_with("\\X4\\"))
                    {
                        const std::size_t w = rest[2] == '2' ? 4 : 8;
                        std::size_t j = 4;
                        std::uint32_t pending_high = 0;
                        while (j < rest.size() && rest[j] != '\\')
                        {
                            if (j + w > rest.size())
                                return false;
                            std::uint32_t cp = 0;
                            for (std::size_t k = 0; k < w; ++k)
                            {
                                const int d = hex(rest[j + k]);
                                if (d < 0)
                                    return false;
                                cp = cp * 16 + std::uint32_t(d);
                            }
                            j += w;
                            if (w == 4 && cp >= 0xD800 && cp < 0xDC00)
                                pending_high = cp;
                            else if (w == 4 && cp >= 0xDC00 && cp < 0xE000 && pending_high)
                            {
                                append_utf8(out, 0x10000 + ((pending_high - 0xD800) << 10) + (cp - 0xDC00));
                                pending_high = 0;
                            }
                            else
                                append_utf8(out, cp);
                        }
                        if (!rest.substr(j).starts_with("\\X0\\"))
                            return false;
                        i += j + 4;
                    }
                    else if (rest.size() >= 4 && rest[1] == 'P' && rest[3] == '\\')
                        i += 4; // \PA\ ... code page switch: ignored, \S\ maps to the upper half of ISO 8859-1
                    else
                        return false;
                }
            }
            return true;
        }

        class Lexer
        {
            std::string_view s_;
            std::size_t pos_{0}, line_{1}, line_start_{0};

        public:
            Tok tok{Tok::End};
            std::string_view text; ///< token text (keyword, number, string body, enum name, binary digits)
            std::size_t tok_line{1}, tok_col{1};

            explicit Lexer(std::string_view s) : s_{s} {}

            [[nodiscard]] P21Error error(std::string msg) const { return {tok_line, tok_col, std::move(msg)}; }

            /// Skips raw text up to and including the next "ENDSEC" keyword (sections whose
            /// content is not made of entity instances, e.g. ANCHOR with its <name> tokens).
            bool skip_to_endsec()
            {
                const auto e = s_.find("ENDSEC", pos_);
                if (e == std::string_view::npos)
                    return false;
                for (auto k = pos_; k < e; ++k)
                    if (s_[k] == '\n')
                        ++line_, line_start_ = k + 1;
                pos_ = e + 6;
                return true;
            }

            /// Advances to the next token; returns an error for an invalid character or an unterminated token.
            std::optional<P21Error> next()
            {
                // whitespace and comments
                for (;;)
                {
                    while (pos_ < s_.size() && std::isspace(static_cast<unsigned char>(s_[pos_])))
                    {
                        if (s_[pos_] == '\n')
                            ++line_, line_start_ = pos_ + 1;
                        ++pos_;
                    }
                    if (pos_ + 1 < s_.size() && s_[pos_] == '/' && s_[pos_ + 1] == '*')
                    {
                        const auto e = s_.find("*/", pos_ + 2);
                        tok_line = line_, tok_col = pos_ - line_start_ + 1;
                        if (e == std::string_view::npos)
                            return error("unterminated comment");
                        for (auto k = pos_; k < e; ++k)
                            if (s_[k] == '\n')
                                ++line_, line_start_ = k + 1;
                        pos_ = e + 2;
                        continue;
                    }
                    break;
                }
                tok_line = line_, tok_col = pos_ - line_start_ + 1;
                if (pos_ >= s_.size())
                {
                    tok = Tok::End;
                    return std::nullopt;
                }
                const char c = s_[pos_];
                const auto start = pos_;
                auto single = [&](Tok t) { tok = t; text = s_.substr(pos_, 1); ++pos_; return std::nullopt; };
                switch (c)
                {
                case '(': return single(Tok::LParen);
                case ')': return single(Tok::RParen);
                case ',': return single(Tok::Comma);
                case ';': return single(Tok::Semicolon);
                case '=': return single(Tok::Equal);
                case '$': return single(Tok::Dollar);
                case '*': return single(Tok::Star);
                default: break;
                }
                if (c == '#')
                {
                    ++pos_;
                    while (pos_ < s_.size() && std::isdigit(static_cast<unsigned char>(s_[pos_])))
                        ++pos_;
                    if (pos_ == start + 1)
                        return error("'#' not followed by an instance number");
                    tok = Tok::Ref, text = s_.substr(start + 1, pos_ - start - 1);
                    return std::nullopt;
                }
                if (c == '\'')
                {
                    ++pos_;
                    for (;;)
                    {
                        if (pos_ >= s_.size())
                            return error("unterminated string");
                        if (s_[pos_] == '\'')
                        {
                            if (pos_ + 1 < s_.size() && s_[pos_ + 1] == '\'')
                            {
                                pos_ += 2;
                                continue;
                            }
                            break;
                        }
                        if (s_[pos_] == '\n')
                            ++line_, line_start_ = pos_ + 1;
                        ++pos_;
                    }
                    tok = Tok::String, text = s_.substr(start + 1, pos_ - start - 1);
                    ++pos_;
                    return std::nullopt;
                }
                if (c == '"')
                {
                    const auto e = s_.find('"', pos_ + 1);
                    if (e == std::string_view::npos)
                        return error("unterminated binary");
                    tok = Tok::Binary, text = s_.substr(start + 1, e - start - 1);
                    pos_ = e + 1;
                    return std::nullopt;
                }
                if (c == '.' && pos_ + 1 < s_.size() && (std::isalpha(static_cast<unsigned char>(s_[pos_ + 1])) || s_[pos_ + 1] == '_'))
                {
                    const auto e = s_.find('.', pos_ + 1);
                    if (e == std::string_view::npos)
                        return error("unterminated enumeration");
                    tok = Tok::Enum, text = s_.substr(start + 1, e - start - 1);
                    pos_ = e + 1;
                    return std::nullopt;
                }
                if (c == '+' || c == '-' || c == '.' || std::isdigit(static_cast<unsigned char>(c)))
                {
                    bool real = false;
                    ++pos_;
                    while (pos_ < s_.size())
                    {
                        const char d = s_[pos_];
                        if (std::isdigit(static_cast<unsigned char>(d)))
                            ++pos_;
                        else if (d == '.')
                            real = true, ++pos_;
                        else if ((d == 'E' || d == 'e'))
                        {
                            real = true, ++pos_;
                            if (pos_ < s_.size() && (s_[pos_] == '+' || s_[pos_] == '-'))
                                ++pos_;
                        }
                        else
                            break;
                    }
                    tok = real || c == '.' ? Tok::Real : Tok::Integer;
                    text = s_.substr(start, pos_ - start);
                    return std::nullopt;
                }
                if (std::isalpha(static_cast<unsigned char>(c)) || c == '_' || c == '!')
                {
                    ++pos_;
                    while (pos_ < s_.size() && (std::isalnum(static_cast<unsigned char>(s_[pos_])) || s_[pos_] == '_' || s_[pos_] == '-'))
                        ++pos_;
                    tok = Tok::Keyword, text = s_.substr(start, pos_ - start);
                    return std::nullopt;
                }
                return error(std::string("unexpected character '") + c + "'");
            }
        };
    } // namespace detail

    // =========================================================================
    // Parser
    // =========================================================================

    class P21Parser
    {
        detail::Lexer lx_;
        P21File f_;
        std::vector<std::uint32_t> scratch_;
        std::vector<std::pair<std::uint64_t, std::uint32_t>> ids_; // (#id, instance index)

        using Err = std::unexpected<P21Error>;

        std::optional<P21Error> advance() { return lx_.next(); }

        std::optional<P21Error> expect(detail::Tok t, const char *what)
        {
            if (lx_.tok != t)
                return lx_.error(std::string("expected ") + what);
            return advance();
        }

        std::uint32_t intern(std::string_view name)
        {
            std::string up{name};
            std::ranges::transform(up, up.begin(), [](unsigned char c) { return char(std::toupper(c)); });
            auto [it, inserted] = f_.type_ids_.try_emplace(up, std::uint32_t(f_.types_.size()));
            if (inserted)
            {
                f_.types_.push_back(up);
                f_.by_type_.emplace_back();
            }
            return it->second;
        }

        std::uint32_t add_node(detail::Node n)
        {
            f_.nodes_.push_back(n);
            return std::uint32_t(f_.nodes_.size() - 1);
        }

        std::uint32_t add_chars(std::string_view s)
        {
            const auto off = std::uint32_t(f_.chars_.size());
            f_.chars_ += s;
            return off;
        }

        /// '(' params ')' -> List node; the current token is '('.
        std::expected<std::uint32_t, P21Error> parse_list()
        {
            if (auto e = expect(detail::Tok::LParen, "'('"))
                return Err(*e);
            const auto mark = scratch_.size();
            if (lx_.tok != detail::Tok::RParen)
            {
                for (;;)
                {
                    auto v = parse_value();
                    if (!v)
                        return Err(v.error());
                    scratch_.push_back(*v);
                    if (lx_.tok == detail::Tok::Comma)
                    {
                        if (auto e = advance())
                            return Err(*e);
                        continue;
                    }
                    break;
                }
            }
            if (auto e = expect(detail::Tok::RParen, "',' or ')'"))
                return Err(*e);
            detail::Node n;
            n.kind = ValueKind::List;
            n.a = std::uint32_t(f_.items_.size());
            n.n = std::uint32_t(scratch_.size() - mark);
            f_.items_.insert(f_.items_.end(), scratch_.begin() + std::ptrdiff_t(mark), scratch_.end());
            scratch_.resize(mark);
            return add_node(n);
        }

        std::expected<std::uint32_t, P21Error> parse_value()
        {
            using detail::Tok;
            detail::Node n;
            switch (lx_.tok)
            {
            case Tok::Ref:
            {
                n.kind = ValueKind::Reference;
                std::from_chars(lx_.text.data(), lx_.text.data() + lx_.text.size(), n.ref);
                break;
            }
            case Tok::Integer:
            {
                n.kind = ValueKind::Integer;
                const auto *b = lx_.text.data() + (lx_.text.front() == '+' ? 1 : 0);
                auto [p, ec] = std::from_chars(b, lx_.text.data() + lx_.text.size(), n.i);
                if (ec != std::errc{} || p != lx_.text.data() + lx_.text.size())
                    return Err(lx_.error("invalid integer '" + std::string(lx_.text) + "'"));
                break;
            }
            case Tok::Real:
            {
                n.kind = ValueKind::Real;
                if (!detail::parse_real(lx_.text, n.r))
                    return Err(lx_.error("invalid real '" + std::string(lx_.text) + "'"));
                break;
            }
            case Tok::String:
            {
                n.kind = ValueKind::String;
                std::string decoded;
                if (!detail::decode_string(lx_.text, decoded))
                    return Err(lx_.error("malformed escape in string"));
                n.a = add_chars(decoded);
                n.n = std::uint32_t(decoded.size());
                break;
            }
            case Tok::Enum:
            {
                n.kind = ValueKind::Enum;
                std::string up{lx_.text};
                std::ranges::transform(up, up.begin(), [](unsigned char c) { return char(std::toupper(c)); });
                n.a = add_chars(up);
                n.n = std::uint32_t(up.size());
                break;
            }
            case Tok::Binary:
                n.kind = ValueKind::Binary;
                n.a = add_chars(lx_.text);
                n.n = std::uint32_t(lx_.text.size());
                break;
            case Tok::Dollar:
                n.kind = ValueKind::Omitted;
                break;
            case Tok::Star:
                n.kind = ValueKind::Derived;
                break;
            case Tok::LParen:
                return parse_list();
            case Tok::Keyword:
            {
                // typed parameter: KEYWORD '(' value ')'
                const auto type = intern(lx_.text);
                if (auto e = advance())
                    return Err(*e);
                if (auto e = expect(Tok::LParen, "'(' after a type name"))
                    return Err(*e);
                auto child = parse_value();
                if (!child)
                    return Err(child.error());
                if (auto e = expect(Tok::RParen, "')' closing a typed parameter"))
                    return Err(*e);
                n.kind = ValueKind::Typed;
                n.a = type;
                n.n = *child;
                return add_node(n);
            }
            default:
                return Err(lx_.error("expected a parameter value"));
            }
            const auto id = add_node(n);
            if (auto e = advance())
                return Err(*e);
            return id;
        }

        /// KEYWORD '(' params ')' -> one part.
        std::expected<detail::Part, P21Error> parse_part()
        {
            if (lx_.tok != detail::Tok::Keyword)
                return Err(lx_.error("expected an entity name"));
            const auto type = intern(lx_.text);
            if (auto e = advance())
                return Err(*e);
            if (lx_.tok != detail::Tok::LParen)
                return Err(lx_.error("expected '(' after the entity name"));
            auto args = parse_list();
            if (!args)
                return Err(args.error());
            return detail::Part{type, *args};
        }

        /// Simple or complex entity body, then ';'. Appends the instance to `into`.
        std::optional<P21Error> parse_entity(std::uint64_t id, std::vector<detail::Instance> &into)
        {
            detail::Instance inst;
            inst.id = id;
            inst.line = std::uint32_t(lx_.tok_line);
            inst.first_part = std::uint32_t(f_.parts_.size());
            if (lx_.tok == detail::Tok::LParen)
            {
                inst.complex = true;
                if (auto e = advance())
                    return e;
                while (lx_.tok == detail::Tok::Keyword)
                {
                    auto p = parse_part();
                    if (!p)
                        return p.error();
                    f_.parts_.push_back(*p);
                }
                if (auto e = expect(detail::Tok::RParen, "')' closing a complex instance"))
                    return e;
            }
            else
            {
                auto p = parse_part();
                if (!p)
                    return p.error();
                f_.parts_.push_back(*p);
            }
            inst.n_parts = std::uint32_t(f_.parts_.size() - inst.first_part);
            if (inst.n_parts == 0)
                return lx_.error("empty complex instance");
            if (auto e = expect(detail::Tok::Semicolon, "';' ending the instance"))
                return e;
            into.push_back(inst);
            return std::nullopt;
        }

        std::optional<P21Error> parse_header()
        {
            if (auto e = expect(detail::Tok::Semicolon, "';' after HEADER"))
                return e;
            while (!(lx_.tok == detail::Tok::Keyword && lx_.text == "ENDSEC"))
            {
                if (lx_.tok == detail::Tok::End)
                    return lx_.error("missing ENDSEC of the HEADER section");
                if (auto e = parse_entity(0, f_.header_))
                    return e;
            }
            if (auto e = advance())
                return e;
            return expect(detail::Tok::Semicolon, "';' after ENDSEC");
        }

        std::optional<P21Error> parse_data()
        {
            if (lx_.tok == detail::Tok::LParen) // DATA('name', (schema)) of edition 3: ignored
            {
                auto l = parse_list();
                if (!l)
                    return l.error();
            }
            if (auto e = expect(detail::Tok::Semicolon, "';' after DATA"))
                return e;
            while (!(lx_.tok == detail::Tok::Keyword && lx_.text == "ENDSEC"))
            {
                if (lx_.tok != detail::Tok::Ref)
                    return lx_.error("expected an instance '#n=' or ENDSEC");
                std::uint64_t id{};
                std::from_chars(lx_.text.data(), lx_.text.data() + lx_.text.size(), id);
                if (auto e = advance())
                    return e;
                if (auto e = expect(detail::Tok::Equal, "'=' after the instance number"))
                    return e;
                if (auto e = parse_entity(id, f_.instances_))
                    return e;
                ids_.emplace_back(id, std::uint32_t(f_.instances_.size() - 1));
            }
            if (auto e = advance())
                return e;
            return expect(detail::Tok::Semicolon, "';' after ENDSEC");
        }

        /// Sections of edition 3 that do not carry entity instances (ANCHOR, REFERENCE, SIGNATURE): skipped.
        std::optional<P21Error> skip_section()
        {
            // the current token is the one after the section keyword (usually ';'): it may already
            // be ENDSEC for an empty section
            if (!(lx_.tok == detail::Tok::Keyword && lx_.text == "ENDSEC") && !lx_.skip_to_endsec())
                return lx_.error("unterminated section");
            if (auto e = advance())
                return e;
            if (lx_.tok == detail::Tok::Keyword && lx_.text == "ENDSEC") // empty section: ENDSEC was the current token
                if (auto e = advance())
                    return e;
            return expect(detail::Tok::Semicolon, "';' after ENDSEC");
        }

        std::optional<P21Error> finish()
        {
            // #id index: dense table when ids are reasonably compact, hash map otherwise
            std::uint64_t max_id = 0;
            for (const auto &[id, i] : ids_)
                max_id = std::max(max_id, id);
            if (max_id <= 4 * ids_.size() + 1024)
            {
                f_.dense_.assign(std::size_t(max_id) + 1, P21File::npos);
                for (const auto &[id, i] : ids_)
                {
                    if (f_.dense_[std::size_t(id)] != P21File::npos)
                        return P21Error{f_.instances_[i].line, 1, "instance #" + std::to_string(id) + " defined twice"};
                    f_.dense_[std::size_t(id)] = i;
                }
            }
            else
                for (const auto &[id, i] : ids_)
                    if (!f_.sparse_.try_emplace(id, i).second)
                        return P21Error{f_.instances_[i].line, 1, "instance #" + std::to_string(id) + " defined twice"};
            // type index, every partial type of a complex instance
            for (std::uint32_t i{}; i < f_.instances_.size(); ++i)
            {
                const auto &in = f_.instances_[i];
                for (std::uint32_t k{}; k < in.n_parts; ++k)
                    f_.by_type_[f_.parts_[in.first_part + k].type].push_back(i);
            }
            // schemas of FILE_SCHEMA(('...'))
            for (const auto &h : f_.header_)
                if (f_.types_[f_.parts_[h.first_part].type] == "FILE_SCHEMA")
                {
                    InstanceView v{&f_, P21File::npos - 1 - std::uint32_t(&h - f_.header_.data())};
                    const auto args = v.args();
                    if (!args.empty() && args[0].is(ValueKind::List))
                        for (auto s : args[0].as_list())
                            if (s.is(ValueKind::String))
                                f_.schemas_.emplace_back(s.as_string());
                }
            return std::nullopt;
        }

    public:
        explicit P21Parser(std::string_view text) : lx_{text} {}

        std::expected<P21File, P21Error> run()
        {
            auto fail = [](P21Error e) { return std::expected<P21File, P21Error>(std::unexpect, std::move(e)); };
            if (auto e = advance())
                return fail(*e);
            if (!(lx_.tok == detail::Tok::Keyword && lx_.text == "ISO-10303-21"))
                return fail(lx_.error("not a STEP Part 21 file: expected ISO-10303-21"));
            if (auto e = advance())
                return fail(*e);
            if (auto e = expect(detail::Tok::Semicolon, "';' after ISO-10303-21"))
                return fail(*e);
            bool header_seen = false;
            for (;;)
            {
                if (lx_.tok == detail::Tok::End)
                    return fail(lx_.error("missing END-ISO-10303-21"));
                if (lx_.tok != detail::Tok::Keyword)
                    return fail(lx_.error("expected a section keyword"));
                const auto kw = lx_.text;
                if (kw == "END-ISO-10303-21")
                    break;
                if (auto e = advance())
                    return fail(*e);
                std::optional<P21Error> err;
                if (kw == "HEADER")
                    err = parse_header(), header_seen = true;
                else if (kw == "DATA")
                    err = parse_data();
                else if (kw == "ANCHOR" || kw == "REFERENCE" || kw == "SIGNATURE")
                    err = skip_section();
                else
                    return fail(lx_.error("unknown section " + std::string(kw)));
                if (err)
                    return fail(*err);
            }
            if (!header_seen)
                return fail(lx_.error("missing HEADER section"));
            if (auto e = finish())
                return fail(*e);
            return std::move(f_);
        }
    };

    /// Parses the text of a Part 21 exchange structure.
    [[nodiscard]] inline auto parse_p21(std::string_view text) -> std::expected<P21File, P21Error>
    {
        return P21Parser{text}.run();
    }

    /// Reads and parses a Part 21 file.
    [[nodiscard]] inline auto read_p21(const std::filesystem::path &file) -> std::expected<P21File, P21Error>
    {
        std::ifstream in(file, std::ios::binary);
        if (!in)
            return std::unexpected(P21Error{0, 0, "cannot open " + file.string()});
        std::string text{std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
        return parse_p21(text);
    }

} // namespace gbs::step
