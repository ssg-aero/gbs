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
 * into a P21File and stay valid as long as it lives, and so do their iterators.
 *
 * The lexer only uses ASCII rules: the global C locale has no effect.
 */

#include <algorithm>
#include <array>
#include <bit>
#include <charconv>
#include <compare>
#include <cstdint>
#include <cstdlib>
#include <expected>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iterator>
#include <limits>
#include <optional>
#include <ranges>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <unordered_map>
#include <utility>
#include <vector>

#if defined(_LIBCPP_VERSION)
#include <cerrno>
#include <clocale>
#include <cmath>
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

    /// Raised by the views when a value is not of the requested kind, a list index is out of range or a
    /// reference dangles (content error, not syntax).
    class P21AccessError : public std::runtime_error
    {
    public:
        using std::runtime_error::runtime_error;
    };

    /// Deepest nesting of lists and typed parameters the parser accepts (real files nest a few levels);
    /// deeper is a syntax error rather than a stack overflow.
    inline constexpr std::size_t max_nesting = 256;

    class P21File;

    namespace detail
    {
        /// One parameter value: 12 bytes, the payload of every kind fits in two 32-bit words.
        struct Node
        {
            /// String/Enum/Binary: offset in chars_ and length; List: first item and count;
            /// Typed: type id and child node; Integer/Real/Reference: the 64-bit value (get/set).
            std::uint32_t a{}, n{};
            ValueKind kind{ValueKind::Omitted};

            template <class T>
            [[nodiscard]] T get() const noexcept { return std::bit_cast<T>(std::array<std::uint32_t, 2>{a, n}); }
            template <class T>
            void set(T v) noexcept
            {
                const auto w = std::bit_cast<std::array<std::uint32_t, 2>>(v);
                a = w[0], n = w[1];
            }
        };
        static_assert(sizeof(Node) == 12);

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

        // ASCII character classes: unlike <cctype>, independent of the global C locale
        constexpr bool is_digit(char c) noexcept { return c >= '0' && c <= '9'; }
        constexpr bool is_lower(char c) noexcept { return c >= 'a' && c <= 'z'; }
        constexpr bool is_alpha(char c) noexcept { return (c >= 'A' && c <= 'Z') || is_lower(c); }
        constexpr bool is_space(char c) noexcept { return c == ' ' || (c >= '\t' && c <= '\r'); }
        constexpr char to_upper(char c) noexcept { return is_lower(c) ? char(c - 'a' + 'A') : c; }
        inline std::string to_upper(std::string_view s)
        {
            std::string up(s.size(), '\0');
            std::ranges::transform(s, up.begin(), [](char c) { return to_upper(c); });
            return up;
        }

        /// Hash of the type name table, for lookups by std::string_view without a copy.
        struct TypeHash
        {
            using is_transparent = void;
            std::size_t operator()(std::string_view s) const noexcept { return std::hash<std::string_view>{}(s); }
        };
    } // namespace detail

    class ListView;
    class InstanceView;

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

    namespace detail
    {
        /// Value at a position of the item table (list children).
        struct ListItem
        {
            ValueView operator()(const P21File *f, std::uint32_t pos) const;
        };

        /**
         * Random-access iterator over consecutive indices of a P21File, dereferenced by Get to a view.
         * It holds the file and the index, not a pointer to the range that made it: an iterator of a
         * temporary range stays valid as long as the file lives (the ranges are borrowed ranges).
         */
        template <class Get>
        class IndexIterator
        {
            const P21File *f_{};
            std::uint32_t i_{};

        public:
            using value_type = decltype(Get{}(nullptr, 0u));
            using reference = value_type;
            using pointer = void;
            using difference_type = std::ptrdiff_t;
            using iterator_concept = std::random_access_iterator_tag;
            using iterator_category = std::input_iterator_tag; // operator* returns a value

            IndexIterator() = default;
            IndexIterator(const P21File *f, std::uint32_t i) : f_{f}, i_{i} {}

            value_type operator*() const { return Get{}(f_, i_); }
            value_type operator[](difference_type k) const { return *(*this + k); }
            IndexIterator &operator++() { ++i_; return *this; }
            IndexIterator operator++(int) { auto t = *this; ++i_; return t; }
            IndexIterator &operator--() { --i_; return *this; }
            IndexIterator operator--(int) { auto t = *this; --i_; return t; }
            IndexIterator &operator+=(difference_type k) { i_ = std::uint32_t(difference_type(i_) + k); return *this; }
            IndexIterator &operator-=(difference_type k) { return *this += -k; }
            friend IndexIterator operator+(IndexIterator it, difference_type k) { return it += k; }
            friend IndexIterator operator+(difference_type k, IndexIterator it) { return it += k; }
            friend IndexIterator operator-(IndexIterator it, difference_type k) { return it -= k; }
            friend difference_type operator-(const IndexIterator &a, const IndexIterator &b) { return difference_type(a.i_) - difference_type(b.i_); }
            bool operator==(const IndexIterator &o) const { return i_ == o.i_; }
            auto operator<=>(const IndexIterator &o) const { return i_ <=> o.i_; }
        };
    } // namespace detail

    /// Read-only view of a list of values (a borrowed random-access range).
    class ListView : public std::ranges::view_interface<ListView>
    {
        const P21File *f_{};
        std::uint32_t first_{}, n_{};

    public:
        using iterator = detail::IndexIterator<detail::ListItem>;

        ListView() = default;
        ListView(const P21File *f, std::uint32_t first, std::uint32_t n) : f_{f}, first_{first}, n_{n} {}
        [[nodiscard]] std::size_t size() const noexcept { return n_; }
        [[nodiscard]] bool empty() const noexcept { return n_ == 0; }
        /// Item i; throws P21AccessError if i >= size().
        [[nodiscard]] ValueView operator[](std::size_t i) const;
        /// Same as operator[].
        [[nodiscard]] ValueView at(std::size_t i) const { return (*this)[i]; }

        [[nodiscard]] iterator begin() const { return {f_, first_}; }
        [[nodiscard]] iterator end() const { return {f_, first_ + n_}; }
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

    namespace detail
    {
        /// Data instance of a given index.
        struct InstanceAt
        {
            InstanceView operator()(const P21File *f, std::uint32_t i) const { return {f, i}; }
        };
    } // namespace detail

    /// Data instances of a P21File in file order (a borrowed random-access range).
    class InstanceRange : public std::ranges::view_interface<InstanceRange>
    {
        const P21File *f_{};
        std::uint32_t n_{};

    public:
        using iterator = detail::IndexIterator<detail::InstanceAt>;

        InstanceRange() = default;
        InstanceRange(const P21File *f, std::uint32_t n) : f_{f}, n_{n} {}
        [[nodiscard]] std::size_t size() const noexcept { return n_; }
        [[nodiscard]] iterator begin() const { return {f_, 0}; }
        [[nodiscard]] iterator end() const { return {f_, n_}; }
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
        friend struct detail::ListItem;

        std::vector<detail::Node> nodes_;
        std::vector<std::uint32_t> items_; ///< list children, contiguous per list
        std::string chars_;                ///< decoded strings, enums, binaries
        std::vector<detail::Part> parts_;
        std::vector<detail::Instance> instances_; ///< data section instances
        std::vector<detail::Instance> header_;    ///< header entities (id 0)
        std::vector<std::string> types_;          ///< interned type names (uppercase)
        std::unordered_map<std::string, std::uint32_t, detail::TypeHash, std::equal_to<>> type_ids_;
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
            // names are interned in uppercase: exact lookup, an uppercase copy only for a name with lowercase letters
            auto it = type_ids_.find(type);
            if (it == type_ids_.end() && std::ranges::any_of(type, detail::is_lower))
                it = type_ids_.find(detail::to_upper(type));
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
        [[nodiscard]] InstanceRange instances() const { return {this, std::uint32_t(instances_.size())}; }

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
} // namespace gbs::step

// iterators of a temporary ListView or InstanceRange stay valid: they only refer to the P21File
template <>
inline constexpr bool std::ranges::enable_borrowed_range<gbs::step::ListView> = true;
template <>
inline constexpr bool std::ranges::enable_borrowed_range<gbs::step::InstanceRange> = true;

namespace gbs::step
{
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

    inline std::uint64_t ValueView::as_ref() const { detail::expect(kind(), ValueKind::Reference); return f_->nodes_[node_].get<std::uint64_t>(); }
    inline std::int64_t ValueView::as_int() const { detail::expect(kind(), ValueKind::Integer); return f_->nodes_[node_].get<std::int64_t>(); }
    inline double ValueView::as_real() const
    {
        if (kind() == ValueKind::Integer)
            return double(f_->nodes_[node_].get<std::int64_t>());
        detail::expect(kind(), ValueKind::Real);
        return f_->nodes_[node_].get<double>();
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

    namespace detail
    {
        [[noreturn]] inline void throw_list_index(std::size_t i, std::size_t n)
        {
            throw P21AccessError("list index " + std::to_string(i) + " out of range (size " + std::to_string(n) + ")");
        }

        inline ValueView ListItem::operator()(const P21File *f, std::uint32_t pos) const { return {f, f->items_[pos]}; }
    } // namespace detail

    inline ValueView ListView::operator[](std::size_t i) const
    {
        if (i >= n_) [[unlikely]]
            detail::throw_list_index(i, n_);
        return {f_, f_->items_[first_ + i]};
    }

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

        /**
         * For a real token whose value is beyond the range of double: true if it is tiny (underflow),
         * false if it is huge (overflow). Read from the digits, as std::from_chars leaves the value unset:
         * the sign of the decimal exponent of the first significant digit.
         */
        inline bool real_underflows(std::string_view s)
        {
            if (s.starts_with('+') || s.starts_with('-'))
                s.remove_prefix(1);
            const auto e = std::min(s.find_first_of("Ee"), s.size());
            const auto mantissa = s.substr(0, e);
            std::int64_t exponent = 0; // saturated: beyond 10^9, only its sign matters
            if (e < s.size())
            {
                auto x = s.substr(e + 1);
                const bool negative = x.starts_with('-');
                if (negative || x.starts_with('+'))
                    x.remove_prefix(1);
                for (const char c : x)
                    if (is_digit(c) && exponent < 1'000'000'000)
                        exponent = exponent * 10 + (c - '0');
                if (negative)
                    exponent = -exponent;
            }
            const auto first = mantissa.find_first_of("123456789");
            if (first == std::string_view::npos)
                return true; // zero
            const auto dot = std::min(mantissa.find('.'), mantissa.size());
            const auto magnitude = first < dot ? std::int64_t(dot - first) - 1 : -std::int64_t(first - dot);
            return magnitude + exponent < 0;
        }

        /**
         * Locale-independent conversion of a decimal real, alike on every platform: a leading '+' is
         * accepted, an overflow gives result_out_of_range, an underflow gives 0 (or the denormal), the
         * length is not limited. Returns invalid_argument for a malformed token.
         */
        inline std::errc parse_real(std::string_view s, double &out)
        {
            if (s.starts_with('+'))
            {
                s.remove_prefix(1);
                if (s.starts_with('+') || s.starts_with('-'))
                    return std::errc::invalid_argument;
            }
#if defined(_LIBCPP_VERSION)
            // libc++ (conda-forge on macOS) marks floating-point from_chars unavailable:
            // use strtod_l with the "C" locale, as exact and locale independent.
            static const locale_t c_locale = newlocale(LC_ALL_MASK, "C", locale_t(0));
            char small[64];
            std::string large;
            const char *text = small;
            if (s.size() < sizeof small)
            {
                std::ranges::copy(s, small);
                small[s.size()] = '\0';
            }
            else
                text = (large = s).c_str();
            char *end = nullptr;
            errno = 0;
            out = strtod_l(text, &end, c_locale);
            if (end != text + s.size())
                return std::errc::invalid_argument;
            return errno == ERANGE && std::isinf(out) ? std::errc::result_out_of_range : std::errc{};
#else
            const auto [p, ec] = std::from_chars(s.data(), s.data() + s.size(), out);
            if (p != s.data() + s.size())
                return std::errc::invalid_argument;
            if (ec == std::errc::result_out_of_range && real_underflows(s))
            {
                out = s.starts_with('-') ? -0. : 0.;
                return std::errc{};
            }
            return ec;
#endif
        }

        /// Appends a code point in UTF-8; a surrogate or a value beyond U+10FFFF becomes U+FFFD.
        inline void append_utf8(std::string &out, std::uint32_t cp)
        {
            if ((cp >= 0xD800 && cp < 0xE000) || cp > 0x10FFFF)
                cp = 0xFFFD;
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

        /// Length of the valid UTF-8 sequence that starts s (with a byte >= 0x80), 0 if there is none.
        inline std::size_t utf8_length(std::string_view s)
        {
            const auto byte = [&](std::size_t k) { return std::uint32_t(static_cast<unsigned char>(s[k])); };
            const auto c = byte(0);
            const std::size_t n = c >= 0xC2 && c <= 0xDF ? 2 : c >= 0xE0 && c <= 0xEF ? 3 : c >= 0xF0 && c <= 0xF4 ? 4 : 0;
            if (n == 0 || s.size() < n)
                return 0;
            std::uint32_t cp = c & (0x7Fu >> n);
            for (std::size_t k = 1; k < n; ++k)
            {
                if ((byte(k) & 0xC0) != 0x80)
                    return 0;
                cp = cp << 6 | (byte(k) & 0x3F);
            }
            constexpr std::uint32_t shortest[] = {0, 0, 0x80, 0x800, 0x10000}; // no overlong form
            return cp < shortest[n] || (cp >= 0xD800 && cp < 0xE000) || cp > 0x10FFFF ? 0 : n;
        }

        inline int hex(char c)
        {
            if (c >= '0' && c <= '9') return c - '0';
            if (c >= 'A' && c <= 'F') return c - 'A' + 10;
            if (c >= 'a' && c <= 'f') return c - 'a' + 10;
            return -1;
        }

        /**
         * Decodes the escape that starts `rest` (with '\'), giving its code points to put(cp, utf16);
         * returns the number of characters used, 0 if `rest` does not start with a valid escape.
         */
        template <class Put>
        std::size_t decode_escape(std::string_view rest, Put &&put)
        {
            if (rest.starts_with("\\\\"))
            {
                put('\\', false);
                return 2;
            }
            if (rest.starts_with("\\S\\") && rest.size() >= 4 && rest[3] >= ' ' && rest[3] <= '~')
            {
                put(std::uint32_t(rest[3]) + 0x80, false);
                return rest[3] == '\'' && rest.size() >= 5 ? 5 : 4; // an apostrophe is written ''
            }
            if (rest.starts_with("\\X\\") && rest.size() >= 5 && hex(rest[3]) >= 0 && hex(rest[4]) >= 0)
            {
                put(std::uint32_t(hex(rest[3]) * 16 + hex(rest[4])), false);
                return 5;
            }
            if (rest.starts_with("\\X2\\") || rest.starts_with("\\X4\\"))
            {
                const std::size_t w = rest[2] == '2' ? 4 : 8;
                const auto end = rest.find("\\X0\\", 4);
                if (end == std::string_view::npos || (end - 4) % w != 0 ||
                    !std::ranges::all_of(rest.substr(4, end - 4), [](char c) { return hex(c) >= 0; }))
                    return 0;
                for (std::size_t j = 4; j < end; j += w)
                {
                    std::uint32_t cp = 0;
                    for (std::size_t k = 0; k < w; ++k)
                        cp = cp * 16 + std::uint32_t(hex(rest[j + k]));
                    put(cp, w == 4);
                }
                return end + 4;
            }
            if (rest.size() >= 4 && rest[1] == 'P' && rest[2] >= 'A' && rest[2] <= 'I' && rest[3] == '\\')
                return 4; // \PA\ ... code page switch: ignored, \S\ maps to the upper half of ISO 8859-1
            return 0;
        }

        /**
         * Appends the body of a Part 21 string (between the quotes), decoded, to `out` in UTF-8:
         * '' -> ', \\ -> \, \S\c -> ISO 8859-1 upper half, \X\hh -> 8-bit code,
         * \X2\hhhh...\X0\ -> UTF-16 (surrogate pairs combined, also across adjacent blocks),
         * \X4\hhhhhhhh...\X0\ -> UTF-32, \Pc\ code page switches ignored, end-of-line characters ignored.
         * Never fails, and the result is valid UTF-8: a malformed or unknown escape is kept as text,
         * raw bytes that are not UTF-8 are read as ISO 8859-1 (the Part 21 default), and a lone
         * surrogate or a code point beyond U+10FFFF becomes U+FFFD.
         */
        inline void decode_string(std::string_view raw, std::string &out)
        {
            const auto plain = [](char c) { return c != '\\' && c != '\'' && c != '\n' && c != '\r' && static_cast<unsigned char>(c) < 0x80; };
            if (std::ranges::all_of(raw, plain))
            {
                out += raw;
                return;
            }
            if (raw.find_first_of("\r\n") != std::string_view::npos) // a long string may be wrapped anywhere
            {
                std::string joined;
                std::ranges::copy_if(raw, std::back_inserter(joined), [](char c) { return c != '\n' && c != '\r'; });
                decode_string(joined, out);
                return;
            }
            std::uint32_t high = 0; // UTF-16 high surrogate waiting for its low half, possibly in the next \X2\ block
            const auto flush = [&] {
                if (high)
                    append_utf8(out, 0xFFFD), high = 0;
            };
            const auto put = [&](std::uint32_t cp, bool utf16) {
                if (utf16 && high && cp >= 0xDC00 && cp < 0xE000)
                {
                    append_utf8(out, 0x10000 + ((high - 0xD800) << 10) + (cp - 0xDC00));
                    high = 0;
                    return;
                }
                flush();
                if (utf16 && cp >= 0xD800 && cp < 0xDC00)
                    high = cp;
                else
                    append_utf8(out, cp);
            };
            for (std::size_t i = 0; i < raw.size();)
            {
                const char c = raw[i];
                if (c == '\\')
                    if (const auto used = decode_escape(raw.substr(i), put))
                    {
                        i += used;
                        continue;
                    }
                if (c == '\'' && i + 1 < raw.size() && raw[i + 1] == '\'')
                    put('\'', false), i += 2;
                else if (static_cast<unsigned char>(c) < 0x80)
                    put(std::uint32_t(c), false), ++i; // a '\' that starts no valid escape is kept
                else if (const auto n = utf8_length(raw.substr(i)))
                {
                    flush();
                    out += raw.substr(i, n);
                    i += n;
                }
                else
                    put(static_cast<unsigned char>(c), false), ++i; // ISO 8859-1
            }
            flush();
        }

        class Lexer
        {
            std::string_view s_;
            std::size_t pos_{0}, line_{1}, line_start_{0};

        public:
            Tok tok{Tok::End};
            std::string_view text; ///< token text (keyword, number, string body, enum name, binary digits)
            std::size_t tok_line{1}, tok_col{1};

            explicit Lexer(std::string_view s) : s_{s}
            {
                if (s_.starts_with("\xEF\xBB\xBF")) // UTF-8 byte order mark
                    s_.remove_prefix(3);
            }

            /// Error at the current token.
            [[nodiscard]] P21Error error(std::string msg) const { return {tok_line, tok_col, std::move(msg)}; }
            /// Error at a position of the current line.
            [[nodiscard]] P21Error error_at(std::size_t pos, std::string msg) const { return {line_, pos - line_start_ + 1, std::move(msg)}; }

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
                    while (pos_ < s_.size() && is_space(s_[pos_]))
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
                    while (pos_ < s_.size() && is_digit(s_[pos_]))
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
                    // BINARY = '"' ("0" | "1" | "2" | "3") { HEX } '"' (HEX: 0-9, A-F); a long binary may be
                    // wrapped: its line breaks are skipped
                    std::size_t digits = 0;
                    for (++pos_;; ++pos_)
                    {
                        if (pos_ >= s_.size())
                            return error("unterminated binary");
                        const char d = s_[pos_];
                        if (d == '"')
                            break;
                        if (d == '\n')
                            ++line_, line_start_ = pos_ + 1;
                        else if (d != '\r' && (digits++ == 0 ? d < '0' || d > '3' : hex(d) < 0 || is_lower(d)))
                            return error_at(pos_, "malformed binary: the first digit is 0 to 3, then hexadecimal digits 0-9 A-F");
                    }
                    if (digits == 0)
                        return error("malformed binary: empty");
                    tok = Tok::Binary, text = s_.substr(start + 1, pos_ - start - 1);
                    ++pos_;
                    return std::nullopt;
                }
                if (c == '.' && pos_ + 1 < s_.size() && (is_alpha(s_[pos_ + 1]) || s_[pos_ + 1] == '_'))
                {
                    // ENUMERATION = "." UPPER { UPPER | DIGIT } "." with UPPER = A-Z and '_' (lowercase
                    // letters accepted and uppercased, as in keywords)
                    for (pos_ += 2; pos_ < s_.size() && (is_alpha(s_[pos_]) || is_digit(s_[pos_]) || s_[pos_] == '_');)
                        ++pos_;
                    if (pos_ >= s_.size() || s_[pos_] != '.')
                        return error("malformed enumeration: letters, digits and '_' between two dots");
                    tok = Tok::Enum, text = s_.substr(start + 1, pos_ - start - 1);
                    ++pos_;
                    return std::nullopt;
                }
                if (c == '+' || c == '-' || c == '.' || is_digit(c))
                {
                    bool real = false;
                    ++pos_;
                    while (pos_ < s_.size())
                    {
                        const char d = s_[pos_];
                        if (is_digit(d))
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
                if (is_alpha(c) || c == '_' || c == '!')
                {
                    ++pos_;
                    while (pos_ < s_.size() && (is_alpha(s_[pos_]) || is_digit(s_[pos_]) || s_[pos_] == '_' || s_[pos_] == '-'))
                        ++pos_;
                    tok = Tok::Keyword, text = s_.substr(start, pos_ - start);
                    return std::nullopt;
                }
                if (c >= ' ' && c <= '~')
                    return error(std::string("unexpected character '") + c + "'");
                constexpr char digits[] = "0123456789ABCDEF";
                const auto b = static_cast<unsigned char>(c);
                return error(std::string("unexpected character (byte 0x") + digits[b >> 4] + digits[b & 15] + ")");
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

        /// List or typed parameter being read (parse_value keeps them on this stack, not on the call stack).
        struct Open
        {
            std::uint32_t start{}; ///< list: its first item in scratch_; typed parameter: type id
            bool typed{false};
        };
        std::vector<Open> open_;

        using Err = std::unexpected<P21Error>;

        /// Number of the current '#123' token.
        std::expected<std::uint64_t, P21Error> instance_number() const
        {
            std::uint64_t id{};
            if (std::from_chars(lx_.text.data(), lx_.text.data() + lx_.text.size(), id).ec != std::errc{}) // only digits: out of range
                return Err(lx_.error("instance number too large: #" + std::string(lx_.text.substr(0, 40))));
            return id;
        }

        std::optional<P21Error> advance() { return lx_.next(); }

        std::optional<P21Error> expect(detail::Tok t, const char *what)
        {
            if (lx_.tok != t)
                return lx_.error(std::string("expected ") + what);
            return advance();
        }

        std::uint32_t intern(std::string_view name)
        {
            if (const auto it = f_.type_ids_.find(name); it != f_.type_ids_.end()) // usual case: known, uppercase
                return it->second;
            auto up = detail::to_upper(name);
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

        /// List node of the items of scratch_ from `mark`, which are moved to items_.
        std::uint32_t add_list(std::size_t mark)
        {
            detail::Node n;
            n.kind = ValueKind::List;
            n.a = std::uint32_t(f_.items_.size());
            n.n = std::uint32_t(scratch_.size() - mark);
            f_.items_.insert(f_.items_.end(), scratch_.begin() + std::ptrdiff_t(mark), scratch_.end());
            scratch_.resize(mark);
            return add_node(n);
        }

        /**
         * Parameter value -> node, the current token being its first: '(' params ')', KEYWORD '(' value ')'
         * or a simple value. Iterative: the lists and typed parameters being read are kept on open_, so
         * that the call stack does not grow with the nesting, which is limited to max_nesting levels.
         */
        std::expected<std::uint32_t, P21Error> parse_value()
        {
            using detail::Tok;
            open_.clear();
            for (;;)
            {
                // down to a simple value or an empty list
                std::uint32_t node{};
                if (lx_.tok == Tok::LParen || lx_.tok == Tok::Keyword)
                {
                    if (open_.size() == max_nesting)
                        return Err(lx_.error("nesting too deep: more than " + std::to_string(max_nesting) + " levels of lists or typed parameters"));
                    if (lx_.tok == Tok::Keyword) // typed parameter
                    {
                        open_.push_back({intern(lx_.text), true});
                        if (auto e = advance())
                            return Err(*e);
                        if (auto e = expect(Tok::LParen, "'(' after a type name"))
                            return Err(*e);
                        continue;
                    }
                    open_.push_back({std::uint32_t(scratch_.size()), false});
                    if (auto e = advance())
                        return Err(*e);
                    if (lx_.tok != Tok::RParen)
                        continue;
                    if (auto e = advance())
                        return Err(*e);
                    node = add_list(open_.back().start);
                    open_.pop_back();
                }
                else // simple value
                {
                    detail::Node n;
                    switch (lx_.tok)
                    {
                    case Tok::Ref:
                    {
                        n.kind = ValueKind::Reference;
                        const auto id = instance_number();
                        if (!id)
                            return Err(id.error());
                        n.set(*id);
                        break;
                    }
                    case Tok::Integer:
                    {
                        n.kind = ValueKind::Integer;
                        const auto *b = lx_.text.data() + (lx_.text.front() == '+' ? 1 : 0);
                        std::int64_t i{};
                        auto [p, ec] = std::from_chars(b, lx_.text.data() + lx_.text.size(), i);
                        if (ec != std::errc{} || p != lx_.text.data() + lx_.text.size())
                            return Err(lx_.error("invalid integer '" + std::string(lx_.text) + "'"));
                        n.set(i);
                        break;
                    }
                    case Tok::Real:
                    {
                        n.kind = ValueKind::Real;
                        double r{};
                        if (const auto ec = detail::parse_real(lx_.text, r); ec != std::errc{})
                            return Err(lx_.error((ec == std::errc::result_out_of_range ? "real out of range '" : "invalid real '") +
                                                 std::string(lx_.text.substr(0, 40)) + "'"));
                        n.set(r);
                        break;
                    }
                    case Tok::String:
                    {
                        n.kind = ValueKind::String;
                        const auto off = f_.chars_.size();
                        detail::decode_string(lx_.text, f_.chars_);
                        n.a = std::uint32_t(off);
                        n.n = std::uint32_t(f_.chars_.size() - off);
                        break;
                    }
                    case Tok::Enum:
                    {
                        n.kind = ValueKind::Enum;
                        n.a = add_chars(lx_.text);
                        n.n = std::uint32_t(lx_.text.size());
                        std::ranges::transform(f_.chars_.begin() + n.a, f_.chars_.end(), f_.chars_.begin() + n.a, [](char c) { return detail::to_upper(c); });
                        break;
                    }
                    case Tok::Binary:
                    {
                        n.kind = ValueKind::Binary;
                        const auto off = f_.chars_.size();
                        std::ranges::copy_if(lx_.text, std::back_inserter(f_.chars_), [](char c) { return c != '\n' && c != '\r'; });
                        n.a = std::uint32_t(off);
                        n.n = std::uint32_t(f_.chars_.size() - off);
                        break;
                    }
                    case Tok::Dollar:
                        n.kind = ValueKind::Omitted;
                        break;
                    case Tok::Star:
                        n.kind = ValueKind::Derived;
                        break;
                    default:
                        return Err(lx_.error("expected a parameter value"));
                    }
                    node = add_node(n);
                    if (auto e = advance())
                        return Err(*e);
                }
                // up: the value completes typed parameters and closed lists, until a list goes on after ','
                for (;;)
                {
                    if (open_.empty())
                        return node;
                    const auto top = open_.back();
                    if (top.typed)
                    {
                        if (auto e = expect(Tok::RParen, "')' closing a typed parameter"))
                            return Err(*e);
                        detail::Node n;
                        n.kind = ValueKind::Typed;
                        n.a = top.start;
                        n.n = node;
                        node = add_node(n);
                        open_.pop_back();
                        continue;
                    }
                    scratch_.push_back(node);
                    if (lx_.tok == Tok::Comma)
                    {
                        if (auto e = advance())
                            return Err(*e);
                        break;
                    }
                    if (auto e = expect(Tok::RParen, "',' or ')'"))
                        return Err(*e);
                    node = add_list(top.start);
                    open_.pop_back();
                }
            }
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
            auto args = parse_value(); // the list of arguments
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
                auto l = parse_value();
                if (!l)
                    return l.error();
            }
            if (auto e = expect(detail::Tok::Semicolon, "';' after DATA"))
                return e;
            while (!(lx_.tok == detail::Tok::Keyword && lx_.text == "ENDSEC"))
            {
                if (lx_.tok != detail::Tok::Ref)
                    return lx_.error("expected an instance '#n=' or ENDSEC");
                const auto id = instance_number();
                if (!id)
                    return id.error();
                if (auto e = advance())
                    return e;
                if (auto e = expect(detail::Tok::Equal, "'=' after the instance number"))
                    return e;
                if (auto e = parse_entity(*id, f_.instances_))
                    return e;
                ids_.emplace_back(*id, std::uint32_t(f_.instances_.size() - 1));
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
            // type index, every partial type of a complex instance (counted first: no slack)
            std::vector<std::uint32_t> count(f_.types_.size());
            for (const auto &in : f_.instances_)
                for (std::uint32_t k{}; k < in.n_parts; ++k)
                    ++count[f_.parts_[in.first_part + k].type];
            for (std::size_t t{}; t < count.size(); ++t)
                f_.by_type_[t].reserve(count[t]);
            for (std::uint32_t i{}; i < f_.instances_.size(); ++i)
            {
                const auto &in = f_.instances_[i];
                for (std::uint32_t k{}; k < in.n_parts; ++k)
                    f_.by_type_[f_.parts_[in.first_part + k].type].push_back(i);
            }
            // capacity reserved from the text size and left untouched costs no physical memory; give back a
            // large excess only (a text sparser than estimated), as shrinking copies the table
            const auto fit = [](auto &v) {
                if (v.capacity() > v.size() + v.size() / 2)
                    v.shrink_to_fit();
            };
            fit(f_.nodes_), fit(f_.items_), fit(f_.parts_), fit(f_.instances_);
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
        explicit P21Parser(std::string_view text) : lx_{text}
        {
            // tables sized from the text (measured on real files: 7 to 9 bytes of text per value, 9 to 12 per
            // list item, 35 to 45 per partial entity, 37 and more per instance): no growth for a usual file
            const auto n = text.size();
            f_.nodes_.reserve(n / 7 + 16);
            f_.items_.reserve(n / 8 + 16);
            f_.parts_.reserve(n / 32 + 4);
            f_.instances_.reserve(n / 36 + 4);
            ids_.reserve(n / 36 + 4);
        }

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

    namespace detail
    {
        /// Whole content of a file, read in one block sized from the end position (a stream read when the
        /// file cannot seek); nullopt if it cannot be opened or read.
        inline auto read_text(const std::filesystem::path &file) -> std::optional<std::string>
        {
            if (std::error_code ec; std::filesystem::is_directory(file, ec))
                return std::nullopt;
            std::ifstream in(file, std::ios::binary);
            if (!in)
                return std::nullopt;
            std::string text;
            in.seekg(0, std::ios::end);
            if (const auto end = in.tellg(); end != std::streampos(-1) && in.seekg(0, std::ios::beg))
            {
                text.resize(static_cast<std::size_t>(end));
                in.read(text.data(), static_cast<std::streamsize>(text.size()));
                text.resize(static_cast<std::size_t>(in.gcount())); // in case the file shrank meanwhile
            }
            else // pipe or special file
            {
                in.clear();
                text.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
            }
            if (in.bad())
                return std::nullopt;
            return text;
        }
    } // namespace detail

    /// Reads and parses a Part 21 file.
    [[nodiscard]] inline auto read_p21(const std::filesystem::path &file) -> std::expected<P21File, P21Error>
    {
        const auto text = detail::read_text(file);
        if (!text)
            return std::unexpected(P21Error{0, 0, "cannot open " + file.string()});
        return parse_p21(*text);
    }

} // namespace gbs::step
