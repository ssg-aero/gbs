#pragma once

/**
 * @file units.h
 * @brief Units and uncertainty of a STEP representation context.
 *
 * Design: docs/sources/design/step_reader.md, section 6; architecture note
 * docs/sources/design/step_pr04_geometry.md.
 *
 * A geometric_representation_context is a complex instance carrying
 * GLOBAL_UNIT_ASSIGNED_CONTEXT (length, plane angle and solid angle units) and
 * GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT (the distance accuracy). Units are SI
 * units with an optional prefix, or conversion based units (inch, degree…)
 * defined by a measure in another unit.
 */

#include <gbs-io/step/p21.h>

#include <cmath>
#include <optional>
#include <string>
#include <string_view>
#include <utility>

namespace gbs::step
{
    /// Error in the content of a STEP file met while interpreting it (unknown unit, invalid geometry…).
    class StepError : public std::runtime_error
    {
        std::uint64_t id_;

    public:
        StepError(std::uint64_t id, const std::string &what) : std::runtime_error{"#" + std::to_string(id) + ": " + what}, id_{id} {}
        /// #id of the offending instance.
        [[nodiscard]] std::uint64_t id() const noexcept { return id_; }
    };

    /// An entity that is valid STEP but not supported by the reader (hyperbola, offset curve…).
    class StepUnsupported : public StepError
    {
        std::string type_;

    public:
        StepUnsupported(std::uint64_t id, std::string_view type)
            : StepError{id, "unsupported entity " + std::string(type)}, type_{type} {}
        [[nodiscard]] const std::string &type() const noexcept { return type_; }
    };

    /// Units of a representation context, as factors from the file units.
    struct Units
    {
        double length_mm{1.};             ///< millimetres per file length unit
        double length{1.};                ///< factor applied to lengths: file unit -> target unit
        double angle{1.};                 ///< factor applied to plane angles: file unit -> radian
        std::string length_name{};        ///< e.g. "MILLI METRE", "METRE", "INCH"; empty if the context declares none
        std::string angle_name{};         ///< e.g. "RADIAN", "DEGREE"
        std::optional<double> uncertainty; ///< distance accuracy of the file, in target units
    };

    namespace detail
    {
        inline double si_prefix(std::string_view p)
        {
            constexpr std::pair<std::string_view, double> table[]{
                {"EXA", 1e18}, {"PETA", 1e15}, {"TERA", 1e12}, {"GIGA", 1e9}, {"MEGA", 1e6}, {"KILO", 1e3},
                {"HECTO", 1e2}, {"DECA", 1e1}, {"DECI", 1e-1}, {"CENTI", 1e-2}, {"MILLI", 1e-3}, {"MICRO", 1e-6},
                {"NANO", 1e-9}, {"PICO", 1e-12}, {"FEMTO", 1e-15}, {"ATTO", 1e-18}};
            for (auto [name, f] : table)
                if (name == p)
                    return f;
            throw std::runtime_error("unknown SI prefix " + std::string(p));
        }

        struct UnitValue
        {
            double factor; ///< length: mm per unit ; angle: radian per unit
            std::string name;
        };

        // Value and unit of a measure_with_unit, simple (LENGTH_MEASURE_WITH_UNIT(value, unit))
        // or complex ((LENGTH_MEASURE_WITH_UNIT() MEASURE_WITH_UNIT(value, unit))).
        inline std::pair<double, std::uint64_t> measure_with_unit(const InstanceView &m)
        {
            auto a = m.part("MEASURE_WITH_UNIT") ? m.args("MEASURE_WITH_UNIT") : m.args();
            if (a.size() < 2)
                throw StepError(m.id(), "measure with unit without value or unit");
            return {a[0].as_measure(), a[1].as_ref()};
        }

        // Own (prefix, name) arguments of an SI_UNIT, or (name, conversion_factor) of a
        // CONVERSION_BASED_UNIT: a simple instance has the inherited dimensions first.
        inline std::pair<ListView, std::size_t> own_unit_args(const InstanceView &u, std::string_view type)
        {
            auto a = u.part(type)->args();
            if (a.size() < 2)
                throw StepError(u.id(), "missing arguments of " + std::string(type));
            return {a, u.is_complex() ? 0 : a.size() - 2};
        }

        inline UnitValue unit_value(const P21File &f, std::uint64_t id, int depth = 0)
        {
            if (depth > 8)
                throw StepError(id, "circular unit definition");
            const auto u = f.instance(id);
            if (u.part("SI_UNIT"))
            {
                const auto [a, k] = own_unit_args(u, "SI_UNIT");
                const double prefix = a[k].is_null() ? 1. : si_prefix(a[k].as_enum());
                const auto name = a[k + 1].as_enum();
                const std::string full = a[k].is_null() ? std::string(name) : std::string(a[k].as_enum()) + " " + std::string(name);
                if (name == "METRE")
                    return {1000. * prefix, full};
                if (name == "RADIAN")
                    return {prefix, full};
                if (name == "STERADIAN")
                    return {prefix, full};
                throw StepError(id, "unsupported SI unit " + std::string(name));
            }
            if (u.part("CONVERSION_BASED_UNIT"))
            {
                const auto [a, k] = own_unit_args(u, "CONVERSION_BASED_UNIT");
                const auto [value, base] = measure_with_unit(f.instance(a[k + 1].as_ref()));
                return {value * unit_value(f, base, depth + 1).factor, std::string(a[k].as_string())};
            }
            throw StepUnsupported(id, u.type());
        }

        enum class UnitKind : std::uint8_t
        {
            Length,
            PlaneAngle,
            Other, ///< solid angle, time, mass…
        };

        /**
         * What a unit measures: from LENGTH_UNIT or PLANE_ANGLE_UNIT in a complex instance. A simple
         * instance lists no such type: an SI_UNIT is classified by its name (METRE, RADIAN), a
         * CONVERSION_BASED_UNIT by the unit of its conversion factor.
         */
        inline UnitKind unit_kind(const P21File &f, std::uint64_t id, int depth = 0)
        {
            if (depth > 8)
                throw StepError(id, "circular unit definition");
            const auto u = f.instance(id);
            if (u.has_type("LENGTH_UNIT"))
                return UnitKind::Length;
            if (u.has_type("PLANE_ANGLE_UNIT"))
                return UnitKind::PlaneAngle;
            if (u.part("SI_UNIT"))
            {
                const auto [a, k] = own_unit_args(u, "SI_UNIT");
                const auto name = a[k + 1].as_enum();
                return name == "METRE" ? UnitKind::Length : name == "RADIAN" ? UnitKind::PlaneAngle : UnitKind::Other;
            }
            if (u.part("CONVERSION_BASED_UNIT"))
            {
                const auto [a, k] = own_unit_args(u, "CONVERSION_BASED_UNIT");
                return unit_kind(f, measure_with_unit(f.instance(a[k + 1].as_ref())).second, depth + 1);
            }
            return UnitKind::Other;
        }
    }

    /**
     * @brief Units and distance uncertainty of a geometric representation context.
     *
     * @param context        the GEOMETRIC_REPRESENTATION_CONTEXT complex instance
     * @param target_unit_mm length of the target unit in millimetres (1: millimetre, 25.4: inch)
     * A context without units gives millimetres and radians, with empty unit names.
     */
    inline Units read_units(const P21File &f, std::uint64_t context, double target_unit_mm = 1.)
    {
        if (!(target_unit_mm > 0.) || !std::isfinite(target_unit_mm))
            throw std::invalid_argument("read_units: target unit must be positive");
        Units u;
        const auto ctx = f.instance(context);
        if (auto ua = ctx.part("GLOBAL_UNIT_ASSIGNED_CONTEXT"))
            for (auto ref : ua->args()[0].as_list())
            {
                const auto kind = detail::unit_kind(f, ref.as_ref());
                if (kind == detail::UnitKind::Length)
                {
                    auto v = detail::unit_value(f, ref.as_ref());
                    u.length_mm = v.factor;
                    u.length_name = v.name;
                }
                else if (kind == detail::UnitKind::PlaneAngle)
                {
                    auto v = detail::unit_value(f, ref.as_ref());
                    u.angle = v.factor;
                    u.angle_name = v.name;
                }
            }
        u.length = u.length_mm / target_unit_mm;
        if (auto ug = ctx.part("GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT"))
            for (auto ref : ug->args()[0].as_list())
            {
                const auto m = f.instance(ref.as_ref());
                const auto [value, unit] = detail::measure_with_unit(m);
                if (detail::unit_kind(f, unit) != detail::UnitKind::Length)
                    continue;
                const double v = value * detail::unit_value(f, unit).factor / target_unit_mm;
                // the distance accuracy if named so, otherwise the first length uncertainty
                auto a = m.part("UNCERTAINTY_MEASURE_WITH_UNIT") ? m.args("UNCERTAINTY_MEASURE_WITH_UNIT") : m.args();
                const bool named = m.is_complex() ? (a.size() > 0 && a[0].as_string() == "distance_accuracy_value")
                                                  : (a.size() > 2 && a[2].as_string() == "distance_accuracy_value");
                if (named || !u.uncertainty)
                    u.uncertainty = v;
            }
        return u;
    }
}
