#pragma once
#include <gbs/bscurve.h>

namespace gbs{

    template <typename T, size_t dim>
    class CurveReparametrized : public Curve<T,dim>
    {
        const std::shared_ptr<Curve<T,dim>>   m_p_crv;
        const BSCfunction<T> m_f_u;
    public:
        CurveReparametrized(const std::shared_ptr<Curve<T,dim>> &crv, const BSCfunction<T> &f_u) : m_p_crv{crv}, m_f_u{f_u} {}
        /**
         * @brief Curve evaluation at parameter u
         *
         * @param u : parameter on curve
         * @param d : derivative order
         * @return std::array<T, dim>
         */
        virtual auto value(T u, size_t d = 0) const -> std::array<T, dim>  override
        {
            switch (d)
            {
            case 0:
                return m_p_crv->value(m_f_u(u));
                break;
            case 1:
                return m_p_crv->value(m_f_u(u),1)*m_f_u(u,1);
                break;
            case 2:
            {
                auto fd0 = m_f_u(u);
                auto fd1 = m_f_u(u, 1);
                return m_p_crv->value(fd0, 1) * m_f_u(u, 2) + m_p_crv->value(fd0, 2) * fd1 * fd1;
            }
                break;
            default:
                throw std::runtime_error("Not implemented yet.");
                break;
            }

        }
        auto bounds() const -> std::array<T, 2> override
        {
            return m_f_u.bounds();
        }
        auto basisCurve() const -> const std::shared_ptr<Curve<T,dim>>&
        {
            return m_p_crv;
        }
        auto paramFunction() const -> const BSCfunction<T>&
        {
            return m_f_u;
        }
        auto setParamFunction(const BSCfunction<T> &f_u)
        {
            m_f_u = f_u;
        }
    };
    }