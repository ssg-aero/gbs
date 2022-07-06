#pragma once
#include <algorithm>
#include <vector>
#include <list>
#include <array>
#include <exception>
#include <TopTools_ListOfShape.hxx>
#include <GeomAPI.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS.hxx>
#include <gp_Pln.hxx>
#include <math_Vector.hxx>
#include <NCollection_Sequence.hxx>
#include <NCollection_Array2.hxx>
#include <gbs/occt/geomprim.h>
#include <gbs/occt/brepbuilders.h>


namespace occt_utils
{
    template <typename container>
    inline auto to_sh_list(const container &curves, TopTools_ListOfShape &L) -> void
    {
        to_sh_list(curves.begin(), curves.end(), L);
    }

    template <typename container>
    inline auto to_sh_list(const container &curves) -> TopTools_ListOfShape
    {
        TopTools_ListOfShape L;
        to_sh_list(curves,L);
        return L;
    }

    inline auto to_math_vec(const std::vector<double> &v )
    {
        return math_Vector(v.data(),1,v.size());
    }

    inline auto from_math_vec( const math_Vector &v)
    {
        std::vector<double> v_( v.Length() );
        //Pas d'accès au ptr de base de math_Vector
        for(auto i = v.Lower() ; i <= v.Upper() ; i++)
        {
            v_[i - v.Lower() ] = v(i);
        }

        return v_;
    }

    template <class T>
    NCollection_Array1<T> col_from_vec(const std::vector<T> &v, int iStart = 1)
    {
        auto length = Standard_Integer(v.size());
        return NCollection_Array1<T>(v[0], iStart, iStart + length - 1);
    }

    template <class T>
    NCollection_Array1<T> col_from_vec(const std::vector<std::array<double, 2>> &v, int iStart = 1)
    {
        auto length = Standard_Integer(v.size());
        NCollection_Array1<T> arr(iStart, iStart + length - 1);
        std::transform(v.begin(), v.end(), arr.begin(),
                       [&](const std::array<double, 2> &X_) { return T(gp_XY(X_[0], X_[1])); });
        return arr;
    }

    template <class T>
    NCollection_Array1<T> col_from_vec(const std::vector<std::array<double, 3>> &v, int iStart = 1)
    {
        auto length = Standard_Integer(v.size());
        NCollection_Array1<T> arr(iStart, iStart + length - 1);
        std::transform(v.begin(), v.end(), arr.begin(),
                       [&](const std::array<double, 3> &X_) { return T(gp_XYZ(X_[0], X_[1], X_[2])); });
        return arr;
    }

    template <class T>
    NCollection_Array2<T> col2_from_vec(const std::vector<std::array<double, 3>> &v, size_t nj, int iStart = 1)
    {
        auto length = Standard_Integer(v.size());
        if (length % nj)
        {
            throw std::length_error(std::string("Bad shaped array"));
        }
        auto ni = length / nj;

        // auto col =  col_from_vec<T>(v);

        // NCollection_Array2<T> arr(col.ChangeValue(iStart), iStart,ni+iStart-1,iStart,nj+iStart-1); // Not the most efficent but more compact

        NCollection_Array2<T> arr(iStart, ni + iStart - 1, iStart, nj + iStart - 1);
        for (int i = 0; i < ni; i++)
        {
            for (int j = 0; j < nj; j++)
            {
                auto X_ = v[i + ni * j];
                arr.SetValue(i + iStart, j + iStart, T(gp_XYZ(X_[0], X_[1], X_[2])));
            }
        }

        return arr;
    }

    template <typename T,typename _InIt>
    inline NCollection_Sequence<T> col_from_list(const _InIt &_First, const _InIt &_Last)
    {
        NCollection_Sequence<T> S;
        std::for_each(_First,_Last,[&](const T & t_){S.Append(t_);});
        return S;
    }
    template <class T>
    inline NCollection_Sequence<T> col_from_list(const std::vector<T> &lst)
    {
        return col_from_list<T>(lst.begin(),lst.end() );
    }

    template <class T>
    inline NCollection_Sequence<T> col_from_list(const std::list<T> &lst)
    {
        return col_from_list<T>(lst.begin(),lst.end());
    }

} // namespace occt_utils