#include <gbs/bscapprox.h>
#include <charconv>
#include <cctype>
#include <fstream>
namespace gbs
{
    void split(const std::string &s, char c,
               std::vector<std::string> &v)
    {
        std::string::size_type i = 0;
        std::string::size_type j = s.find(c);

        while (j != std::string::npos)
        {
            v.push_back(s.substr(i, j - i));
            i = ++j;
            j = s.find(c, j);

            if (j == std::string::npos)
                v.push_back(s.substr(i, s.length()));
        }

        // remove withe spaces
        for( auto &v_ : v)
        {
            v_.erase(
                std::remove_if(
                    v_.begin(), v_.end(),
                    [](const auto &c) { return std::isspace(c); }
                ),
                v_.end()
            );
        }
        // remove empty strings
        v.erase(
            std::remove_if(
                v.begin(), v.end(),
                [](const auto &s) { return s.size() == 0; }
            ),
            v.end()
        );

    }

    template<typename T, size_t dim>
    auto bscurve_approx_from_points(const std::string &fname, size_t p, gbs::KnotsCalcMode mode = gbs::KnotsCalcMode::CHORD_LENGTH, size_t n_skip_line = 0, char split_char = ' ') -> BSCurve<T,dim>
    {
        std::string line;
        std::ifstream myfile(fname);
        if(!myfile.is_open())
        {
            throw std::runtime_error("Could not open file:" + fname);
        }
        points_vector<T, dim> pts;
        while(n_skip_line)
        {
            getline(myfile, line);
            n_skip_line--;
        }
        while (getline(myfile, line))
        {
            std::vector<std::string> v;
            split(line, split_char, v);
            point<T, dim> pt;
            if( v.size() != dim)
            {
                throw std::runtime_error("Incorrect file format");
            }
            std::transform(v.begin(),v.end(),pt.begin(),
                [](const auto &s_)
                {
                    T result{};
                    auto [ptr, ec] { std::from_chars(s_.data(), s_.data() + s_.size(), result) };
                    if (ec != std::errc{})
                    {
                        throw std::runtime_error("Incorrect file format");
                    }
                    return result;
                }
            );
            pts.push_back(pt);
        }
        myfile.close();

        return gbs::approx(pts, p, mode,true);
    }

    template<typename T, size_t dim>
    auto bscurve_approx_from_points(const char* fname, size_t p, gbs::KnotsCalcMode mode, size_t n_skip_line) -> BSCurve<T,dim>
    {
        return bscurve_approx_from_points<T,dim>(std::string(fname), p, mode, n_skip_line);
    }
}