 #include <topology\edge.h>
namespace gbs
{
    // Explicit instantiation of the Edge class for required dimensions and types.
    template class gbs::Edge<float, 2>;
    template class gbs::Edge<float, 3>;
    template class gbs::Edge<double, 2>;
    template class gbs::Edge<double, 3>;
}