 #include <topology\edge.h>
namespace gbs
{
    // Explicit instantiation of the Edge class for required dimensions and types.
    template class Edge<float, 2>;
    template class Edge<float, 3>;
    template class Edge<double, 2>;
    template class Edge<double, 3>;
}