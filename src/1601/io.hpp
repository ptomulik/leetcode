#ifndef LC_1601_IO_HPP
#define LC_1601_IO_HPP

#include <ostream>

constexpr const char* sep(size_t i, const char* s2 = ", ", const char *s1 = "") {
    return 0 == i ? s1 : s2;
}

template<typename CharT, class Traits, typename T>
std::basic_ostream<CharT,Traits>& operator<< (std::basic_ostream<CharT, Traits>& os, std::vector<T> const& vector) {
    os << "[";
    for(size_t i = 0; i < vector.size(); i++) {
        os << sep(i) << vector[i];
    }
    os << "]";
    return os;
}

namespace Templates {

//template<typename CharT, class Traits, typename T, size_t N>
//std::basic_ostream<CharT,Traits>& operator<< (std::basic_ostream<CharT, Traits>& os, jector<T, N> const& vector) {
//    os << "[";
//    for(size_t i = 0; i < vector.size(); i++) {
//        os << sep(i) << vector(i);
//    }
//    os << "]";
//    return os;
//}

template<typename CharT, class Traits, size_t N, typename T>
std::basic_ostream<CharT,Traits>& operator<< (std::basic_ostream<CharT, Traits>& os, Matrix<T, N> const& matrix) {
    for(size_t i = 0; i < matrix.size(); i++) {
        os << "[";
        for (size_t j = 0; j < matrix.size(); j++) {
            os << sep(j) << matrix(i, j);
        }
        os << "]" << std::endl;
    }
    return os;
}

template<typename CharT, class Traits, typename T, size_t N>
std::basic_ostream<CharT,Traits>& operator<< (std::basic_ostream<CharT, Traits>& os, SortSeq<T, N> const& seq) {
    size_t k = 0;
    for (auto v: seq) {
        os << sep(k++, ", ") << (int)v;
    }
    return os;
}

template<typename CharT, class Traits, typename T, size_t M>
std::basic_ostream<CharT,Traits>& operator<< (std::basic_ostream<CharT, Traits>& os, Edge<T, M> const& edge) {
    os << "(" << (int) edge.i << ", " << (int) edge.j << ")[" << (int) edge.m << "]";
    return os;
}

template<typename CharT, class Traits, typename T, size_t M>
std::basic_ostream<CharT,Traits>& operator<< (std::basic_ostream<CharT, Traits>& os, Endpoint<T, M> const& endpoint) {
    os << (int) endpoint.node << "[" << (int) endpoint.edge << "]";
    return os;
}

template<typename CharT, class Traits, typename ObjectT, size_t N>
std::basic_ostream<CharT,Traits>& operator<< (std::basic_ostream<CharT, Traits>& os, EncodedObjectSortSeq<ObjectT, N> const& seq) {
    size_t k = 0;
    for (auto obj: seq) {
        os << sep(k++, ", ") << (ObjectT)obj;
    }
    return os;
}

template<typename CharT, class Traits, typename EdgeT, size_t N>
std::basic_ostream<CharT, Traits>& operator<< (std::basic_ostream<CharT, Traits>& os, AdjLst<EdgeT, N> const& list) {
    size_t k = 0;
    for (auto const& ref: list) {
        os << sep(k++, "\n") << (int)ref.i() << " -> {" <<  ref.list() << "}";
    }
    return os;
}

//template<typename CharT, class Traits>
//std::basic_ostream<CharT, Traits>& operator<< (std::basic_ostream<CharT, Traits>& os, ShortestPaths const& paths) {
//    os << std::endl
//       << "-- distances:"
//       << std::endl
//       << paths.distances()
//       << "-- predecessors:"
//       << std::endl
//       << paths.predecessors();
//    return os;
//}
//
template<typename CharT, class Traits, typename T, size_t N, size_t M>
std::basic_ostream<CharT, Traits>& operator<< (std::basic_ostream<CharT, Traits>& os, Graph<T, N, M> const& graph) {
    typedef typename Graph<T, N, M>::edge_type edge_type;

    os << "--- nodes(" << graph.nodes().size() << "): " << std::endl << "[" << graph.nodes() << "]" << std::endl;
    os << "--- outbound(" << graph.outbound().size() << "): " << std::endl << graph.outbound() << std::endl;
    os << "--- inbound(" << graph.inbound().size() << "): " << std::endl << graph.inbound() << std::endl;
    os << "--- edges: " << std::endl;

    size_t k = 0;
    for (auto ref: graph.outbound()) {
        for (auto head: ref.list()) {
            auto e = edge_type(ref.i(), head);
            os << sep(k++, "\n") << e << " = { c: " << (int)graph.cost(e) << ", u:" << (int)graph.caps(e) << ", x: " << (int)graph.flow(e) << " }";
        }
    }
    return os;
}

} /* namespace Templates */

#endif /* LC_1601_IO_HPP */
