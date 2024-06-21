#ifndef LC_1601_IO_CPP
#define LC_1601_IO_CPP

#include <ostream>

constexpr const char* sep(unum_t i, const char* s2 = ", ", const char *s1 = "") {
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

template<typename CharT, class Traits, size_t N, typename T>
std::basic_ostream<CharT,Traits>& operator<< (std::basic_ostream<CharT, Traits>& os, Vector<N, T> const& vector) {
    os << "[";
    for(size_t i = 0; i < vector.size(); i++) {
        os << sep(i) << vector(i);
    }
    os << "]";
    return os;
}

template<typename CharT, class Traits, size_t N, typename T>
std::basic_ostream<CharT,Traits>& operator<< (std::basic_ostream<CharT, Traits>& os, Matrix<N, T> const& matrix) {
    for(unum_t i = 0; i < matrix.size(); i++) {
        os << "[";
        for (unum_t j = 0; j < matrix.size(); j++) {
            os << sep(j) << matrix(i, j);
        }
        os << "]" << std::endl;
    }
    return os;
}

template <typename CharT, class Traits, size_t N>
std::basic_ostream<CharT, Traits>& operator<< (std::basic_ostream<CharT, Traits>& os, NodeSeq<N> const& nodes) {
    for(unum_t k = 0; k < nodes.size(); ++k) {
        std::cout << sep(k) << nodes(k);
    }
    return os;
}

template<typename CharT, class Traits>
std::basic_ostream<CharT, Traits>& operator<< (std::basic_ostream<CharT, Traits>& os, AdjLst const& list) {
    unum_t k = 0;
    for (auto const& ref: list) {
        os << sep(k++, "\n") << ref.i() << " -> {" <<  ref.list() << "}";
    }
    os << std::endl;
    return os;
}

template<typename CharT, class Traits>
std::basic_ostream<CharT, Traits>& operator<< (std::basic_ostream<CharT, Traits>& os, Pair const& pair) {
    os << "(" << pair.i() << ", " << pair.j() << ")";
    return os;
}

template<typename CharT, class Traits>
std::basic_ostream<CharT, Traits>& operator<< (std::basic_ostream<CharT, Traits>& os, ShortestPaths const& paths) {
    os << std::endl
       << "-- distances:"
       << std::endl
       << paths.distances()
       << "-- predecessors:"
       << std::endl
       << paths.predecessors();
    return os;
}

template<typename CharT, class Traits>
std::basic_ostream<CharT, Traits>& operator<< (std::basic_ostream<CharT, Traits>& os, Graph const& graph) {
    os << "--- flow(" << graph.flow().size() << "): " << std::endl << graph.flow();
    os << "--- inflow: " << std::endl; for (auto i: graph.nodes()) { os << i << ": " << graph.inflow(i) << std::endl; }
    os << "--- outflow: " << std::endl; for (auto i: graph.nodes()) { os << i << ": " << graph.outflow(i) << std::endl; }
    os << "--- outbound("<< graph.outbound().size() << "): " << std::endl << graph.outbound();
    os << "--- inbound("<< graph.inbound().size() << "): " << std::endl << graph.inbound();
    return os;
}

#endif /* LC_1601_IO_CPP */
