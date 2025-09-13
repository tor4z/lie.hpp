#include <iostream>
#define LIE_IMPLEMENTATION
#include "lie.hpp"

int main()
{
    lie::Matrix<float, 2, 2> m;
    m << 1, 2,
         4, 5;

    auto m2 = m * m;
    std::cout << m << "\n";

    lie::Matrix<float, 3, 3> m33;
    m33 << 1, 2, 3,
           2, 4, 6,
           7, 8, 108;
    std::cout << "m3x3: " << m33 << "\n";
    std::cout << "det of m3x3: " << m33.det() << "\n";
    std::cout << "det of eye3: " << lie::Matrix3f::eye().det() << "\n";

    m33.null_space();
    return 0;
}
