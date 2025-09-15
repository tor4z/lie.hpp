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

    auto m33_null{m33.null_space()};
    std::cout << "null space of m33: " << m33_null << "\n";

    lie::Matrix<float, 3, 3> m33_inv;
    m33_inv << 1, 2, 3,
               2, 50, 6,
               7, 8, 108;
    std::cout << "another m33: " << m33_inv << "\n";
    std::cout << "m33.inv: " << m33_inv.inv() << "\n";
    std::cout << "m33 * m33.inv: " << m33_inv * m33_inv.inv() << "\n";
    return 0;
}
