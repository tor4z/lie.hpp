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

    lie::SO3f so3;
    std::cout << so3 << "\n";
    return 0;
}
