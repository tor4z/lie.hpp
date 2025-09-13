#include <iostream>
#define LIE_IMPLEMENTATION
#include "lie.hpp"

int main()
{
    lie::Vector<float, 3> v;
    v << 1, 2, 4;
    std::cout << v << "\n";
    std::cout << "norm2 of v: " << v.norm2() << "\n";
    std::cout << "hat: " << v.hat() << "\n";
    std::cout << "is_skew_sym: " << v.hat().is_skew_sym() << "\n";
    return 0;
}
