#define LIE_IMPLEMENTATION
#define LIE_PRINT_SETPRECISION 5
#include "lie.hpp"
#include <iostream>

int main()
{
    std::cout << "Hello, localization\n";
    auto curr{lie::SE3f::eye()};

    int cnt{0};
    while (cnt++ < 10) {
        lie::Vector6f v{};
        v << 0.1f, 0, 0, 0.0, 0, 0;
        curr = curr.rplus(v);
        std::cout << curr.Log() << "\n";
    }
    return 0;
}
