#include <iostream>
#define LIE_IMPLEMENTATION
#include "lie.hpp"

int main()
{
    lie::SE3f se3;
    lie::Vector3f v;
    se3 << 1, 2, 3, 4,
           5, 6, 7, 8,
           9, 10, 11, 12,
           13, 14, 15, 16;
    v << 1, 2, 3;

    std::cout << "se3: " << se3 << "\n";
    std::cout << "se3 rot: " << se3.rot() << "\n";
    std::cout << "se3 offset: " << se3.offset() << "\n";
    std::cout << "act se3 on v: " << se3 * v << "\n";

    {
        lie::Vectorf<6> v;
        v << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6;
        std::cout << "Exp: " << v.SE_Exp() << "\n";
        std::cout << "Exp.Log: " << v.SE_Exp().Log() << "\n";
    }
    return 0;
}
