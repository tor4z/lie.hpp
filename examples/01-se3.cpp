#include <iostream>
#define LIE_IMPLEMENTATION
#include "lie.hpp"

int main()
{
    lie::SE3f se3;
    lie::Vector3f v;
    se3 << 0.902, -0.376, 0.213, 4,
           0.396, 0.916, -0.057, 8,
           -0.173, 0.136, 0.975, 12,
           0, 0, 0, 1;
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
        std::cout << "Exp.Log.hat: " << v.SE_Exp().Log().hat() << "\n";
    }

    std::cout << "=========== SE inv ==========\n";
    {   
        lie::Vectorf<6> v;
        v << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6;
        auto valid_se3{v.SE_Exp()};
        std::cout << "valid_se3: " << valid_se3 << "\n";
        std::cout << "se3 * se3.inv(): " << valid_se3 * valid_se3.inv() << "\n";
    }

    std::cout << "=========== Adj ==============\n";
    {
        lie::Vectorf<6> v;
        v << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6;
        lie::Vectorf<6> v2;
        v2 << 0.1, 0.9, 0.1, 0.5, 0.5, 0.9;
        auto se32{v2.SE_Exp()};
    
        std::cout << "v: " << v << "\n";
        std::cout << "adj: " << se32.adj(v) << "\n";
        std::cout << "adj2: " << se32.Adj() * v << "\n";
    }
    return 0;
}
