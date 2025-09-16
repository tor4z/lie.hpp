#include <iostream>
#define LIE_IMPLEMENTATION
#include "lie.hpp"

int main()
{
    lie::Vector<float, 3> v;
    v << 0.1, 0.2, 0.4;
    lie::Vector<float, 3> v2;
    v2 << 0.7, 0.5, 0.8;
    auto so32{v2.SO_Exp()};
    std::cout << "============ exp ==========\n";
    {
        auto m{v.hat()};
        auto so3{m.SO_exp()};
        std::cout << "so3: " << so3 << "\n";
        std::cout << "so3^T: " << so3.t() << "\n";
        std::cout << "so3 * so3.t(): " << so3 * so3.t() << "\n";
    }

    std::cout << "============ Exp ==========\n";
    {
        auto so3{v.SO_Exp()};
        std::cout << "so3: " << so3 << "\n";
        std::cout << "so3^T: " << so3.t() << "\n";
        std::cout << "so3 * so3.t(): " << so3 * so3.t() << "\n";
    }

    std::cout << "=========== Log ==============\n";
    {
        auto so3{v.SO_Exp()};
        std::cout << "v: " << v << "\n";
        std::cout << "so3: " << so3 << "\n";
        std::cout << "Log of so3: " << so3.Log() << "\n";
        std::cout << "log of so3: " << so3.log() << "\n";
    }

    std::cout << "=========== Exp.Log ==============\n";
    std::cout << "v: " << v << "\n";
    std::cout << "v.Exp().Log(): " << v.SO_Exp().Log() << "\n";

    std::cout << "=========== Adj ==============\n";
    std::cout << "v: " << v << "\n";
    std::cout << "adj: " << so32.adj(v) << "\n";
    std::cout << "adj2: " << so32.Adj() * v << "\n";
    return 0;
}
