#include <iostream>
#define LIE_IMPLEMENTATION
#include "lie.hpp"

int main()
{
    lie::Vector<float, 3> v;
    v << 0.1, 0.2, 0.4;
    std::cout << "============ exp ==========\n";
    {
        auto m{v.hat()};
        auto so3{m.exp()};
        std::cout << "so3: " << so3 << "\n";
        std::cout << "so3^T: " << so3.t() << "\n";
        std::cout << "so3 * so3.t(): " << so3 * so3.t() << "\n";
    }

    std::cout << "============ Exp ==========\n";
    {
        auto so3{v.Exp()};
        std::cout << "so3: " << so3 << "\n";
        std::cout << "so3^T: " << so3.t() << "\n";
        std::cout << "so3 * so3.t(): " << so3 * so3.t() << "\n";
    }

    std::cout << "=========== Log ==============\n";
    {
        auto so3{v.Exp()};
        std::cout << "v: " << v << "\n";
        std::cout << "so3: " << so3 << "\n";
        std::cout << "Log of so3: " << so3.Log() << "\n";
        std::cout << "log of so3: " << so3.log() << "\n";
    }
    return 0;
}
