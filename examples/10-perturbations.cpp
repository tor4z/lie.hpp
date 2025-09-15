#include <iostream>
#define LIE_IMPLEMENTATION
#include "lie.hpp"

int main()
{
    auto id{lie::SO3f::eye()};
    lie::Vector3f v;
    v << 0.1, 0, 0;

    std::cout << "Identity so3: " << id << "\n";
    std::cout << "Identity so3 inv: " << id.inv() << "\n";
    std::cout << "pert v: " << v << "\n";
    std::cout << "pert v Exp: " << v.Exp() << "\n";

    auto id_with_perb{id.plus(v)};
    std::cout << "Identity so3 with pert: " << id_with_perb << "\n";
    std::cout << "The pert: " << id_with_perb.minus(id) << "\n";

    return 0;
}
