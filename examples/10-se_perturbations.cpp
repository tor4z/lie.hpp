#include <iostream>
#define LIE_IMPLEMENTATION
#include "lie.hpp"

int main()
{
    auto id{lie::SE3f::eye()};
    lie::Vector6f v;
    v << 0.1, 0, 0, 0.2, 0, 0;

    std::cout << "Identity se3: " << id << "\n";
    std::cout << "Identity se3 inv: " << id.inv() << "\n";
    std::cout << "pert v: " << v << "\n";
    std::cout << "pert v Exp: " << v.SE_Exp() << "\n";

    std::cout << "======== right plus/minus ===========\n";
    auto id_with_perb{id.plus(v)};
    std::cout << "id_with_perb: " << id_with_perb << "\n";
    std::cout << "The pert: " << id_with_perb.minus(id) << "\n";

    return 0;
}
