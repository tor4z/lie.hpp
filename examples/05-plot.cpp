#define GNUPLOT_IMPLEMENTATION
#include "gnuplot.hpp"

int main()
{
    gp::Plotter p;
    p.set_title("hello");
    std::vector<float> x{1, 2, 3};
    std::vector<float> y{1, 2, 3};
    p.line(x, y);
    p.hold(true);
    p.function("sin(x)");
    return 0;
}
