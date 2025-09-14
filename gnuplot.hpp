#ifndef GNUPLOT_HPP_
#define GNUPLOT_HPP_

#include <cstddef>
#include <cstdio>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>

namespace gp {

struct Plotter
{
    explicit Plotter(bool persist=true);
    ~Plotter();
    inline void hold(bool h) { hold_ = h; }

    bool set_title(const char* title);
    bool line(std::vector<float> x, std::vector<float> y);
    bool line(float* x, float* y, size_t len);
    bool function(const char* fun);
    bool plot(const char* plot_cmd);
    bool cmd(const char* c);
private:
    const char* plot_cmd() const;

    FILE *gp_fp_;
    bool hold_;
    bool has_plot_;
}; // struct Plotter

}; // namespace gp

#endif // GNUPLOT_HPP_


#define GNUPLOT_IMPLEMENTATION // delete me


#ifdef GNUPLOT_IMPLEMENTATION
#ifndef GNUPLOT_CPP_
#define GNUPLOT_CPP_

namespace gp {

#define GP_MAX_PLOT_CMD_LEN 128

std::string unique_name()
{
    static int cnt{0};
    static std::mutex locker;
    std::stringstream ss;

    {
        std::lock_guard<std::mutex> guard{locker};
        ss << "gp_" <<getpid() << "_" << cnt++ << + "_dat";
    }
    return ss.str();
}

Plotter::Plotter(bool persist)
    : gp_fp_(nullptr)
    , hold_(false)
    , has_plot_(false)
{
#define GP_MAX_POPEN_CMD_LEN 64
    char buff[GP_MAX_POPEN_CMD_LEN];
    snprintf(buff, GP_MAX_POPEN_CMD_LEN, "gnuplot %s", (persist ? "-persist" : ""));
    gp_fp_ = popen(buff, "w");
#undef GP_MAX_POPEN_CMD_LEN
}

Plotter::~Plotter()
{
    if (gp_fp_) {
        fclose(gp_fp_);
        gp_fp_ = nullptr;
    }
}

bool Plotter::set_title(const char* title)
{
    if (!title) return false;

    char c[GP_MAX_PLOT_CMD_LEN];
    snprintf(c, GP_MAX_PLOT_CMD_LEN, "set title \"%s\"\n", title);
    return cmd(c);
}

bool Plotter::line(std::vector<float> x, std::vector<float> y)
{
    if (x.size() != y.size()) return false;
    return line(x.data(), y.data(), x.size());
}

bool Plotter::line(float* x, float* y, size_t len)
{
    if (!gp_fp_) return false;
    if (!x || !y || len == 0) return false;

    auto data_name{unique_name()};
    fprintf(gp_fp_, "$%s<<EOD\n", data_name.c_str());
    for (size_t i = 0; i < len; ++i) {
        fprintf(gp_fp_, "%f %f\n", x[i], y[i]);
    }
    fprintf(gp_fp_, "EOD\n");

    char c[GP_MAX_PLOT_CMD_LEN];
    snprintf(c, GP_MAX_PLOT_CMD_LEN, "%s $%s with line\n", plot_cmd(), data_name.c_str());
    return plot(c);
}

bool Plotter::function(const char* fun)
{
    if (!fun) return false;

    char c[GP_MAX_PLOT_CMD_LEN];
    snprintf(c, GP_MAX_PLOT_CMD_LEN, "%s %s\n", plot_cmd(), fun);
    return plot(c);
}

bool Plotter::plot(const char* plot_cmd)
{
    if (cmd(plot_cmd)) {
        has_plot_ = true;
        return true;
    } else {
        return false;
    }
}

bool Plotter::cmd(const char* c)
{
    if (!gp_fp_) return false;
    if (!c) return false;

    fprintf(gp_fp_, "%s\n", c);
    return true;
}

const char* Plotter::plot_cmd() const
{
    static const char* plot{"plot"};
    static const char* replot{"replot"};

    if (hold_ && has_plot_) {
        return replot;
    } else {
        return plot;
    }
}

} // namespace gp

#endif // GNUPLOT_CPP_
#endif // GNUPLOT_IMPLEMENTATION
