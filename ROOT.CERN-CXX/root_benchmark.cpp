// root_benchmark.cpp
// Benchmark of ROOT::Math::BrentRootFinder vs ModABRootFinder on 92 test problems.

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <cstdarg>
#include <cstring>
#include <functional>
#include <string>
#include <vector>
#include <limits>
#include <algorithm>
#include <memory>
#include <chrono>
#include <numeric>

#include "Math/BrentRootFinder.h"
#include "Math/ModABRootFinder.h"
#include "Math/Functor.h"
#include "Math/IFunction.h"

// -------------------------------------------------------------------
// Function-call counter
// -------------------------------------------------------------------
static int g_call_count = 0;

template <typename F>
double counted(F f, double x) {
    ++g_call_count;
    return f(x);
}

// -------------------------------------------------------------------
// Problem definition
// -------------------------------------------------------------------
struct Problem {
    std::string name;
    std::function<double(double)> f;
    double a, b;
    double target; // find x such that f(x) = target
};

// Helper used by f27
static double P(double x) { return x + 1.11111; }

// Sign function
static double sign(double x) {
    return (x > 0) ? 1.0 : (x < 0 ? -1.0 : 0.0);
}

// -------------------------------------------------------------------
// Build all 92 test problems
// -------------------------------------------------------------------
std::vector<Problem> build_problems() {
    std::vector<Problem> ps;
    auto add = [&](const char* name, std::function<double(double)> f,
        double a, double b, double target = 0.0) {
            ps.push_back({ name, f, a, b, target });
        };

    // -- Galdino set (f01–f33) --
    add("f01", [](double x) { return std::pow(x, 3) - 1; }, 0.5, 1.5);
    add("f02", [](double x) { return x * x * (x * x / 3 + std::sqrt(2.0) * std::sin(x)) - std::sqrt(3.0) / 18; }, 0.1, 1);
    add("f03", [](double x) { return 11 * std::pow(x, 11) - 1; }, 0.1, 1);
    add("f04", [](double x) { return std::pow(x, 3) + 1; }, -1.8, 0);
    add("f05", [](double x) { return std::pow(x, 3) - 2 * x - 5; }, 2, 3);
    add("f06", [](double x) { return 2 * x * std::exp(-5.0) + 1 - 2 * std::exp(-5 * x); }, 0, 1);
    add("f07", [](double x) { return 2 * x * std::exp(-10.0) + 1 - 2 * std::exp(-10 * x); }, 0, 1);
    add("f08", [](double x) { return 2 * x * std::exp(-20.0) + 1 - 2 * std::exp(-20 * x); }, 0, 1);
    add("f09", [](double x) { double n = 5;  return (1 + std::pow(1 - n, 2)) * x * x - std::pow(1 - n * x, 2); }, 0, 1);
    add("f10", [](double x) { double n = 10; return (1 + std::pow(1 - n, 2)) * x * x - std::pow(1 - n * x, 2); }, 0, 1);
    add("f11", [](double x) { double n = 20; return (1 + std::pow(1 - n, 2)) * x * x - std::pow(1 - n * x, 2); }, 0, 1);
    add("f12", [](double x) { return x * x - std::pow(1 - x, 5); }, 0, 1);
    add("f13", [](double x) { return x * x - std::pow(1 - x, 10); }, 0, 1);
    add("f14", [](double x) { return x * x - std::pow(1 - x, 20); }, 0, 1);
    add("f15", [](double x) { double n = 5;  return (1 + std::pow(1 - n, 4)) * x - std::pow(1 - n * x, 4); }, 0, 1);
    add("f16", [](double x) { double n = 10; return (1 + std::pow(1 - n, 4)) * x - std::pow(1 - n * x, 4); }, 0, 1);
    add("f17", [](double x) { double n = 20; return (1 + std::pow(1 - n, 4)) * x - std::pow(1 - n * x, 4); }, 0, 1);
    add("f18", [](double x) { return std::exp(-5 * x) * (x - 1) + std::pow(x, 5); }, 0, 1);
    add("f19", [](double x) { return std::exp(-10 * x) * (x - 1) + std::pow(x, 10); }, 0, 1);
    add("f20", [](double x) { return std::exp(-20 * x) * (x - 1) + std::pow(x, 20); }, 0, 1);
    add("f21", [](double x) { return x * x + std::sin(x / 5.0) - 0.25; }, 0, 1);
    add("f22", [](double x) { return x * x + std::sin(x / 10.0) - 0.25; }, 0, 1);
    add("f23", [](double x) { return x * x + std::sin(x / 20.0) - 0.25; }, 0, 1);
    add("f24", [](double x) { return (x + 2) * (x + 1) * std::pow(x - 3, 3); }, 2.6, 4.6);
    add("f25", [](double x) { return std::pow(x - 4, 5) * std::log(x); }, 3.6, 5.6);
    add("f26", [](double x) { return std::pow(std::sin(x) - x / 4.0, 3); }, 2, 4);
    add("f27", [](double x) {
        double p = P(x);
        double val = 81 - p * (108 - p * (54 - p * (12 - p)));
        double s = (p < 3) ? 1.0 : (p > 3 ? -1.0 : 0.0);
        return val * s;
        }, 1, 3);
    add("f28", [](double x) { return std::sin(std::pow(x - 7.143, 3)); }, 7, 8);
    add("f29", [](double x) { return std::exp(std::pow(x - 3, 5)) - 1; }, 2.6, 4.6);
    add("f30", [](double x) { return std::exp(std::pow(x - 3, 5)) - std::exp(x - 1); }, 4, 5);
    add("f31", [](double x) { return M_PI - 1.0 / x; }, 0.05, 5);
    add("f32", [](double x) { return 4 - std::tan(x); }, 0, 1.5);
    add("f33", [](double x) { return std::cos(x) - std::pow(x, 3); }, 0, 4);

    // -- Stage set (f34–f40) --
    add("f34", [](double x) { return std::cos(x) - x; }, -11, 9);
    add("f35", [](double x) { return std::sqrt(std::abs(x - 2.0 / 3)) * (x <= 2.0 / 3 ? 1 : -1) - 0.1; }, -11, 9);
    add("f36", [](double x) { return std::pow(std::abs(x - 2.0 / 3), 0.2) * (x <= 2.0 / 3 ? 1 : -1); }, -11, 9);
    add("f37", [](double x) { return std::pow(x - 7.0 / 9, 3) + (x - 7.0 / 9) * 1e-3; }, -11, 9);
    add("f38", [](double x) { return (x <= 1.0 / 3) ? -0.5 : 0.5; }, -11, 9);
    add("f39", [](double x) { return (x <= 1.0 / 3) ? -1e-3 : 1.0 - 1e-3; }, -11, 9);
    add("f40", [](double x) { return (x == 0) ? 0.0 : 1.0 / (x - 2.0 / 3); }, -11, 9);

    // -- Swift & Lindfield set (f41–f62) --
    add("f41", [](double x) { return 2 * x * std::exp(-5.0) - 2 * std::exp(-5 * x) + 1; }, 0, 10);
    add("f42", [](double x) { return (x * x - x - 6) * (x * x - 3 * x + 2); }, 0, M_PI);
    add("f43", [](double x) { return std::pow(x, 3); }, -1, 1.5);
    add("f44", [](double x) { return std::pow(x, 5); }, -1, 1.5);
    add("f45", [](double x) { return std::pow(x, 7); }, -1, 1.5);
    add("f46", [](double x) { return (std::exp(-5 * x) - x - 0.5) / std::pow(x, 5); }, 0.09, 0.7);
    add("f47", [](double x) { return 1 / std::sqrt(x) - 2 * std::log(5e3 * std::sqrt(x)) + 0.8; }, 0.0005, 0.5);
    add("f48", [](double x) { return 1 / std::sqrt(x) - 2 * std::log(5e7 * std::sqrt(x)) + 0.8; }, 0.0005, 0.5);
    add("f49", [](double x) { return (x <= 0) ? (-x * x * x - x - 1) : (std::cbrt(x) - x - 1); }, -1, 1);
    add("f50", [](double x) { return x * x * x - 3 * x + 3; }, -3, 2);
    add("f51", [](double x) { return std::log(x); }, 0.5, 5);
    add("f52", [](double x) { return (10 - x) * std::exp(-10 * x) - std::pow(x, 10) + 1; }, 0.5, 8);
    add("f53", [](double x) { return std::exp(std::sin(x)) - x - 1; }, 1.0, 4);
    add("f54", [](double x) { return 2 * std::sin(x) - 1; }, 0.1, M_PI / 3);
    add("f55", [](double x) { return (x - 1) * std::exp(-x); }, 0.0, 1.5);
    add("f56", [](double x) { return std::pow(x - 1, 3) - 1; }, 1.5, 3);
    add("f57", [](double x) { return std::exp(x * x + 7 * x - 30) - 1; }, 2.6, 3.5);
    add("f58", [](double x) { return std::atan(x) - 1; }, 1.0, 8);
    add("f59", [](double x) { return std::exp(x) - 2 * x - 1; }, 0.2, 3);
    add("f60", [](double x) { return std::exp(-x) - x - std::sin(x); }, 0.0, 2);
    add("f61", [](double x) { return x * x - std::pow(std::sin(x), 2) - 1; }, -1, 2);
    add("f62", [](double x) { return std::sin(x) - x / 2.0; }, M_PI / 2, M_PI);

    // -- Oliveira & Takahashi set (f63–f83) --
    add("f63", [](double x) { return x * std::exp(x) - 1; }, -1, 1);
    add("f64", [](double x) { return std::tan(x - 0.1); }, -1, 1);
    add("f65", [](double x) { return std::sin(x) + 0.5; }, -1, 1);
    add("f66", [](double x) { return 4 * std::pow(x, 5) + x * x + 1; }, -1, 1);
    add("f67", [](double x) { return x + std::pow(x, 10) - 1; }, -1, 1);
    add("f68", [](double x) { return std::pow(M_PI, x) - M_E; }, -1, 1);
    add("f69", [](double x) { return std::log(std::abs(x - 10.0 / 9)); }, -1, 1);
    add("f70", [](double x) {
        return 1.0 / 3 + sign(x) * std::pow(std::abs(x), 1.0 / 3) + x * x * x;
        }, -1, 1);
    add("f71", [](double x) { return (x + 2.0 / 3) / (x + 101.0 / 100); }, -1, 1);
    add("f72", [](double x) { return std::pow(x * 1e6 - 1, 3); }, -1, 1);
    add("f73", [](double x) { return std::exp(x) * std::pow(x * 1e6 - 1, 3); }, -1, 1);
    add("f74", [](double x) { return std::pow(x - 1.0 / 3, 2) * std::atan(x - 1.0 / 3); }, -1, 1);
    add("f75", [](double x) {
        double t = 3 * x - 1;
        return sign(t) * (1 - std::sqrt(std::max(0.0, 1 - t * t / 81.0)));
        }, -1, 1);
    add("f76", [](double x) { return (x > (1.0 - 1e6) / 1e6) ? (1.0 + 1e6) / 1e6 : -1.0; }, -1, 1);
    add("f77", [](double x) { return (x != 1.0 / 21) ? 1.0 / (21 * x - 1) : 0.0; }, -1, 1);
    add("f78", [](double x) { return x * x / 4.0 + std::ceil(x / 2.0) - 0.5; }, -1, 1);
    add("f79", [](double x) { return std::ceil(10 * x - 1) + 0.5; }, -1, 1);
    add("f80", [](double x) { return x + std::sin(x * 1e6) / 10.0 + 1e-3; }, -1, 1);
    add("f81", [](double x) { return (x > -1) ? 1 + std::sin(1.0 / (x + 1)) : -1.0; }, -1, 1);
    add("f82", [](double x) { return 202 * x - 2 * std::floor((2 * x + 1e-2) / (2e-2)) - 0.1; }, -1, 1);
    add("f83", [](double x) { return std::pow(202 * x - 2 * std::floor((2 * x + 1e-2) / (2e-2)) - 0.1, 3); }, -1, 1);

    // -- SciML set (f84–f92) --
    add("f84", [](double x) { return (x - 1) * (x - 2) * (x - 3) * (x - 4) * (x - 5) - 0.05; }, 0.5, 5.5);
    add("f85", [](double x) { return std::sin(x) - 0.5 * x - 0.3; }, -10.0, 10.0);
    add("f86", [](double x) { return std::exp(x) - 1 - x - x * x / 2 - 0.005; }, -2.0, 2.0);
    add("f87", [](double x) { return 1.0 / (x - 0.5) - 2 - 0.05; }, 0.6, 2.0);
    add("f88", [](double x) { return std::log(x) - x + 2 - 0.05; }, 0.1, 3.0);
    add("f89", [](double x) { return std::sin(20 * x) + 0.1 * x - 0.1; }, -4.0, 5.0);
    add("f90", [](double x) { return x * x * x - 2 * x * x + x - 0.025; }, -1.0, 2.0);
    add("f91", [](double x) { return x * std::sin(1.0 / x) - 0.1 - 0.01; }, 0.01, 1.0);
    add("f92", [](double x) { return x * x * x - 0.001; }, -10, 10);

    return ps;
}

// -------------------------------------------------------------------
// Solver abstraction
// -------------------------------------------------------------------

struct SolverResult {
    double root;
    bool   ok;
};

struct SolverBase {
    std::string name;
    virtual ~SolverBase() = default;
    virtual SolverResult solve(const std::function<double(double)>& f,
        double a, double b, double target, double eps) = 0;
};
const int MAX_ITER = 200;
// Wrapper for ROOT::Math::BrentRootFinder (no GSL needed)
struct BrentSolver : SolverBase {
    BrentSolver() { name = "Brent"; }

    SolverResult solve(const std::function<double(double)>& f,
        double a, double b, double target, double eps) override {
        auto g = [&](double x) -> double { return counted(f, x) - target; };
        double lo = std::min(a, b), hi = std::max(a, b);

        ROOT::Math::BrentRootFinder finder;
        ROOT::Math::Functor1D func(g);
        if (!finder.SetFunction(func, lo, hi))
            return { std::numeric_limits<double>::quiet_NaN(), false };
        if (!finder.Solve(MAX_ITER, eps, 4 * std::numeric_limits<double>::epsilon()))
            return { std::numeric_limits<double>::quiet_NaN(), false };
        return { finder.Root(), true };
    }
};

// Wrapper for ROOT::Math::ModABRootFinder (no GSL needed)
struct ModABSolver : SolverBase {
    ModABSolver() { name = "ModAB"; }

    SolverResult solve(const std::function<double(double)>& f,
        double a, double b, double target, double eps) override {
        auto g = [&](double x) -> double { return counted(f, x) - target; };
        double lo = std::min(a, b), hi = std::max(a, b);

        ROOT::Math::ModABRootFinder finder;
        ROOT::Math::Functor1D func(g);
        if (!finder.SetFunction(func, lo, hi))
            return { std::numeric_limits<double>::quiet_NaN(), false };
        if (!finder.Solve(MAX_ITER, eps, 4 * std::numeric_limits<double>::epsilon()))
            return { std::numeric_limits<double>::quiet_NaN(), false };
        return { finder.Root(), true };
    }
};

// -------------------------------------------------------------------
// Dual output: prints to both stdout and a results file
// -------------------------------------------------------------------
static FILE* g_outfile = nullptr;

void tee(const char* fmt, ...) {
    va_list args;
    va_start(args, fmt);
    std::vprintf(fmt, args);
    va_end(args);

    if (g_outfile) {
        va_start(args, fmt);
        std::vfprintf(g_outfile, fmt, args);
        va_end(args);
    }
}

// -------------------------------------------------------------------
// Main benchmark
// -------------------------------------------------------------------
int main() {
    g_outfile = std::fopen("root_benchmark_results.txt", "w");
    if (!g_outfile) {
        std::fprintf(stderr, "Error: cannot open root_benchmark_results.txt for writing\n");
        return 1;
    }

    auto problems = build_problems();

    std::vector<std::unique_ptr<SolverBase>> solvers;
    solvers.emplace_back(new BrentSolver());
    solvers.emplace_back(new ModABSolver());

    const double eps = 1e-14;
    const int col_w = 22;
    const int name_w = 4;
    const int cnt_w = 10;
    const size_t nSolvers = solvers.size();

    tee("=====================================================================\n");
    tee("  ROOT-FINDING BENCHMARK  (%zu problems, tol = %.0e, maxIter = %d)\n", problems.size(), eps, MAX_ITER);
    tee("  Solvers: Brent (MathCore), ModAB\n");
    tee("=====================================================================\n\n");

    // -- Roots found --
    tee("-- Roots found -----------------------------------------------------\n");
    tee("%*s", name_w, "Func");
    for (auto& s : solvers) tee(" | %*s", col_w, s->name.c_str());
    tee("\n");
    for (auto& p : problems) {
        tee("%*s", name_w, p.name.c_str());
        for (auto& s : solvers) {
            g_call_count = 0;
            auto res = s->solve(p.f, p.a, p.b, p.target, eps);
            if (!res.ok)
                tee(" | %*s", col_w, "FAIL");
            else
                tee(" | %*.15g", col_w, res.root);
        }
        tee("\n");
    }

    // -- Function values at roots --
    tee("\n-- f(root) - target ------------------------------------------------\n");
    tee("%*s", name_w, "Func");
    for (auto& s : solvers) tee(" | %*s", col_w, s->name.c_str());
    tee("\n");
    for (auto& p : problems) {
        tee("%*s", name_w, p.name.c_str());
        for (auto& s : solvers) {
            g_call_count = 0;
            auto res = s->solve(p.f, p.a, p.b, p.target, eps);
            if (!res.ok) {
                tee(" | %*s", col_w, "FAIL");
            }
            else {
                double fval = p.f(res.root) - p.target;
                tee(" | %*.6e", col_w, fval);
            }
        }
        tee("\n");
    }

    // -- Evaluation counts --
    tee("\n-- Function evaluations --------------------------------------------\n");
    tee("%*s", name_w, "Func");
    for (auto& s : solvers) tee(" | %*s", cnt_w, s->name.c_str());
    tee("\n");

    std::vector<int> totals(nSolvers, 0);
    std::vector<int> failures(nSolvers, 0);

    for (auto& p : problems) {
        tee("%*s", name_w, p.name.c_str());
        for (size_t j = 0; j < nSolvers; ++j) {
            g_call_count = 0;
            auto res = solvers[j]->solve(p.f, p.a, p.b, p.target, eps);
            if (!res.ok) {
                tee(" | %*s", cnt_w, "FAIL");
                ++failures[j];
                totals[j] += MAX_ITER;
            }
            else {
                tee(" | %*d", cnt_w, g_call_count);
                totals[j] += g_call_count;
            }
        }
        tee("\n");
    }

    // Summary for evaluation counts
    tee("--------------------------------------------------------------------\n");
    tee("%*s", name_w, "SUM");
    for (size_t j = 0; j < nSolvers; ++j)
        tee(" | %*d", cnt_w, totals[j]);
    tee("\n");

    tee("%*s", name_w, "FAIL");
    for (size_t j = 0; j < nSolvers; ++j)
        tee(" | %*d", cnt_w, failures[j]);
    tee("\n");

    tee("%*s", name_w, "AVG");
    for (size_t j = 0; j < nSolvers; ++j) {
        int solved = (int)problems.size() - failures[j];
        double avg = solved > 0 ? (double)totals[j] / solved : 0;
        tee(" | %*.1f", cnt_w, avg);
    }
    tee("\n");

    // -- Per-problem timing --
    // Each problem is solved REPEAT times to get stable measurements.
    const int REPEAT = 100000;
    const int time_w = 12;

    tee("\n-- Time per problem (microseconds, avg over %d runs) -------------\n", REPEAT);
    tee("%*s", name_w, "Func");
    for (auto& s : solvers) tee(" | %*s", time_w, s->name.c_str());
    tee("\n");

    // Accumulate total time per solver (in microseconds)
    std::vector<double> total_us(nSolvers, 0.0);
    std::vector<int>    time_failures(nSolvers, 0);

    for (auto& p : problems) {
        tee("%*s", name_w, p.name.c_str());
        for (size_t j = 0; j < nSolvers; ++j) {
            // Quick check: can this solver handle this problem?
            g_call_count = 0;
            auto test = solvers[j]->solve(p.f, p.a, p.b, p.target, eps);
            if (!test.ok) {
                tee(" | %*s", time_w, "FAIL");
                ++time_failures[j];
                continue;
            }

            // Timed loop
            auto t0 = std::chrono::high_resolution_clock::now();
            for (int r = 0; r < REPEAT; ++r) {
                g_call_count = 0;
                solvers[j]->solve(p.f, p.a, p.b, p.target, eps);
            }
            auto t1 = std::chrono::high_resolution_clock::now();

            double us = std::chrono::duration<double, std::micro>(t1 - t0).count() / REPEAT;
            total_us[j] += us;
            tee(" | %*.3f", time_w, us);
        }
        tee("\n");
    }

    // Timing summary
    tee("--------------------------------------------------------------------\n");
    tee("%*s", name_w, "TOT");
    for (size_t j = 0; j < nSolvers; ++j)
        tee(" | %*.3f", time_w, total_us[j]);
    tee("  (total us across all solved problems)\n");

    tee("%*s", name_w, "AVG");
    for (size_t j = 0; j < nSolvers; ++j) {
        int solved = (int)problems.size() - time_failures[j];
        double avg = solved > 0 ? total_us[j] / solved : 0;
        tee(" | %*.3f", time_w, avg);
    }
    tee("  (avg us per problem)\n");

    // -- Aggregate timing: total wall-clock for all problems --
    tee("\n-- Total wall-clock time (all %zu problems x %d repeats) ----------\n",
        problems.size(), REPEAT);
    for (size_t j = 0; j < nSolvers; ++j) {
        // Time the entire suite in one shot
        auto t0 = std::chrono::high_resolution_clock::now();
        for (int r = 0; r < REPEAT; ++r) {
            for (auto& p : problems) {
                g_call_count = 0;
                solvers[j]->solve(p.f, p.a, p.b, p.target, eps);
            }
        }
        auto t1 = std::chrono::high_resolution_clock::now();

        double total_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        double per_call_us = (total_ms * 1000.0) / (REPEAT * problems.size());
        tee("  %-12s: %10.2f ms total  |  %.3f us/problem avg\n",
            solvers[j]->name.c_str(), total_ms, per_call_us);
    }
    tee("\n");

    std::fclose(g_outfile);
    std::printf("Results written to root_benchmark_results.txt\n");

    return 0;
}