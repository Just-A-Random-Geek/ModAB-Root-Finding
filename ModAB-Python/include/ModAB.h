#ifndef MODAB_H
#define MODAB_H

#ifdef _WIN32
    #ifdef MODAB_EXPORTS
        #define MODAB_API __declspec(dllexport)
    #else
        #define MODAB_API __declspec(dllimport)
    #endif
#else
    #define MODAB_API __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C" {
#endif

// Get the number of function evaluations from the last root-finding call
MODAB_API int get_evaluation_count(void);

// Finds the root of "F(x) = 0" within the interval [x1, x2]
// with the specified precisions - absolute: aTol and relative: rTol,
// using an improved version of the modified Anderson Bjork's method.
// F(x) must be continuous and sign(F(x1)) != sign(F(x2))
// Returns NAN if no root is found or inputs are invalid
MODAB_API double modAB_find_root(double (*f)(double), double x1, double x2,
                                  double aTol, double rTol, int maxIter);

// Convenience wrapper with default tolerances (aTol=1e-14, rTol=1e-14, maxIter=200)
MODAB_API double modAB_default(double (*f)(double), double x1, double x2);

#ifdef __cplusplus
}
#endif

#endif // MODAB_H
