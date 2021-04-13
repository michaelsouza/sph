#pragma once

#include <math.h>
#include <stdio.h>
#include "vec.h"

/// set the values of x
inline void vec_set(double x[3], double x0, double x1, double x2) {
    x[0] = x0;
    x[1] = x1;
    x[2] = x2;
}

/// calculates the euclidean norm of v
inline double vec_norm(const double v[3]) {
    return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

/// y = a*x+b*y
inline void vec_axpby(double a, const double x[3], double b, double y[3]) {
    y[0] = a * x[0] + b * y[0];
    y[1] = a * x[1] + b * y[1];
    y[2] = a * x[2] + b * y[2];
}

/// z = a*x+b*y
inline void vec_axpby(double a, const double x[3], double b, const double y[3], double z[3]) {
    z[0] = a * x[0] + b * y[0];
    z[1] = a * x[1] + b * y[1];
    z[2] = a * x[2] + b * y[2];
}

/// y += x
inline void vec_add(const double x[3], double y[3]) {
    y[0] += x[0];
    y[1] += x[1];
    y[2] += x[2];
}

/// y -= x
inline void vec_sub(const double x[3], double y[3]) {
    y[0] -= x[0];
    y[1] -= x[1];
    y[2] -= x[2];
}

/// z = x - y
inline void vec_sub(const double x[3], const double y[3], double z[3]){
    z[0] = x[0] - y[0];
    z[1] = x[1] - y[1];
    z[2] = x[2] - y[2];
}

/// y = a*x+y
inline void vec_axpy(double a, const double x[3], double y[3]) {
    y[0] += a * x[0];
    y[1] += a * x[1];
    y[2] += a * x[2];
}

/// y = x+b*y
inline void vec_xpby(const double x[3], double b, double y[3]) {
    y[0] = x[0] + b * y[0];
    y[1] = x[1] + b * y[1];
    y[2] = x[2] + b * y[2];
}

/// \brief v = u
inline void vec_copy(const double u[3], double v[3]) {
    v[0] = u[0];
    v[1] = u[1];
    v[2] = u[2];
}

/// inner product
inline double vec_dot(const double u[3], const double v[3]) {
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
}

/// cross product w = u x v
inline void vec_cross(const double u[3], const double v[3], double w[3]) {
    w[0] = u[1] * v[2] - u[2] * v[1];
    w[1] = u[2] * v[0] - u[0] * v[2];
    w[2] = u[0] * v[1] - u[1] * v[0];
}

/// scale vector x = a * x
inline void vec_scale(double a, double x[3]) {
    x[0] *= a;
    x[1] *= a;
    x[2] *= a;
}

/// scale vector y = a * x
inline void vec_scale(double a, const double x[3], double y[3]) {
    y[0] = a * x[0];
    y[1] = a * x[1];
    y[2] = a * x[2];
}

/// print vector
inline void vec_print(const char *name, double v[3]) {
    printf("%s = {% 3.6f % 3.6f % 3.6f };\n", name, v[0], v[1], v[2]);
}

/// normalizes the vector v
inline void vec_normalize(double v[3]) {
    double norm_v = vec_norm(v);
    v[0] /= norm_v;
    v[1] /= norm_v;
    v[2] /= norm_v;
}

/// rotates "x" around the axis "d = v - u" anchored on "u" by an angle theta.
void vec_rotate(double cos_theta, double sin_theta, const double u[3],
                const double v[3], double x[3]);

/// calculates the euclidean distance from x to y
inline double vec_dist(const double x[3], const double y[3]) {
    double u[3] = { x[0] - y[0], x[1] - y[1], x[2] - y[2] };
    return vec_norm(u);
}

/// normal to the plane (a,b,c), ie, n = (b - a) x (c - b)
inline void vec_plane_normal(const double a[3], const double b[3],
                             const double c[3], double n[3]) {
    // u = c - b
    const double u0 = c[0] - b[0];
    const double u1 = c[1] - b[1];
    const double u2 = c[2] - b[2];

    // v = a - b
    const double v0 = a[0] - b[0];
    const double v1 = a[1] - b[1];
    const double v2 = a[2] - b[2];

    // n = u x v
    n[0] = u1 * v2 - u2 * v1;
    n[1] = u2 * v0 - u0 * v2;
    n[2] = u0 * v1 - u1 * v0;

    vec_normalize(n);
}
