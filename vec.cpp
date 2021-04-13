#include <math.h>
#include "vec.h"

/// rotates "x" around the axis "d = v - u" anchored on "u" by an angle theta.
void vec_rotate(double cos_theta, double sin_theta, const double u[3],
                const double v[3], double x[3])
{
    // See https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
    // d : normal vector in the rotation direction u->v
    double d[3] = {v[0], v[1], v[2]}; // d = v
    bool has_u = (u && vec_norm(u) > 1E-8);
    if (has_u)
    {
        vec_sub(u, d); // d = v - u
    }
    vec_normalize(d);

    // translates to the origin
    if (has_u)
    {
        vec_sub(u, x);
    }

    // w = cross(d, x)
    double w[3];
    vec_cross(d, x, w);

    // w = sin_theta * w + alpha * d
    double alpha = (1.0 - cos_theta) * vec_dot(x, d);
    vec_axpby(alpha, d, sin_theta, w);

    // x = w + cos_theta * x
    vec_xpby(w, cos_theta, x);

    // translates to the correct position
    if (has_u)
    {
        vec_add(u, x);
    }
}
