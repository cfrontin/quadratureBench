
#ifndef MATH_H
#define MATH_H

#include "types.h"

namespace quadmath
{

void invert3x3(const std::vector<std::vector<quad>> &A,
               std::vector<std::vector<quad>> &Ainv);

void invert4x4(const std::vector<std::vector<quad>> &A,
               std::vector<std::vector<quad>> &Ainv);

void invert3x3(
    const quad &A_11, const quad &A_12, const quad &A_13,
    const quad &A_21, const quad &A_22, const quad &A_23,
    const quad &A_31, const quad &A_32, const quad &A_33,
    quad &Ainv_11, quad &Ainv_12, quad &Ainv_13,
    quad &Ainv_21, quad &Ainv_22, quad &Ainv_23,
    quad &Ainv_31, quad &Ainv_32, quad &Ainv_33);

void invert4x4(
    const quad &A_11, const quad &A_12, const quad &A_13, const quad &A_14,
    const quad &A_21, const quad &A_22, const quad &A_23, const quad &A_24,
    const quad &A_31, const quad &A_32, const quad &A_33, const quad &A_34,
    const quad &A_41, const quad &A_42, const quad &A_43, const quad &A_44,
    quad &Ainv_11, quad &Ainv_12, quad &Ainv_13, quad &Ainv_14,
    quad &Ainv_21, quad &Ainv_22, quad &Ainv_23, quad &Ainv_24,
    quad &Ainv_31, quad &Ainv_32, quad &Ainv_33, quad &Ainv_34,
    quad &Ainv_41, quad &Ainv_42, quad &Ainv_43, quad &Ainv_44);

void phys2bary3d(const std::vector<quad> &ref,
                 const std::vector<quad> &r1,
                 const std::vector<quad> &r2,
                 const std::vector<quad> &r3,
                 const std::vector<quad> &r4,
                 std::vector<quad> &r);

void bary2phys3d(const std::vector<quad> &r,
                 const std::vector<quad> &r1,
                 const std::vector<quad> &r2,
                 const std::vector<quad> &r3,
                 const std::vector<quad> &r4,
                 std::vector<quad> &ref);

void phys2ref3d(const std::vector<quad> &r,
                const std::vector<quad> &r1,
                const std::vector<quad> &r2,
                const std::vector<quad> &r3,
                const std::vector<quad> &r4,
                std::vector<quad> &ref);

void ref2phys3d(const std::vector<quad> &ref,
                const std::vector<quad> &r1,
                const std::vector<quad> &r2,
                const std::vector<quad> &r3,
                const std::vector<quad> &r4,
                std::vector<quad> &r);

void phys2bary4d(const std::vector<quad> &lambda,
                 const std::vector<quad> &r1,
                 const std::vector<quad> &r2,
                 const std::vector<quad> &r3,
                 const std::vector<quad> &r4,
                 const std::vector<quad> &r5,
                 std::vector<quad> &r);

void bary2phys4d(const std::vector<quad> &r,
                 const std::vector<quad> &r1,
                 const std::vector<quad> &r2,
                 const std::vector<quad> &r3,
                 const std::vector<quad> &r4,
                 const std::vector<quad> &r5,
                 std::vector<quad> &lambda);

void phys2ref4d(const std::vector<quad> &r,
                const std::vector<quad> &r1,
                const std::vector<quad> &r2,
                const std::vector<quad> &r3,
                const std::vector<quad> &r4,
                const std::vector<quad> &r5,
                std::vector<quad> &ref);

void ref2phys4d(const std::vector<quad> &ref,
                const std::vector<quad> &r1,
                const std::vector<quad> &r2,
                const std::vector<quad> &r3,
                const std::vector<quad> &r4,
                const std::vector<quad> &r5,
                std::vector<quad> &r);

} // namespace quadmath

#endif // MATH_H
