
#include "math.h"

namespace quadmath
{

void invert3x3(const std::vector<std::vector<quad>> &A,
               std::vector<std::vector<quad>> &Ainv)
{
  assert(A.size() == 3);
  for (auto Ai : A)
    assert(Ai.size() == 3);
  assert(Ainv.size() == 3);
  for (auto Ainv_i : Ainv)
    assert(Ainv_i.size() == 3);

  invert3x3(
      A[0][0], A[0][1], A[0][2],
      A[1][0], A[1][1], A[1][2],
      A[2][0], A[2][1], A[2][2],
      Ainv[0][0], Ainv[0][1], Ainv[0][2],
      Ainv[1][0], Ainv[1][1], Ainv[1][2],
      Ainv[2][0], Ainv[2][1], Ainv[2][2]);
};

void invert4x4(const std::vector<std::vector<quad>> &A,
               std::vector<std::vector<quad>> &Ainv)
{
  assert(A.size() == 4);
  for (auto Ai : A)
    assert(Ai.size() == 4);
  assert(Ainv.size() == 4);
  for (auto Ainv_i : Ainv)
    assert(Ainv_i.size() == 4);

  invert4x4(
      A[0][0], A[0][1], A[0][2], A[0][3],
      A[1][0], A[1][1], A[1][2], A[1][3],
      A[2][0], A[2][1], A[2][2], A[2][3],
      A[3][0], A[3][1], A[3][2], A[3][3],
      Ainv[0][0], Ainv[0][1], Ainv[0][2], Ainv[0][3],
      Ainv[1][0], Ainv[1][1], Ainv[1][2], Ainv[1][3],
      Ainv[2][0], Ainv[2][1], Ainv[2][2], Ainv[2][3],
      Ainv[3][0], Ainv[3][1], Ainv[3][2], Ainv[3][3]);
};

void invert3x3(
    const quad &a_1, const quad &a_2, const quad &a_3,
    const quad &b_1, const quad &b_2, const quad &b_3,
    const quad &c_1, const quad &c_2, const quad &c_3,
    quad &T_11, quad &T_12, quad &T_13,
    quad &T_21, quad &T_22, quad &T_23,
    quad &T_31, quad &T_32, quad &T_33)
{

  quad detT= -(a_3*b_2*c_1) + a_2*b_3*c_1 + a_3*b_1*c_2 - a_1*b_3*c_2
      - a_2*b_1*c_3 + a_1*b_2*c_3;

  T_11= (-(b_3*c_2) + b_2*c_3)/detT;
  T_12= (a_3*c_2 - a_2*c_3)/detT;
  T_13= (-(a_3*b_2) + a_2*b_3)/detT;

  T_21= (b_3*c_1 - b_1*c_3)/detT;
  T_22= (-(a_3*c_1) + a_1*c_3)/detT;
  T_23= (a_3*b_1 - a_1*b_3)/detT;

  T_31= (-(b_2*c_1) + b_1*c_2)/detT;
  T_32= (a_2*c_1 - a_1*c_2)/detT;
  T_33= (-(a_2*b_1) + a_1*b_2)/detT;

};

void invert4x4(
    const quad &a_1, const quad &a_2, const quad &a_3, const quad &a_4,
    const quad &b_1, const quad &b_2, const quad &b_3, const quad &b_4,
    const quad &c_1, const quad &c_2, const quad &c_3, const quad &c_4,
    const quad &d_1, const quad &d_2, const quad &d_3, const quad &d_4,
    quad &T_11, quad &T_12, quad &T_13, quad &T_14,
    quad &T_21, quad &T_22, quad &T_23, quad &T_24,
    quad &T_31, quad &T_32, quad &T_33, quad &T_34,
    quad &T_41, quad &T_42, quad &T_43, quad &T_44)
{

  quad detT=
      a_4*b_3*c_2*d_1 - a_3*b_4*c_2*d_1 - a_4*b_2*c_3*d_1 + a_2*b_4*c_3*d_1
      + a_3*b_2*c_4*d_1 - a_2*b_3*c_4*d_1 - a_4*b_3*c_1*d_2 + a_3*b_4*c_1*d_2
      + a_4*b_1*c_3*d_2 - a_1*b_4*c_3*d_2 - a_3*b_1*c_4*d_2 + a_1*b_3*c_4*d_2
      + a_4*b_2*c_1*d_3 - a_2*b_4*c_1*d_3 - a_4*b_1*c_2*d_3 + a_1*b_4*c_2*d_3
      + a_2*b_1*c_4*d_3 - a_1*b_2*c_4*d_3 - a_3*b_2*c_1*d_4 + a_2*b_3*c_1*d_4
      + a_3*b_1*c_2*d_4 - a_1*b_3*c_2*d_4 - a_2*b_1*c_3*d_4 + a_1*b_2*c_3*d_4;

  T_11= (-(b_4*c_3*d_2) + b_3*c_4*d_2 + b_4*c_2*d_3 - b_2*c_4*d_3 - b_3*c_2*d_4 + b_2*c_3*d_4)/detT;
  T_12= (a_4*c_3*d_2 - a_3*c_4*d_2 - a_4*c_2*d_3 + a_2*c_4*d_3 + a_3*c_2*d_4 - a_2*c_3*d_4)/detT;
  T_13= (-(a_4*b_3*d_2) + a_3*b_4*d_2 + a_4*b_2*d_3 - a_2*b_4*d_3 - a_3*b_2*d_4 + a_2*b_3*d_4)/detT;
  T_14= (a_4*b_3*c_2 - a_3*b_4*c_2 - a_4*b_2*c_3 + a_2*b_4*c_3 + a_3*b_2*c_4 - a_2*b_3*c_4)/detT;

  T_21= (b_4*c_3*d_1 - b_3*c_4*d_1 - b_4*c_1*d_3 + b_1*c_4*d_3 + b_3*c_1*d_4 - b_1*c_3*d_4)/detT;
  T_22= (-(a_4*c_3*d_1) + a_3*c_4*d_1 + a_4*c_1*d_3 - a_1*c_4*d_3 - a_3*c_1*d_4 + a_1*c_3*d_4)/detT;
  T_23= (a_4*b_3*d_1 - a_3*b_4*d_1 - a_4*b_1*d_3 + a_1*b_4*d_3 + a_3*b_1*d_4 - a_1*b_3*d_4)/detT;
  T_24= (-(a_4*b_3*c_1) + a_3*b_4*c_1 + a_4*b_1*c_3 - a_1*b_4*c_3 - a_3*b_1*c_4 + a_1*b_3*c_4)/detT;

  T_31= (-(b_4*c_2*d_1) + b_2*c_4*d_1 + b_4*c_1*d_2 - b_1*c_4*d_2 - b_2*c_1*d_4 + b_1*c_2*d_4)/detT;
  T_32= (a_4*c_2*d_1 - a_2*c_4*d_1 - a_4*c_1*d_2 + a_1*c_4*d_2 + a_2*c_1*d_4 - a_1*c_2*d_4)/detT;
  T_33= (-(a_4*b_2*d_1) + a_2*b_4*d_1 + a_4*b_1*d_2 - a_1*b_4*d_2 - a_2*b_1*d_4 + a_1*b_2*d_4)/detT;
  T_34= (a_4*b_2*c_1 - a_2*b_4*c_1 - a_4*b_1*c_2 + a_1*b_4*c_2 + a_2*b_1*c_4 - a_1*b_2*c_4)/detT;

  T_41= (b_3*c_2*d_1 - b_2*c_3*d_1 - b_3*c_1*d_2 + b_1*c_3*d_2 + b_2*c_1*d_3 - b_1*c_2*d_3)/detT;
  T_42= (-(a_3*c_2*d_1) + a_2*c_3*d_1 + a_3*c_1*d_2 - a_1*c_3*d_2 - a_2*c_1*d_3 + a_1*c_2*d_3)/detT;
  T_43= (a_3*b_2*d_1 - a_2*b_3*d_1 - a_3*b_1*d_2 + a_1*b_3*d_2 + a_2*b_1*d_3 - a_1*b_2*d_3)/detT;
  T_44= (-(a_3*b_2*c_1) + a_2*b_3*c_1 + a_3*b_1*c_2 - a_1*b_3*c_2 - a_2*b_1*c_3 + a_1*b_2*c_3)/detT;

};

void phys2bary3d(const std::vector<quad> &r,
                 const std::vector<quad> &r1,
                 const std::vector<quad> &r2,
                 const std::vector<quad> &r3,
                 const std::vector<quad> &r4,
                 std::vector<quad> &lambda)
{
  assert(r.size() == 3);
  assert(r1.size() == 3);
  assert(r2.size() == 3);
  assert(r3.size() == 3);
  assert(r4.size() == 3);
  assert(lambda.size() == 4);

  std::vector<std::vector<quad>> T= {{0.0Q, 0.0Q, 0.0Q}, {0.0Q, 0.0Q, 0.0Q}, {0.0Q, 0.0Q, 0.0Q}};
  std::vector<std::vector<quad>> Tinv= {{0.0Q, 0.0Q, 0.0Q}, {0.0Q, 0.0Q, 0.0Q}, {0.0Q, 0.0Q, 0.0Q}};

  T[0][0]= r1[0] - r4[0];
  T[0][1]= r2[0] - r4[0];
  T[0][2]= r3[0] - r4[0];
  T[1][0]= r1[1] - r4[1];
  T[1][1]= r2[1] - r4[1];
  T[1][2]= r3[1] - r4[1];
  T[2][0]= r1[2] - r4[2];
  T[2][1]= r2[2] - r4[2];
  T[2][2]= r3[2] - r4[2];

  invert3x3(T, Tinv);

  std::vector<quad> b= {r[0] - r4[0], r[1] - r4[1], r[2] - r4[2]};

  for (auto val : lambda)
    val= 0.0Q;

  for (int j= 0; j < 3; j++)
    for (int i= 0; i < 3; i++)
    {
      lambda[j] += Tinv[j][i]*b[i];
    }
  lambda[3]= 1.0Q - lambda[0] - lambda[1] - lambda[2];

}

void bary2phys3d(const std::vector<quad> &lambda,
                 const std::vector<quad> &r1,
                 const std::vector<quad> &r2,
                 const std::vector<quad> &r3,
                 const std::vector<quad> &r4,
                 std::vector<quad> &r)
{
  assert(r.size() == 3);
  assert(r1.size() == 3);
  assert(r2.size() == 3);
  assert(r3.size() == 3);
  assert(r4.size() == 3);
  assert(lambda.size() == 4);

  std::vector<std::vector<quad>> T= {{0.0Q, 0.0Q, 0.0Q}, {0.0Q, 0.0Q, 0.0Q}, {0.0Q, 0.0Q, 0.0Q}};
  // std::vector<std::vector<quad>> Tinv= {{0.0Q, 0.0Q, 0.0Q}, {0.0Q, 0.0Q, 0.0Q}, {0.0Q, 0.0Q, 0.0Q}};

  T[0][0]= r1[0] - r4[0];
  T[0][1]= r2[0] - r4[0];
  T[0][2]= r3[0] - r4[0];
  T[1][0]= r1[1] - r4[1];
  T[1][1]= r2[1] - r4[1];
  T[1][2]= r3[1] - r4[1];
  T[2][0]= r1[2] - r4[2];
  T[2][1]= r2[2] - r4[2];
  T[2][2]= r3[2] - r4[2];

  for (int j= 0; j < 3; j++)
    r[j]= lambda[0]*r1[j] + lambda[1]*r2[j] + lambda[2]*r3[j] + lambda[3]*r4[j];

}

void ref2phys3d(const std::vector<quad> &ref,
                const std::vector<quad> &r1,
                const std::vector<quad> &r2,
                const std::vector<quad> &r3,
                const std::vector<quad> &r4,
                std::vector<quad> &r)
{

  const std::vector<quad> ref0= {-1.0Q, -1.0Q, -1.0Q};
  const std::vector<quad> ref1= {+1.0Q, -1.0Q, -1.0Q};
  const std::vector<quad> ref2= {-1.0Q, +1.0Q, -1.0Q};
  const std::vector<quad> ref3= {-1.0Q, -1.0Q, +1.0Q};

  // get the barycentric coordinates on the reference element
  std::vector<quad> lambda(4, 0.0);

  phys2bary3d(ref, ref0, ref1, ref2, ref3, lambda);

  // now convert to physical element position from the barycentrics
  bary2phys3d(lambda, r1, r2, r3, r4, r);

}

void phys2ref3d(const std::vector<quad> &r,
                const std::vector<quad> &r1,
                const std::vector<quad> &r2,
                const std::vector<quad> &r3,
                const std::vector<quad> &r4,
                std::vector<quad> &ref)
{

  const std::vector<quad> ref0= {-1.0Q, -1.0Q, -1.0Q};
  const std::vector<quad> ref1= {+1.0Q, -1.0Q, -1.0Q};
  const std::vector<quad> ref2= {-1.0Q, +1.0Q, -1.0Q};
  const std::vector<quad> ref3= {-1.0Q, -1.0Q, +1.0Q};

  // get the barycentric coordinates on the physical element
  std::vector<quad> lambda(4, 0.0);

  phys2bary3d(r, r1, r2, r3, r4, lambda);

  // now convert to reference element position from the barycentrics
  bary2phys3d(lambda, ref0, ref1, ref2, ref3, ref);

}

void phys2bary4d(const std::vector<quad> &r,
                 const std::vector<quad> &r1,
                 const std::vector<quad> &r2,
                 const std::vector<quad> &r3,
                 const std::vector<quad> &r4,
                 const std::vector<quad> &r5,
                 std::vector<quad> &lambda)
{
  assert(r.size() == 4);
  assert(r1.size() == 4);
  assert(r2.size() == 4);
  assert(r3.size() == 4);
  assert(r4.size() == 4);
  assert(r5.size() == 4);
  assert(lambda.size() == 5);

  std::vector<std::vector<quad>> T= {{0.0Q, 0.0Q, 0.0Q, 0.0Q}, {0.0Q, 0.0Q, 0.0Q, 0.0Q}, {0.0Q, 0.0Q, 0.0Q, 0.0Q}, {0.0Q, 0.0Q, 0.0Q, 0.0Q}};
  std::vector<std::vector<quad>> Tinv= {{0.0Q, 0.0Q, 0.0Q, 0.0Q}, {0.0Q, 0.0Q, 0.0Q, 0.0Q}, {0.0Q, 0.0Q, 0.0Q, 0.0Q}, {0.0Q, 0.0Q, 0.0Q, 0.0Q}};

  T[0][0]= r1[0] - r5[0];
  T[0][1]= r2[0] - r5[0];
  T[0][2]= r3[0] - r5[0];
  T[0][3]= r4[0] - r5[0];
  T[1][0]= r1[1] - r5[1];
  T[1][1]= r2[1] - r5[1];
  T[1][2]= r3[1] - r5[1];
  T[1][3]= r4[1] - r5[1];
  T[2][0]= r1[2] - r5[2];
  T[2][1]= r2[2] - r5[2];
  T[2][2]= r3[2] - r5[2];
  T[2][3]= r4[2] - r5[2];
  T[3][0]= r1[3] - r5[3];
  T[3][1]= r2[3] - r5[3];
  T[3][2]= r3[3] - r5[3];
  T[3][3]= r4[3] - r5[3];

  invert4x4(T, Tinv);

  std::vector<quad> b= {r[0] - r5[0], r[1] - r5[1], r[2] - r5[2], r[3] - r5[3]};

  for (auto val : lambda)
    val= 0.0Q;

  for (int j= 0; j < 4; j++)
    for (int i= 0; i < 4; i++)
    {
      lambda[j] += Tinv[j][i]*b[i];
    }
  lambda[4]= 1.0Q - lambda[0] - lambda[1] - lambda[2] - lambda[3];

}

void bary2phys4d(const std::vector<quad> &lambda,
                 const std::vector<quad> &r1,
                 const std::vector<quad> &r2,
                 const std::vector<quad> &r3,
                 const std::vector<quad> &r4,
                 const std::vector<quad> &r5,
                 std::vector<quad> &r)
{
  assert(r.size() == 4);
  assert(r1.size() == 4);
  assert(r2.size() == 4);
  assert(r3.size() == 4);
  assert(r4.size() == 4);
  assert(r5.size() == 4);
  assert(lambda.size() == 5);

  std::vector<std::vector<quad>> T= {{0.0Q, 0.0Q, 0.0Q, 0.0Q}, {0.0Q, 0.0Q, 0.0Q, 0.0Q}, {0.0Q, 0.0Q, 0.0Q, 0.0Q}, {0.0Q, 0.0Q, 0.0Q, 0.0Q}};
  // std::vector<std::vector<quad>> Tinv= {{0.0Q, 0.0Q, 0.0Q, 0.0Q}, {0.0Q, 0.0Q, 0.0Q, 0.0Q}, {0.0Q, 0.0Q, 0.0Q, 0.0Q}, {0.0Q, 0.0Q, 0.0Q, 0.0Q}};

  T[0][0]= r1[0] - r5[0];
  T[0][1]= r2[0] - r5[0];
  T[0][2]= r3[0] - r5[0];
  T[0][3]= r4[0] - r5[0];
  T[1][0]= r1[1] - r5[1];
  T[1][1]= r2[1] - r5[1];
  T[1][2]= r3[1] - r5[1];
  T[1][3]= r4[1] - r5[1];
  T[2][0]= r1[2] - r5[2];
  T[2][1]= r2[2] - r5[2];
  T[2][2]= r3[2] - r5[2];
  T[2][3]= r4[2] - r5[2];
  T[3][0]= r1[3] - r5[3];
  T[3][1]= r2[3] - r5[3];
  T[3][2]= r3[3] - r5[3];
  T[3][3]= r4[3] - r5[3];

  for (int j= 0; j < 4; j++)
    r[j]= lambda[0]*r1[j] + lambda[1]*r2[j] + lambda[2]*r3[j] + lambda[3]*r4[j] + lambda[4]*r5[j];

}

void ref2phys4d(const std::vector<quad> &ref,
                const std::vector<quad> &r1,
                const std::vector<quad> &r2,
                const std::vector<quad> &r3,
                const std::vector<quad> &r4,
                const std::vector<quad> &r5,
                std::vector<quad> &r)
{

  const std::vector<quad> ref0= {-1.0Q, -1.0Q, -1.0Q, -1.0Q};
  const std::vector<quad> ref1= {+1.0Q, -1.0Q, -1.0Q, -1.0Q};
  const std::vector<quad> ref2= {-1.0Q, +1.0Q, -1.0Q, -1.0Q};
  const std::vector<quad> ref3= {-1.0Q, -1.0Q, +1.0Q, -1.0Q};
  const std::vector<quad> ref4= {-1.0Q, -1.0Q, -1.0Q, +1.0Q};

  // get the barycentric coordinates on the reference element
  std::vector<quad> lambda(5, 0.0);

  phys2bary4d(ref, ref0, ref1, ref2, ref3, ref4, lambda);

  // now convert to physical element position from the barycentrics
  bary2phys4d(lambda, r1, r2, r3, r4, r5, r);

}

void phys2ref4d(const std::vector<quad> &r,
                const std::vector<quad> &r1,
                const std::vector<quad> &r2,
                const std::vector<quad> &r3,
                const std::vector<quad> &r4,
                const std::vector<quad> &r5,
                std::vector<quad> &ref)
{

  const std::vector<quad> ref0= {-1.0Q, -1.0Q, -1.0Q, -1.0Q};
  const std::vector<quad> ref1= {+1.0Q, -1.0Q, -1.0Q, -1.0Q};
  const std::vector<quad> ref2= {-1.0Q, +1.0Q, -1.0Q, -1.0Q};
  const std::vector<quad> ref3= {-1.0Q, -1.0Q, +1.0Q, -1.0Q};
  const std::vector<quad> ref4= {-1.0Q, -1.0Q, -1.0Q, +1.0Q};

  // get the barycentric coordinates on the physical element
  std::vector<quad> lambda(5, 0.0);

  phys2bary4d(r, r1, r2, r3, r4, r5, lambda);

  // now convert to reference element position from the barycentrics
  bary2phys4d(lambda, ref0, ref1, ref2, ref3, ref4, ref);

}

} // namespace quadmath
