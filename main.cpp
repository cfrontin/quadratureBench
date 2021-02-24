
#include "main.h"

#define DO_ELEMENT

// #define DO_BONUS
// #define DO_ODD
// #define DO_EXPONENTIAL

#ifdef DO_BONUS

quad f_test(const quad &x, const quad &y, const quad &z, const quad &t)
{
  return exp(4.0Q*pow(x, 5) + 3.0Q*pow(y, 4) + 2.0Q*pow(z, 3) + 1.0Q*pow(t, 2));
}

const quad Jtrue= 34.606338088755952564789867693078275939690544824360Q;

#elif defined(DO_ODD)

quad f_test(const quad &x, const quad &y, const quad &z, const quad &t)
{
  return sin(1.0Q*pow(x, 2) + 1.0Q*pow(y, 2) + 1.0Q*pow(z, 2) + 1.0Q*pow(t, 2));
}

const quad Jtrue= 0.8103884703464141536417680253761570893076Q;

#elif defined(DO_EXPONENTIAL)

quad f_test(const quad &x, const quad &y, const quad &z, const quad &t)
{
  return exp(1.0Q*pow(x, 2) + 2.0Q*pow(y, 3) + 3.0Q*pow(z, 4) + 4.0Q*pow(t, 5));
}

const quad Jtrue= 34.60633808875595256478986769307827593969Q;

#else

quad f_test(const quad &x, const quad &y, const quad &z, const quad &t)
{
  quad rSquared= (1.0Q*pow(x, 2) + 2.0Q*pow(y, 3) + 3.0Q*pow(z, 4) + 4.0Q*pow(t, 5));
  return sin(rSquared);
}

const quad Jtrue= 0.4104949988512827755216606231419224883815Q;

#endif

void get_expos(const int m, std::list<std::vector<int>> &expo_list)
{
  expo_list.clear();

  for (int i= 0; i <= m; i++)
  {
    for (int j= 0; j <= (m - i); j++)
    {
      for (int k= 0; k <= (m - i - j); k++)
      {
        for (int l= 0; l <= (m - i - j - k); l++)
        {
          std::vector<int> expos= {i, j, k, l};
          expo_list.push_back(expos);
        }
      }
    }
  }

}

quad factorial2quad(int x)
{
  quad fac= 1.0Q;

  for (int i= x; i > 0; i--)
  {
    fac *= static_cast<quad>(i);
  }

  return fac;
}

quad function_monomial(quad coeff, std::vector<int> expos, std::vector<quad> xVec)
{
  quad x= xVec[0];
  quad y= xVec[1];
  quad z= xVec[2];
  quad t= xVec[3];

  quad a= static_cast<quad>(expos[0]);
  quad b= static_cast<quad>(expos[1]);
  quad c= static_cast<quad>(expos[2]);
  quad d= static_cast<quad>(expos[3]);

  quad f= coeff*pow(x, a)*pow(y, b)*pow(z, c)*pow(t, d);

  return f;
}

quad true_integral_monomial_simplex(quad coeff, std::vector<int> expos)
{
  quad J= coeff;

  int a= expos[0];
  int b= expos[1];
  int c= expos[2];
  int d= expos[3];

  J *= factorial2quad(a);
  J *= factorial2quad(b);
  J *= factorial2quad(c);
  J *= factorial2quad(d);
  J /= factorial2quad(a + b + c + d + 4);

  return J;

}

quad true_integral_monomial_prism(quad coeff, std::vector<int> expos)
{
  quad J= coeff;

  int a= expos[0];
  int b= expos[1];
  int c= expos[2];
  int d= expos[3];

  J *= factorial2quad(a);
  J *= factorial2quad(b);
  J *= factorial2quad(c);
  J /= factorial2quad(a + b + c + 3);
  J /= (1.0Q + static_cast<quad>(d));

  return J;

}

quad true_integral_monomial_tesseract(quad coeff, std::vector<int> expos)
{
  quad J= coeff;

  int a= expos[0];
  int b= expos[1];
  int c= expos[2];
  int d= expos[3];

  J /= (1.0Q + static_cast<quad>(a));
  J /= (1.0Q + static_cast<quad>(b));
  J /= (1.0Q + static_cast<quad>(c));
  J /= (1.0Q + static_cast<quad>(d));

  return J;

}

#ifndef DO_ELEMENT

int main(int argc, char **argv)
{
  using namespace boost::multiprecision;
  using namespace quadquad;
  using namespace quadmath;

  std::ofstream outfile;
#ifdef DO_BONUS
  outfile.open("data_bonus.txt", std::ios::app);
#elif defined(DO_ODD)
  outfile.open("data_odd.txt", std::ios::app);
#elif defined(DO_EXPONENTIAL)
  outfile.open("data_exponential.txt", std::ios::app);
#else
  outfile.open("data_sinusoid.txt", std::ios::app);
#endif
  // outfile.open("data.txt", std::ios::app);

  // define computational box

  const int d= 4;

  int Nelem_1d;
  int pReq;
  int shapeArg;
  quadShape shape;

  // read from file
  if (argc <= 1)
  {
    Nelem_1d= 4;
    pReq= 3;
    shapeArg= 0;
  }
  else if (argc <= 2)
  {
    Nelem_1d= atoi(argv[1]);
    pReq= 3;
    shapeArg= 0;
  }
  else
  {
    Nelem_1d= atoi(argv[1]);
    pReq= atoi(argv[2]);
    shapeArg= atoi(argv[3]);
  }

  if (shapeArg == 0)
    shape= quadShape::pentatope;
  else if (shapeArg == 1)
    shape= quadShape::prism;
  else if (shapeArg == 2)
    shape= quadShape::tesseract;
  else
    throw std::runtime_error("argument not found");

  // number of gridpoints
  int Nx= Nelem_1d + 1;
  int Ny= Nelem_1d + 1;
  int Nz= Nelem_1d + 1;
  int Nt= Nelem_1d + 1;

  // limits of pentatope
  const quad x0= 0.0;
  const quad y0= 0.0;
  const quad z0= 0.0;
  const quad t0= 0.0;

  const quad x1= 1.0;
  const quad y1= 1.0;
  const quad z1= 1.0;
  const quad t1= 1.0;

  // evenly spaced divisions
  quad dx= (x1 - x0)/static_cast<quad>(Nx - 1);
  quad dy= (y1 - y0)/static_cast<quad>(Ny - 1);
  quad dz= (z1 - z0)/static_cast<quad>(Nz - 1);
  quad dt= (t1 - t0)/static_cast<quad>(Nt - 1);

  // variables for points
  quad xx0;
  quad xx1;
  quad yy0;
  quad yy1;
  quad zz0;
  quad zz1;
  quad tt0;
  quad tt1;

  std::cout << "current shape: " << shape << std::endl;

  // integral
  quad J= 0.0Q;
  int nElem= 0;
  int pAct= -1;

  // loop over box subdivisions
  for(int i= 0; i < Nx - 1; i++)
  {
    for (int j= 0; j < Ny - 1; j++)
    {
      for (int k= 0; k < Nz - 1; k++)
      {
        for (int m= 0; m < Nt - 1; m++)
        {

          // limits of the sub-box
          xx0= static_cast<quad>(i)*dx;
          yy0= static_cast<quad>(j)*dy;
          zz0= static_cast<quad>(k)*dz;
          tt0= static_cast<quad>(m)*dt;

          xx1= static_cast<quad>(i + 1)*dx;
          yy1= static_cast<quad>(j + 1)*dy;
          zz1= static_cast<quad>(k + 1)*dz;
          tt1= static_cast<quad>(m + 1)*dt;

          // std::cout << "sub-box (i, j, k, m)= (" << i << ", ";
          // std::cout << j << ", " << k << ", " << m << ")" << std::endl;
          // std::cout << std::endl;
          //
          // std::cout << "(xx0, yy0, zz0, tt0)= (" << xx0 << ", ";
          // std::cout << yy0 << ", " << zz0 << ", " << tt0 << ")" << std::endl;
          // std::cout << "(xx1, yy1, zz1, tt1)= (" << xx1 << ", ";
          // std::cout << yy1 << ", " << zz1 << ", " << tt1 << ")" << std::endl;
          // std::cout << "f= " << f_test(xx0, yy0, zz0, tt0) << std::endl;
          // std::cout << std::endl;

          MasterBox box(xx0, xx1, yy0, yy1, zz0, zz1, tt0, tt1, shape, pReq);

          box.integrate(f_test, J);
          nElem += box.nElem();
          pAct= box.quad_order();

          // std::cout << "nElem in this subbox: " << box.nElem() << std::endl;
          // std::cout << "integral after this subbox: " << J << std::endl;

        }
      }
    }
  }

  switch (shape)
  {
    case quadShape::pentatope:
      outfile << "PENTA" << "\t";
      break;
    case quadShape::prism:
      outfile << "PRISM" << "\t";
      break;
    case quadShape::tesseract:
      outfile << "TESSA" << "\t";
      break;
  }

  std::cout << "nElem: " << nElem << std::endl;
  std::cout << "pReq: " << pReq << std::endl;
  std::cout << "pAct: " << pAct << std::endl;
  outfile << nElem << "\t" << pReq << "\t" << pAct << "\t";

  std::cout << std::setprecision(std::numeric_limits<float128>::max_digits10);
  std::cout << "true integral: " << Jtrue << std::endl;
  std::cout << "final integral: " << J << std::endl;
  std::cout << "final integral error: " << J - Jtrue << std::endl;

  outfile << std::setprecision(std::numeric_limits<float128>::max_digits10);
  outfile << J << "\t";
  outfile << Jtrue << "\t";
  outfile << J - Jtrue << std::endl;

  // Quadrature quadrat(pReq, shape);
  //
  // int Nquad= quadrat.Nquad();
  //
  // std::cout << "Nquad= " << Nquad << std::endl;
  //
  // for (int iQ= 0; iQ < Nquad; iQ++)
  // {
  //   quad xq;
  //   quad yq;
  //   quad zq;
  //   quad tq;
  //   quad wq;
  //
  //   quadrat.coordinate(iQ, xq, yq, zq, tq);
  //   quadrat.weight(iQ, wq);
  //
  //   std::cout << "iQ= " << iQ << ":\t(x, y, z, t)= (" << xq << "," << yq << ","
  //       << wq << "," << tq << ")\tw= " << wq << std::endl;
  // }
  // std::cout << std::endl;
  //
  // std::cout << std::endl << std::endl;
  //
  // {
  //   quad T11= 5.0Q;
  //   quad T12= 10.0Q;
  //   quad T13= 5.0Q;
  //   quad T21= 7.0Q;
  //   quad T22= 9.0Q;
  //   quad T23= 4.0Q;
  //   quad T31= 2.0Q;
  //   quad T32= 7.0Q;
  //   quad T33= 2.0Q;
  //
  //   quad Tinv11e= -2.0Q/9.0Q;
  //   quad Tinv12e= 1.0Q/3.0Q;
  //   quad Tinv13e= -1.0Q/9.0Q;
  //   quad Tinv21e= -2.0Q/15.0Q;
  //   quad Tinv22e= 0.0Q;
  //   quad Tinv23e= 1.0Q/3.0Q;
  //   quad Tinv31e= 31.0Q/45.0Q;
  //   quad Tinv32e= -1.0Q/3.0Q;
  //   quad Tinv33e= -5.0Q/9.0Q;
  //
  //   quad Tinv11c;
  //   quad Tinv12c;
  //   quad Tinv13c;
  //   quad Tinv21c;
  //   quad Tinv22c;
  //   quad Tinv23c;
  //   quad Tinv31c;
  //   quad Tinv32c;
  //   quad Tinv33c;
  //
  //   std::vector<std::vector<quad>> Tmat= {
  //       {T11, T12, T13},
  //       {T21, T22, T23},
  //       {T31, T32, T33}};
  //
  //   std::vector<std::vector<quad>> Tinv_mat= {
  //       {Tinv11c, Tinv12c, Tinv13c},
  //       {Tinv21c, Tinv22c, Tinv23c},
  //       {Tinv31c, Tinv32c, Tinv33c}};
  //
  //   // invert3x3(T11, T12, T13,
  //   //   T21, T22, T23,
  //   //   T31, T32, T33,
  //   //   Tinv11c, Tinv12c, Tinv13c,
  //   //   Tinv21c, Tinv22c, Tinv23c,
  //   //   Tinv31c, Tinv32c, Tinv33c);
  //
  //   invert3x3(Tmat, Tinv_mat);
  //
  //   Tinv11c= Tinv_mat[0][0];
  //   Tinv12c= Tinv_mat[0][1];
  //   Tinv13c= Tinv_mat[0][2];
  //   Tinv21c= Tinv_mat[1][0];
  //   Tinv22c= Tinv_mat[1][1];
  //   Tinv23c= Tinv_mat[1][2];
  //   Tinv31c= Tinv_mat[2][0];
  //   Tinv32c= Tinv_mat[2][1];
  //   Tinv33c= Tinv_mat[2][2];
  //
  //   std::cout << "3x3 matrix inversion check:" << std::endl;
  //
  //   std::cout << "\terror in Tinv11: " << (Tinv11c - Tinv11e) << "\tTinv11e: " << Tinv11e << std::endl;
  //   std::cout << "\terror in Tinv12: " << (Tinv12c - Tinv12e) << "\tTinv12e: " << Tinv12e << std::endl;
  //   std::cout << "\terror in Tinv13: " << (Tinv13c - Tinv13e) << "\tTinv13e: " << Tinv13e << std::endl;
  //   std::cout << "\terror in Tinv21: " << (Tinv21c - Tinv21e) << "\tTinv21e: " << Tinv21e << std::endl;
  //   std::cout << "\terror in Tinv22: " << (Tinv22c - Tinv22e) << "\tTinv22e: " << Tinv22e << std::endl;
  //   std::cout << "\terror in Tinv23: " << (Tinv23c - Tinv23e) << "\tTinv23e: " << Tinv23e << std::endl;
  //   std::cout << "\terror in Tinv31: " << (Tinv31c - Tinv31e) << "\tTinv31e: " << Tinv31e << std::endl;
  //   std::cout << "\terror in Tinv32: " << (Tinv32c - Tinv32e) << "\tTinv32e: " << Tinv32e << std::endl;
  //   std::cout << "\terror in Tinv33: " << (Tinv33c - Tinv33e) << "\tTinv33e: " << Tinv33e << std::endl;
  // }
  //
  // std::cout << std::endl;
  //
  // {
  //   quad T11= 10.0Q;
  //   quad T12= 9.0Q;
  //   quad T13= 3.0Q;
  //   quad T14= 5.0Q;
  //   quad T21= 4.0Q;
  //   quad T22= 3.0Q;
  //   quad T23= 0.0Q;
  //   quad T24= 8.0Q;
  //   quad T31= 0.0Q;
  //   quad T32= 6.0Q;
  //   quad T33= 0.0Q;
  //   quad T34= 7.0Q;
  //   quad T41= 8.0Q;
  //   quad T42= 3.0Q;
  //   quad T43= 2.0Q;
  //   quad T44= 4.0Q;
  //
  //   quad Tinv11e= 9.0Q/32.0Q;
  //   quad Tinv12e= 25.0Q/64.0Q;
  //   quad Tinv13e= -13.0Q/32.0Q;
  //   quad Tinv14e= -27.0Q/64.0Q;
  //   quad Tinv21e= 7.0Q/24.0Q;
  //   quad Tinv22e= 7.0Q/48.0Q;
  //   quad Tinv23e= -1.0Q/8.0Q;
  //   quad Tinv24e= -7.0Q/16.0Q;
  //   quad Tinv31e= -17.0Q/16.0Q;
  //   quad Tinv32e= -49.0Q/32.0Q;
  //   quad Tinv33e= 21.0Q/16.0Q;
  //   quad Tinv34e= 67.0Q/32.0Q;
  //   quad Tinv41e= -1.0Q/4.0Q;
  //   quad Tinv42e= -1.0Q/8.0Q;
  //   quad Tinv43e= 1.0Q/4.0Q;
  //   quad Tinv44e= 3.0Q/8.0Q;
  //
  //   quad Tinv11c;
  //   quad Tinv12c;
  //   quad Tinv13c;
  //   quad Tinv14c;
  //   quad Tinv21c;
  //   quad Tinv22c;
  //   quad Tinv23c;
  //   quad Tinv24c;
  //   quad Tinv31c;
  //   quad Tinv32c;
  //   quad Tinv33c;
  //   quad Tinv34c;
  //   quad Tinv41c;
  //   quad Tinv42c;
  //   quad Tinv43c;
  //   quad Tinv44c;
  //
  //   std::vector<std::vector<quad>> Tmat= {
  //       {T11, T12, T13, T14},
  //       {T21, T22, T23, T24},
  //       {T31, T32, T33, T34},
  //       {T41, T42, T43, T44}};
  //
  //   std::vector<std::vector<quad>> Tinv_mat= {
  //       {Tinv11c, Tinv12c, Tinv13c, Tinv14c},
  //       {Tinv21c, Tinv22c, Tinv23c, Tinv24c},
  //       {Tinv31c, Tinv32c, Tinv33c, Tinv34c},
  //       {Tinv41c, Tinv42c, Tinv43c, Tinv44c}};
  //
  //   // invert4x4(T11, T12, T13, T14,
  //   //   T21, T22, T23, T24,
  //   //   T31, T32, T33, T34,
  //   //   T41, T42, T43, T44,
  //   //   Tinv11c, Tinv12c, Tinv13c, Tinv14c,
  //   //   Tinv21c, Tinv22c, Tinv23c, Tinv24c,
  //   //   Tinv31c, Tinv32c, Tinv33c, Tinv34c,
  //   //   Tinv41c, Tinv42c, Tinv43c, Tinv44c);
  //
  //   invert4x4(Tmat, Tinv_mat);
  //
  //   Tinv11c= Tinv_mat[0][0];
  //   Tinv12c= Tinv_mat[0][1];
  //   Tinv13c= Tinv_mat[0][2];
  //   Tinv14c= Tinv_mat[0][3];
  //   Tinv21c= Tinv_mat[1][0];
  //   Tinv22c= Tinv_mat[1][1];
  //   Tinv23c= Tinv_mat[1][2];
  //   Tinv24c= Tinv_mat[1][3];
  //   Tinv31c= Tinv_mat[2][0];
  //   Tinv32c= Tinv_mat[2][1];
  //   Tinv33c= Tinv_mat[2][2];
  //   Tinv34c= Tinv_mat[2][3];
  //   Tinv41c= Tinv_mat[3][0];
  //   Tinv42c= Tinv_mat[3][1];
  //   Tinv43c= Tinv_mat[3][2];
  //   Tinv44c= Tinv_mat[3][3];
  //
  //   std::cout << "4x4 matrix inversion check:" << std::endl;
  //
  //   std::cout << "\terror in Tinv11: " << Tinv11c - Tinv11e << "\tTinv11e: " << Tinv11e << std::endl;
  //   std::cout << "\terror in Tinv12: " << Tinv12c - Tinv12e << "\tTinv12e: " << Tinv12e << std::endl;
  //   std::cout << "\terror in Tinv13: " << Tinv13c - Tinv13e << "\tTinv13e: " << Tinv13e << std::endl;
  //   std::cout << "\terror in Tinv14: " << Tinv14c - Tinv14e << "\tTinv14e: " << Tinv14e << std::endl;
  //   std::cout << "\terror in Tinv21: " << Tinv21c - Tinv21e << "\tTinv21e: " << Tinv21e << std::endl;
  //   std::cout << "\terror in Tinv22: " << Tinv22c - Tinv22e << "\tTinv22e: " << Tinv22e << std::endl;
  //   std::cout << "\terror in Tinv23: " << Tinv23c - Tinv23e << "\tTinv23e: " << Tinv23e << std::endl;
  //   std::cout << "\terror in Tinv24: " << Tinv24c - Tinv24e << "\tTinv24e: " << Tinv24e << std::endl;
  //   std::cout << "\terror in Tinv31: " << Tinv31c - Tinv31e << "\tTinv31e: " << Tinv31e << std::endl;
  //   std::cout << "\terror in Tinv32: " << Tinv32c - Tinv32e << "\tTinv32e: " << Tinv32e << std::endl;
  //   std::cout << "\terror in Tinv33: " << Tinv33c - Tinv33e << "\tTinv33e: " << Tinv33e << std::endl;
  //   std::cout << "\terror in Tinv34: " << Tinv34c - Tinv34e << "\tTinv34e: " << Tinv34e << std::endl;
  //   std::cout << "\terror in Tinv41: " << Tinv41c - Tinv41e << "\tTinv41e: " << Tinv41e << std::endl;
  //   std::cout << "\terror in Tinv42: " << Tinv42c - Tinv42e << "\tTinv42e: " << Tinv42e << std::endl;
  //   std::cout << "\terror in Tinv43: " << Tinv43c - Tinv43e << "\tTinv43e: " << Tinv43e << std::endl;
  //   std::cout << "\terror in Tinv44: " << Tinv44c - Tinv44e << "\tTinv44e: " << Tinv44e << std::endl;
  // }
  //
  // std::cout << std::endl;
  //
  // std::vector<quad> r= {1.0Q/4.0Q, 1.0Q/4.0Q, 1.0Q/4.0Q};
  // std::vector<quad> r1= {0.0Q, 0.0Q, 0.0Q};
  // std::vector<quad> r2= {1.0Q, 0.0Q, 0.0Q};
  // std::vector<quad> r3= {0.0Q, 1.0Q, 0.0Q};
  // std::vector<quad> r4= {0.0Q, 0.0Q, 1.0Q};
  //
  // std::vector<quad> lambda(4, 0.0Q);
  // // lambda[0]= 0.0Q;
  // // lambda[1]= 0.0Q;
  // // lambda[2]= 0.0Q;
  // // lambda[3]= 0.0Q;
  //
  // phys2bary3d(r, r1, r2, r3, r4, lambda);
  // // bary2phys3d(lambda, r1, r2, r3, r4, r);
  //
  // std::cout << "lambda[0]= " << lambda[0] << std::endl;
  // std::cout << "lambda[1]= " << lambda[1] << std::endl;
  // std::cout << "lambda[2]= " << lambda[2] << std::endl;
  // std::cout << "lambda[3]= " << lambda[3] << std::endl << std::endl;
  //
  // std::cout << "r= (" << r[0] << ", " << r[1] << ", " << r[2] << ")" << std::endl;
  //
  // std::cout << std::endl;

  outfile.close();

}

#else // DO_ELEMENT

int main(int argc, char **argv)
{
  using namespace boost::multiprecision;
  using namespace quadquad;
  using namespace quadmath;

  std::ofstream outfile;
  outfile.open("data_monomial.txt", std::ios::app);

  for (int m= 4; m <= 20; m++)
  {

    std::list<std::vector<int>> expo_list= {};
    get_expos(m, expo_list);

    quad coeff= 1.0Q;

    for (std::vector<int> expo_set : expo_list)
    {
      quad a= expo_set[0];
      quad b= expo_set[1];
      quad c= expo_set[2];
      quad d= expo_set[3];

      std::cout << "(a, b, c, d)= (" << a << ", " << b << ", " << c << ", " << d << ")" << std::endl;

      quad J_true_simplex= true_integral_monomial_simplex(coeff, expo_set);
      quad J_true_prism= true_integral_monomial_prism(coeff, expo_set);
      quad J_true_tesseract= true_integral_monomial_tesseract(coeff, expo_set);

      std::cout << "J_true_simplex: " << J_true_simplex << std::endl;
      std::cout << "J_true_prism: " << J_true_prism << std::endl;
      std::cout << "J_true_tesseract: " << J_true_tesseract << std::endl;

      for (int qOrder= 6; qOrder <= 12; qOrder++)
      {

        std::vector<quad> node0= {0.0Q, 0.0Q, 0.0Q, 0.0Q};
        std::vector<quad> node1= {1.0Q, 0.0Q, 0.0Q, 0.0Q};
        std::vector<quad> node2= {0.0Q, 1.0Q, 0.0Q, 0.0Q};
        std::vector<quad> node3= {0.0Q, 0.0Q, 1.0Q, 0.0Q};
        std::vector<quad> node4= {0.0Q, 0.0Q, 0.0Q, 1.0Q};

        std::vector<quad> node3d0= {0.0Q, 0.0Q, 0.0Q};
        std::vector<quad> node3d1= {1.0Q, 0.0Q, 0.0Q};
        std::vector<quad> node3d2= {0.0Q, 1.0Q, 0.0Q};
        std::vector<quad> node3d3= {0.0Q, 0.0Q, 1.0Q};

        quadShape shape;

        quad J_q_simplex= 0.0Q;
        quad J_q_prism= 0.0Q;
        quad J_q_tesseract= 0.0Q;

        shape= quadShape::pentatope;
        Quadrature quadrature_pentatope(qOrder, shape);

        for (int iQ= 0; iQ < quadrature_pentatope.Nquad(); iQ++)
        {
          quad xi0;
          quad xi1;
          quad xi2;
          quad xi3;
          quad w;

          // get coordinate and weight
          quadrature_pentatope.coordinate(iQ, xi0, xi1, xi2, xi3);
          quadrature_pentatope.weight(iQ, w);

          // pack it up in a vector
          std::vector<quad> xi= {xi0, xi1, xi2, xi3};
          // vector for holding the physical coords
          std::vector<quad> xVec(4, 0.0Q);

          // translate to physical coordinates
          ref2phys4d(xi, node0, node1, node2, node3, node4, xVec);

          // adjust weight assuming even pentatopal subdivisions of
          // hyperrectangular box
          w *= 1.0Q/(2.0Q*2.0Q*2.0Q*2.0Q);

          quad f= function_monomial(1.0Q, expo_set, xVec);

          J_q_simplex += w*f;

        }

        shape= quadShape::prism;
        Quadrature quadrature_prism(qOrder, shape);

        for (int iQ= 0; iQ < quadrature_prism.Nquad(); iQ++)
        {
          quad xi0;
          quad xi1;
          quad xi2;
          quad xi3;
          quad w;

          // get coordinate and weight
          quadrature_prism.coordinate(iQ, xi0, xi1, xi2, xi3);
          quadrature_prism.weight(iQ, w);

          // pack it up in a vector
          std::vector<quad> xi3D= {xi0, xi1, xi2};
          // vector for holding the physical coords
          std::vector<quad> xVec3D(3, 0.0Q);
          std::vector<quad> xVec(4, 0.0Q);

          // translate to physical coordinates
          ref2phys3d(xi3D, node3d0, node3d1, node3d2, node3d3, xVec3D);

          // adjust weight assuming even pentatopal subdivisions of
          // hyperrectangular box
          w *= 1.0Q/(2.0Q*2.0Q*2.0Q*2.0Q);

          xVec[0]= xVec3D[0];
          xVec[1]= xVec3D[1];
          xVec[2]= xVec3D[2];
          xVec[3]= 0.5Q*(xi3 + 1.0Q);

          quad f= function_monomial(1.0Q, expo_set, xVec);

          J_q_prism += w*f;

        }

        std::cout << "\tqOrder= " << qOrder << ";\tJq= " << J_q_prism << std::endl;

        shape= quadShape::tesseract;
        Quadrature quadrature_tesseract(qOrder, shape);

        for (int iQ= 0; iQ < quadrature_tesseract.Nquad(); iQ++)
        {
          quad xi0;
          quad xi1;
          quad xi2;
          quad xi3;
          quad w;

          // get coordinate and weight
          quadrature_tesseract.coordinate(iQ, xi0, xi1, xi2, xi3);
          quadrature_tesseract.weight(iQ, w);

          // pack it up in a vector for holding the physical coords
          std::vector<quad> xVec= {
            0.5Q*(xi0 + 1.0Q),
            0.5Q*(xi1 + 1.0Q),
            0.5Q*(xi2 + 1.0Q),
            0.5Q*(xi3 + 1.0Q)
          };

          // adjust weight assuming even pentatopal subdivisions of
          // hyperrectangular box
          w *= 1.0Q/(2.0Q*2.0Q*2.0Q*2.0Q);

          quad f= function_monomial(1.0Q, expo_set, xVec);

          J_q_tesseract += w*f;

        }

        std::cout << "\tqOrder= " << qOrder << ";\tJq= " << J_q_tesseract << std::endl;

        outfile << std::setprecision(std::numeric_limits<float128>::max_digits10);
        outfile << a << "\t" << b << "\t" << c << "\t" << d << "\t";
        outfile << qOrder << "\t";
        outfile << J_true_simplex << "\t" << J_q_simplex << "\t" << (J_q_simplex - J_true_simplex) << "\t";
        outfile << J_true_prism << "\t" << J_q_prism << "\t" << (J_q_prism - J_true_prism) << "\t";
        outfile << J_true_tesseract << "\t" << J_q_tesseract << "\t" << (J_q_tesseract - J_true_tesseract) << std::endl;

      }

    }

  }

  outfile.close();

}

#endif // DO_ELEMENT
