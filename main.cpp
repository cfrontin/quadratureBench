
#include "main.h"

quad f_test(const quad &x, const quad &y, const quad &z, const quad &t)
{
  return exp(1.0Q*pow(x, 2) + 2.0Q*pow(y, 3) + 3.0Q*pow(y, 4) + 4.0Q*pow(t, 5));
}

quad Jtrue= 65.94533037500062508665035943680133764992;

int main(int argc, char **argv)
{
  using namespace boost::multiprecision;
  using namespace quadquad;
  using namespace quadmath;

  std::ofstream outfile;
  outfile.open("data.txt", std::ios::app);

  // define computational box

  const qint d= 4;

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
  qint Nx= static_cast<qint>(Nelem_1d + 1);
  qint Ny= static_cast<qint>(Nelem_1d + 1);
  qint Nz= static_cast<qint>(Nelem_1d + 1);
  qint Nt= static_cast<qint>(Nelem_1d + 1);

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
  for(qint i= 0; i < Nx - 1; i++)
  {
    for (qint j= 0; j < Ny - 1; j++)
    {
      for (qint k= 0; k < Nz - 1; k++)
      {
        for (qint m= 0; m < Nt - 1; m++)
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
