
#include "master_box.h"

namespace quadquad {

using namespace boost::multiprecision;
using namespace quadmath;

// DEGREES OF FREEDOM OF A 3-BOX
const quad hex_dof[8][3]= {
  {0.0Q,  0.0Q,  0.0Q},
  {1.0Q,  0.0Q,  0.0Q},
  {0.0Q,  1.0Q,  0.0Q},
  {1.0Q,  1.0Q,  0.0Q},
  {0.0Q,  0.0Q,  1.0Q},
  {1.0Q,  0.0Q,  1.0Q},
  {0.0Q,  1.0Q,  1.0Q},
  {1.0Q,  1.0Q,  1.0Q}
};

// MAP [i][j] TO ith TET OF A HEX DECOMP, EACH DEFINED BY 4 POINTS
const int tet2hex[6][4]= {
  {0, 1, 2, 4},
  {1, 4, 3, 2},
  {6, 2, 3, 4},

  {1, 4, 5, 3},
  {4, 6, 5, 3},
  {7, 6, 3, 5}
};

// DEGREES OF FREEDOM OF THE TESSERACT
const quad tesseract_dof[16][4]= {
  {0.0Q, 0.0Q, 0.0Q, 0.0Q},
  {0.0Q, 0.0Q, 0.0Q, 1.0Q},
  {0.0Q, 0.0Q, 1.0Q, 0.0Q},
  {0.0Q, 0.0Q, 1.0Q, 1.0Q},
  {0.0Q, 1.0Q, 0.0Q, 0.0Q},
  {0.0Q, 1.0Q, 0.0Q, 1.0Q},
  {0.0Q, 1.0Q, 1.0Q, 0.0Q},
  {0.0Q, 1.0Q, 1.0Q, 1.0Q},
  {1.0Q, 0.0Q, 0.0Q, 0.0Q},
  {1.0Q, 0.0Q, 0.0Q, 1.0Q},
  {1.0Q, 0.0Q, 1.0Q, 0.0Q},
  {1.0Q, 0.0Q, 1.0Q, 1.0Q},
  {1.0Q, 1.0Q, 0.0Q, 0.0Q},
  {1.0Q, 1.0Q, 0.0Q, 1.0Q},
  {1.0Q, 1.0Q, 1.0Q, 0.0Q},
  {1.0Q, 1.0Q, 1.0Q, 1.0Q}
};

// MAP [i][j] TO ith PENTATOPE OF A KFL TESSERACT DECOMP, EACH DEFINED BY 5 PTS
const int pentatope2tesseract[24][4 + 1]= {
  {0, 8, 10, 15, 14},
  {0, 8, 10, 11, 15},
  {0, 8, 9, 13, 15},
  {0, 8, 9, 15, 11},
  {0, 2, 10, 14, 15},
  {0, 2, 10, 15, 11},
  {0, 8, 12, 14, 15},
  {0, 8, 12, 15, 13},
  {0, 4, 12, 15, 14},
  {0, 4, 12, 13, 15},
  {0, 4, 6, 14, 15},
  {0, 4, 6, 15, 7},
  {0, 4, 5, 15, 13},
  {0, 4, 5, 7, 15},
  {0, 2, 6, 15, 14},
  {0, 2, 6, 7, 15},
  {0, 2, 3, 11, 15},
  {0, 2, 3, 15, 7},
  {0, 1, 9, 15, 13},
  {0, 1, 9, 11, 15},
  {0, 1, 5, 13, 15},
  {0, 1, 5, 15, 7},
  {0, 1, 3, 15, 11},
  {0, 1, 3, 7, 15}
};

MasterBox::MasterBox(const quad x0, const quad x1, const quad y0, const quad y1,
                     const quad z0, const quad z1, const quad t0, const quad t1,
                     const quadShape shape, const int pReq) :
                         x0_(x0), y0_(y0), z0_(z0), t0_(t0),
                         x1_(x1), y1_(y1), z1_(z1), t1_(t1),
                         shape_(shape), pReq_(pReq), quad_(pReq, shape)
{
  // std::cout << "in masterbox shape is " << shape_ << std::endl;
}

int MasterBox::quad_order()
{
  return quad_.order();
}

int MasterBox::nElem()
{
  switch (shape_)
  {
    case quadShape::pentatope:
      return nElem_pentatope;
    case quadShape::prism:
      return nElem_prism;
    case quadShape::tesseract:
      return nElem_tesseract;
  }
}

void MasterBox::integrate(quad (*func)(const quad &x, const quad &y,
    const quad &z, const quad &t), quad &J)
{
  switch (shape_)
  {
    case quadShape::tesseract:
      {
        for (int iq= 0; iq < quad_.Nquad(); iq++)
        {
          // get quadrature point and weight
          quad xi0;
          quad xi1;
          quad xi2;
          quad xi3;
          quad w;

          quad_.coordinate(iq, xi0, xi1, xi2, xi3);
          quad_.weight(iq, w);

          // transform, assume rectangular shape
          // from (-1, 1) to (x0, x1), etc...
          quad x= 0.5Q*(x0_ + x1_) + 0.5Q*xi0*(x1_ - x0_);
          quad y= 0.5Q*(y0_ + y1_) + 0.5Q*xi1*(y1_ - y0_);
          quad z= 0.5Q*(z0_ + z1_) + 0.5Q*xi2*(z1_ - z0_);
          quad t= 0.5Q*(t0_ + t1_) + 0.5Q*xi3*(t1_ - t0_);

          // std::cout << "(x0_, y0_, z0_, t0_)= (" << x0_ << ", " << x0_ << ", "
          //     << z0_ << ", " << t0_ << ")" << std::endl;
          //
          // std::cout << "(x1_, y1_, z1_, t1_)= (" << x1_ << ", " << x1_ << ", "
          //     << z1_ << ", " << t1_ << ")" << std::endl;
          //
          // std::cout << "xi= (" << xi0 << ", " << xi1 << ", " << xi2 << ", "
          //     << xi3 << ")" << std::endl;
          //
          // std::cout << "(x, y, z, t)= (" << x << ", " << y << ", " << z
          //     << ", " << t << ")" << std::endl;

          // correct weight
          w *= L_x()*L_y()*L_z()*L_t()/(2.0Q*2.0Q*2.0Q*2.0Q);

          // grab function value at point
          quad f= (*func)(x, y, z, t);

          // add to the integral
          J += w*f;

        }

        break;
      }
    case quadShape::prism:
      {
        for (int idx_elem= 0; idx_elem < nElem(); idx_elem++)
        {
          std::vector<quad> node0(3);
          std::vector<quad> node1(3);
          std::vector<quad> node2(3);
          std::vector<quad> node3(3);

          int d;
          d= 0; node0[d]= x0_ + (x1_ - x0_)*hex_dof[tet2hex[idx_elem][0]][d];
          d= 1; node0[d]= y0_ + (y1_ - y0_)*hex_dof[tet2hex[idx_elem][0]][d];
          d= 2; node0[d]= z0_ + (z1_ - z0_)*hex_dof[tet2hex[idx_elem][0]][d];

          d= 0; node1[d]= x0_ + (x1_ - x0_)*hex_dof[tet2hex[idx_elem][1]][d];
          d= 1; node1[d]= y0_ + (y1_ - y0_)*hex_dof[tet2hex[idx_elem][1]][d];
          d= 2; node1[d]= z0_ + (z1_ - z0_)*hex_dof[tet2hex[idx_elem][1]][d];

          d= 0; node2[d]= x0_ + (x1_ - x0_)*hex_dof[tet2hex[idx_elem][2]][d];
          d= 1; node2[d]= y0_ + (y1_ - y0_)*hex_dof[tet2hex[idx_elem][2]][d];
          d= 2; node2[d]= z0_ + (z1_ - z0_)*hex_dof[tet2hex[idx_elem][2]][d];

          d= 0; node3[d]= x0_ + (x1_ - x0_)*hex_dof[tet2hex[idx_elem][3]][d];
          d= 1; node3[d]= y0_ + (y1_ - y0_)*hex_dof[tet2hex[idx_elem][3]][d];
          d= 2; node3[d]= z0_ + (z1_ - z0_)*hex_dof[tet2hex[idx_elem][3]][d];

          for (int iq= 0; iq < quad_.Nquad(); iq++)
          {
            // get quadrature point and weight
            quad xi0;
            quad xi1;
            quad xi2;
            quad xi3;
            quad w;

            quad_.coordinate(iq, xi0, xi1, xi2, xi3);
            quad_.weight(iq, w);

            // get it into the physical frame
            std::vector<quad> xi= {xi0, xi1, xi2};
            std::vector<quad> x3Vec(3, 0.0Q);

            // translate to physical coordinates
            ref2phys3d(xi, node0, node1, node2, node3, x3Vec);
            quad t= 0.5Q*(t0_ + t1_) + 0.5Q*xi3*(t1_ - t0_);

            // adjust weight assuming even pentatopal subdivisions of
            // hyperrectangular box
            w *= (L_x()*L_y()*L_z()*L_t())/(2.0Q*2.0Q*2.0Q*2.0Q);

            // grab function value at point
            quad f= func(x3Vec[0], x3Vec[1], x3Vec[2], t);

            // add to the integral
            J += w*f;

          }
        }
      }
      break;
    case quadShape::pentatope:
      {
        for (int idx_elem= 0; idx_elem < nElem(); idx_elem++)
        {
          std::vector<quad> node0(4);
          std::vector<quad> node1(4);
          std::vector<quad> node2(4);
          std::vector<quad> node3(4);
          std::vector<quad> node4(4);

          int d;
          d= 0; node0[d]= x0_ + (x1_ - x0_)*tesseract_dof[pentatope2tesseract[idx_elem][0]][d];
          d= 1; node0[d]= y0_ + (y1_ - y0_)*tesseract_dof[pentatope2tesseract[idx_elem][0]][d];
          d= 2; node0[d]= z0_ + (z1_ - z0_)*tesseract_dof[pentatope2tesseract[idx_elem][0]][d];
          d= 3; node0[d]= t0_ + (t1_ - t0_)*tesseract_dof[pentatope2tesseract[idx_elem][0]][d];

          d= 0; node1[d]= x0_ + (x1_ - x0_)*tesseract_dof[pentatope2tesseract[idx_elem][1]][d];
          d= 1; node1[d]= y0_ + (y1_ - y0_)*tesseract_dof[pentatope2tesseract[idx_elem][1]][d];
          d= 2; node1[d]= z0_ + (z1_ - z0_)*tesseract_dof[pentatope2tesseract[idx_elem][1]][d];
          d= 3; node1[d]= t0_ + (t1_ - t0_)*tesseract_dof[pentatope2tesseract[idx_elem][1]][d];

          d= 0; node2[d]= x0_ + (x1_ - x0_)*tesseract_dof[pentatope2tesseract[idx_elem][2]][d];
          d= 1; node2[d]= y0_ + (y1_ - y0_)*tesseract_dof[pentatope2tesseract[idx_elem][2]][d];
          d= 2; node2[d]= z0_ + (z1_ - z0_)*tesseract_dof[pentatope2tesseract[idx_elem][2]][d];
          d= 3; node2[d]= t0_ + (t1_ - t0_)*tesseract_dof[pentatope2tesseract[idx_elem][2]][d];

          d= 0; node3[d]= x0_ + (x1_ - x0_)*tesseract_dof[pentatope2tesseract[idx_elem][3]][d];
          d= 1; node3[d]= y0_ + (y1_ - y0_)*tesseract_dof[pentatope2tesseract[idx_elem][3]][d];
          d= 2; node3[d]= z0_ + (z1_ - z0_)*tesseract_dof[pentatope2tesseract[idx_elem][3]][d];
          d= 3; node3[d]= t0_ + (t1_ - t0_)*tesseract_dof[pentatope2tesseract[idx_elem][3]][d];

          d= 0; node4[d]= x0_ + (x1_ - x0_)*tesseract_dof[pentatope2tesseract[idx_elem][4]][d];
          d= 1; node4[d]= y0_ + (y1_ - y0_)*tesseract_dof[pentatope2tesseract[idx_elem][4]][d];
          d= 2; node4[d]= z0_ + (z1_ - z0_)*tesseract_dof[pentatope2tesseract[idx_elem][4]][d];
          d= 3; node4[d]= t0_ + (t1_ - t0_)*tesseract_dof[pentatope2tesseract[idx_elem][4]][d];

          for (int iq= 0; iq < quad_.Nquad(); iq++)
          {
            // get quadrature point and weight
            quad xi0;
            quad xi1;
            quad xi2;
            quad xi3;
            quad w;

            quad_.coordinate(iq, xi0, xi1, xi2, xi3);
            quad_.weight(iq, w);

            // get it into the physical frame
            std::vector<quad> xi= {xi0, xi1, xi2, xi3};
            std::vector<quad> xVec(4, 0.0Q);

            // translate to physical coordinates
            ref2phys4d(xi, node0, node1, node2, node3, node4, xVec);

            // adjust weight assuming even pentatopal subdivisions of
            // hyperrectangular box
            w *= (L_x()*L_y()*L_z()*L_t())/(2.0Q*2.0Q*2.0Q*2.0Q);

            // grab function value at point
            quad f= func(xVec[0], xVec[1], xVec[2], xVec[3]);

            // add to the integral
            J += w*f;

          }
        }
      }
      break;
    default:
      throw std::runtime_error("should not have gotten here.");
  }
}

}

//
