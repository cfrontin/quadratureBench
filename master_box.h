
#ifndef MASTER_BOX_H
#define MASTER_BOX_H

#include "types.h"
#include "quad.h"
#include "math.h"

namespace quadquad
{

using namespace boost::multiprecision;

class MasterBox
{
public:
  explicit MasterBox(const quad x0, const quad x1, const quad y0, const quad y1,
      const quad z0, const quad z1, const quad t0, const quad t1,
      const quadShape shape, const int pReq);
  MasterBox(const MasterBox&)= delete; // no copy constructor

  int nElem();
  int quad_order();

  void integrate(quad (*func)(const quad &x, const quad &y, const quad &z,
      const quad &t), quad &J);

  quad L_x() { return x1_ - x0_; };
  quad L_y() { return y1_ - y0_; };
  quad L_z() { return z1_ - z0_; };
  quad L_t() { return t1_ - t0_; };

private:

  int nElem_pentatope= 24;
  int nElem_prism= 6;
  int nElem_tesseract= 1;

  // quadrature
  Quadrature quad_;

  // order requested
  int pReq_;

  // type of quadrature on box
  quadShape shape_;

  // limits
  quad x0_;
  quad x1_;
  quad y0_;
  quad y1_;
  quad z0_;
  quad z1_;
  quad t0_;
  quad t1_;

};









} // namespace quadquad

#endif // MASTER_BOX_H
