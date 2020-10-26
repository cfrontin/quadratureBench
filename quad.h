
#ifndef QUAD_H
#define QUAD_H

#include "types.h"

namespace quadquad
{

enum quadShape{pentatope, prism, tesseract};

class Quadrature
{
public:
  explicit Quadrature(const int pReq, const quadShape shape);
  Quadrature(const Quadrature&)= delete; // no copy constructor
  ~Quadrature();

  int order() const {return order_;}
  int Nquad() const {return Nquad_;}
  void coordinate(int n, quad &x, quad &y, quad &z, quad &t);
  void weight(int n, quad &w);

private:
  quadShape shape_;
  int order_;
  int orderReq_;
  int Nquad_;
  quad *x_;
  quad *y_;
  quad *z_;
  quad *t_;
  quad *w_;
};

} // namespace quadquad

#endif // QUAD_H
