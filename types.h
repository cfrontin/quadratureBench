
#ifndef TYPES_H
#define TYPES_H

#include <boost/multiprecision/float128.hpp>
#include <boost/multiprecision/cpp_int.hpp>

// typedefs for quad-precision numbers
typedef boost::multiprecision::float128 quad;
typedef boost::multiprecision::int128_t qint;

// define Q pi
static const quad pi= 3.1415926535897932384626433832795028841971693993751058Q;

#endif // !defined(TYPES_H)
