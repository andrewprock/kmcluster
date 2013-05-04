#ifndef CLUSTER_TRIANGLE_H_
#define CLUSTER_TRIANGLE_H_

#include <string>
#include <boost/format.hpp>
#include <boost/tuple/tuple.hpp>

namespace rps
{
  // use a generic boost interface for now
  typedef boost::tuple<double,double> bpoint2_t;
  typedef boost::tuple<double,double,double> bpoint3_t;

  class Triangle
  {

  public:
    Triangle (bpoint2_t p1,
              bpoint2_t p2,
              bpoint2_t p3)
      : _p1(p1)
      , _p2(p2)
      , _p3(p3)
      , _det ((X(_p1)-X(_p3)) * (Y(_p2)-Y(_p3)) -
              (X(_p2)-X(_p3)) * (Y(_p1)-Y(_p3)))
    { }

    // see: 
    // http://en.wikipedia.org/wiki/Barycentric_coordinates_(mathematics)
    bpoint3_t getBarycentricCoordinates (const bpoint2_t& P) const
    {
#if 0
      double lambda1 = ( (Y(_p2)-Y(_p3))*(X(P)-X(_p3)) - (X(_p2)-X(_p3))*(Y(P)-Y(_p3)))/_det;
      double lambda2 = (-(Y(_p1)-Y(_p3))*(X(P)-X(_p3)) + (X(_p1)-X(_p3))*(Y(P)-Y(_p3)))/_det;
      double lambda3 = 1 - lambda1 - lambda2;
#else
      double lambda1 = ( (Y(_p2)-Y(_p3))*(X(P)-X(_p3)) - (X(_p2)-X(_p3))*(Y(P)-Y(_p3)));
      double lambda2 = (-(Y(_p1)-Y(_p3))*(X(P)-X(_p3)) + (X(_p1)-X(_p3))*(Y(P)-Y(_p3)));
      double lambda3 = 1 - (lambda1 + lambda2)/_det;
      // the division is delayed to make the computation of lambda3 more stable
      lambda1 /= _det;
      lambda2 /= _det;

#endif
      return bpoint3_t (lambda1, lambda2, lambda3);
    }

    bool contains (const bpoint2_t& point) const
    {
      bpoint3_t bary = getBarycentricCoordinates (point);
      return (bary.get<0>() >= 0 &&
              bary.get<1>() >= 0 &&
              bary.get<2>() >= 0);
    }

    bool near_contains (const bpoint2_t& point, double epsilon) const
    {
      bpoint3_t bary = getBarycentricCoordinates (point);
      return (bary.get<0>() >= -epsilon &&
              bary.get<1>() >= -epsilon &&
              bary.get<2>() >= -epsilon);
    }

    std::string str() const
    {
      return (boost::format ("(%f,%f) (%f,%f) (%f,%f)")
              % _p1.get<0>()
              % _p1.get<1>()
              % _p2.get<0>()
              % _p2.get<1>()
              % _p3.get<0>()
              % _p3.get<1>()
              ).str();
      
    }

    bpoint2_t getCenter () const
    {
      return bpoint2_t (
                        (X(_p1)+X(_p2)+X(_p3))/3,
                        (Y(_p1)+Y(_p2)+Y(_p3))/3
                        );
    }
    
    // the operator here works by comparing the smallest x value the
    // triangle centers
    bool operator<(const Triangle& other) const
    {
      return X(getCenter()) < X(other.getCenter());
    }

  private:
    // vertices
    bpoint2_t _p1;
    bpoint2_t _p2;
    bpoint2_t _p3;
    double    _det;

    // to make indexing more mathmatically convenient
    double X(const bpoint2_t& p) const { return p.get<0>(); }
    double Y(const bpoint2_t& p) const { return p.get<1>(); }

  };

}

#endif  // CLUSTER_TRIANGLE_H_
