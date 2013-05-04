#ifndef CLUSTER_TRIANGULATION_H_
#define CLUSTER_TRIANGULATION_H_

#include <utility>
#include <boost/format.hpp>
#include "Triangle.h"

namespace kmcluster
{
  typedef boost::tuple<int,int,int> triad_t;

  const int TRIAD_SIZE = 3;

  class Triangulation
  {
  public:
    /**
     *
     */
    Triangulation (const std::vector<Triangle>& tlist)
      : _faces(tlist)
      , _searchList(tlist.size())
      , _containsOrigin(true)
    {
      for (size_t i=0; i<_faces.size(); i++)
        _searchList[i] = make_pair(_faces[i].getCenter().get<0>(), i);
      std::sort (_searchList.begin(), _searchList.end());
    }

    /**
     * return the number of triangles 
     */
    size_t size() const
    {
      return _faces.size ();
    }

    /**
     * Trinangles can be retrieved by their id
     */
    Triangle getTriangle (int id) const
    {
      return _faces[id];
    }
    
    /**
     * Find the mesh points and coordinates.
     *
     * The return value is a pair.
     * The triad_t is the list of vertices of the face
     * The bpoint3_t is the barrycdntric coordinates in terms of the
     * vertices
     *
     * Each face has it's own set of vertices, even if points that
     * coorespond to those vertices are shared by another face.
     * In a mesh with two faces that share a side, as shown in the
     * primitive asci diagram:
     *
     * .-----.-----.
     *  \1  2|4  5/
     *   \   |   /
     *    \  |  /
     *     \3|6/
     *      \|/
     *       .
     *
     * There are six vertex labels even though there are only four points
     * in the mesh.
     */
    std::pair<triad_t,bpoint3_t> getBarycentricCoordinates (const bpoint2_t& point) const
    {
#if 0
      cout << "--- start debug output ---\n";
      cout << boost::format ("pt:  %.5f %.5f\n") % point.get<0>() % point.get<1>();
      pair<triad_t,bpoint3_t> fpt = getBarycentricCoordinatesFAST (point);
      pair<triad_t,bpoint3_t> spt = getBarycentricCoordinatesSLOW (point);
      if (fpt != spt && point.get<0>() != 1.0 && point.get<1>() != 1.0)
        cout << "fail\n";
#endif
      return getBarycentricCoordinatesFAST (point);
    }

  private:
    std::pair<triad_t,bpoint3_t> getBarycentricCoordinatesSLOW (const bpoint2_t& point) const
    {
      // TODO: this is slow, make fast
      for (size_t i=0; i<_faces.size(); i++)
        {
          if (_faces[i].contains (point))
            {
//               cout << "sfound: " << i << " " << _faces[i].str() << endl;
              return std::make_pair (triad_t(TRIAD_SIZE*i,TRIAD_SIZE*i+1,TRIAD_SIZE*i+2),
                                     _faces[i].getBarycentricCoordinates (point));
            }
        }

      if (point.get<0>() == 0 && point.get<1>() == 0)
        return std::make_pair (triad_t (0,0,0), bpoint3_t(0,0,0));

      cout << "not found\n";
      for (size_t i=0; i<_faces.size(); i++)
        if (_faces[i].near_contains (point, .000001))
            return std::make_pair (triad_t(TRIAD_SIZE*i,TRIAD_SIZE*i+1,TRIAD_SIZE*i+2),
                                   _faces[i].getBarycentricCoordinates (point));

      cout << "REALLY not found\n";

#if 0
      cout << "--\n";
      cout << (boost::format ("error, could not find face for: %f %f\n") % point.get<0>() % point.get<1>());
      for (size_t i=0; i<_faces.size(); i++)
        {
          bpoint3_t bc = _faces[i].getBarycentricCoordinates (point);
          cout << _faces[i].str() << " " << bc.get<0>() << " " << bc.get<1>() << " " << bc.get<2>() << endl;
        }
      cout << "--\n";
#endif
      throw (std::runtime_error ((boost::format ("error, could not find face for point: %f %f") % point.get<0>() % point.get<1>()).str()));
      return std::make_pair (triad_t (0,0,0), bpoint3_t(0,0,0));
    }
    
    std::pair<triad_t,bpoint3_t> getBarycentricCoordinatesFAST (const bpoint2_t& point) const
    {
      using namespace std;


      if (point.get<0>() == 0 && point.get<1>() == 0 &&
          _containsOrigin == false)
        return std::make_pair (triad_t (0,0,0), bpoint3_t(0,0,0));

      vector<pair<double,int> >::const_iterator forward = lower_bound (_searchList.begin(), _searchList.end(),
                                                                       make_pair(point.get<0>(),0));
      vector<pair<double,int> >::const_iterator backward = forward;

      //cout << boost::format ("fit: %.3f %d\n") % forward->first % forward->second;
      //cout << boost::format ("bit: %.3f %d\n") % backward->first % backward->second;


      while (forward  != _searchList.end() ||
             backward != _searchList.begin())
        {
          //cout << ".";
          if (forward != _searchList.end())
            {
              if (_faces[forward->second].contains(point))
                {
                  int i=forward->second;
                  //cout << "\nffound: " << i << " " << _faces[i].str() << endl;
                  return make_pair (triad_t(TRIAD_SIZE*i,TRIAD_SIZE*i+1,TRIAD_SIZE*i+2),
                                    _faces[i].getBarycentricCoordinates (point));
                }
              forward++;
            }
          if (backward != _searchList.begin())
            {
              backward--;
              if (_faces[backward->second].contains(point))
                {
                  int i=backward->second;
//                   cout << "\nbfound: " << i << " " << _faces[i].str() << endl;
//                   cout << _faces[backward->second].contains(point) << endl;
//                   cout << boost::format ("pt:  %.5f %.5f\n") % point.get<0>() % point.get<1>();
//                   cout << "--\n";
                  return make_pair (triad_t(TRIAD_SIZE*i,TRIAD_SIZE*i+1,TRIAD_SIZE*i+2),
                                    _faces[i].getBarycentricCoordinates (point));
                }
            }
        }

      if (point.get<0>() == 0 && point.get<1>() == 0)
        {
          _containsOrigin = false;
          return std::make_pair (triad_t (0,0,0), bpoint3_t(0,0,0));
        }

      cout << "\nusing slow\n";
      return getBarycentricCoordinatesSLOW (point);
    }

    std::vector<Triangle> _faces;
    std::vector<std::pair<double,int> > _searchList;
    mutable bool _containsOrigin;
  };
}

#endif  // CLUSTER_TRIANGULATION_H_
