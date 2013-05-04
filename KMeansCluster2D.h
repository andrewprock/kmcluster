#ifndef CLUSTER_KMEANSCLUSTER_2D_H_
#define CLUSTER_KMEANSCLUSTER_2D_H_

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <cmath>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <rps/random/aprand.h>

using namespace std;

namespace rps
{
  /**
   * use KMeans++ to create K cluster for the input
   * data
   */
  struct Point2D
  {
    double x;
    double y;

    Point2D ()
      : x(0)
      , y(0)
    { }

    Point2D (double xin, double yin)
      : x(xin)
      , y(yin)
    { }

    double distanceSquared (const Point2D& p)
    {
      return pow((p.x-x)*(p.x-x) + (p.y-y)*(p.y-y),1.0/64.0);
      //return sqrt((p.x-x)*(p.x-x) + (p.y-y)*(p.y-y));
    }

    std::string str () const
    {
      return (boost::format("%.6f %.6f") % x % y).str();
    }

    bool operator<(const Point2D& p) const { return x+y < p.x+p.y; }
  };

  class KMeansCluster2D
  {
  public:
    KMeansCluster2D (const std::vector< std::pair<double,double> >& inputData, size_t nClusters)
      : _data(inputData.size())
      , _clusters()
      , _nClusters (nClusters)
    {
      for (size_t i=0; i<inputData.size(); i++)
        _data[i] = PointData(Point2D(inputData[i].first,inputData[i].second));
    }

    std::vector<Point2D> cluster () 
    {
      // select initial seeds for clusters
      for (size_t i=0; i<_nClusters; i++)
        {
          cerr << "+";
          weightDataPoints ();
          selectClusterCenter ();
        }
      cerr << "done initializing\n";
      return KMeansCluster ();
    }

    std::string str() const
    {
      std::string out;
      for (size_t i=0; i<_clusters.size(); i++)
        out += _clusters[i].str();
      return out;
    }

  private:

    void weightDataPoints ()
    {
      size_t newestClusterIndex = _clusters.size()-1;
      for (size_t i=0; i<_data.size(); i++)
        {
          if (_clusters.size() == 0)
            {
              _data[i].weight = 1;
            }
          else
            {
              double dmetric = _data[i].point.distanceSquared (_clusters[newestClusterIndex].getCenter());
              if (_data[i].weight > dmetric)
                _data[i].weight = dmetric;
            }
        }
    }

    struct PointData
    {
      Point2D   point;
      int       clusterid;
      double    weight;

      PointData (const Point2D& p)
        : point(p)
        , clusterid(-1)
        , weight (0.0)
      { }

      PointData ()
        : point()
        , clusterid(-1)
        , weight (0.0)
      { }

      std::string str() const
      {
        return (boost::format("%s %d %.3f\n") 
                % point.str()
                % clusterid
                % weight).str();
      }
    };

    class Cluster
    {
    public:
      Cluster (const Point2D& p)
      {
        _points.insert (p);
        _center = p;
      }

      Cluster ()
      {
      }

      Point2D getCenter () const
      {
        return _center;
      }

      size_t size () const
      {
        return _points.size();
      }

      void remove (const Point2D& p)
      {
        _points.erase (p);
      }

      void insert (const Point2D& p)
      {
        _points.insert (p);
      }

      void calculateCentroid ()
      {
        if (_points.size() == 0)
          {
            _center.x = _center.y = .5;
            return;
          }

        double x=0.0, y=0.0;
        for (std::set<Point2D>::iterator 
               i  = _points.begin();
               i != _points.end();
             ++i)
          {
            x += i->x;
            y += i->y;
          }
        x /= _points.size();
        y /= _points.size();
        _center.x = x;
        _center.y = y;
      }

      std::string str() const
      {
        return _center.str() + "\n";
      }
      
    private:
      std::set<Point2D>      _points;
      Point2D                _center;
    };

    std::vector<PointData>   _data;
    std::vector<Cluster>     _clusters;
    size_t                   _nClusters;

    std::vector<Point2D> KMeansCluster ()
    {
      std::vector<Point2D> centers;
      bool changed = true;
      while (changed)
        {
          changed = false;
          assignAllPoints (changed);
          calculateCentriods ();
        }
      return centers;
    }

    void assignAllPoints (bool & changed)
    {
      static int count = 0;
      int del = 0;
      for (size_t i=0; i<_data.size(); i++)
        {
          Point2D p = _data[i].point;
          size_t c = getNearestCluster (p);
          if (c != _data[i].clusterid)
            {
              if (_data[i].clusterid >= 0)
                _clusters[_data[i].clusterid].remove (p);
              _clusters[c].insert (p);
              _data[i].clusterid = c;
              changed = true;
              del++;
            }
        }
      cerr << "del: " << del << " of: " << _data.size() << " " << ++count << endl;
      if (count > 200 && count > del/2)
        {
          cerr << "non-convergent clustring, inspect cluster visually to verify\n";
          changed = false;
        }
    }

    size_t getNearestCluster (const Point2D& p)
    {
      size_t closest = 0;
      double distance = _clusters[0].getCenter().distanceSquared (p);
      for (size_t i=1; i<_clusters.size(); i++)
        {
          double d = _clusters[i].getCenter().distanceSquared (p);
          if (d < distance)
            {
              distance = d;
              closest = i;
            }
        }
      return closest;
    }

    void calculateCentriods ()
    {
      for (size_t c=0; c<_clusters.size(); c++)
        {
          _clusters[c].calculateCentroid ();
        }
    }

    void selectClusterCenter ()
    {
      double pick = aprand (getTotalPointWeight ());
      double running = 0.0;
      for (size_t i=0; i<_data.size(); i++)
        {
          if (running + _data[i].weight > pick)
            {
              _clusters.push_back (Cluster(_data[i].point));
              return;
            }
          running += _data[i].weight;
        }
    }

    double getTotalPointWeight ()
    {
      double tot=0.0;
      for (size_t i=0; i<_data.size(); i++)
        tot += _data[i].weight;
      return tot;
    }

  };

}

#endif  // CLUSTER_KMEANSCLUSTER_2D_H_
