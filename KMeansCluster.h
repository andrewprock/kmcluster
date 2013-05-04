#ifndef CLUSTER_KMEANSCLUSTER_H_
#define CLUSTER_KMEANSCLUSTER_H_

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <functional>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <boost/bind.hpp>
#include <rps/random/aprand.h>

using namespace std;

namespace rps
{
  struct PointND
  {
    std::string label;
    std::vector<double> x;


    PointND ()
      : label()
      , x()
    { }

    // simple shorthand for 2D points
    PointND (double a, double b)
      : label()
      , x(2)
    { 
      x[0] = a;
      x[1] = b;
    }

    PointND (const std::vector<double> input)
      : label()
      , x(input.begin(), input.end())
    { 
    }

    PointND (const std::string& inlabel, const std::vector<double> input)
      : label(inlabel)
      , x(input.begin(), input.end())
    { }

    double distance (const PointND& p) const
    {
      if (p.x.size () != x.size())
        cout << boost::format ("size mismatch: %d != %d\n %s\n %s\n") % p.x.size() % x.size () % p.str() % this->str();
      assert (p.x.size() == x.size());
      double dist = 0;
      for (size_t i=0; i<x.size(); i++)
        dist += pow (x[i]-p.x[i], 2);
      //dist += x[i]*p.x[i];
      return sqrt(dist);
      //return (dist);
    }
    
    std::string str() const
    {
      std::string out;
      if (label.size() > 0)
        out += label + ",";
      for (size_t i=0; i<x.size(); i++)
        out += boost::str (boost::format("%.12f,") % x[i]);
      return out;
    }

    bool operator<(const PointND& p) const
    { 
      double myval = std::accumulate(x.begin(), x.end(),0.0);
      double pval  = std::accumulate(p.x.begin(), p.x.end(),0.0);
      return myval < pval;
    }

    void clear()
    {
      x.clear ();
    }

    PointND& operator+=(const PointND& p)
    {
      if (x.size() == 0)
        x.resize (p.x.size(), 0);
      std::transform (x.begin(), x.end(), p.x.begin(), x.begin(), std::plus<double>());
      return *this;
    }

    PointND& operator/=(double v)
    {
      std::transform (x.begin(), x.end(), x.begin(), boost::bind(std::divides<double>(),_1,v));
      return *this;
    }
  };

  class KMeansClusterND
  {
  public:
    KMeansClusterND (size_t nClusters)
      : _data()
      , _clusters()
      , _nClusters (nClusters)
    { }

    void add (const PointND& p)
    {
      _data.push_back (p);
    }

    std::vector<PointND> cluster () 
    {
      // select initial seeds for clusters
      for (size_t i=0; i<_nClusters; i++)
        {
          //std::cerr << "+";
          weightDataPoints ();
          selectClusterCenter ();
        }
      //std::cerr << "done initializing\n";
      return KMeansCluster ();
    }

    std::string str() const
    {
      std::string out;
      for (size_t i=0; i<_clusters.size(); i++)
        out += _clusters[i].str();
      return out;
    }

    std::string clusterSets () const
    {
      std::vector<std::string> clusterStrings(_nClusters);

      double totalSpread = 0.0;
      for (size_t i=0; i<_clusters.size(); i++)
        {
          totalSpread += _clusters[i].spread();
          clusterStrings[i] = (boost::format("cluster %d spread %f\n") % i % _clusters[i].spread()).str();
        }
      cout << "total spread: " << totalSpread << endl;
      
      for (size_t i=0; i<_data.size(); i++)
        {
          PointND p = _data[i].point;
          size_t c = getNearestCluster (p);
          clusterStrings[c] += p.str () + "\n";
        }

      std::string ret;

      for (size_t i=0; i<_nClusters; i++)
        ret += clusterStrings[i] + "\n";

      return ret;
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
              //cout << _data[i].point.str() << endl;
              double dmetric = _data[i].point.distance (_clusters[newestClusterIndex].getCenter());
              if (_data[i].weight > dmetric)
                _data[i].weight = dmetric;
            }
        }
    }

    //
    // private struct
    //
    
    struct PointData
    {
      PointND   point;
      int       clusterid;
      double    weight;

      PointData (const PointND& p)
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

    //
    // private class
    //

    class Cluster
    {
    public:
      Cluster (const PointND& p)
      {
        _points.insert (p);
        _center = p;
      }

      Cluster ()
      {
      }

      PointND getCenter () const
      {
        return _center;
      }

      size_t size () const
      {
        return _points.size();
      }

      void remove (const PointND& p)
      {
        _points.erase (p);
      }

      void insert (const PointND& p)
      {
        _points.insert (p);
      }

      double spread () const
      {
        double spread = 0.0;
        for (std::set<PointND>::iterator 
               i  = _points.begin();
               i != _points.end();
             ++i)
          {
            spread += i->distance(_center);
          }
        //cout << spread << endl;
        return spread/_points.size();
      }

      void calculateCentroid ()
      {
        if (_points.size() == 0)
          {
            fill (_center.x.begin(), _center.x.end(), 0.5);
            return;
          }

        _center.clear();
        for (std::set<PointND>::iterator 
               i  = _points.begin();
               i != _points.end();
             ++i)
          {
            _center += *i;
          }
        _center /= _points.size();
      }

      std::string str() const
      {
        return _center.str() + "\n";
      }
      
    private:
      std::set<PointND>      _points;
      PointND                _center;
    };


    //
    // private data
    //

    std::vector<PointData>   _data;
    std::vector<Cluster>     _clusters;
    size_t                   _nClusters;

    std::vector<PointND> KMeansCluster ()
    {
      std::vector<PointND> centers;
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
          PointND p = _data[i].point;
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
      //std::cerr << "del: " << del << " of: " << _data.size() << " " << ++count << std::endl;
      if (count > 200 && count > del/2)
        {
          std::cerr << "non-convergent clustring, inspect cluster visually to verify\n";
          changed = false;
        }
    }
    
    size_t getNearestCluster (const PointND& p) const
    {
      size_t closest = 0;
      double distance = _clusters[0].getCenter().distance (p);
      for (size_t i=1; i<_clusters.size(); i++)
        {
          double d = _clusters[i].getCenter().distance (p);
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

#endif  // CLUSTER_KMEANSCLUSTER_H_
