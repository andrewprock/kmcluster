#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <kmcluster/KMeansCluster.h>

// compile:  g++ testcluster.cpp ../random/rand_isaac.cpp -I../.. -o cluster

using namespace std;
using namespace boost;

int main (int argc, char ** argv)
{
  string fname     = argv[1];
  int    nClusters = atoi(argv[2]);
  
  rps::KMeansClusterND clusters (nClusters);

  // we support two different file types
  // .txt is just a flat file with one 2D point per line
  ifstream fin (fname.c_str());
  if (!fin)
    {
      cerr << "could not open file: " << fname << endl;
      exit(-1);
    }

  if (boost::ends_with (fname, ".txt"))
    {
      do
        {
          double x,y;
          fin >> x >> y;
          if (!fin.eof ())
            {
              rps::PointND pt(x,y);
              clusters.add (pt);
            }
        }
      while (!fin.eof());
    }

  // we also support a csv format where each line is a
  // labeled points, with the label stored in the first
  // field of the csv
  else if (boost::ends_with (fname, ".csv"))
    {
      //cout << "doing csv\n";
      string line;
      while (!fin.eof())
        {
          getline (fin, line);

          // filer out empty lines
          if (line.size() == 0)
            continue;

          vector<string> fields;
          boost::split (fields, line, boost::is_any_of(","));

          vector<double> pt;
          for (size_t i=1; i<fields.size(); i++)
            pt.push_back (lexical_cast<double>(fields[i]));

          clusters.add (rps::PointND(fields[0], pt));
        }
    }
  
  //srand ();
  
  clusters.cluster ();
  //cout << clusters.str () << endl;

  cout << clusters.clusterSets () << endl;
}
