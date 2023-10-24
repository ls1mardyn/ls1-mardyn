#ifndef DOMAIN_H
#define DOMAIN_H

#include <map>
#include <set>


class Domain
{
public:
   Domain(int ingrid);

   unsigned size() { return cavities.size(); }
   // bool contains(int x, int y, int z) { return (cavities.count(encode(x,y,z)) != 0); }
   void insert(int x, int y, int z) { cavities.insert(encode(x,y,z)); }

   long encode(int xgrid, int ygrid, int zgrid);
   void decode(int* qgrid, long code);

   bool neighbours(long code, std::set<long>* vicinity);
   void build_verlet();

   void detectClusters();
   void attach(long vertex, long cluster, bool detach_previous);
   unsigned countClusters(std::map<unsigned, unsigned>* thresholds);

private:
   int grid;
   std::set<long> cavities;

   std::map< long, std::set<long> > verlet;  // for all cavities, contains the neighbour cavities

   std::map<long, long> clusterID;  // contains cluster ID of the cavities (i.e. smallest cavity ID from the cluster)
   std::map< long, std::set<long> > clusterVertices;  // all cavities within the cluster
};

#endif
