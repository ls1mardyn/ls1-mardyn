#include "Domain.h"

#include <cstdlib>
#include <iostream>
#include <stack>

Domain::Domain(int ingrid)
{
   this->grid = ingrid;
   this->cavities = set<long>();
   this->verlet = map< long, set<long> >();
   this->clusterID = map<long, long>();
   this->clusterVertices = map< long, set<long> >();
}

long Domain::encode(int xgrid, int ygrid, int zgrid)
{
   if((xgrid < 0) || (ygrid < 0) || (zgrid < 0) || (xgrid >= grid) || (ygrid >= grid) || (zgrid >= grid)) exit(21);
   
   return xgrid*grid*grid + ygrid*grid + zgrid;
}

void Domain::decode(int* qgrid, long code)
{
   long rest;
   if((code < 0) || (code >= grid*grid*grid)) { cout << "invalid code " << code << ".\n"; exit(22); }
   
   rest = code % grid;
   qgrid[2] = rest;
   code = (code - rest) / grid;
   rest = code % grid;
   qgrid[1] = rest;
   code = (code - rest) / grid;
   rest = code % grid;
   qgrid[0] = rest;
   
   return;
}

bool Domain::neighbours(long code, set<long>* vicinity)
{
   // cout << "\t{";
   vicinity->clear();
   if(this->cavities.count(code) > 0)
   {
      // cout << code;
      int x[3];
      this->decode(x, code);
      // cout << " ; (" << x[0] << "/" << x[1] << "/" << x[2] << ")";
      
      int c[3];
      int y[3];
      long ycode;
      for(c[0] = -1; 1 >= c[0]; c[0]++)
      {
         for(c[1] = -1; 1 >= c[1]; c[1]++)
         {
            for(c[2] = -1; 1 >= c[2]; c[2]++)
            {
               if((c[0] == 0) || (c[1] == 0) || (c[2] == 0))
               {
                  for(int d = 0; d < 3; d++)
                  {
                     y[d] = (x[d] + c[d]);
                     if(y[d] < 0) y[d] += grid;
                     else if(y[d] >= grid) y[d] -= grid;
                  }
                  // cout << " -?-> (" << y[0] << "/" << y[1] << "/" << y[2] << ")";
                  ycode = this->encode(y[0], y[1], y[2]);
                  if(this->cavities.count(ycode) > 0) vicinity->insert(ycode);
               }
            }
         }
      }
   }
   // cout << "}";
   return !vicinity->empty();
}

void Domain::build_verlet()
{
   set<long>::iterator cavit;
   for(cavit = this->cavities.begin(); cavit != this->cavities.end(); cavit++)
   {
      this->verlet[*cavit] = set<long>();
      this->neighbours(*cavit, &(this->verlet[*cavit]));
   }
}

void Domain::detectClusters()
{
   set<long>::iterator acti, actj;
   
   long lowlink, tnode;
   set<long> processed_nodes;
   set<long> present_nodes;
   set<long> unprocessed_nodes;
   
   stack<long> dfs_stack;
   map<long, set<long>::iterator> edgeit;
   
   for(acti = this->cavities.begin(); acti != this->cavities.end(); acti++)
   {
      // cout << "\t{" << *acti << "}";
      this->attach(*acti, *acti, false);
      unprocessed_nodes.insert(*acti);
   }
   
   // cout << "\nBuilding Verlet list:";
   this->build_verlet();
   // cout << " done.\n";
   
   while(!unprocessed_nodes.empty())
   {
      present_nodes.clear();
      dfs_stack.push( *(unprocessed_nodes.begin()) );
      lowlink = dfs_stack.top();
      
      while(!dfs_stack.empty())
      {
         tnode = dfs_stack.top();
         if(unprocessed_nodes.count(tnode) > 0)
         {
            if(this->clusterID[tnode] < lowlink)
            {
               lowlink = this->clusterID[tnode];
            }
            unprocessed_nodes.erase(tnode);
            present_nodes.insert(tnode);
            edgeit[tnode] = this->verlet[tnode].begin();
         }
         
         while((edgeit[tnode] != this->verlet[tnode].end()) && (unprocessed_nodes.count(*edgeit[tnode]) == 0))
         {
            edgeit[tnode] ++;
         }
         if(edgeit[tnode] != this->verlet[tnode].end())
         {
            tnode = *edgeit[tnode];
            dfs_stack.push(tnode);
         }
         else
         {
            dfs_stack.pop();
         }
      }
      for(acti = present_nodes.begin(); acti != present_nodes.end(); acti++)
      {
         processed_nodes.insert(*acti);
         attach(*acti, lowlink, true);
      }
   }
}

void Domain::attach(long vertex, long cluster, bool detach_previous)
{
   if(cluster > vertex)
   {
      cout << "\nCavity " << vertex << " cannot be attached to cluster " << cluster << ".\n";
      exit(23);
   }
   
   if(detach_previous && (clusterID.count(vertex) > 0))
   {
      if(clusterVertices.count(clusterID[vertex]) > 0)
      {
         this->clusterVertices[clusterID[vertex]].erase(vertex);
         if(this->clusterVertices[clusterID[vertex]].empty())
         {
            this->clusterVertices.erase(clusterID[vertex]);
         }
      }
   }
   
   this->clusterID[vertex] = cluster;
   if(this->clusterVertices.count(cluster) == 0)
   {
      this->clusterVertices[cluster] = set<long>();
   }
   this->clusterVertices[cluster].insert(vertex);
}

unsigned Domain::countClusters(map<unsigned, unsigned>* thresholds)
{
   map<unsigned, unsigned> exactPopulation;
   map< long, set<long> >::iterator cluvit;
   for(cluvit = this->clusterVertices.begin(); cluvit != this->clusterVertices.end(); cluvit++)
   {
      unsigned cavsize = cluvit->second.size();
      if(exactPopulation.count(cavsize) == 0) exactPopulation[cavsize] = 1;
      else exactPopulation[cavsize] ++;
   }
   
   unsigned largest = 0;
   map<unsigned, unsigned>::iterator threshit, popit;
   for(threshit = thresholds->begin(); threshit != thresholds->end(); *threshit++)
   {
      (*thresholds)[threshit->first] = 0;
   }
   for(popit = exactPopulation.begin(); popit != exactPopulation.end(); popit++)
   {
      if(popit->first > largest) largest = popit->first;
      for(threshit = thresholds->begin(); threshit != thresholds->end(); *threshit++)
      {
         if(popit->first >= threshit->first) (*thresholds)[threshit->first] += popit->second;
      }
   }
   return largest;
}
