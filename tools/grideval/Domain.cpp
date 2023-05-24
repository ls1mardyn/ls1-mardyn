#include "Domain.h"

#include <cstdlib>
#include <iostream>
#include <stack>

Domain::Domain(int ingrid)
{
   this->grid = ingrid;
   this->cavities = std::set<long>();
   this->verlet = std::map< long, std::set<long> >();
   this->clusterID = std::map<long, long>();
   this->clusterVertices = std::map< long, std::set<long> >();
}

long Domain::encode(int xgrid, int ygrid, int zgrid)
{
   if((xgrid < 0) || (ygrid < 0) || (zgrid < 0) || (xgrid >= grid) || (ygrid >= grid) || (zgrid >= grid)) exit(21);
   
   return xgrid*grid*grid + ygrid*grid + zgrid;
}

void Domain::decode(int* qgrid, long code)
{
   long rest;
   if((code < 0) || (code >= grid*grid*grid)) { std::cout << "invalid code " << code << ".\n"; exit(22); }
   
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

bool Domain::neighbours(long code, std::set<long>* vicinity)
{
   // std::cout << "\t{";
   vicinity->clear();
   if(this->cavities.count(code) > 0)
   {
      // std::cout << code;
      int x[3];
      this->decode(x, code);
      // std::cout << " ; (" << x[0] << "/" << x[1] << "/" << x[2] << ")";
      
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
                  // std::cout << " -?-> (" << y[0] << "/" << y[1] << "/" << y[2] << ")";
                  ycode = this->encode(y[0], y[1], y[2]);
                  if(this->cavities.count(ycode) > 0) vicinity->insert(ycode);
               }
            }
         }
      }
   }
   // std::cout << "}";
   return !vicinity->empty();
}

void Domain::build_verlet()
{
   std::set<long>::iterator cavit;
   for(cavit = this->cavities.begin(); cavit != this->cavities.end(); cavit++)
   {
      this->verlet[*cavit] = std::set<long>();
      this->neighbours(*cavit, &(this->verlet[*cavit]));
   }
}

void Domain::detectClusters()
{
   std::set<long>::iterator acti, actj;
   
   long lowlink, tnode;
   std::set<long> processed_nodes;
   std::set<long> present_nodes;
   std::set<long> unprocessed_nodes;
   
   std::stack<long> dfs_stack;
   std::map<long, std::set<long>::iterator> edgeit;
   
   for(acti = this->cavities.begin(); acti != this->cavities.end(); acti++)
   {
      // std::cout << "\t{" << *acti << "}";
      this->attach(*acti, *acti, false);
      unprocessed_nodes.insert(*acti);
   }
   
   // std::cout << "\nBuilding Verlet list:";
   this->build_verlet();
   // std::cout << " done.\n";
   
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
      std::cout << "\nCavity " << vertex << " cannot be attached to cluster " << cluster << ".\n";
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
      this->clusterVertices[cluster] = std::set<long>();
   }
   this->clusterVertices[cluster].insert(vertex);
}

unsigned Domain::countClusters(std::map<unsigned, unsigned>* thresholds)
{
   std::map<unsigned, unsigned> exactPopulation;
   std::map< long, std::set<long> >::iterator cluvit;
   for(cluvit = this->clusterVertices.begin(); cluvit != this->clusterVertices.end(); cluvit++)
   {
      unsigned cavsize = cluvit->second.size();
      if(exactPopulation.count(cavsize) == 0) exactPopulation[cavsize] = 1;
      else exactPopulation[cavsize] ++;
   }
   
   unsigned largest = 0;
   std::map<unsigned, unsigned>::iterator threshit, popit;
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
