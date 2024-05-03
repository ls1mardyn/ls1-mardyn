/*
 * Created on Wed Apr 24 2024
 * Author: Jose A. Pinzon Escobar
 * Email to: jose.escobar@hsu-hh.de
 * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
 */

#pragma once
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "Maps.h"
#include<iostream>
#include "Subset.h"
//Set array type
using IdxArray=std::array<int, 3>;
using DataArray = std::array<double, 3>;

//TODO: all these need iterators and that kind of stuff
//Use vector only where necessary, stop using tuple its stupid
class Element{
    private:

    struct Properties{

        int index;
        double x_width;//Not set
        double y_width;//Not set
        double z_width;//Not set
        double volume;//Not set
        IdxArray local_indeces;

    };

    Properties properties;

    public:
    Properties& GetProperties(){
        return properties;
    }

};


class Node{
    public:

    DataArray& GetPosition(){
        return this->position;
    }
    void SetPosition(double x, double y, double z){
        this->position={x,y,z};
    }
    
    private:

    DataArray position;
  
};

//Should hanlde both other classes entirely
// Is blind wrt the Sampler or its actual use
class Grid{
    private:

    struct NodeInfo{
        IdxArray nodes_per_dimension{1,1,1};
        int total_nodes = 1;
    };

    struct ElementInfo{

        IdxArray elements_per_dimension{1,1,1};
        DataArray element_width_per_dimension{0,0,0};
        int total_elements;
        int largest_index;
        double volume;
        
    };

    private:
    NodeInfo nodal_info;
    ElementInfo element_info;
    std::vector<Node> nodes;
    std::vector<Element> elements;
    std::array<double,3> lower_corner;
    std::array<double,3> upper_corner;
    ParticleContainer* particle_container;
    std::vector<Subset> subsets;


    public:

    Grid();
    Grid(std::array<double, 3> lower, std::array<double, 3> upper, int x, int y, int z);
    ~Grid(){};
    void SetMeshLimits(DataArray lower, DataArray upper);
    void MeshAllDomain();
    void StartGrid(int els_x, int els_y, int els_z);
    void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain);
    std::array<double, 3>& GetUpperCorner();
    std::array<double, 3>& GetLowerCorner();
    ElementInfo& GetElementInfo(){
        return element_info;
    }
    NodeInfo& GetNodeInfo(){
        return nodal_info;
    }
    std::vector<Node>& GetNodes(){
        return this->nodes;
    }

    std::array<int, 8> GetElementGlobalNodeIndeces(int el);
    void AddSubset(std::vector<int> list, std::string Name);
    Subset& GetSubset(std::string& Name);
    std::vector<Subset>& GetSubsets();
    friend std::ostream& operator<<(std::ostream& out, Grid& grid);


    private:
    void InitNodePositions();

};

