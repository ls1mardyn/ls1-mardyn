/*
 * Created on Wed Feb 14 2024
 * Author: Jose A. Pinzon Escobar
 * Email to: jose.escobar@hsu-hh.de
 * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
 */

#pragma once

#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include<string>
#include<numeric>
struct Node{
    std::vector<double> position{0,0,0};
    void SetPosition(double x, double y, double z){
        position={x,y,z};
    }
    std::vector<double>& GetPosition(){
        return position;
    }
};


struct NodeInformation{
    std::vector<int> nodes_per_dimension{1,1,1};
    std::vector<Node> nodes;
    int total_nodes = 1;

    std::vector<Node>& GetNodes(){
        return this->nodes;
    }

    std::tuple<int, int, int> MapGlobalToLocal(int idx){

        int x, y, z;

        std::tuple<int, int, int> local;

        x= std::fmod(idx,nodes_per_dimension[0]);
        y= std::fmod((idx/nodes_per_dimension[0]),nodes_per_dimension[1]);
        z=idx/(nodes_per_dimension[0]*nodes_per_dimension[1]);
        local = std::make_tuple(x,y,z);

        return local;
    }

    int MapLocalToGlobalIndex(std::tuple<int, int, int> index_vector){
        return MapLocalToGlobalIndex(std::get<0>(index_vector),std::get<1>(index_vector),std::get<2>(index_vector));
    }

    int MapLocalToGlobalIndex(int i, int j, int k){

        int global_index = i + j*nodes_per_dimension[0] + k*nodes_per_dimension[0]*nodes_per_dimension[1];

        return global_index;
    }
};


struct ElementInfo{

    double length_x, length_y, length_z;
    int elements_in_x, elements_in_y, elements_in_z;
    std::vector<int> elements_per_dimension{1,1,1};
    std::vector<double> element_width_per_dimension{0,0,0};
    int index;
    double x_width;
    double y_width;
    double z_width;
    double volume;


    std::tuple<int, int, int> local_indeces;

    std::array<int, 3> TupleToArray(std::tuple<int,int,int> indeces){
        std::array<int, 3> indeces_array;
        indeces_array[0]=std::get<0>(indeces);
        indeces_array[1]=std::get<1>(indeces);
        indeces_array[2]=std::get<2>(indeces);

        return indeces_array;
    }

    std::tuple<int, int, int> GlobalToLocalIndex(int idx){

        int x, y, z;

        std::tuple<int, int, int> local;

        x= std::fmod(idx,elements_in_x);
        y= std::fmod((idx/elements_in_x),elements_in_y);
        z=idx/(elements_in_x*elements_in_y);
        local = std::make_tuple(x,y,z);

        return local;
    }

    int GetNeighbor(int idx, int direction, int dim){
        int neighbor;
        std::array<int, 3> local_current = TupleToArray(GlobalToLocalIndex(idx));
        //check if out of bounds, handle at discretization (grid) level
        neighbor = local_current[dim]+direction;
        local_current[dim] += direction;
        if(IsBorder(neighbor,dim)){
            return idx;
        }
        neighbor=LocalToGlobalIndex(local_current[0],local_current[1],local_current[2]);

        return neighbor;


    }

    bool IsBorder(int neighbor, int dim){
        if(neighbor<0 || neighbor>=elements_per_dimension[dim]){
            return true;
        }
        return false;
    }

    int GetGlobalRightNeighbor(int idx){
        int right_neighbor;
        std::tuple<int,int,int> local_current = GlobalToLocalIndex(idx);
        //check if out of bounds, hanlde at discretization (grid) level
        right_neighbor=std::get<0>(local_current)+1;
        if(right_neighbor<0 || right_neighbor>=elements_in_x){
            return -1;//is boundary to the right
        }
        //Get the global index
        right_neighbor=LocalToGlobalIndex(right_neighbor,std::get<1>(local_current),std::get<2>(local_current));

    }

    int LocalToGlobalIndex(std::tuple<int, int, int> index_vector){
        return LocalToGlobalIndex(std::get<0>(index_vector),std::get<1>(index_vector),std::get<2>(index_vector));
    }

    int LocalToGlobalIndex(int i, int j, int k){

        int global_index = i + j*elements_in_x + k*elements_in_x*elements_in_y;

        return global_index;
    }

    std::vector<double> GetElementLowerCorner(int localx, int localy, int localz){

        std::vector<double> lower_corner(3);
        lower_corner[0]=localx*x_width;
        lower_corner[1]=localy*y_width;
        lower_corner[2]=localz*z_width;
        
        return lower_corner;
    }

    std::vector<double> GetElementUpperCorner(int localx, int localy, int localz){

        std::vector<double> upper_corner(3);
        upper_corner[0]=(localx+1)*x_width;
        upper_corner[1]=(localy+1)*y_width;
        upper_corner[2]=(localz+1)*z_width;
        
        return upper_corner;
    }

    std::vector<double> GetElementLowerCorner(int global_index){
        std::tuple<int,int,int> local_index = GlobalToLocalIndex(global_index);
        return GetElementLowerCorner(std::get<0>(local_index),std::get<1>(local_index),std::get<2>(local_index));
    }

    std::vector<double> GetElementUpperCorner(int global_index){
        std::tuple<int,int,int> local_index = GlobalToLocalIndex(global_index);
        return GetElementUpperCorner(std::get<0>(local_index),std::get<1>(local_index),std::get<2>(local_index));
    }

    std::vector<double> GetElementCenter(int global_index){
        std::vector<double> center(3);
        std::vector<double> lc = GetElementLowerCorner(global_index);
        std::vector<double> uc = GetElementUpperCorner(global_index);

        center[0]= 0.5*(lc[0]+uc[0]);
        center[1]= 0.5*(lc[1]+uc[1]);
        center[2]= 0.5*(lc[2]+uc[2]);

        return center;
    }


    void MeshTraversal(){
        for(int idx=0;idx<=index;idx++){
            local_indeces = GlobalToLocalIndex(idx);
            int indx_x, indx_y, indx_z;
            indx_x=std::get<0>(local_indeces);
            indx_y=std::get<1>(local_indeces);
            indx_z=std::get<2>(local_indeces);

            std::vector<double> ucorner, lcorner;

            ucorner = GetElementUpperCorner(indx_x, indx_y, indx_z);
            lcorner = GetElementLowerCorner(indx_x, indx_y, indx_z);

            std::cout<<"\tGlobal element: "<<idx<<" has local indeces: \n";
            std::cout<<"\t("<<indx_x<<","<<indx_y<<","<<indx_z<<")\n";
            std::cout<<"\tThe bounding box of global element: "<<LocalToGlobalIndex(indx_x, indx_y, indx_z)<<" is: \n"
            <<"\t("<<lcorner[0]<<","<<lcorner[1]<<","<<lcorner[2]<<")x("<<ucorner[0]<<","<<ucorner[1]<<","<<ucorner[2]<<")\n";
        }
    }




};

class PropertySampler{
    
    private:
        std::vector<int> particles_per_cell;
        std::vector<double> material_density;
        ElementInfo* info;
        NodeInformation* node_info;
        double measure_radius;
        std::vector<int> particles_per_node; 
    public:

    PropertySampler();

    void ParticlePerCellCount(ParticleContainer* pc);
    void ComputeMaterialDensityPerCell(ParticleContainer* pc);
    std::vector<int>& GetParticlesPerCell();
    std::vector<double>& GetMaterialDensityPerCell();
    std::vector<int>& GetParticlesPerNode();
    void init(ElementInfo* info, NodeInformation* node_info);
    //Iterate all nodes, create a virtual sphere centered at node position, check which particles are inside which sphere (can be repeated)
    //Or iterate all particles?
    //Use Euclidean norm?
    void SampleAtNodes(ParticleContainer* pc);
    double ComputeMaterialDensityAtPosition(ParticleContainer* pc, std::array<double, 3>& pos);
    bool ParticleInsideMeasuringSpace(std::vector<double> nodal_pos, std::array<double, 3> par_pos);

    private:

    std::tuple<int, int, int> GetParticleLocalCellIndices(ParticleIterator it);
    //need access to the nodes in grid
};

class GridGenerator{

    private:
    ElementInfo element_information;
    std::vector<int> elements_per_dimension{1,1,1};
    std::vector<double> element_width_per_dimension{0,0,0};
    std::vector<double> lower_corner{3};
    std::vector<double> upper_corner{3};
    int total_elements;
    PropertySampler sampler;
    ParticleContainer* particle_container;
    NodeInformation node_information;
    //std::vector<Node> nodes;//initialized once we know how many elements per dimension


    public:
        GridGenerator();
        ~GridGenerator(){};
        void SetGridGenerator(int els_x, int els_y, int els_z);
        //void readXML(XMLfileUnits& xmlconfig) override;
        void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain);
        //void beforeEventNewTimestep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override;
        //void beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override{}
        //void afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override{}
        //void endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep) override{}
        //void finish(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override{}

        int GetTotalElements(){
            return this->total_elements;
        }
        ElementInfo& GetElementInfo(){
            return this->element_information;
        }
        PropertySampler& GetPropertySampler();
        void Output();
        void OutputMeshInformation();
        void OutputNodalInformation();
        void OutputElementConnectivity();
        void OutputNodalDensityValues(unsigned long step);
        void OutputPropertyPerCell(unsigned long step);
        void OutputMaterialDensityPerCell(unsigned long cell);

        //Node stuff
        std::array<int, 8> GetElementGlobalNodeIndeces(int el);
        NodeInformation& GetNodeInfo(){
            return this->node_information;
        }

    private:

    void MeshEntireDomain();
    void SetElementInfo();
    void SetTotalElements();
    void SetNodeInfo();

    void InitNodePositions();

};


