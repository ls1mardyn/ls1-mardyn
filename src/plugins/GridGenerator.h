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
#include "plugins/PluginBase.h"

struct ElementInfo{

    double length_x, length_y, length_z;
    int elements_in_x, elements_in_y, elements_in_z;
    int index;
    double x_width;
    double y_width;
    double z_width;
    double volume;

    std::tuple<int, int, int> local_indeces;

    std::tuple<int, int, int> GlobalToLocalIndex(int idx){

        int x, y, z;

        std::tuple<int, int, int> local;

        x= std::fmod(idx,elements_in_x);
        y= std::fmod((idx/elements_in_x),elements_in_y);
        z=idx/(elements_in_x*elements_in_y);
        local = std::make_tuple(x,y,z);

        return local;
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
        ElementInfo* info;
    public:

    PropertySampler();

    void ParticlePerCellCount(ParticleContainer* pc);
    std::vector<int>& GetParticlesPerCell();
    void init(ElementInfo* info);

    private:


    std::tuple<int, int, int> GetParticleLocalCellIndices(ParticleIterator it);


};

class GridGenerator: public PluginBase{

    private:
    ElementInfo element_information;
    std::vector<int> elements_per_dimension{1,1,1};
    std::vector<double> element_width_per_dimension{0,0,0};
    std::vector<double> lower_corner{3};
    std::vector<double> upper_corner{3};
    int total_elements;
    PropertySampler sampler;
    ParticleContainer* particle_container;


    public:
        GridGenerator();
        ~GridGenerator(){};
        void readXML(XMLfileUnits& xmlconfig) override;
        void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override;
        void beforeEventNewTimestep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override;
        void beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override{}
        void afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override{}
        void endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep) override{}
        void finish(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override{}

        void Output(std::string prefix, long unsigned accumulated);

        std::string getPluginName()  {return "GridGenerator";}
    static PluginBase* createInstance() {return new GridGenerator(); }
    private:

    void MeshEntireDomain();
    void SetElementInfo();
    void SetTotalElements();

};


