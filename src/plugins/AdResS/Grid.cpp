/*
 * Created on Wed Apr 24 2024
 * Author: Jose A. Pinzon Escobar
 * Email to: jose.escobar@hsu-hh.de
 * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
 */

#include "Grid.h"

Grid::Grid(){

}

Grid::Grid(DataArray low, DataArray up, int x, int y, int z):lower_corner(low),upper_corner(up){
    StartGrid(x,y,z);
}

void Grid::MeshAllDomain(){
    DataArray upper;
    upper[0] = _simulation.getDomain()->getGlobalLength(0);
    upper[1] = _simulation.getDomain()->getGlobalLength(1);
    upper[2] = _simulation.getDomain()->getGlobalLength(2);

    DataArray lower={0.0,0.0,0.0};
    this->SetMeshLimits(lower, upper);
}


void Grid::SetMeshLimits(DataArray lower, DataArray upper){

    this->lower_corner=lower;
    this->upper_corner=upper;

}

//TODO:called from another class or in constructor? or in init? by class that reads xml
void Grid::StartGrid(int x, int y, int z){

    element_info.elements_per_dimension[0]=x;
    element_info.elements_per_dimension[1]=y;
    element_info.elements_per_dimension[2]=z;

    element_info.total_elements=element_info.elements_per_dimension[0]*element_info.elements_per_dimension[1]*element_info.elements_per_dimension[2];

    //TODO: only if limits already set
    for(int d=0;d<3;d++){
        element_info.element_width_per_dimension[d]= (upper_corner[d]-lower_corner[d])/(double)element_info.elements_per_dimension[d];
    }

    element_info.volume=element_info.element_width_per_dimension[0]*element_info.element_width_per_dimension[1]*element_info.element_width_per_dimension[2];

    //Start all elements
    this->elements.resize(element_info.total_elements);
    element_info.largest_index=element_info.total_elements-1;
    for(int idx=0;idx<elements.size();idx++){
        elements[idx].GetProperties().index=idx;
        elements[idx].GetProperties().local_indeces = MapGlobalToLocal<IdxArray>(idx,element_info.elements_per_dimension);
    }

    //Set node info
    nodal_info.total_nodes=(element_info.elements_per_dimension[0]+1)*(element_info.elements_per_dimension[1]+1)*(element_info.elements_per_dimension[2]+1);
    nodes.resize(nodal_info.total_nodes);
    nodal_info.nodes_per_dimension[0]=element_info.elements_per_dimension[0]+1;
    nodal_info.nodes_per_dimension[1]=element_info.elements_per_dimension[1]+1;
    nodal_info.nodes_per_dimension[2]=element_info.elements_per_dimension[2]+1;

    InitNodePositions();


}

void Grid::init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain){

    DataArray lower = {0.0,0.0,0.0};
    DataArray upper = {domain->getGlobalLength(0),domain->getGlobalLength(1),domain->getGlobalLength(2)};
    SetMeshLimits(lower, upper);
    this->InitNodePositions();

}

void Grid::InitNodePositions(){
    for(int idx=0;idx<nodes.size();idx++){
        IdxArray local_indeces =  MapGlobalToLocal<IdxArray>(idx,nodal_info.nodes_per_dimension);
        //IdxArray local_indeces = this->node_information.MapGlobalToLocal(idx);
        double x = (double)local_indeces[0]*element_info.element_width_per_dimension[0];
        double y = (double)local_indeces[1]*element_info.element_width_per_dimension[1];
        double z = (double)local_indeces[2]*element_info.element_width_per_dimension[2];
        nodes[idx].SetPosition(x,y,z);
        
    }
}

std::array<int, 8> Grid::GetElementGlobalNodeIndeces(int elIdx){
    std::array<int, 8> element_nodes_idx_global;

    IdxArray el_local_indcs = MapGlobalToLocal<IdxArray>(elIdx, element_info.elements_per_dimension);
    int counter =0;
    for(int z=0;z<=1;z++){
        for(int y=0;y<=1;y++){
            for(int x=0;x<=1;x++){
                element_nodes_idx_global[counter] = MapLocalToGlobal<IdxArray>(el_local_indcs[0]+x,el_local_indcs[1]+y,el_local_indcs[2]+z,nodal_info.nodes_per_dimension);
                counter++;
            }
        }
    }


    return element_nodes_idx_global;
}



std::ostream& operator<<(std::ostream& out, Grid& grid){
    std::string prefix = "//[Grid]: ";
    out<<prefix+"Meshed region is: ("<<grid.lower_corner[0]<<","<<grid.lower_corner[1]<<","<<grid.lower_corner[2]<<")x("<<grid.upper_corner[0]<<","<<grid.upper_corner[1]<<","<<grid.upper_corner[2]<<")\n";
    out<<prefix+"Total elements: "<<grid.element_info.total_elements<<"\n"
       <<prefix+"Length per direction: ["<<grid.element_info.element_width_per_dimension[0]<<","
       <<grid.element_info.element_width_per_dimension[1]<<","
       <<grid.element_info.element_width_per_dimension[2]<<"]\n"
       <<prefix+"Elements per direction: ["<<grid.element_info.element_width_per_dimension[0]
       <<","
       <<grid.element_info.element_width_per_dimension[1]
       <<","
       <<grid.element_info.element_width_per_dimension[2]<<"]\n"
       <<prefix+"Volume per element: "<<grid.element_info.volume<<"\n"
       <<prefix+"Largest Element Index: "<<grid.element_info.largest_index<<"\n"
       <<prefix+"Total nodes: "<<grid.nodal_info.total_nodes<<"\n"
       <<prefix+"Nodes per direction: ("<<grid.nodal_info.nodes_per_dimension[0]<<","<<grid.nodal_info.nodes_per_dimension[1]<<","<<grid.nodal_info.nodes_per_dimension[2]<<")\n";
    out<<" \n\n";
    out<<"idx\t x\t y\t z \n";
    int total_nodes = grid.nodal_info.total_nodes;
    for(int i=0;i<total_nodes;i++){
        out<<i<<"\t"<<grid.GetNodes()[i].GetPosition()[0]<<"\t"<<grid.GetNodes()[i].GetPosition()[1]<<"\t"<<grid.GetNodes()[i].GetPosition()[2]<<"\n";
    }
    out<<"\n\n";
    out<<"idx\t 0\t 1\t 2\t 3\t 4\t 5\t 6\t 7\n";
    int total_elements = grid.element_info.total_elements;
    for(int i=0;i<total_elements;i++){
        std::array<int, 8> nodes = grid.GetElementGlobalNodeIndeces(i);
        out<<i<<"\t";
        for(int j=0;j<8;j++){
            out<<nodes[j]<<"\t";
        }
        out<<"\n";

    }

    out<<prefix+"Total subsets: "<<grid.subsets.size()<<"\n";
    if(grid.subsets.size()>0){
        for(int i=0;i<grid.subsets.size();i++){
            int number_nodes = grid.subsets[i].GetNodes().size();
            out<<prefix+"Subset "+grid.subsets[i].GetName()<<" total nodes: "<<number_nodes<<"\n";
            for(int j=0;j<number_nodes;j++){
                out<<grid.subsets[i].GetNodes()[j]<<"\t";
            }
            out<<"\n";
        }
    }

    return out;
}

void Grid::AddSubset(std::vector<int> list, std::string Name){
    this->subsets.push_back(Subset{Name, list});
}

Subset& Grid::GetSubset(std::string& Name){
    auto found_subset = std::find_if(subsets.begin(),subsets.end(), 
    [Name](Subset& ss){
        return ss.GetName()==Name;
    });

    return *found_subset;
}

std::array<double, 3>& Grid::GetUpperCorner(){
    return upper_corner;
}


std::array<double, 3>& Grid::GetLowerCorner(){
    return lower_corner;
}

std::vector<Subset>& Grid::GetSubsets(){
    return this->subsets;
}