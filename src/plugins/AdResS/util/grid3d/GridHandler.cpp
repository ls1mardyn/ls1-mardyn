#include "GridHandler.h"

std::vector<int> GridHandler::GetNodesOnPlane(Grid* grid, double p){
    std::vector<int> list_nodes;//global indeces
    //Assume n=(0,0,1)
    //Look for all points that r_z = pz
    std::vector<Node>& nodes = grid->GetNodes();
    for(int i=0;i<nodes.size();i++){
        double rz=nodes[i].GetPosition()[2];
        if(rz == p){
            list_nodes.push_back(i);
        }
    }

    return list_nodes;

}

std::vector<int> GridHandler::GetNodesOnPlanebyComponent(Grid* grid, double p, int cmp){
    std::vector<int> list_nodes;//global indeces
    //Assume n=(0,0,1)
    //Look for all points that r_z = pz
    std::vector<Node>& nodes = grid->GetNodes();
    for(int i=0;i<nodes.size();i++){
        double rz=nodes[i].GetPosition()[cmp];
        if(rz == p){
            list_nodes.push_back(i);
        }
    }

    return list_nodes;
}

void GridHandler::SetGridBoundarySubsets(Grid* grid){
    double ux, uy, uz, lx, ly, lz;
        ux = grid->GetUpperCorner()[0];
        uy = grid->GetUpperCorner()[1];
        uz = grid->GetUpperCorner()[2];

        lx = grid->GetLowerCorner()[0];
        ly = grid->GetLowerCorner()[1];
        lz = grid->GetLowerCorner()[2];
    //start at upper boundary
    //plane is (0,0,max_z)
    std::vector<int> nodes = this->GetNodesOnPlane(grid, uz);
    grid->AddSubset(nodes,"upper");

    //lower boundary 
    nodes.clear();
    nodes = this->GetNodesOnPlane(grid, lz);
    grid->AddSubset(nodes,"lower");
    
    //inlet
    nodes.clear();
    nodes = this->GetNodesOnPlanebyComponent(grid, lx, 0);
    grid->AddSubset(nodes,"inlet");

    //outlet
    nodes.clear();
    nodes = this->GetNodesOnPlanebyComponent(grid, ux, 0);
    grid->AddSubset(nodes,"outlet");

}