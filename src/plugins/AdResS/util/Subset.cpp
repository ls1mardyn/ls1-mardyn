#include "Subset.h"

Subset::Subset(std::string Name, std::vector<int> nodes):name(Name), local_nodes(nodes){

}

void Subset::SetNodes(std::vector<int> list){
    this->local_nodes=list;
}

std::vector<int>& Subset::GetNodes(){
    return local_nodes;
}

std::string& Subset::GetName(){
    return name;
}

void Subset::SetName(std::string nme){
    this->name=nme;
}