/*
 * Created on Wed Feb 14 2024
 * Author: Jose A. Pinzon Escobar
 * Email to: jose.escobar@hsu-hh.de
 * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
 */


#include "GridGenerator.h"

std::string GridGenerator::getPluginName(){
    return std::string{"[GridGenerator]"};
}

void GridGenerator::readXML(XMLfileUnits& xmlconfig){
    Log::global_log->info()<<"[GridGenerator] enabled"<<std::endl;

    xmlconfig.getNodeValue("elementsX",elements_per_dimension[0]);
    xmlconfig.getNodeValue("elementsY",elements_per_dimension[1]);
    xmlconfig.getNodeValue("elementsZ",elements_per_dimension[2]);

    //SET LOWER AND UPPER CORNERS OR USE DEFAULT
}

void GridGenerator::init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain){
    
    
    MeshEntireDomain();
    SetTotalElements();
    SetElementInfo();
    element_information.index=this->total_elements-1;
    
}

void GridGenerator::beforeEventNewTimestep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep){
    Log::global_log->info()<<this->getPluginName()<<" information\n";
    Log::global_log->info()<<"Total elements: "<<this->total_elements<<"\n"
                           <<"Elemens per direction: ["<<this->elements_per_dimension[0]<<","<<this->elements_per_dimension[1]<<","
                           <<this->elements_per_dimension[2]<<"]"<<std::endl;
                           
}

void GridGenerator::MeshEntireDomain(){
    for(int d=0;d<3;d++){
        this->lower_corner[d]=0;
        this->upper_corner[d]=_simulation.getDomain()->getGlobalLength(d);
    }
}

void GridGenerator::SetElementInfo(){
    //Compute widths
    for(int d=0;d<3;d++){
        element_width_per_dimension[d]= (upper_corner[d]-lower_corner[d])/elements_per_dimension[d];
    }

    this->element_information.volume=element_width_per_dimension[0]*element_width_per_dimension[1]*element_width_per_dimension[2];

}

void GridGenerator::SetTotalElements(){
    for(int d=0;d<3;d++){
        this->total_elements += elements_per_dimension[d];
    }
}