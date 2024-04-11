// /*
//  * Created on Wed Feb 14 2024
//  * Author: Jose A. Pinzon Escobar
//  * Email to: jose.escobar@hsu-hh.de
//  * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
//  */


// #include "GridGenerator.h"

// GridGenerator::GridGenerator(){

// }

// void GridGenerator::readXML(XMLfileUnits& xmlconfig){
//     //Log::global_log->info()<<"[GridGenerator] enabled"<<std::endl;

//     xmlconfig.getNodeValue("elementsX",elements_per_dimension[0]);
//     xmlconfig.getNodeValue("elementsY",elements_per_dimension[1]);
//     xmlconfig.getNodeValue("elementsZ",elements_per_dimension[2]);

//     //SET LOWER AND UPPER CORNERS OR USE DEFAULT
// }

// void GridGenerator::init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain){
    
//     MeshEntireDomain();
//     Log::global_log->info()<<"["<<this->getPluginName()<<"]\n";
//     Log::global_log->info()<<"Meshed region is: ("<<lower_corner[0]<<","<<lower_corner[1]<<","<<lower_corner[2]<<")x("<<upper_corner[0]<<","<<upper_corner[1]<<","<<upper_corner[2]<<")\n";
//     SetTotalElements();
//     SetElementInfo();
//     element_information.index=this->total_elements-1;
//     //this->element_information.MeshTraversal();
//     Log::global_log->info()<<"Total elements: "<<this->total_elements<<"\n"
//                            <<"Length of domain: "<<this->element_information.length_x<<","
//                            <<this->element_information.length_y<<","
//                            <<this->element_information.length_z<<"\n"
//                            <<"Elements per direction: ["<<this->elements_per_dimension[0]<<","<<this->elements_per_dimension[1]<<","
//                            <<this->elements_per_dimension[2]<<"]\n"
//                            <<"Volume per element: "<<element_information.volume
//                            <<"Largest Index: "<<this->element_information.index<<std::endl;
//     this->sampler.init(&element_information);
//     this->Output("not-in-use",1);
// }

// void GridGenerator::beforeEventNewTimestep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep){

//     sampler.ParticlePerCellCount(particleContainer);   
//     this->Output("not-in-use",1);                       
// }

// void GridGenerator::Output(std::string prefix, long unsigned accumulated){
//     Log::global_log->info()<<"[GridGenerator] output"<<std::endl;

//     std::ofstream outfile("mesh-information.txt");
//     outfile.precision(6);

//     std::ofstream centerfile("centers.txt");
//     centerfile.precision(6);

//     outfile <<"//GridGenerator plugin mesh information \n"
//             <<"//Total number of elements: "<<this->total_elements<<"\n";
    
//     //Table header
//     outfile <<"idx \t lx \t ly \t lz \t ux \t uy \t uz \n";
//     centerfile <<"idx \t cx \t cy \t cz \t rho \n";
//     std::vector<double> lc(3),uc(3),ec(3);
//     for(int i=0; i<= element_information.index;i++){
//         lc = this->element_information.GetElementLowerCorner(i);
//         uc = this->element_information.GetElementUpperCorner(i);
//         ec = this->element_information.GetElementCenter(i);
//         outfile << i <<"\t"<<lc[0]<<"\t"<<lc[1]<<"\t"<<lc[2]<<"\t"
//                 <<uc[0]<<"\t"<<uc[1]<<"\t"<<uc[2]<<"\n";
//         centerfile << i <<"\t"<<ec[0]<<"\t"<<ec[1]<<"\t"<<ec[2]<<"\t"
//                    <<sampler.GetParticlesPerCell()[i]<<"\n";
//     }

//     outfile.close();

//     centerfile.close();

// }

// void GridGenerator::MeshEntireDomain(){
//     for(int d=0;d<3;d++){
//         this->lower_corner[d]=0;
//         this->upper_corner[d]=_simulation.getDomain()->getGlobalLength(d);
//     }
// }

// void GridGenerator::SetElementInfo(){
//     //Compute widths
//     for(int d=0;d<3;d++){
//         element_width_per_dimension[d]= (upper_corner[d]-lower_corner[d])/elements_per_dimension[d];
//     }
//     this->element_information.x_width=element_width_per_dimension[0];
//     this->element_information.y_width=element_width_per_dimension[1];
//     this->element_information.z_width=element_width_per_dimension[2];

//     this->element_information.volume=element_width_per_dimension[0]*element_width_per_dimension[1]*element_width_per_dimension[2];

//     this->element_information.length_x=(upper_corner[0]-lower_corner[0]);
//     this->element_information.length_y=(upper_corner[1]-lower_corner[1]);
//     this->element_information.length_z=(upper_corner[2]-lower_corner[2]);

//     this->element_information.elements_in_x=elements_per_dimension[0];
//     this->element_information.elements_in_y=elements_per_dimension[1];
//     this->element_information.elements_in_z=elements_per_dimension[2];

// }

// void GridGenerator::SetTotalElements(){
//     this->total_elements=1;
//     for(int d=0;d<3;d++){
//         this->total_elements *= elements_per_dimension[d];
//     }
// }

// PropertySampler::PropertySampler(){
// }

// void PropertySampler::init(ElementInfo* info){
//     this->info = info;
//     particles_per_cell.resize(this->info->index);
//     std::fill(particles_per_cell.begin(), particles_per_cell.end(),0);
//     std::cout<<"\nContainer size is: "<<particles_per_cell.size()<<"\n";

// }

// void PropertySampler::ParticlePerCellCount(ParticleContainer* particle_container){
//     //Iterate all non-halo cells
//     ParticleIterator it = particle_container->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
//     std::fill(particles_per_cell.begin(), particles_per_cell.end()+1, 0);//last element not inclusive
//     for(it;it.isValid();++it){

//         std::tuple<int, int, int> indeces= GetParticleLocalCellIndices(it);
//         int particle_global_index = this->info->LocalToGlobalIndex(indeces);
//         particles_per_cell[particle_global_index]= particles_per_cell[particle_global_index]+1;
//     }
// } 

// std::tuple<int, int, int> PropertySampler::GetParticleLocalCellIndices(ParticleIterator it)
// {
    
//     std::array<double, 3> position = it->r_arr();
//     std::tuple<int, int, int> local_inds;

//     int x, y, z;

//     x = std::floor(position[0]/info->x_width);
//     y = std::floor(position[1]/info->y_width);
//     z = std::floor(position[0]/info->z_width);

//     local_inds= std::make_tuple(x,y,z);

//     return local_inds;


// }

// std::vector<int>& PropertySampler::GetParticlesPerCell(){
//     return this->particles_per_cell;
// }