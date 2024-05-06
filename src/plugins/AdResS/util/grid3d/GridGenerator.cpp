/*
 * Created on Wed Feb 14 2024
 * Author: Jose A. Pinzon Escobar
 * Email to: jose.escobar@hsu-hh.de
 * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
 */
#include "GridGenerator.h"

void Grid3D::GridGenerator::SetGridGenerator(int x, int y, int z){
    elements_per_dimension[0]=x;
    elements_per_dimension[1]=y;
    elements_per_dimension[2]=z;
 }

void Grid3D::GridGenerator::init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain){
    
    MeshEntireDomain();
    Log::global_log->info()<<"Meshed region is: ("<<lower_corner[0]<<","<<lower_corner[1]<<","<<lower_corner[2]<<")x("<<upper_corner[0]<<","<<upper_corner[1]<<","<<upper_corner[2]<<")\n";
    SetTotalElements();
    SetElementInfo();
    this->node_information.GetNodes().resize(node_information.total_nodes);
    element_information.index=this->total_elements-1;
    //this->element_information.MeshTraversal();
    Log::global_log->info()<<"Total elements: "<<this->total_elements<<"\n"
                           <<"Length of domain: "<<this->element_information.length_x<<","
                           <<this->element_information.length_y<<","
                           <<this->element_information.length_z<<"\n"
                           <<"Elements per direction: ["<<this->elements_per_dimension[0]<<","<<this->elements_per_dimension[1]<<","
                           <<this->elements_per_dimension[2]<<"]\n"
                           <<"Volume per element: "<<element_information.volume<<"\n"
                           <<"Largest Element Index: "<<this->element_information.index<<std::endl;
    Log::global_log->info()<<"Node information: total "<<this->node_information.total_nodes<<"\n"
                           <<"Nodes per dimensions: ("<<this->node_information.nodes_per_dimension[0]<<","<<this->node_information.nodes_per_dimension[1]<<","<<this->node_information.nodes_per_dimension[2]<<")\n"
                           <<"Container size: "<<this->node_information.GetNodes().size()<<"\n";
    this->sampler.init(&element_information, &node_information);
    this->InitNodePositions();
    this->OutputMeshInformation();
    this->OutputNodalInformation();
    OutputElementConnectivity();
}

void Grid3D::GridGenerator::Output(){
    Log::global_log->info()<<"[AdResS] output"<<std::endl;

    std::ofstream outfile("mesh-information.txt");
    outfile.precision(6);

    std::ofstream centerfile("centers.txt");
    centerfile.precision(6);

    outfile <<"//GridGenerator plugin mesh information \n"
            <<"//Total number of elements: "<<this->total_elements<<"\n";
    
    //Table header
    outfile <<"idx \t lx \t ly \t lz \t ux \t uy \t uz \n";
    centerfile <<"idx \t cx \t cy \t cz \t rho \n";
    std::vector<double> lc(3),uc(3),ec(3);
    for(int i=0; i<= element_information.index;i++){
        lc = this->element_information.GetElementLowerCorner(i);
        uc = this->element_information.GetElementUpperCorner(i);
        ec = this->element_information.GetElementCenter(i);
        outfile << i <<"\t"<<lc[0]<<"\t"<<lc[1]<<"\t"<<lc[2]<<"\t"
                <<uc[0]<<"\t"<<uc[1]<<"\t"<<uc[2]<<"\n";
        centerfile << i <<"\t"<<ec[0]<<"\t"<<ec[1]<<"\t"<<ec[2]<<"\t"
                   <<sampler.GetParticlesPerCell()[i]<<"\n";
    }

    outfile.close();

    centerfile.close();

}

void Grid3D::GridGenerator::OutputMeshInformation(){
    std::ofstream mesh_outfile("mesh-information.txt");
    mesh_outfile.precision(6);

    mesh_outfile <<"//Grid information \n"
                 <<"//Total number of elements: "<<this->total_elements<<"\n";
    
    mesh_outfile<<"idx \t lx \t ly \t lz \t ux \t uy \t uz \n";
    std::vector<double> lc(3), uc(3), ec(3);
    for(int i=0;i<=element_information.index;i++){
        lc=this->element_information.GetElementLowerCorner(i);
        uc=this->element_information.GetElementUpperCorner(i);

        mesh_outfile << i <<"\t"<<lc[0]<<"\t"<<lc[1]<<"\t"<<lc[2]<<"\t"
                     << uc[0]<<"\t"<<uc[1]<<"\t"<<uc[2]<<"\n";
    }

    mesh_outfile.close();    
}

void Grid3D::GridGenerator::OutputNodalInformation(){
    std::ofstream nodes_file("node-positions.txt");
    nodes_file.precision(8);

    nodes_file<<"//Node Information\n"
              <<"//Total number of nodes: "<<this->node_information.GetNodes().size()<<"\n";
    nodes_file<<"idx\t x \t y \t z \n";
    DataArray position { };
    for(int idx=0;idx<node_information.GetNodes().size();idx++){
        position = node_information.GetNodes()[idx].GetPosition();
        nodes_file<<idx <<"\t"<<position[0]<<"\t"<<position[1]<<"\t"<<position[2]<<"\n";
    }
    nodes_file.close();
}

void Grid3D::GridGenerator::OutputElementConnectivity(){
    std::ofstream connectivity_file("connectivity.txt");
    connectivity_file.precision(2);
    
    std::array<int,8> idcs;
    for(int idx=0;idx<total_elements;idx++){
        idcs = GetElementGlobalNodeIndeces(idx);
        connectivity_file<<idx<<"\t"<<idcs[0]<<"\t"<<idcs[1]<<"\t"<<idcs[2]<<"\t"<<idcs[3]<<"\t"<<idcs[4]<<"\t"<<idcs[5]<<"\t"<<idcs[6]<<"\t"<<idcs[7]<<"\n";
    }
}

void Grid3D::GridGenerator::OutputPropertyPerCell(unsigned long step){
    //output measured property at center
    std::string file_name = "property_"+ std::to_string(step)+".txt";
    std::ofstream property_output(file_name);

    property_output<<"//Output particle number density\n"
                   <<"//The values correspond to the element center\n";
    property_output <<"//idx \t cx \t cy \t cz \t N_idx \n";


    std::vector<double> ec(3);
    for(int i=0; i<= element_information.index;i++){
        ec = this->element_information.GetElementCenter(i);
        property_output << i <<"\t"<<ec[0]<<"\t"<<ec[1]<<"\t"<<ec[2]<<"\t"
                   <<sampler.GetParticlesPerCell()[i]<<"\n";
    }

    property_output.close();
}

void Grid3D::GridGenerator::OutputMaterialDensityPerCell(unsigned long step){
    std::string file_name="density_"+std::to_string(step)+".txt";
    std::ofstream density_output(file_name);

    density_output<<"//Output material density per local volume\n"
                  <<"//The output values are placed at the volume center\n"
                  <<"//The mass density is : rho=n_parts*mass/local_volume\n";
    density_output<<"//idx \t cx \t cy \t cz \t rho \n";

    std::vector<double> ec(3);
    for(int i=0; i<= element_information.index;i++){
        ec = this->element_information.GetElementCenter(i);
        density_output << i <<"\t"<<ec[0]<<"\t"<<ec[1]<<"\t"<<ec[2]<<"\t"
                   <<sampler.GetMaterialDensityPerCell()[i]<<"\n";
    }

    density_output.close();


}

void Grid3D::GridGenerator::OutputNodalDensityValues(unsigned long step){
    std::string file_name="nodal_density_"+std::to_string(step)+".txt";
    std::ofstream density_output(file_name);
    int total_particles_sampled = std::reduce(sampler.GetParticlesPerNode().begin(),sampler.GetParticlesPerNode().end());
    density_output<<"//Output material density per node position\n"
                  <<"//The output values are placed at r_arr of the particle\n"
                  <<"//The mass density is : rho=n_parts*mass/local_volume\n"
                  <<"//The total number of particles sampled was: "<<total_particles_sampled<<"\n";
    density_output<<"//idx \t cx \t cy \t cz \t rho \n";

    DataArray ec { };
    for(int i=0; i< node_information.GetNodes().size();i++){
        ec = this->node_information.GetNodes()[i].GetPosition();
        density_output << i <<"\t"<<ec[0]<<"\t"<<ec[1]<<"\t"<<ec[2]<<"\t"//<<"\n";
                   <<sampler.GetParticlesPerNode()[i]<<"\n";
    }

    density_output.close();
}

void Grid3D::GridGenerator::MeshEntireDomain(){
    for(int d=0;d<3;d++){
        this->lower_corner[d]=0;
        this->upper_corner[d]=_simulation.getDomain()->getGlobalLength(d);
    }
}

void Grid3D::GridGenerator::SetElementInfo(){
    //Compute widths
    for(int d=0;d<3;d++){
        element_width_per_dimension[d]= (upper_corner[d]-lower_corner[d])/elements_per_dimension[d];
    }
    this->element_information.x_width=element_width_per_dimension[0];
    this->element_information.y_width=element_width_per_dimension[1];
    this->element_information.z_width=element_width_per_dimension[2];

    this->element_information.volume=element_width_per_dimension[0]*element_width_per_dimension[1]*element_width_per_dimension[2];

    this->element_information.length_x=(upper_corner[0]-lower_corner[0]);
    this->element_information.length_y=(upper_corner[1]-lower_corner[1]);
    this->element_information.length_z=(upper_corner[2]-lower_corner[2]);

    this->element_information.elements_in_x=elements_per_dimension[0];
    this->element_information.elements_in_y=elements_per_dimension[1];
    this->element_information.elements_in_z=elements_per_dimension[2];
    this->element_information.elements_per_dimension=this->elements_per_dimension;
    this->element_information.element_width_per_dimension=this->element_width_per_dimension;
    //After setting element info, set node information
    this->SetNodeInfo();

}

void Grid3D::GridGenerator::SetNodeInfo(){
    node_information.total_nodes=(elements_per_dimension[0]+1)*(elements_per_dimension[1]+1)*(elements_per_dimension[2]+1);
    node_information.nodes_per_dimension[0]=elements_per_dimension[0]+1;
    node_information.nodes_per_dimension[1]=elements_per_dimension[1]+1;
    node_information.nodes_per_dimension[2]=elements_per_dimension[2]+1;
}

void Grid3D::GridGenerator::SetTotalElements(){
    this->total_elements=1;
    for(int d=0;d<3;d++){
        this->total_elements *= elements_per_dimension[d];
    }
}

void Grid3D::GridGenerator::InitNodePositions(){
    for(int idx=0;idx<node_information.GetNodes().size();idx++){
        
        std::tuple<int,int,int> local_indeces = this->node_information.MapGlobalToLocal(idx);
        double x = (double)std::get<0>(local_indeces)*element_width_per_dimension[0];
        double y = (double)std::get<1>(local_indeces)*element_width_per_dimension[1];
        double z = (double)std::get<2>(local_indeces)*element_width_per_dimension[2];
        node_information.GetNodes()[idx].SetPosition(x,y,z);
    }
}

Grid3D::PropertySampler& Grid3D::GridGenerator::GetPropertySampler(){
    return this->sampler;
}

std::array<int, 8> Grid3D::GridGenerator::GetElementGlobalNodeIndeces(int el){
    std::array<int, 8> element_nodes_idx_global;
    //Get local index of element
    std::tuple<int,int,int> element_local_idx = element_information.GlobalToLocalIndex(el);
    //std::cout<<"Index is: "<<element_local_idx
    //generate a list of neighbors
    //Navigate using element corner closer to origin, loop specific to cube shaped elements
    int counter=0;
    for(int z=0;z<=1;z++){
        for(int y=0;y<=1;y++){
            for(int x=0;x<=1;x++){
                element_nodes_idx_global[counter]=node_information.MapLocalToGlobalIndex(std::get<0>(element_local_idx)+x,std::get<1>(element_local_idx)+y,std::get<2>(element_local_idx)+z);
                counter++;
            }
        }
    }

    return element_nodes_idx_global;
}

void Grid3D::PropertySampler::init(ElementInfo* info, NodeInformation* node_info){
    measure_radius = _simulation.getcutoffRadius()*4.0;
    this->info = info;
    this->node_info=node_info;
    int total_els = info->elements_in_x*info->elements_in_y*info->elements_in_z;
    particles_per_cell.resize(total_els);
    material_density.resize(total_els);
    particles_per_node.resize(node_info->total_nodes);
    std::fill(particles_per_cell.begin(), particles_per_cell.end(),0);
    std::fill(material_density.begin(), material_density.end(),0);
	Log::global_log->info() <<"Cell container size is: " << particles_per_cell.size() << "\n";
	Log::global_log->info() <<"Node container size is: " << particles_per_node.size() << "\n";
	Log::global_log->info() <<"The measure radius is: " << measure_radius << "\n";

}

void Grid3D::PropertySampler::ParticlePerCellCount(ParticleContainer* particle_container){
    //Iterate all non-halo cells
    ParticleIterator it = particle_container->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
    std::fill(particles_per_cell.begin(), particles_per_cell.end(), 0);//last element not inclusive
    for(it;it.isValid();++it){

        std::tuple<int, int, int> indeces= GetParticleLocalCellIndices(it);
        int particle_global_index = this->info->LocalToGlobalIndex(indeces);
        particles_per_cell[particle_global_index]= particles_per_cell[particle_global_index]+1;
    }
}

void Grid3D::PropertySampler::ComputeMaterialDensityPerCell(ParticleContainer* particle_container){
    //Iterate all non-halo cells
    ParticleIterator it = particle_container->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
    std::fill(material_density.begin(), material_density.end(), 0);//last element not inclusive
    for(it;it.isValid();++it){
        double molecule_mass = (*it).mass();
        std::tuple<int, int, int> indeces= GetParticleLocalCellIndices(it);
        int particle_global_index = this->info->LocalToGlobalIndex(indeces);
        material_density[particle_global_index]= material_density[particle_global_index]+molecule_mass;
    }
    for(int i=0;i<material_density.size();++i){
        material_density[i]= material_density[i]/info->volume;
    }
}

double Grid3D::PropertySampler::ComputeMaterialDensityAtPosition(ParticleContainer* pc, std::array<double, 3>& pos){
    ParticleIterator it = pc->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
    int particles_inside_sphere =0;

    DataArray pos_as_vector { };
    pos_as_vector[0]=pos[0];
    pos_as_vector[1]=pos[1];
    pos_as_vector[2]=pos[2];

    for(it;it.isValid();++it){
        std::array<double, 3> particle_position = it->r_arr();
        //std::cout<<"Particle "<<it->getID()<<" has mass of "<<it->mass()<<"\n";
        //Assume current position is node (virtual, since it does not exist)
        if(ParticleInsideMeasuringSpace(pos_as_vector,particle_position)){
            particles_inside_sphere +=it->mass();
        }
    }//end for
    double volume = (4.0/3.0)*M_PI*measure_radius*measure_radius*measure_radius;
    double density = (double)particles_inside_sphere/volume;

    return density;
}

void Grid3D::PropertySampler::SampleAtNodes(ParticleContainer* pc){
    ParticleIterator it = pc->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
    std::fill(particles_per_node.begin(), particles_per_node.end(), 0);
    std::vector<Node>& all_nodes = node_info->GetNodes();
    for(it;it.isValid();++it){
        std::array<double, 3> part_pos = it->r_arr();
        //Check for every node, if a particle is inside of its "bubble"
        for(int nidx=0;nidx<particles_per_node.size();nidx++){
            
            auto nodal_pos=all_nodes[nidx].GetPosition();
            //check euclidean distance between particle as center of sphere/bubble
            if(ParticleInsideMeasuringSpace(nodal_pos, part_pos)){
                particles_per_node[nidx] +=1;
            }
        }
    }
}

bool Grid3D::PropertySampler::ParticleInsideMeasuringSpace(std::array<double, 3> nodal_pos, std::array<double, 3> par_pos){
    bool is_inside=false;

    DataArray distance { };
    distance[0] = par_pos[0]-nodal_pos[0];
    distance[1] = par_pos[1]-nodal_pos[1];
    distance[2] = par_pos[2]-nodal_pos[2];

    double norm = std::sqrt(distance[0]*distance[0]+distance[1]*distance[1]+distance[2]*distance[2]);

    if(norm<measure_radius){
        is_inside=true;
    }
    return is_inside;
 }

std::tuple<int, int, int> Grid3D::PropertySampler::GetParticleLocalCellIndices(ParticleIterator it)
{
    std::array<double, 3> position = it->r_arr();
    std::tuple<int, int, int> local_inds;

    int x, y, z;
    x = std::floor(position[0]/info->x_width);
    y = std::floor(position[1]/info->y_width);
    z = std::floor(position[2]/info->z_width);
    local_inds= std::make_tuple(x,y,z);

    return local_inds;
}

std::vector<int>& Grid3D::PropertySampler::GetParticlesPerCell(){
    return this->particles_per_cell;
}

std::vector<double>& Grid3D::PropertySampler::GetMaterialDensityPerCell(){
    return this->material_density;
}

std::vector<int>& Grid3D::PropertySampler::GetParticlesPerNode(){
    return this->particles_per_node;
}
