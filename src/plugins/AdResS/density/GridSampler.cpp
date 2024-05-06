#include"GridSampler.h"

GridSampler::GridSampler(){

}

void GridSampler::SetMeasureRadius(double r){
    this->measure_radius=r;
}

void GridSampler::init(Grid* grid){
    //measure_radius = _simulation.getcutoffRadius()*this->measire;
    this->grid=grid;
    int total_els = grid->GetElementInfo().total_elements;
    int total_nodes = grid->GetNodeInfo().total_nodes;
    samples.particles_per_cell.resize(total_els);
    samples.material_density.resize(total_nodes);
    samples.particles_per_node.resize(total_nodes);
    std::fill(samples.particles_per_cell.begin(), samples.particles_per_cell.end(),0);
    std::fill(samples.material_density.begin(), samples.material_density.end(),0);
    std::cout<<"Cell container size is: "<<samples.particles_per_cell.size()<<"\n";
    std::cout<<"Node container size is: "<<samples.particles_per_node.size()<<"\n";
    std::cout<<"The measure radius is: "<<measure_radius<<"\n";

    this->SetTargetValue();

}

IdxArray GridSampler::GetParticleLocalCellIndices(ParticleIterator it){
    std::array<double, 3> position = it->r_arr();
    IdxArray local_inds;

    int x, y, z;
    x = std::floor(position[0]/grid->GetElementInfo().element_width_per_dimension[0]);
    y = std::floor(position[1]/grid->GetElementInfo().element_width_per_dimension[1]);
    z = std::floor(position[2]/grid->GetElementInfo().element_width_per_dimension[2]);

    return local_inds;
}

void GridSampler::SampleAtNodes(ParticleContainer* pc){

    ParticleIterator it = pc->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
    std::fill(samples.particles_per_node.begin(), samples.particles_per_node.end(), 0);
    std::vector<Node>& all_nodes = grid->GetNodes();
    for(it;it.isValid();++it){
        std::array<double, 3> part_pos = it->r_arr();
        //Check for every node, if a particle is inside of its "bubble"
        for(int nidx=0;nidx<samples.particles_per_node.size();nidx++){
            
            auto nodal_pos=all_nodes[nidx].GetPosition();
            //check euclidean distance between particle as center of sphere/bubble
            if(ParticleInsideMeasuringSpace(nodal_pos, part_pos)){
                samples.particles_per_node[nidx] +=1;
            }
        }
    }
    //Convert to material density values
    double sphere_volume = 4.0/3.0 * M_PI * measure_radius*measure_radius*measure_radius;
    for(int nidx=0;nidx<samples.particles_per_node.size();nidx++){
        
        samples.material_density[nidx] = (double)samples.particles_per_node[nidx]/sphere_volume;
    }

}

bool GridSampler::ParticleInsideMeasuringSpace(std::array<double, 3> nodal_pos, std::array<double, 3> par_pos){
    bool is_inside=false;

    std::array<double, 3> distance;

    distance[0] = par_pos[0]-nodal_pos[0];
    distance[1] = par_pos[1]-nodal_pos[1];
    distance[2] = par_pos[2]-nodal_pos[2];

    double norm = std::sqrt(distance[0]*distance[0]+distance[1]*distance[1]+distance[2]*distance[2]);

    if(norm<measure_radius){
        is_inside=true;
    }
    return is_inside;

 }

std::ostream& GridSampler::WriteSample(std::ostream& out, std::vector<double>& smpl){
    std::string prefix = "//[Sampler]: ";
    out<<prefix+" target value: "<<this->samples.target_number_density<<"\n";
    out<<prefix+"total sampled nodes: "<<smpl.size()<<"\n";
    for(int i=0;i<smpl.size();i++){
        out<<i<<"\t"<<smpl[i]<<"\n";
    }
    return out;
}

std::ostream& GridSampler::WriteSample(std::ostream& out, std::vector<int>& smpl){
    std::string prefix = "//[Sampler]: ";
    out<<prefix+" target value: "<<this->samples.target_number_density<<"\n";
    out<<prefix+"total sampled nodes: "<<smpl.size()<<"\n";
    for(int i=0;i<smpl.size();i++){
        out<<i<<"\t"<<smpl[i]<<"\n";
    }
    return out;
}

std::ostream& GridSampler::WriteInfo(std::ostream& out){
    std::string prefix ="//[Sampler]: ";
    
}

GridHandler& GridSampler::GetGridHandler(){
    return handler;
}

void GridSampler::SetSubsetMaterialDensityValues(){
    //assign target number density to all subsets
    int nu_ss = grid->GetSubsets().size();
    for(int i=0;i<nu_ss;i++){
        auto ss=grid->GetSubsets();
        auto nodes = ss[i].GetNodes();
        int nu_nodes = ss[i].GetNodes().size();
        for(int j=0;j<nu_nodes;j++){
            this->samples.material_density[nodes[j]]=samples.target_number_density;
        }
    }

}

void GridSampler::SetTargetValue(){
    //this->samples.target_number_density=(double)_simulation.getMoleculeContainer()->getNumberOfParticles(ParticleIterator::ONLY_INNER_AND_BOUNDARY)*_simulation.getEnsemble()->getComponent(0)->m()/_simulation.getEnsemble()->domain()->V();
    this->samples.target_number_density=(double)_simulation.getMoleculeContainer()->getNumberOfParticles(ParticleIterator::ONLY_INNER_AND_BOUNDARY)/_simulation.getEnsemble()->domain()->V();
}

