#include"GridSampler.h"

void Grid3D::GridSampler::SetMeasureRadius(double r){
    this->measure_radius=r;
}

void Grid3D::GridSampler::init(Grid* grid){
    //measure_radius = _simulation.getcutoffRadius()*this->measire;
    this->grid=grid;
    int total_els = grid->GetElementInfo().total_elements;
    int total_nodes = grid->GetNodeInfo().total_nodes;
    samples.particles_per_cell.resize(total_els);
    samples.material_density.resize(total_nodes);
    samples.particles_per_node.resize(total_nodes);
    std::fill(samples.particles_per_cell.begin(), samples.particles_per_cell.end(),0);
    std::fill(samples.material_density.begin(), samples.material_density.end(),0);
    Log::global_log->info() << "Cell container size is: " << samples.particles_per_cell.size() << "\n";
	Log::global_log->info() << "Node container size is: " << samples.particles_per_node.size() << "\n";
	Log::global_log->info() << "The measure radius is: " << measure_radius << "\n";
    this->SetTargetValue();
}

Grid3D::i3 Grid3D::GridSampler::GetParticleLocalCellIndices(ParticleIterator it){
    std::array<double, 3> position = it->r_arr();
    i3 local_inds;

    int x, y, z;
    x = std::floor(position[0]/grid->GetElementInfo().element_width_per_dimension[0]);
    y = std::floor(position[1]/grid->GetElementInfo().element_width_per_dimension[1]);
    z = std::floor(position[2]/grid->GetElementInfo().element_width_per_dimension[2]);

    return local_inds;
}

void Grid3D::GridSampler::SampleAtNodes(const std::array<std::vector<double>,3>& positions){
    std::fill(samples.particles_per_node.begin(), samples.particles_per_node.end(), 0);
    std::vector<Node>& all_nodes = grid->GetNodes();

	// put mol pos into grid of spacing 2x measure_radius
	const auto box_dim = 2 * measure_radius;
	const std::array<int, 3> bins_per_dim = { static_cast<int>(std::ceil(_simulation.getDomain()->getGlobalLength(0) / box_dim)),
											  static_cast<int>(std::ceil(_simulation.getDomain()->getGlobalLength(1) / box_dim)),
											  static_cast<int>(std::ceil(_simulation.getDomain()->getGlobalLength(2) / box_dim))};
	const auto bins = bins_per_dim[0]*bins_per_dim[1]*bins_per_dim[2];
	std::vector<std::vector<std::array<double, 3>>> mol_pos_binned;
	mol_pos_binned.resize(bins);
	for(std::size_t idx = 0; idx < positions[0].size(); ++idx) {
		std::array<double, 3> part_pos {positions[0][idx], positions[1][idx], positions[2][idx]};
		int bin_x = std::min(static_cast<int>(part_pos[0] / box_dim), bins_per_dim[0]-1);
		int bin_y = std::min(static_cast<int>(part_pos[1] / box_dim), bins_per_dim[1]-1);
		int bin_z = std::min(static_cast<int>(part_pos[2] / box_dim), bins_per_dim[2]-1);

		mol_pos_binned[bins_per_dim[0]*bins_per_dim[1]*bin_z + bins_per_dim[0]*bin_y + bin_x].push_back(part_pos);
	}

	#if defined(_OPENMP)
	#pragma omp parallel for
	#endif
	for(int nidx = 0; nidx < samples.particles_per_node.size(); nidx++) {
		auto nodal_pos=all_nodes[nidx].GetPosition();
		int bin_x = static_cast<int>(nodal_pos[0] / box_dim);
		int bin_y = static_cast<int>(nodal_pos[1] / box_dim);
		int bin_z = static_cast<int>(nodal_pos[2] / box_dim);

		for(int offset_z = -1; offset_z <= 1; offset_z++) {
			int bz = offset_z + bin_z;
			if(bz < 0 || bz >= bins_per_dim[2]) continue;
			for(int offset_y = -1; offset_y <= 1; offset_y++) {
				int by = offset_y + bin_y;
				if(by < 0 || by >= bins_per_dim[1]) continue;
				for(int offset_x = -1; offset_x <= 1; offset_x++) {
					int bx = offset_x + bin_x;
					if(bx < 0 || bx >= bins_per_dim[0]) continue;

					for(auto& mol_pos : mol_pos_binned[bins_per_dim[0]*bins_per_dim[1]*bz + bins_per_dim[0]*by + bx]) {
						if(ParticleInsideMeasuringSpace(nodal_pos, mol_pos)){
							samples.particles_per_node[nidx] +=1;
						}
					}
				}
			}
		}
    }
    //Convert to material density values
    double sphere_volume = 4.0/3.0 * M_PI * measure_radius*measure_radius*measure_radius;
	#if defined(_OPENMP)
	#pragma omp parallel for
	#endif
    for(int nidx = 0; nidx < samples.particles_per_node.size(); nidx++){
        
        samples.material_density[nidx] = (double)samples.particles_per_node[nidx]/sphere_volume;
    }
}

bool Grid3D::GridSampler::ParticleInsideMeasuringSpace(std::array<double, 3> nodal_pos, std::array<double, 3> par_pos){
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

std::ostream& Grid3D::GridSampler::WriteSample(std::ostream& out, std::vector<double>& smpl){
    std::string prefix = "//[Sampler]: ";
    out<<prefix+" target value: "<<this->samples.target_number_density<<"\n";
    out<<prefix+"total sampled nodes: "<<smpl.size()<<"\n";
    for(int i=0;i<smpl.size();i++){
        out<<i<<"\t"<<smpl[i]<<"\n";
    }
    return out;
}

std::ostream& Grid3D::GridSampler::WriteSample(std::ostream& out, std::vector<int>& smpl){
    std::string prefix = "//[Sampler]: ";
    out<<prefix+" target value: "<<this->samples.target_number_density<<"\n";
    out<<prefix+"total sampled nodes: "<<smpl.size()<<"\n";
    for(int i=0;i<smpl.size();i++){
        out<<i<<"\t"<<smpl[i]<<"\n";
    }
    return out;
}

std::ostream& Grid3D::GridSampler::WriteInfo(std::ostream& out){
    std::string prefix ="//[Sampler]: ";
    return out;
}

Grid3D::GridHandler& Grid3D::GridSampler::GetGridHandler(){
    return handler;
}

void Grid3D::GridSampler::SetSubsetMaterialDensityValues(){
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

void Grid3D::GridSampler::SetTargetValue(){
    this->samples.target_number_density=(double)_simulation.getDomain()->getglobalNumMolecules()
			/ _simulation.getEnsemble()->domain()->V();
}
