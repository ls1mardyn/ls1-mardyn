/*
 * Created on Thu Apr 25 2024
 * Author: Jose A. Pinzon Escobar
 * Email to: jose.escobar@hsu-hh.de
 * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
 */


#pragma once

class ForceComputation{
    public:

    private:



};

// struct ForceComputation{

//     //One per cell
//     std::vector<Vec3> gradient_per_cell;
//     //Store old property for convergence check
//     std::vector<double> old_property;
//     //the current property is stored inside the sampler

//     double target_density;
//     //need boundary values for systems A--HY--B a and b
 


//     //Resize vectors, other initializations?
//     void init(GridGenerator& grid){
//         gradient_per_cell.resize(grid.GetTotalElements());
//         old_property.resize(grid.GetTotalElements());
//     }
    
//     void computeGradients(GridGenerator& grid){
//         //Central difference for first derivative
//         //use local indeces of grid to find neighbors, one data point per cell
//         //start with global index, request local, simple summation for neighbors
//         //borders need to be handled
//         //gridgenerator uses tuples

//         //get density vector
//         std::vector<double>& density = grid.GetPropertySampler().GetMaterialDensityPerCell();
//         ElementInfo& elements=grid.GetElementInfo();
//         for(int i=0;i<gradient_per_cell.size();i++){


//             for(int d=0;d<3;d++){//spatial dimensions for loop
//             //Get the two neighbors on the current spatial dimensions
//                 int backward_neighbor, forward_neighbor;
//                 backward_neighbor = elements.GetNeighbor(i,-1,d);
//                 forward_neighbor = elements.GetNeighbor(i,1,d);
//                 double central_difference = density[forward_neighbor]-density[backward_neighbor];
//                 //std::cout<<"Central difference "<<d<<" at "<<i<<" is "<<central_difference<<"\n";
//                 //std::cout<<"The forwards and backward values are: "<<density[forward_neighbor]<<" and "<<density[backward_neighbor]<<"\n";
//                 gradient_per_cell[i].at(d)=central_difference/elements.element_width_per_dimension[d];
//             }
//             //std::cout<<"The gradient at cell "<<i<<" is :["<<gradient_per_cell[i][0]<<","<<gradient_per_cell[i][1]<<","<<gradient_per_cell[i][2]<<"]\n";
//         }
//     }


//     bool convergence_check(GridGenerator& grid){
//         std::vector<double>& density = grid.GetPropertySampler().GetMaterialDensityPerCell();
//         std::vector<double> density_diff(density.size());
//         //make lambda instead
//         for(int i=0;i<density.size();i++){
//             density_diff[i]=std::abs(density[i]-target_density)/target_density;
//         }

//         auto max_diff = std::max_element(density_diff.begin(),density_diff.end());

//         if(*max_diff>0.02){
//             return false;
//         }
//         else{        
//             return true;
//         }
//     }

//     void SampleMassDensities(GridGenerator& grid, ParticleContainer* pc){
//         grid.GetPropertySampler().ComputeMaterialDensityPerCell(pc);
//     }


// };

// struct ParticleValues{
//     double density;
//     Vec3 f_x;
//     Vec3 f_inc;
// };

// struct NaiveForceComputation{
//     //Convergence parameters
//     double convergence_tolerance=0.001;//set to a small number
//     double step_size=0.001;//constant c that multiplies f_inc
//     double max_steps=300;//number of updates experienced by f_th
//     double step_stride=1;//how many steps between each update

//     std::unordered_map<unsigned long, ParticleValues> particle_values_map;//need not be sized
//     std::vector<std::array<double, 3>> particle_forces;
//     std::vector<double> particle_densities;
//     double max_density_diff=1.0;
//     double target_density;
//     std::array<double, 3> deltas = {_simulation.getcutoffRadius()/10.0,_simulation.getcutoffRadius()/10.0,_simulation.getcutoffRadius()/10.0};

//     //resizing and value setting
//     void init(double den){
//         particle_forces.resize(_simulation.getMoleculeContainer()->getNumberOfParticles(ParticleIterator::ONLY_INNER_AND_BOUNDARY));
//         particle_densities.resize(_simulation.getMoleculeContainer()->getNumberOfParticles(ParticleIterator::ONLY_INNER_AND_BOUNDARY));
//         target_density = den;

//         //fill up with zeroes
//         //std::fill(particle_forces.begin(),particle_forces.end(),0.0);
//         std::fill(particle_densities.begin(),particle_densities.end(),0.0);

//         std::cout<<"[Force Computer] The target density is: "<<target_density<<"\n"
//                  <<"[Force Computer] The force(x), rho(x) containers contain: "<<particle_forces.size()<<","<<particle_densities.size()<<" entries\n"
//                  <<"[Force Computer] Convergence properties set to: "<<"max steps="<<max_steps<<"  convergence tolerance="<<convergence_tolerance<<"   stride="<<step_stride<<"   step size c="<<step_size<<"\n";
                 
//     }

//     void SetMaxDiffDensity(double value){
//         if(value>max_density_diff){
//             max_density_diff=value;
//         }
//     }
//     //Access PropertySampler with GridGenerator
//     std::array<double, 3> computeGradients(GridGenerator& grid, ParticleContainer* pc){

//         ParticleIterator it = pc->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
//         PropertySampler& ps = grid.GetPropertySampler();

//         std::array<double, 3> part_pos;
//         std::array<double, 3> forward_pos;
//         //iterate all particles
//         for(it;it.isValid();++it){
//             part_pos = it->r_arr();
//             forward_pos[0] = forward_pos[0]+part_pos[0];
//             forward_pos[1] = forward_pos[1]+part_pos[1];
//             forward_pos[2] = forward_pos[2]+part_pos[2];



//             double rho_x = ps.ComputeMaterialDensityAtPosition(pc, part_pos);
//             max_density_diff = (std::abs(target_density-rho_x))/target_density;

//             //Compute forward difference density?
//             double rho_deltaX = ps.ComputeMaterialDensityAtPosition(pc, forward_pos);

//             //Compute forward difference
            

//         }


//     }

//     bool convergenceCheck(){

//         if(max_density_diff <= convergence_tolerance){
//             return true;
//         }
//         return false;

//     }

//     double ComputeForwardDifference1D(double fx, double fx_plus_delta, double delta_x){
//         double diff =0.0;
//         diff = (fx_plus_delta-fx)/delta_x;
//         return diff;
//     }

//     //Must go over each particle, compute density gradient and apply correction
//     void apply_FTH(GridGenerator& grid, ParticleContainer* pc){
//         ParticleIterator it = pc->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
//         PropertySampler& ps = grid.GetPropertySampler();

//         std::array<double, 3> part_pos;
//         std::array<double, 3> forward_pos;

//         //iterate all particles

//         for(it;it.isValid();++it){
//             unsigned long idx = it->getID();
//             part_pos = it->r_arr();
//             forward_pos[0] = forward_pos[0]+part_pos[0];
//             forward_pos[1] = forward_pos[1]+part_pos[1];
//             forward_pos[2] = forward_pos[2]+part_pos[2];

//             double rho_x = ps.ComputeMaterialDensityAtPosition(pc, part_pos);
//             ////store density at particle
//             particle_densities[idx] = rho_x;
//             double density_diff = (std::abs(target_density-rho_x))/target_density;

//             SetMaxDiffDensity(rho_x);
//             std::cout<<"The density is: "<<rho_x<<"\n";
//             std::cout<<"The difference in density is: "<<density_diff<<"\n";      
//             ////Compute forward difference density?
//             double rho_deltaX = ps.ComputeMaterialDensityAtPosition(pc, forward_pos);
//     //
//             ////Apply only difference on x-dimension
//             double force_value = ComputeForwardDifference1D(rho_x, rho_deltaX, deltas[0]);
//             std::cout<<"The density at x+dx is: "<<rho_deltaX<<"\n";
//             std::cout<<"The forward difference is: "<<force_value<<"\n";
//             ////Update force of particle
//             std::array<double, 3> force_at_x; 
//             force_at_x[0]=force_value;
//             force_at_x[1]=0.0;
//             force_at_x[2]=0.0;
//             particle_forces[idx][0] += force_at_x[0];
//             particle_forces[idx][1] += force_at_x[1];
//             particle_forces[idx][2] += force_at_x[2];
//             ////it->Fadd(particle_forces[idx].data());
//         }
//     }

// };