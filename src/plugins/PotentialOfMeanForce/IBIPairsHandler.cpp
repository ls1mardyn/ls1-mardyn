#include "IBIPairsHandler.h"
#include "molecules/potforce.h"
#include "IBI.h"

IBIPairsHandler::IBIPairsHandler() {
    const int number_threads = mardyn_get_max_threads();
    Log::global_log->info()<<"[InteractionForceAdapter]: allocate data for "<<number_threads<<" threads."<<std::endl;

    thread_data.resize(number_threads);
    #if defined(_OPENMP)
    #pragma omp parallel
    #endif
    {
        auto* own_data = new ParticlePairs2PotForceAdapter::PP2PFAThreadData();
        const int own_id = mardyn_get_thread_num();
        thread_data[own_id] = own_data;
    }

}

IBIPairsHandler::~IBIPairsHandler() {
    for (int idx = 0; idx < thread_data.size(); idx++) delete thread_data[idx];
}

void IBIPairsHandler::init(){
    Domain* domain = _simulation.getDomain();
	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{
		const int own_id = mardyn_get_thread_num();
		thread_data[own_id]->initComp2Param(domain->getComp2Params());
		thread_data[own_id]->clear();
	}
}

void IBIPairsHandler::finish(){
    double virial = 0, upot6LJ = 0, upotXpoles = 0, myRF = 0;
    for (int idx = 0; idx < thread_data.size(); idx++) {
        virial += thread_data[idx]->_virial;
        upot6LJ += thread_data[idx]->_upot6LJ;
        upotXpoles += thread_data[idx]->_upotXpoles;
        myRF += thread_data[idx]->_myRF;
    }

    Domain* domain = _simulation.getDomain();
    domain->setLocalUpot(upot6LJ / 6. + upotXpoles + myRF + domain->getLocalUpot());
    domain->setLocalVirial(virial + 3.0 * myRF + domain->getLocalVirial());
}

double IBIPairsHandler::processPair(Molecule& m1, Molecule& m2, double distance[3], PairType pair, double dd, bool CalculateLJ){
    const int tid = mardyn_get_thread_num();
    ParticlePairs2PotForceAdapter::PP2PFAThreadData& data = *thread_data[tid];

    switch(pair){
        double Virial[3];
        double dummy1,dummy2,dummy3,dummy4[3];
        case MOLECULE_MOLECULE:
            PotForceOnlyCG(m1,m2,distance,data._upot6LJ,data._upotXpoles,data._myRF,Virial,CalculateLJ);
            return data._upot6LJ+data._upotXpoles;

        case MOLECULE_HALOMOLECULE:
            PotForceOnlyCG(m1,m2,distance,dummy1,dummy2,dummy3,dummy4,CalculateLJ);
            return 0.0;

        case MOLECULE_MOLECULE_FLUID:
            Log::global_log->info()<<"Not implemented for cg"<<std::endl;
            Simulation::exit(670);
        default:
            Simulation::exit(670);
    }
    return 0.0;
}

void IBIPairsHandler::PotForceOnlyCG(Molecule& m1, Molecule& m2, double* distance, double& Upot6LJ, double& UpotXPoles, double& MyRF, double virial[3], bool calcLJ) {
	virial[0] = 0.0;
	virial[1] = 0.0;
	virial[2] = 0.0;

    //for IBI molecule component must have a single site
    std::array<double, 3> com1 = m1.r_arr();
    std::array<double, 3> com2 = m2.r_arr();

    const std::array<double, 3> dr = {com2[0] - com1[0], com2[1] - com1[1], com2[2] - com1[2]};
    const double r = Distance(com1, com2);
    const double Upot = potential_function.EvaluateAt(r);
    std::array<double,3> force = getActingForce(dr, r);

    m1.Fljcenteradd(0, force.data());
    m2.Fljcentersub(0, force.data());
    Upot6LJ += Upot;

    for (int d = 0; d < 3; ++d) {
        virial[d] += 0.5 * dr[d] * force[d];
    }
}

std::array<double, 3> IBIPairsHandler::getActingForce(const std::array<double, 3> &dr, double r) {
    const auto F = force_function.EvaluateAt(r);

    std::array<double, 3> force = dr;
    Normalize(force);
    for (int dim = 0; dim < 3; dim++) force[dim] *= F;

    return force;
}
