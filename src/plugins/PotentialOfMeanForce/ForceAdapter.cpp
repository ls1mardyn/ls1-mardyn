#include "ForceAdapter.h"
#include "molecules/potforce.h"
#include "PMF.h"

InteractionForceAdapter::InteractionForceAdapter(ResolutionHandler& handle,  PMF* pmf):resolution_handler{handle},adres{pmf}{

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

void InteractionForceAdapter::init(){

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

void InteractionForceAdapter::finish(){
 //What to put here?
}

double InteractionForceAdapter::processPair(Molecule& m1, Molecule& m2, double distance[3], PairType pair, double dd, bool CalculateLJ){
    std::vector<FPRegion>& regions = adres->GetRegions();

    std::array<double, 3> com1 = CenterOfMass(m1);
    std::array<double, 3> com2 = CenterOfMass(m2);

    bool has_hybrid = false;
    InteractionType interaction;
    // if(component_handler.GetMoleculeResolution(m1)==ResolutionType::Hybrid || component_handler.GetMoleculeResolution(m2)==ResolutionType::Hybrid)
    if(resolution_handler.GetCOMResolution(com1,regions)==Hybrid || resolution_handler.GetCOMResolution(com2,regions)==Hybrid)
    {
        has_hybrid = true;
        interaction = InteractionType::mixed;
    }

    // if(component_handler.GetMoleculeResolution(m1)==ResolutionType::CoarseGrain && component_handler.GetMoleculeResolution(m2)==CoarseGrain)
    if(resolution_handler.GetCOMResolution(com1,regions)==CoarseGrain && resolution_handler.GetCOMResolution(com2,regions)==CoarseGrain)
    {
        interaction = InteractionType::onlycg;
    }

    // if(component_handler.GetMoleculeResolution(m1)==ResolutionType::FullParticle && component_handler.GetMoleculeResolution(m2) == FullParticle)
    if(resolution_handler.GetCOMResolution(com1,regions)==FullParticle && resolution_handler.GetCOMResolution(com2,regions)==FullParticle)
    {
        interaction = InteractionType::onlyfp;
    }
    
    return processPairBackend(m1, m2, distance, pair, dd, CalculateLJ, interaction);

}

double InteractionForceAdapter::processPairBackend(Molecule& m1, Molecule& m2, double distance[3], PairType pair, double dd, bool calcLJ, InteractionType interaction){
    const int tid = mardyn_get_thread_num();
    ParticlePairs2PotForceAdapter::PP2PFAThreadData& data = *thread_data[tid];
    ParaStrm params = (*data._comp2Param)(m1.componentid(),m2.componentid());
    ParaStrm paramsInv = (*data._comp2Param)(m2.componentid(),m1.componentid());

    switch(pair){
        double Virial[3];
        double dummy1,dummy2,dummy3,dummy4[3];
        case MOLECULE_MOLECULE:
            this->PotForceType(m1,m2,params,paramsInv,distance,data._upot6LJ,data._upotXpoles,data._myRF,Virial,calcLJ,interaction);
			return data._upot6LJ+data._upotXpoles;

        case MOLECULE_HALOMOLECULE:
            this->PotForceType(m1,m2,params,paramsInv,distance,dummy1,dummy2,dummy3,dummy4,calcLJ,interaction);

            return 0.0;

        case MOLECULE_MOLECULE_FLUID:
            this->FluidPotType(m1,m2,params,paramsInv,distance,data._upot6LJ,data._upotXpoles,data._myRF,Virial,calcLJ,interaction);
            return dummy1 / 6.0 + dummy2 + dummy3;
		default:
		Simulation::exit(670);
    }
    return 0.0;
}

void InteractionForceAdapter::PotForceType(Molecule& m1, Molecule& m2, ParaStrm& params, ParaStrm& paramInv, double* drm, double& Upot6LJ, double& UpotXpoles, double& MyRF, double Virial[3], bool calcLJ, InteractionType interaction){
    //Pure interaction case
    if(interaction == InteractionType::onlyfp){
        // Log::global_log->error()<<"PotForceType cannot be FP"<<std::endl;
        PotForce(m1,m2,params,drm,Upot6LJ, UpotXpoles, MyRF, Virial,calcLJ);
    }

    if(interaction == InteractionType::onlycg){
         Log::global_log->error()<<"PotForceType cannot be CG"<<std::endl;

        PotForceOnlyCG(m1,m2,params,drm,Upot6LJ,UpotXpoles,MyRF,Virial,calcLJ);
    }

    if(interaction == InteractionType::mixed){
         Log::global_log->error()<<"PotForceType cannot be Mixed"<<std::endl;
        /**
         * m1 hy vs m2 cg
         * m1 hy vs m2 at
         * viceversa
		 * m1 hy vs m2 hy
         */
        PotForceHybrid(m1,m2,params,drm,Upot6LJ,UpotXpoles,MyRF,Virial,calcLJ);
        //PotForce(m1,m2,params,drm,Upot6LJ, UpotXpoles, MyRF, Virial,calcLJ);
    }
}

void InteractionForceAdapter::FluidPotType(Molecule& m1, Molecule& m2, ParaStrm& params, ParaStrm& paramInv, double* drm, double& Upot6LJ, double& UpotXpoles, double& MyRF, double Virial[3], bool calcLJ, InteractionType interaction){
    //Pure interaction case
    if(interaction == InteractionType::onlyfp){
        FluidPot(m1,m2,params,drm,Upot6LJ, UpotXpoles, MyRF,calcLJ);
    }

    if(interaction == InteractionType::onlycg){
        Log::global_log->error()<<"Not implemented for cg"<<std::endl;
    }

    if(interaction == InteractionType::mixed){
        Log::global_log->error()<<"Not implemented for mixed"<<std::endl;
    }
}


void InteractionForceAdapter::PotForceHybrid(Molecule& m1, Molecule& m2, ParaStrm& params, double* distance, double& Upot6LJ, double& UpotXPoles, double& MyRF, double virial[3], bool calcLJ){
	//determine which one is hybrid and which fpregion is being used

	//Assume only one fpregion
	FPRegion& region = adres->GetRegions()[0];
	PotForceHybridBackend(m1,m2,params,distance,Upot6LJ,UpotXPoles,MyRF,virial,calcLJ,region);

}

void InteractionForceAdapter::PotForceHybridBackend(Molecule& m1, Molecule& m2, ParaStrm& params, double* distance, double& Upot6LJ, double& UpotXPoles, double& MyRF, double virial[3], bool calcLJ, FPRegion& region){

    std::array<double, 3> com1 = CenterOfMass(m1);
    std::array<double, 3> com2 = CenterOfMass(m2);

    //assume both m1 and m2 are hybrid molecules

    //Do LJ sites first for the F_{at}
    // Force Calculation
	double f[3];double f2[3];
	double u;
	double drs[3], dr2; // site distance vector & length^2
	virial[0]=0.;
	virial[1]=0.;
	virial[2]=0.;

    double w1, w2;
    w1 = adres->WeightValue(m1.r_arr(),region);
    w2 = adres->WeightValue(m2.r_arr(),region);
    //All these are F_{at}
    // LJ centers
	// no LJ interaction between solid atoms of the same component
	const unsigned int nc1 = m1.numLJcenters();
	const unsigned int nc2 = m2.numLJcenters();
	for (unsigned int si = 0; si < nc1; ++si) {
		const std::array<double,3> dii = m1.ljcenter_d_abs(si);
		for (unsigned int sj = 0; sj < nc2; ++sj) {
			const std::array<double,3> djj = m2.ljcenter_d_abs(sj);
			SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
			double eps24;
			params >> eps24;
			double sig2;
			params >> sig2;
			double shift6;
			params >> shift6; // must be 0.0 for full LJ
			if (calcLJ) {
				PotForceLJ(drs, dr2, eps24, sig2, f, u);
				u += shift6;

                //Introduce weight function values
                for(int i=0;i<3;i++){
                    f[i] *= w1*w2;//F^{at}
                }

				m1.Fljcenteradd(si, f);
				m2.Fljcentersub(sj, f);
				Upot6LJ += u;
				for (unsigned short d = 0; d < 3; ++d)
					virial[d] += 0.5*distance[d] * f[d];
			}
		}
	}



    m1.Viadd(virial);
	m2.Viadd(virial);

    // check whether all parameters were used
	mardyn_assert(params.eos());

    //Here we have F_{cg}
    double r_com = std::sqrt(SqrdDistanceBetweenCOMs(com1,com2));
    std::array<double,3> f_com;
    {
        double Upot =0.0;
        for(int i=0;i<f_com.size();i++){
            f_com[i] = com2[i]-com1[i];
        }
    
        //Compute potential and force
        Upot = PotentialOfMeanForce(r_com);
        ForceOfPotentialOfMeanForce(f_com,r_com);

    }

    for(int i=0;i<f_com.size();i++){
        f_com[i] *= (1.0-w1*w2);
    }
    //need to map f_{com} to F_{at}
	MapToAtomistic(m1,m2,f_com);
    
}

void InteractionForceAdapter::PotForceOnlyCG(Molecule& m1, Molecule& m2, ParaStrm& params, double* distance, double& Upot6LJ, double& UpotXPoles, double& MyRF, double virial[3], bool calcLJ){

	virial[0] = 0.0;
	virial[1] = 0.0;
	virial[2] = 0.0;

    std::array<double, 3> com1 = CenterOfMass(m1);
    std::array<double, 3> com2 = CenterOfMass(m2);

    // std::array<double, 3> com1 = resolution_handler.GetMoleculeTrackerPosition(m1.getID());
    // std::array<double, 3> com2 = resolution_handler.GetMoleculeTrackerPosition(m2.getID());


    double r_com = std::sqrt(SqrdDistanceBetweenCOMs(com1,com2));
    //Interact only on com sites
    {
        double Upot =0.0;
        std::array<double,3> f;
        for(int i=0;i<f.size();i++){
            f[i] = com2[i]-com1[i];
        }
    
        //Compute potential and force
        Upot = PotentialOfMeanForce(r_com);
        ForceOfPotentialOfMeanForce(f,r_com);


        MapToAtomistic(m1,m2,f);

		Upot6LJ += Upot;

		for(unsigned short d=0;d<3; ++d){
			virial[d] += 0.5*distance[d]*f[d];
		}

    }

}


void InteractionForceAdapter::FluidPotType(Molecule& m1, Molecule& m2, ParaStrm& params, ParaStrm& paramInv, double* drm, double& Upot6LJ, double& UpotXpoles, double& MyRF, double Virial[3], bool calcLJ, bool hybrid){
    //Pure interaction case
    if(!hybrid){
        FluidPot(m1,m2,params,drm,Upot6LJ,UpotXpoles,MyRF,calcLJ);
    }
}


/****************************************************************
 * ************** Potential and Force Calculations **************
 ***************************************************************/

double InteractionForceAdapter::PotentialOfMeanForce(double r){

    double potential;
    potential = adres->GetPotentialInterpolation().InterpolateAt(r);

    return potential;
}

void InteractionForceAdapter::ForceOfPotentialOfMeanForce(std::array<double,3>& f_com, double r){

    double norm = 0.0;
    for(int i=0;i<f_com.size();++i){
        norm += std::pow(f_com[i],2.0);
    }
    norm = std::sqrt(norm);

    double derivative = 1.0* adres->GetPotentialInterpolation().CentralFiniteDifference(r);
    for(int i=0;i<f_com.size();i++){
        f_com[i] *= derivative/norm;
    }

}

double InteractionForceAdapter::SqrdDistanceBetweenCOMs(std::array<double,3> c1,std::array<double,3> c2){
    double r2 =0.0;

    for(int i=0;i<c1.size();i++){
        r2 += std::abs((c1[i]-c2[i])*(c1[i]-c2[i]));
    }

    return r2;
}

std::array<double, 3> InteractionForceAdapter::CenterOfMass(Molecule& m1){

    std::array<double, 3> com{0,0,0};
    double total_mass = m1.component()->m();

    for(int lj=0;lj<m1.component()->numLJcenters();++lj){
        auto lj_site = m1.ljcenter_d_abs(lj);
        double site_mass = m1.component()->ljcenter(lj).m();
        for(int i=0;i<3;++i){
            com[i] += lj_site[i]*site_mass;
            com[i] = com[i]/total_mass;
        }
    }

    return com;

}

void InteractionForceAdapter::MapToAtomistic(Molecule& m1, Molecule& m2, std::array<double,3>& force){
    double total_mass = m1.component()->m();

    for(int lj=0;lj<m1.component()->numLJcenters();++lj){
        std::array<double,3> site_force{0,0,0};
        double site_mass = m1.component()->ljcenter(lj).m();
        for(int i=0;i<3;++i){
            site_force[i] = force[i] * site_mass/total_mass;
        }

        m1.Fljcenteradd(lj,site_force.data());
        m2.Fljcentersub(lj,site_force.data());

    }

}