#include "ForceAdapter.h"
#include "molecules/potforce.h"
#include "PMF.h"

InteractionForceAdapter::InteractionForceAdapter(ResolutionHandlerBase* handle,  PMF* pmf):resolution_handler{handle},adres{pmf}{

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

    for(auto thread:thread_data){
        values.U_lj += thread->_upot6LJ;
        values.U_poles += thread->_upotXpoles;
        values.virial += thread->_virial;
        values.RF += thread->_myRF;
    }
    values.SetInDomain(_simulation.getDomain());
    values.ClearAll();

}

double InteractionForceAdapter::processPair(Molecule& m1, Molecule& m2, double distance[3], PairType pair, double dd, bool CalculateLJ){

    //TODO:call checkresolution here?
    std::vector<FPRegion>& regions = resolution_handler->GetRegions();


    InteractionType interaction;
    if(resolution_handler->GetMoleculeResolution(m1)==ResolutionType::Hybrid || resolution_handler->GetMoleculeResolution(m2)==ResolutionType::Hybrid)
    {
        interaction = InteractionType::mixed;
    }

    if(resolution_handler->GetMoleculeResolution(m1)==ResolutionType::CoarseGrain && resolution_handler->GetMoleculeResolution(m2)==CoarseGrain)
    {
        interaction = InteractionType::onlycg;
    }

    if(resolution_handler->GetMoleculeResolution(m1)==ResolutionType::FullParticle && resolution_handler->GetMoleculeResolution(m2) == FullParticle)
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
        // Log::global_log->error()<<"PotForceType cannot be CG"<<std::endl;

        PotForceOnlyCG(m1,m2,params,drm,Upot6LJ,UpotXpoles,MyRF,Virial,calcLJ);
    }

    if(interaction == InteractionType::mixed){
        // Log::global_log->error()<<"PotForceType cannot be Mixed"<<std::endl;
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
    ResolutionType type1, type2;
    type1 =resolution_handler->GetMoleculeResolution(m1);
    type2 =resolution_handler->GetMoleculeResolution(m2);

    //Assume only one fpregion
	FPRegion& region = resolution_handler->GetRegions()[0];

    //Pure hybrid
    if(type1 == ResolutionType::Hybrid && type2 == ResolutionType::Hybrid){
        PotForcePureHybridBackend(m1,m2,params,distance,Upot6LJ,UpotXPoles,MyRF,virial,calcLJ,region);
    }

    //One is CG
    if(type1 == ResolutionType::CoarseGrain){
        PotForceHybridCGBackend(m2,m1,params,distance,Upot6LJ,UpotXPoles,MyRF,virial,calcLJ,region);
    }
    if(type2 == ResolutionType::CoarseGrain){
        PotForceHybridCGBackend(m1,m2,params,distance,Upot6LJ,UpotXPoles,MyRF,virial,calcLJ,region);
    }

    //one is FP
    if(type1 == ResolutionType::FullParticle){
        PotForceHybridFPBackend(m2,m1,params,distance,Upot6LJ,UpotXPoles,MyRF,virial,calcLJ,region);
    }
    if(type2 == ResolutionType::FullParticle){
        PotForceHybridFPBackend(m1,m2,params,distance,Upot6LJ,UpotXPoles,MyRF,virial,calcLJ,region);
    }
}

void InteractionForceAdapter::PotForcePureHybridBackend(Molecule& m1, Molecule& m2, ParaStrm& params, double* distance, double& Upot6LJ, double& UpotXPoles, double& MyRF, double virial[3], bool calcLJ, FPRegion& region){


    //Do LJ sites first for the F_{at}
    // Force Calculation
	double f[3];
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
    // skip first LJ site 
	const unsigned int nc1 = m1.numLJcenters();
	const unsigned int nc2 = m2.numLJcenters();
	for (unsigned int si = 0; si < nc1; ++si) {
		const std::array<double,3> dii = m1.ljcenter_d_abs(si);
		for (unsigned int sj = 0; sj < nc2; ++sj) {
            if (si == 0 || sj == 0) {
                double tmp; params >> tmp; params >> tmp; params >> tmp;
                continue;
            }
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
    const std::array<double,3> dcg_1 = m1.ljcenter_d_abs(0);
    const std::array<double,3> dcg_2 = m2.ljcenter_d_abs(0);
    SiteSiteDistanceAbs(dcg_2.data(), dcg_1.data(), drs, dr2);//m2-m1
    // double r_com = std::sqrt(SqrdDistanceBetweenCOMs(com1,com2));
    double r_com = std::sqrt(dr2);
    std::array<double,3> f_com{drs[0],drs[1],drs[2]};
    NormalizeVector(f_com);
    
    double Upot =0.0;
    Upot = PotentialOfMeanForce(r_com);
    ForceOfPotentialOfMeanForce(f_com,r_com);

    Upot6LJ += Upot;

    for(int i=0;i<f_com.size();++i){
        f_com[i] *= (1.0-w1*w2);
    }

    //only at site 0 (CG site)
    m1.Fljcenteradd(0, f_com.data());
    m2.Fljcentersub(0, f_com.data());
    //TODO: what about the virial?
}

void InteractionForceAdapter::PotForceHybridFPBackend(Molecule& m1, Molecule& m2, ParaStrm& params, double* distance, double& Upot6LJ, double& UpotXPoles, double& MyRF, double virial[3], bool calcLJ, FPRegion& region){
    //TODO:we have a mess with the signs within the 
    //m1 is HY
    //m2 is FP

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
    // skip first LJ site 
	const unsigned int nc1 = m1.numLJcenters();//first site is CG, skip
	const unsigned int nc2 = m2.numLJcenters();
	for (unsigned int si = 0; si < nc1; ++si) {
		const std::array<double,3> dii = m1.ljcenter_d_abs(si);
		for (unsigned int sj = 0; sj < nc2; ++sj) {
            if(si == 0 || sj == 0) {
                double tmp; params >> tmp; params >> tmp; params >> tmp;
                continue;
            }
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

    //Compute COM force
    std::array<double,3> fp_com = ComputeCOM(m2);
    const std::array<double,3> hy_com = m1.ljcenter_d_abs(0);
    SiteSiteDistanceAbs(fp_com.data(),hy_com.data(),drs,dr2);
    std::array<double,3> f_com{drs[0],drs[1],drs[2]};
    NormalizeVector(f_com);
    double r_com = std::sqrt(dr2);
    double Upot = PotentialOfMeanForce(r_com);
    ForceOfPotentialOfMeanForce(f_com,r_com);

    //check if forces are finite
    for(int i=0;i<3;++i){
        if(!std::isfinite(f_com[i])){
            Log::global_log->error()<<"[HybridCGBackend] F["<<i<<"] is not finite"<<std::endl;
            mardyn_assert(false);
        }
    }

    for(int i=0;i<3;i++)
    {
        f_com[i] *= (1.0-w1*w2);//F^{com}
    }

    m1.Fljcenteradd(0, f_com.data());
    SubtractAndMapForceToFP(f_com,m2);


}

void InteractionForceAdapter::PotForceHybridCGBackend(Molecule& m1, Molecule& m2, ParaStrm& params, double* distance, double& Upot6LJ, double& UpotXPoles, double& MyRF, double virial[3], bool calcLJ, FPRegion& region){

    //m1 is HY
    //m2 is CG

    double drs[3], dr2; // site distance vector & length^2

    const std::array<double,3> dhy = m1.ljcenter_d_abs(0);//m1
    const std::array<double,3> dcg = m2.ljcenter_d_abs(0);//m2

    SiteSiteDistanceAbs(dcg.data(), dhy.data(), drs, dr2);//m2-m1

    double Upot =0.0;
    std::array<double,3> f{drs[0],drs[1],drs[2]};
    NormalizeVector(f);
    double r_com = std::sqrt(dr2);
    Upot = PotentialOfMeanForce(r_com);
    ForceOfPotentialOfMeanForce(f,r_com);

    //check if forces are finite
    for(int i=0;i<3;++i){
        if(!std::isfinite(f[i])){
            Log::global_log->error()<<"[HybridCGBackend] F["<<i<<"] is not finite"<<std::endl;
            mardyn_assert(false);
        }
    }

    m1.Fljcenteradd(0,f.data());
    m2.Fljcentersub(0,f.data());

    AddAndMapForceToFP(f,m1);

	Upot6LJ += Upot;
	for(unsigned short d=0;d<3; ++d){
		virial[d] += 0.5*distance[d]*f[d];
	}


}

void InteractionForceAdapter::PotForceOnlyCG(Molecule& m1, Molecule& m2, ParaStrm& params, double* distance, double& Upot6LJ, double& UpotXPoles, double& MyRF, double virial[3], bool calcLJ){

	virial[0] = 0.0;
	virial[1] = 0.0;
	virial[2] = 0.0;

    //Interact only on com sites
    {
        double Upot =0.0;
        std::array<double,3> f;
        double r_com = std::sqrt(m1.dist2(m2,f.data()));//m2-m1
        Upot = PotentialOfMeanForce(r_com);
        NormalizeVector(f);
        ForceOfPotentialOfMeanForce(f,r_com);

		Upot6LJ += Upot;

		for(unsigned short d=0;d<3; ++d){
			virial[d] += 0.5*distance[d]*f[d];
		}

        m1.Fljcenteradd(0,f.data());
        m2.Fljcentersub(0,f.data());

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

    // double norm = 0.0;
    // for(int i=0;i<f_com.size();++i){
        // norm += std::pow(f_com[i],2.0);
    // }
    // norm = std::sqrt(norm);

    double derivative = 1.0* adres->GetPotentialInterpolation().CentralFiniteDifference(r);
    for(int i=0;i<f_com.size();i++){
        f_com[i] *= derivative;
    }

}

double InteractionForceAdapter::SqrdDistanceBetweenCOMs(std::array<double,3> c1,std::array<double,3> c2){
    double r2 =0.0;

    for(int i=0;i<c1.size();i++){
        r2 += std::abs((c1[i]-c2[i])*(c1[i]-c2[i]));
    }

    return r2;
}


void InteractionForceAdapter::AddAndMapForceToFP(std::array<double,3>& f_com, Molecule& mol){

    double total_mass = mol.component()->m();
    int start =0;
    //if hybrid then skip first lj site
    ResolutionType type= resolution_handler->GetMoleculeResolution(mol);
    if(type == ResolutionType::Hybrid){
        start =1;
    }

    for(int lj=start;lj<mol.component()->numLJcenters();++lj){

        std::array<double,3> force_i{0,0,0};
        double mass_i = mol.component()->ljcenter(lj).m();
        for(int j=0;j<3;++j){
            force_i[j]= f_com[j]*mass_i/total_mass;
        }

        mol.Fljcenteradd(lj,force_i.data());
    }
}

void InteractionForceAdapter::SubtractAndMapForceToFP(std::array<double,3>& f_com, Molecule& mol){

    double total_mass = mol.component()->m();

    int start =0;
    //if hybrid then skip first lj site
    ResolutionType type= resolution_handler->GetMoleculeResolution(mol);
    if(type == ResolutionType::Hybrid){
        start =1;
    }


    for(int lj=start;lj<mol.component()->numLJcenters();++lj){

        std::array<double,3> force_i{0,0,0};
        double mass_i = mol.component()->ljcenter(lj).m();
        for(int j=0;j<3;++j){
            force_i[j]= f_com[j]*mass_i/total_mass;
        }

        mol.Fljcentersub(lj,force_i.data());
    }
}

