#include "ForceAdapter.h"
#include "molecules/potforce.h"
#include "PMF.h"

InteractionForceAdapter::InteractionForceAdapter(ResolutionHandler& handle, PMF* pmf):resolution_handler{handle},adres{pmf}{

}

void InteractionForceAdapter::init(){

}

void InteractionForceAdapter::finish(){

}

double InteractionForceAdapter::processPair(Molecule& m1, Molecule& m2, double distance[3], PairType pair, double dd, bool CalculateLJ){
    std::vector<FPRegion>& regions = adres->GetRegions();
    //check if any of the 2 molecules is hybrid
    bool has_hybrid = false;
    InteractionType interaction;
    if(adres->GetMoleculeResolution(m1.getID())==Hybrid || adres->GetMoleculeResolution(m2.getID()==Hybrid)){
        has_hybrid = true;
        interaction = InteractionType::mixed;
    }

    //Check if both are 
    if(adres->GetMoleculeResolution(m1.getID())==CoarseGrain && adres->GetMoleculeResolution(m2.getID())==CoarseGrain){
        interaction = InteractionType::onlycg;
    }

    //check if any of the 2 molecules is coarsegrain
    //Check if both are 
    if(adres->GetMoleculeResolution(m1.getID())==CoarseGrain && adres->GetMoleculeResolution(m2.getID())==CoarseGrain){
        interaction = InteractionType::onlyfp;
    }
    return processPairBackend(m1, m2, distance, pair, dd, CalculateLJ, interaction);

}

double InteractionForceAdapter::processPairBackend(Molecule& m1, Molecule& m2, double distance[3], PairType pair, double dd, bool calcLJ, InteractionType interaction){
    const int tid = mardyn_get_thread_num();
    ParticlePairs2PotForceAdapter::PP2PFAThreadData& data = *thread_data[tid];
    ParaStrm params;
    ParaStrm paramsInv;

    switch(pair){
        double Virial[3];
        case MOLECULE_MOLECULE:
            this->PotForceType(m1,m2,params,paramsInv,distance,data._upot6LJ,data._upotXpoles,data._myRF,Virial,calcLJ,interaction);

        case MOLECULE_HALOMOLECULE:

    }
}

void InteractionForceAdapter::PotForceType(Molecule& m1, Molecule& m2, ParaStrm& params, ParaStrm& paramInv, double* drm, double& Upot6LJ, double& UpotXpoles, double& MyRF, double Virial[3], bool calcLJ, InteractionType interaction){
    //Pure interaction case
    if(interaction == InteractionType::onlyfp){
        PotForce(m1,m2,params,drm,Upot6LJ, UpotXpoles, MyRF, Virial,calcLJ);
    }

    if(interaction == InteractionType::onlycg){
        PotForceOnlyCG(m1,m2,params,drm,Upot6LJ,UpotXpoles,MyRF,Virial,calcLJ);
    }

    if(interaction == InteractionType::mixed){
        /**
         * m1 hy vs m2 cg
         * m1 hy vs m2 at
         * viceversa
         */
        PotForceHybrid(m1,m2,params,drm,Upot6LJ,UpotXpoles,MyRF,Virial,calcLJ, region);
    }
}

void InteractionForceAdapter::PotForceHybrid(Molecule& m1, Molecule& m2, ParaStrm& params, double* distance, double& Upot6LJ, double& UpotXPoles, double& MyRF, double virial[3], bool calcLJ, FPRegion& region){

    std::array<double, 3> com1 = adres->GetMoleculeCOMSite(m1.getID()).r();
    std::array<double, 3> com2 = adres->GetMoleculeCOMSite(m2.getID()).r();

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
                    f[i] *= w1*w2;
                }

				m1.Fljcenteradd(si, f);
				m2.Fljcentersub(sj, f);
				Upot6LJ += u;
				for (unsigned short d = 0; d < 3; ++d)
					virial[d] += 0.5*distance[d] * f[d];
			}
		}
	}

    double am1[3], am2[3]; // angular momenta

	const unsigned ne1 = m1.numCharges();
	const unsigned ne2 = m2.numCharges();
	const unsigned int nq1 = m1.numQuadrupoles();
	const unsigned int nq2 = m2.numQuadrupoles();
	const unsigned int nd1 = m1.numDipoles();
	const unsigned int nd2 = m2.numDipoles();
	for (unsigned si = 0; si < ne1; si++) {
		const std::array<double,3> dii = m1.charge_d_abs(si);
		// Charge-Charge
		for (unsigned sj = 0; sj < ne2; sj++) {
			const std::array<double,3> djj = m2.charge_d_abs(sj);
			double q1q2per4pie0; // 4pie0 = 1 in reduced units
			params >> q1q2per4pie0;
			SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
			PotForce2Charge(drs, dr2, q1q2per4pie0, f, u);

            //Introduce weight function values
            for(int i=0;i<3;i++){
                f[i] *= w1*w2;
            }

			m1.Fchargeadd(si, f);
			m2.Fchargesub(sj, f);

			UpotXPoles += u;
			for (unsigned short d = 0; d < 3; d++)
				virial[d] += 0.5*distance[d] * f[d];
		}
		// Charge-Quadrupole
		for (unsigned sj = 0; sj < nq2; sj++) {
			const std::array<double,3> djj = m2.quadrupole_d_abs(sj);
			double qQ05per4pie0; // 4pie0 = 1 in reduced units
			params >> qQ05per4pie0;
			SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
			const std::array<double,3> ejj = m2.quadrupole_e(sj);
			PotForceChargeQuadrupole(drs, dr2, ejj.data(), qQ05per4pie0, f, am2, u);

            //Introduce weight function values
            for(int i=0;i<3;i++){
                f[i] *= w1*w2;
            }

			m1.Fchargeadd(si, f);
			m2.Fquadrupolesub(sj, f);
			m2.Madd(am2);

			UpotXPoles += u;
			for (unsigned short d = 0; d < 3; d++)
				virial[d] += 0.5*distance[d] * f[d];
		}
		// Charge-Dipole
		for (unsigned sj = 0; sj < nd2; sj++) {
			const std::array<double,3> djj = m2.dipole_d_abs(sj);
			double minusqmyper4pie0;
			params >> minusqmyper4pie0;
			SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
			const std::array<double,3> ejj = m2.dipole_e(sj);
			PotForceChargeDipole(drs, dr2, ejj.data(), minusqmyper4pie0, f, am2, u);

            //Introduce weight function values
            for(int i=0;i<3;i++){
                f[i] *= w1*w2;
            }

			m1.Fchargeadd(si, f);
			m2.Fdipolesub(sj, f);
			m2.Madd(am2);

			UpotXPoles += u;
			for (unsigned short d = 0; d < 3; d++)
				virial[d] += 0.5*distance[d] * f[d];
		}
	}
    //quadrupoles
    for (unsigned int si = 0; si < nq1; ++si) {
		const std::array<double,3> dii = m1.quadrupole_d_abs(si);
		const std::array<double,3> eii = m1.quadrupole_e(si);

		// Quadrupole-Charge
		for (unsigned sj = 0; sj < ne2; sj++) {
			const std::array<double,3> djj = m2.charge_d_abs(sj);
			double qQ05per4pie0; // 4pie0 = 1 in reduced units
			params >> qQ05per4pie0;
			minusSiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
			PotForceChargeQuadrupole(drs, dr2, eii.data(), qQ05per4pie0, f, am1, u);

            for(int i=0;i<3;i++){
                f[i] *= w1*w2;
            }

			m1.Fquadrupolesub(si, f);
			m2.Fchargeadd(sj, f);
			m1.Madd(am1);

			UpotXPoles += u;
			for (unsigned short d = 0; d < 3; d++)
				virial[d] -= 0.5*distance[d] * f[d];
		}
		// Quadrupole-Quadrupole -------------------
		for (unsigned int sj = 0; sj < nq2; ++sj) {
			//double drs[3];
			const std::array<double,3> djj = m2.quadrupole_d_abs(sj);
			double q2075;
			params >> q2075;
			SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
			const std::array<double,3> ejj = m2.quadrupole_e(sj);
			PotForce2Quadrupole(drs, dr2, eii.data(), ejj.data(), q2075, f, am1, am2, u);

            for(int i=0;i<3;i++){
                f[i] *= w1*w2;
            }

			m1.Fquadrupoleadd(si, f);
			m2.Fquadrupolesub(sj, f);
			m1.Madd(am1);
			m2.Madd(am2);

			UpotXPoles += u;
			for (unsigned short d = 0; d < 3; ++d)
				virial[d] += 0.5*distance[d] * f[d];
		}
		// Quadrupole-Dipole -----------------------
		for (unsigned int sj = 0; sj < nd2; ++sj) {
			//double drs[3];
			const std::array<double,3> djj = m2.dipole_d_abs(sj);
			double qmy15;
			params >> qmy15;
			minusSiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
			const std::array<double,3> ejj = m2.dipole_e(sj);
			PotForceDiQuadrupole(drs, dr2, ejj.data(), eii.data(), qmy15, f, am2, am1, u);

            for(int i=0;i<3;i++){
                f[i] *= w1*w2;
            }

			m1.Fquadrupolesub(si, f);
			m2.Fdipoleadd(sj, f);
			m1.Madd(am1);
			m2.Madd(am2);
			UpotXPoles += u;
			for (unsigned short d = 0; d < 3; d++)
				virial[d] -= 0.5*distance[d] * f[d];
		}
	}
    //dipoles
    for (unsigned int si = 0; si < nd1; ++si) {
		const std::array<double,3> dii = m1.dipole_d_abs(si);
		const std::array<double,3> eii = m1.dipole_e(si);
		// Dipole-Charge
		for (unsigned sj = 0; sj < ne2; sj++) {
			const std::array<double,3> djj = m2.charge_d_abs(sj);
			double minusqmyper4pie0;
			params >> minusqmyper4pie0;
			minusSiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
			PotForceChargeDipole(drs, dr2, eii.data(), minusqmyper4pie0, f, am1, u);

            for(int i=0;i<3;i++){
                f[i] *= w1*w2;
            }

			m1.Fdipolesub(si, f);
			m2.Fchargeadd(sj, f);
			m1.Madd(am1);

			UpotXPoles += u;
			for (unsigned short d = 0; d < 3; d++)
				virial[d] -= 0.5*distance[d] * f[d];
		}
		// Dipole-Quadrupole -----------------------
		for (unsigned int sj = 0; sj < nq2; ++sj) {
			//double drs[3];
			const std::array<double,3> djj = m2.quadrupole_d_abs(sj);
			double myq15;
			params >> myq15;
			SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
			const std::array<double,3> ejj = m2.quadrupole_e(sj);
			PotForceDiQuadrupole(drs, dr2, eii.data(), ejj.data(), myq15, f, am1, am2, u);

            for(int i=0;i<3;i++){
                f[i] *= w1*w2;
            }

			m1.Fdipoleadd(si, f);
			m2.Fquadrupolesub(sj, f);
			m1.Madd(am1);
			m2.Madd(am2);
			UpotXPoles += u;
			for (unsigned short d = 0; d < 3; ++d)
				virial[d] += 0.5*distance[d] * f[d];
		}
		// Dipole-Dipole ---------------------------
		for (unsigned int sj = 0; sj < nd2; ++sj) {
			const std::array<double,3> djj = m2.dipole_d_abs(sj);
			double my2;
			params >> my2;
			double rffac;
			params >> rffac;
			SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
			const std::array<double,3> ejj = m2.dipole_e(sj);
			PotForce2Dipole(drs, dr2, eii.data(), ejj.data(), my2, rffac, f, am1, am2, u, MyRF);

            for(int i=0;i<3;i++){
                f[i] *= w1*w2;
            }

			m1.Fdipoleadd(si, f);
			m2.Fdipolesub(sj, f);
			m1.Madd(am1);
			m2.Madd(am2);
			UpotXPoles += u;
			for (unsigned short d = 0; d < 3; ++d)
				virial[d] += 0.5*distance[d] * f[d];
		}
	}

    m1.Viadd(virial);
	m2.Viadd(virial);

    //Here we have F_{cg}
    std::array<double, 3> com1 = adres->GetMoleculeCOMSite(m1.getID()).r();
    std::array<double, 3> com2 = adres->GetMoleculeCOMSite(m2.getID()).r();
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
        //Set quantities
        adres->GetMoleculeCOMSite(m1.getID()).AddPotential(Upot);
        adres->GetMoleculeCOMSite(m1.getID()).AddForce(f_com);

    }

    for(int i=0;i<f_com.size();i++){
        f_com[i] *= (1.0-w1*w2);
    }
    

    //need to map f_{com} to F_{at}
    



}

void InteractionForceAdapter::PotForceOnlyCG(Molecule& m1, Molecule& m2, ParaStrm& params, double* distance, double& Upot6LJ, double& UpotXPoles, double& MyRF, double virial[3], bool calcLJ){

    //Get both molecules COMs
    //TODO: Check this assumption: Assume they are updated
    std::array<double, 3> com1 = adres->GetMoleculeCOMSite(m1.getID()).r();
    std::array<double, 3> com2 = adres->GetMoleculeCOMSite(m2.getID()).r();

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
        //Set quantities
        adres->GetMoleculeCOMSite(m1.getID()).AddPotential(Upot);
        adres->GetMoleculeCOMSite(m1.getID()).AddForce(f);

    }

}


void InteractionForceAdapter::HybridFluidPot(Molecule& m1, Molecule& m2, ParaStrm& params, ParaStrm& paramInv, double* drm, double& Upot6LJ, double& UpotXpoles, double& MyRF, double Virial[3], bool calcLJ, bool hybrid){
    //Pure interaction case
    if(!hybrid){
        FluidPot(m1,m2,params,drm,Upot6LJ,UpotXpoles,MyRF,calcLJ);
    }
}


/****************
 * ************** Potential and Force Calculations
 ***************/

double InteractionForceAdapter::PotentialOfMeanForce(double r){
    double dist = adres->GetRDFInterpolation().GetRDFAt(r);
    return -1.0*_simulation.getEnsemble()->T()*std::log(dist);
}

void InteractionForceAdapter::ForceOfPotentialOfMeanForce(std::array<double,3>& f_com, double r){


    double derivative = adres->GetRDFInterpolation().CentralFiniteDifference(r);
    double rdf = adres->GetRDFInterpolation().GetRDFAt(r);
    double f_scalar = -1.0*_simulation.getEnsemble()->T()*derivative/rdf;
    for(int i=0;i<f_com.size();i++){
        f_com[i] *= f_scalar;
    }
}

double InteractionForceAdapter::SqrdDistanceBetweenCOMs(std::array<double,3> c1,std::array<double,3> c2){
    double r2 =0.0;

    for(int i=0;i<c1.size();i++){
        r2 += std::abs((c1[i]-c2[i])*(c1[i]-c2[i]));
    }

    return r2;
}