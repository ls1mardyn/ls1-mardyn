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
}

void InteractionForceAdapter::PotForceOnlyCG(Molecule& m1, Molecule& m2, ParaStrm& params, double* distance, double& Upot6LJ, double& UpotXPoles, double& MyRF, double virial[3], bool calcLJ){

    //Get both molecules COMs
    //Assume they are updated
    std::array<double, 3> com1 = adres->GetMoleculeCOMSite(m1.getID()).r();
    std::array<double, 3> com2 = adres->GetMoleculeCOMSite(m2.getID()).r();

    double temperature = _simulation.getEnsemble()->T();
    double r_com =0.0;
    //Interact only on com sites
    {
        double Upot =0.0;
        std::array<double,3> f;
        //Compute force
        //Compute potential
        //Set quantities
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
    return -1.0*_simulation.getEnsemble()->T()*std::log(r);
}