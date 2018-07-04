//
// Created by kruegener on 6/8/2018.
// moved from /molecules/Wall.cpp by mheinen
//

#include "WallPotential.h"

void WallPotential::readXML(XMLfileUnits &xmlconfig) {
    double density, sigma, epsilon, yoff, ycut;
    xmlconfig.getNodeValue("density", density);
    xmlconfig.getNodeValue("sigma", sigma);
    xmlconfig.getNodeValue("epsilon", epsilon);
    xmlconfig.getNodeValue("yoff", yoff);
    xmlconfig.getNodeValue("ycut", ycut);
    xmlconfig.getNodeValue("width", _dWidth);
    xmlconfig.getNodeValue("delta", _delta);
    _dWidthHalf = _dWidth * 0.5;
    global_log->info() << "[WallPotential] Using plugin with parameters: density=" << density << ", "
                                                                                                 "sigma=" << sigma << ", epsilon=" << epsilon << ", yoff=" << yoff << ", ycut=" << ycut << ", "
                                                                                                                                                                                           "width=" << _dWidth << ", delta=" << _delta << endl;

    int potential;
    xmlconfig.getNodeValue("potential", potential);
    if(potential == 93){
        _potential = LJ9_3;
    }
    else if(potential == 104){
        _potential = LJ10_4;
    }
    else{
        global_log -> info() << "[WallPotential] No valid potential was specified -> Default: LJ9_3" << std::endl;
        _potential = LJ9_3;
        // TODO: is this allowed or should simulation be halted
        // HALT SIM
        //global_simulation -> exit(1);
    }

    XMLfile::Query query = xmlconfig.query("component");
    unsigned int numComponentsConsidered = query.card();
    global_log->info() << "[WallPotential] Setting parameters 'xi', 'eta' for " << numComponentsConsidered << " components." << endl;

    std::vector<Component>* components = global_simulation->getEnsemble()->getComponents();
    unsigned int numComponents = components->size();
    std::vector<double> xi_sf(numComponents);
    std::vector<double> eta_sf(numComponents);

    if(numComponentsConsidered == 0){
        for(auto&& bi : _bConsiderComponent)
            bi = true;
    }
    else {
        _bConsiderComponent.resize(numComponents);
        for (auto &&bi : _bConsiderComponent)
            bi = false;

        string oldpath = xmlconfig.getcurrentnodepath();
        XMLfile::Query::const_iterator componentIter;
        for (componentIter = query.begin(); componentIter; componentIter++) {
            xmlconfig.changecurrentnode(componentIter);
            unsigned int cid;
            xmlconfig.getNodeValue("@id", cid);
            xmlconfig.getNodeValue("xi", xi_sf.at(cid - 1));
            xmlconfig.getNodeValue("eta", eta_sf.at(cid - 1));
            _bConsiderComponent.at(cid - 1) = true;
            global_log->info() << "[WallPotential] " << cid << " " << xi_sf.at(cid - 1) << " " << eta_sf.at(cid - 1)
                               << " " << std::endl;
        }
        xmlconfig.changecurrentnode(oldpath);
    }


    if(_potential == LJ9_3){
        this->initializeLJ93(components, density, sigma, epsilon, xi_sf, eta_sf, yoff, ycut);
        global_log->info() << "[WallPotential] LJ9_3 initialized." << endl;
    }
    else if(_potential == LJ10_4){
        this->initializeLJ1043(components, density, sigma, epsilon, xi_sf, eta_sf, yoff, ycut);
        global_log->info() << "[WallPotential] LJ10_4 initialized." << endl;
    }
    else{
        global_log -> error() << "[WallPotential] UNKNOWN WALL POTENTIAL! EXITING!" << std::endl;
        global_simulation -> exit(11);
    }

}

void WallPotential::initializeLJ93(const std::vector<Component>* components,
                                   double in_rhoWall, double in_sigWall, double in_epsWall, std::vector<double> in_xi, std::vector<double> in_eta,
                                   double in_yOffWall, double in_yWallCut) {

    global_log->info() << "[WallPotential] Initializing the wall function LJ-9-3.\n";
    this->_rhoW = in_rhoWall;
    this->_yc = in_yWallCut;
    this->_yOff = in_yOffWall;

    /*!*** So far: only 1CLJ components allowed ****/
    _nc = components->size();
    _eps_wi = new double[_nc];
    _sig3_wi = new double[_nc];
    _uShift_9_3 = new double[_nc];
    _uPot_9_3 = new double[_nc];

    for (unsigned i = 0; i < _nc; i++) {
        _eps_wi[i] = in_xi[i] * sqrt(in_epsWall * (components->at(i)).ljcenter(0).eps());

        double sig_wi = 0.5 * in_eta[i] * (in_sigWall + (components->at(i)).ljcenter(0).sigma());
        _sig3_wi[i] = sig_wi * sig_wi * sig_wi;
        double sig9_sf = _sig3_wi[i] * _sig3_wi[i] * _sig3_wi[i];

        double y = _yc;
        double y3 = y * y * y;
        double y9 = y3 * y3 * y3;

        _uShift_9_3[i] = 4.0 / 3.0 * M_PI * _rhoW * _eps_wi[i] *_sig3_wi[i] * (sig9_sf / 15.0 / y9 - _sig3_wi[i] / 2.0 / y3);
        _uPot_9_3[i] = 0.0;
    }
}

void WallPotential::initializeLJ1043(const std::vector<Component> *components,
                                     double in_rhoWall, double in_sigWall, double in_epsWall, std::vector<double> in_xi,
                                     std::vector<double> in_eta,
                                     double in_yOffWall, double in_yWallCut) {

    global_log->info() << "[WallPotential] Initializing the wall function LJ-10-4-3 with " << _nc << " components.\n" << endl;

    this->_rhoW = in_rhoWall;
    this->_yc = in_yWallCut;
    this->_yOff = in_yOffWall;

    /*!*** So far: only 1CLJ components allowed ****/
    _nc = components->size();
    _eps_wi = new double[_nc];
    _sig3_wi = new double[_nc];
    _sig2_wi = new double[_nc];
    _sig_wi = new double[_nc];
    _uShift_10_4_3 = new double[_nc];
    _uPot_10_4_3 = new double[_nc];

    for (unsigned i = 0; i < _nc; i++) {
        _eps_wi[i] = in_xi[i] * sqrt(in_epsWall * (components->at(i)).ljcenter(0).eps());

        double sig_wi = 0.5 * in_eta[i] * (in_sigWall + (components->at(i)).ljcenter(0).sigma());
        _sig_wi[i] = sig_wi;
        _sig2_wi[i] = sig_wi * sig_wi;
        _sig3_wi[i] = sig_wi * sig_wi * sig_wi;
        double sig10_sf = _sig3_wi[i] * _sig3_wi[i] * _sig3_wi[i] * _sig_wi[i];
        double sig4_sf = _sig2_wi[i] * _sig2_wi[i];

        double y = _yc;
        double y2 = y * y;
        double y4 = y2 * y2;
        double y10 = y4 * y4 * y2;

        double bracket = y + 0.61 * _delta;
        double bracket3 = bracket * bracket * bracket;

        _uShift_10_4_3[i] = 2 * M_PI * _eps_wi[i] * _rhoW * _sig2_wi[i] * _delta * ((2 / 5) * (sig10_sf / y10) - sig4_sf / y4 - sig4_sf / (3 * _delta * bracket3));
        _uPot_10_4_3[i] = 0.0;
    }
}

void WallPotential::calcTSLJ_9_3(ParticleContainer *partContainer) {
    double regionLowCorner[3], regionHighCorner[3];

    /*! LJ-9-3 potential applied in y-direction */
    if(partContainer->getBoundingBoxMin(1) < _yc+_yOff){ // if linked cell within the potential range (inside the potential's cutoff)
        for(unsigned d = 0; d < 3; d++){
            regionLowCorner[d] = partContainer->getBoundingBoxMin(d);
            regionHighCorner[d] = partContainer->getBoundingBoxMax(d);
        }

        //perform a check if the region is contained by the particleContainer???
        if (partContainer->isRegionInBoundingBox(regionLowCorner, regionHighCorner)){
#if defined (_OPENMP)
#pragma omp parallel shared(regionLowCorner, regionHighCorner)
#endif
            {
                RegionParticleIterator begin = partContainer->regionIterator(regionLowCorner, regionHighCorner);

                double f[3];
                for(unsigned d=0; d<3; d++)
                    f[d] = 0.;

                for(RegionParticleIterator i = begin; i.hasNext(); i.next()){
                    unsigned cid = (*i).componentid();
                    if(false == _bConsiderComponent.at(cid) )
                        continue;  // only add Wall force to molecules of component that should be considered

                    for(unsigned int si=0; si<i->numLJcenters(); ++si) {
                        double y, y3, y9, ry, ryRel;
                        // TODO: ljcenter_d or _d_abs ?
                        const std::array<double, 3> arrSite = i->ljcenter_d_abs(si);
                        const double *posSite = arrSite.data();
                        ry = posSite[1];
                        ryRel = (ry > _yOff) ? (ry - (_yOff + _dWidthHalf)) : (ry - (_yOff - _dWidthHalf));
                        y = abs(ryRel);
                        //y = ryRel;

                        if (y < _yc) {
                            y3 = y * y * y;
                            y9 = y3 * y3 * y3;
                            for (unsigned d = 0; d < 3; d++) {
                                f[d] = 0.0;
                            }

                            double sig9_wi;
                            sig9_wi = _sig3_wi[cid] * _sig3_wi[cid] * _sig3_wi[cid];
                            f[1] = 4.0 * M_PI * _rhoW * _eps_wi[cid] * _sig3_wi[cid] *
                                   (sig9_wi / 5.0 / y9 - _sig3_wi[cid] / 2.0 / y3) / y;
                            if(ryRel < 0){
                                f[1] = -f[1];
                            }
                            _uPot_9_3[cid] += 4.0 * M_PI * _rhoW * _eps_wi[cid] * _sig3_wi[cid] *
                                              (sig9_wi / 45.0 / y9 - _sig3_wi[cid] / 6.0 / y3) - _uShift_9_3[cid];
                            f[0] = 0;
                            f[2] = 0;
                            i->Fljcenteradd(si, f);
                        }
                    }
                }
            }
        }
    }

    double u_pot;
    u_pot = _uPot_9_3[0] + _domain -> getLocalUpotCompSpecific();
    _domain->setLocalUpotCompSpecific(u_pot);
    for(unsigned cid = 0; cid < _nc; cid++) {
        _uPot_9_3[cid] = 0.0;
    }
} // end method calcTSLJ_9_3(...)

void WallPotential::calcTSLJ_10_4(ParticleContainer *partContainer) {

    double regionLowCorner[3], regionHighCorner[3];


    /*! LJ-10-4 potential applied in y-direction */
    //if(partContainer->getBoundingBoxMin(1)< _yc+_yOff ){ // if linked cell within the potential range (inside the potential's cutoff)
    for(unsigned d = 0; d < 3; d++){
        regionLowCorner[d] = partContainer->getBoundingBoxMin(d);
        regionHighCorner[d] = partContainer->getBoundingBoxMax(d);
    }


    //perform a check if the region is contained by the particleContainer???
    if(partContainer->isRegionInBoundingBox(regionLowCorner, regionHighCorner)){

#if defined (_OPENMP)
#pragma omp parallel shared(regionLowCorner, regionHighCorner)
#endif
        {
            RegionParticleIterator begin = partContainer->regionIterator(regionLowCorner, regionHighCorner);

            for(RegionParticleIterator i = begin; i.hasNext() ; i.next()){
                //! so far for 1CLJ only, several 1CLJ-components possible
                double ry, ryRel, y, y2, y4, y5, y10, y11;
                unsigned cid = (*i).componentid();
                if(false == _bConsiderComponent.at(cid) )
                    continue;  // only add Wall force to molecules of component that should be considered

                for(unsigned int si=0; si<i->numLJcenters(); ++si) {
                    const std::array<double,3> arrSite = i->ljcenter_d_abs(si);
                    const double* posSite = arrSite.data();
                    //y = (*i).r(1) - _yOff;
                    ry = posSite[1];
                    ryRel = (ry > _yOff) ? (ry - (_yOff + _dWidthHalf)) : (ry - (_yOff - _dWidthHalf));
                    y = abs(ryRel);
                    if (y < _yc) {
                        y2 = y * y;
                        y4 = y2 * y2;
                        y5 = y4 * y;
                        y10 = y5 * y5;
                        y11 = y10 * y;
                        double f[3];
                        for (unsigned d = 0; d < 3; d++) {
                            f[d] = 0.0;
                        }

                        double sig2_wi = _sig2_wi[cid];
                        double sig4_wi = _sig2_wi[cid] * _sig2_wi[cid];
                        double sig5_wi = sig4_wi * _sig_wi[cid];
                        double sig10_wi = sig5_wi * sig5_wi;
                        double bracket = y + 0.61 * _delta;
                        double bracket3 = bracket * bracket * bracket;
                        double term1 = sig10_wi / y10;
                        double term2 = sig4_wi / y4;
                        double term3 = sig4_wi / (3 * _delta * bracket3);
                        double preFactor = 2 * M_PI * _eps_wi[cid] * _rhoW * sig2_wi * _delta;
                        _uPot_10_4_3[cid] += preFactor * ((2 / 5) * term1 - term2 - term3) - _uShift_10_4_3[cid];
                        f[1] = preFactor * (4 * (sig10_wi / y11) - 4 * (sig4_wi / y5) - term3 * 3 / bracket);
                        if(ryRel < 0){
                            f[1] = -f[1];
                        }
                        f[0] = 0;
                        f[2] = 0;

                        i->Fljcenteradd(si, f);
                        //i->calcFM_site(i->ljcenter_d(si), i->ljcenter_F(si));
                        //i->Fadd(f);
                    }
                }
            }
        }
    }

    double u_pot;
    // TODO: sum up upot for all components
    u_pot = _uPot_10_4_3[0] + _domain -> getLocalUpotCompSpecific();
    _domain->setLocalUpotCompSpecific(u_pot);
    for(unsigned cid = 0; cid < _nc; cid++) {
        _uPot_10_4_3[cid] = 0.0;
    }
}
// end method calcTSLJ_10_4(...)

void WallPotential::beforeForces(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
                                 unsigned long simstep) {

}

void WallPotential::afterForces(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
                                unsigned long simstep) {
    /*
    if(simstep == 0){
        return;
    }
    if(_potential == LJ9_3){
        //global_log->debug() << "[WallPotential] LJ9_3 afterForces." << endl;
        this->calcTSLJ_9_3(particleContainer);
        //global_log->debug() << "[WallPotential] LJ9_3 applied." << endl;
    }
    else if(_potential == LJ10_4){
        //global_log->debug() << "[WallPotential] LJ10_4 afterForces. " << endl;
        this->calcTSLJ_10_4(particleContainer);
        //global_log->debug() << "[WallPotential] LJ10_4 applied." << endl;
    }
    */
}

void WallPotential::forceStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
                              unsigned long simstep) {

    if(simstep == 0){
        return;
    }
    if(_potential == LJ9_3){
        //global_log->debug() << "[WallPotential] LJ9_3 afterForces." << endl;
        this->calcTSLJ_9_3(particleContainer);
        //global_log->debug() << "[WallPotential] LJ9_3 applied." << endl;
    }
    else if(_potential == LJ10_4){
        //global_log->debug() << "[WallPotential] LJ10_4 afterForces. " << endl;
        this->calcTSLJ_10_4(particleContainer);
        //global_log->debug() << "[WallPotential] LJ10_4 applied." << endl;
    }

}