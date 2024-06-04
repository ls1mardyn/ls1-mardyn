/*
 * Copyright (c) 2013      Christoph Niethammer <christoph.niethammer@gmail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#include "utils/generator/Basis.h"
#include "utils/Logger.h"
#include "Simulation.h"
#include "ensemble/EnsembleBase.h"
#include "utils/Coordinate3D.h"



void Basis::readXML(XMLfileUnits& xmlconfig) {
	XMLfile::Query query = xmlconfig.query("site");
	Ensemble* ensemble = _simulation.getEnsemble();
	std::string oldpath = xmlconfig.getcurrentnodepath();
	const size_t numComps = ensemble->getComponents()->size();
	for(auto siteIter = query.begin(); siteIter; siteIter++) {
		Molecule molecule;
		xmlconfig.changecurrentnode(siteIter);
		int componentid;
		if (xmlconfig.getNodeValue("componentid", componentid)) {
			if ((componentid < 1) || (componentid > numComps)) {
				Log::global_log->error() << "[Basis] Specified componentid is invalid. Valid range: 1 <= componentid <= " << numComps << std::endl;
				Simulation::exit(1);
			}
		} else {
			Log::global_log->error() << "[Basis] No componentid specified. Set <componentid>!" << std::endl;
			Simulation::exit(1);
		}
		molecule.setComponent(ensemble->getComponent(componentid - 1));  // Internally stored in vector starting at index 0
		double r[3];
		Coordinate3D sitePosition(xmlconfig, "coordinate");
		sitePosition.get(r);
		molecule.setr(0, r[0]);
		molecule.setr(1, r[1]);
		molecule.setr(2, r[2]);
		Quaternion q(1.0, 0., 0., 0.); /* orientation of molecules has to be set to a value other than 0,0,0,0! */
		molecule.setq(q);
		Log::global_log->info() << "[Basis] Adding molecule cid=" << componentid << ", (x,y,z)=(" << molecule.r(0) << "," << molecule.r(1) << "," << molecule.r(2) << ")" << std::endl;
		addMolecule(molecule);
	}
	xmlconfig.changecurrentnode(oldpath);
}

void Basis::addMolecule(const Molecule& molecule) {
	_molecules.push_back(molecule);
}

Molecule Basis::getMolecule(int i) {
	return _molecules[i];
}


size_t Basis::numMolecules(){
	return _molecules.size();
}
