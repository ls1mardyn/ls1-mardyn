/*
 * RDF.cpp
 *
 * @Date: 11.02.2011
 * @Author: eckhardw
 */

#include "RDF.h"

#include "Domain.h"
#include "molecules/Component.h"
#include "parallel/DomainDecompBase.h"
#include "utils/Logger.h"

#include <sstream>
#include <fstream>

using namespace Log;

RDF::RDF(double intervalLength, unsigned int bins, unsigned int numberOfComponents) :
	_intervalLength(intervalLength),
	_bins(bins),
	_numberOfComponents(numberOfComponents),
	_RDFOutputTimesteps(25000),
	_RDFOutputPrefix("out")
{

	_numberOfRDFTimesteps = 0;
	_accumulatedNumberOfRDFTimesteps = 0;
	_maxDistanceSquare = intervalLength*intervalLength*bins*bins;

	_globalCtr = new unsigned long[_numberOfComponents];
	_globalAccumulatedCtr = new unsigned long[_numberOfComponents];
	_localDistribution = new unsigned long**[_numberOfComponents];
	_globalDistribution = new unsigned long**[_numberOfComponents];
	_globalAccumulatedDistribution = new unsigned long**[_numberOfComponents];
	for(unsigned i = 0; i < _numberOfComponents; i++) {
		this->_globalCtr[i] = 0;
		this->_globalAccumulatedCtr[i] = 0;

		this->_localDistribution[i] = new unsigned long*[_numberOfComponents-i];
		this->_globalDistribution[i] = new unsigned long*[_numberOfComponents-i];
		this->_globalAccumulatedDistribution[i] = new unsigned long*[_numberOfComponents-i];

		for(unsigned k=0; i+k < _numberOfComponents; k++) {
			this->_localDistribution[i][k] = new unsigned long[bins];
			this->_globalDistribution[i][k] = new unsigned long[bins];
			this->_globalAccumulatedDistribution[i][k] = new unsigned long[bins];

			for(unsigned l=0; l < bins; l++) {
				this->_localDistribution[i][k][l] = 0;
				this->_globalDistribution[i][k][l] = 0;
				this->_globalAccumulatedDistribution[i][k][l] = 0;
			}
		}
	}
}


RDF::~RDF() {
	for(unsigned i = 0; i < _numberOfComponents; i++) {
		for(unsigned k=0; i+k < _numberOfComponents; k++) {
			delete[] _localDistribution[i][k];
			delete[] _globalDistribution[i][k];
			delete[] _globalAccumulatedDistribution[i][k];
		}
		delete[] _localDistribution[i];
		delete[] _globalDistribution[i];
		delete[] _globalAccumulatedDistribution[i];
	}

	delete[] _globalCtr;
	delete[] _globalAccumulatedCtr;
	delete[] _localDistribution;
	delete[] _globalDistribution;
	delete[] _globalAccumulatedDistribution;
}

void RDF::accumulateNumberOfMolecules(std::vector<Component>& components) const {
	const int num_components = components.size();
	for (int i = 0; i < num_components; i++) {
		_globalCtr[i] += components[i].getNumMolecules();
	}
}


void RDF::accumulateRDF() {
	if(0 >= _numberOfRDFTimesteps) return;
	_accumulatedNumberOfRDFTimesteps += _numberOfRDFTimesteps;
	for(unsigned i=0; i < _numberOfComponents; i++) {
		this->_globalAccumulatedCtr[i] += this->_globalCtr[i];
		for(unsigned k=0; i+k < _numberOfComponents; k++) {
			for(unsigned l=0; l < _bins; l++) {
				this->_globalAccumulatedDistribution[i][k][l] += this->_globalDistribution[i][k][l];
			}
		}
	}
}

void RDF::collectRDF(DomainDecompBase* dode) {
	dode->collCommInit(_bins * _numberOfComponents * (_numberOfComponents+1)/2);

	for(unsigned i=0; i < _numberOfComponents; i++) {
		for(unsigned k=0; i+k < _numberOfComponents; k++) {
			for(unsigned l=0; l < _bins; l++) {
				dode->collCommAppendUnsLong(_localDistribution[i][k][l]);
			}
		}
	}

	dode->collCommAllreduceSum();
	for(unsigned i=0; i < _numberOfComponents; i++) {
		for(unsigned k=0; i+k < _numberOfComponents; k++) {
			for(unsigned l=0; l < _bins; l++) {
				_globalDistribution[i][k][l] = dode->collCommGetUnsLong();
			}
		}
	}
	dode->collCommFinalize();
}


void RDF::reset() {
	_numberOfRDFTimesteps = 0;
	for(unsigned i=0; i < _numberOfComponents; i++) {
		_globalCtr[i] = 0;
		for(unsigned k=0; i+k < _numberOfComponents; k++) {
			for(unsigned l=0; l < _bins; l++) {
				_localDistribution[i][k][l] = 0;
				_globalDistribution[i][k][l] = 0;
			}
		}
	}
}


void RDF::setOutputTimestep(unsigned int timestep) {
	_RDFOutputTimesteps = timestep;
}

void RDF::setOutputPrefix(std::string& prefix) {
	_RDFOutputPrefix = prefix;
}


void RDF::doOutput(DomainDecompBase* domainDecomposition, const Domain* domain, unsigned long simStep) {
	if (simStep > 0 && simStep % _RDFOutputTimesteps == 0) {
		collectRDF(domainDecomposition);

		if( domain->ownrank() == 0 ) {
			accumulateRDF();
			for (unsigned i = 0; i < _numberOfComponents; i++) {
				for (unsigned j = i; j < _numberOfComponents; j++) {
					std::ostringstream osstrm;
					osstrm << _RDFOutputPrefix << "_" << i << "-" << j << ".";
					osstrm.fill('0');
					osstrm.width(9);
					osstrm << std::right << simStep << ".rdf";
					writeToFile(domain, osstrm.str().c_str(), i, j);
					osstrm.str("");
					osstrm.clear();
				}
			}
		}
		reset();
	}
}


void RDF::writeToFile(const Domain* domain, const char* prefix, unsigned i, unsigned j) const {

	std::ofstream rdfout(prefix);
	if (!rdfout) {
		global_log->warning() << "COULD NOT OPEN FILE" << prefix << std::endl;
		return;
	}

	double V = domain->getGlobalVolume();
	double N_i = _globalCtr[i] / _numberOfRDFTimesteps;
	double N_Ai = _globalAccumulatedCtr[i] / _accumulatedNumberOfRDFTimesteps;
	double N_j = _globalCtr[j] / _numberOfRDFTimesteps;
	double N_Aj = _globalAccumulatedCtr[j] / _accumulatedNumberOfRDFTimesteps;
	double rho_i = N_i / V;
	double rho_Ai = N_Ai / V;
	double rho_j = N_j / V;
	double rho_Aj = N_Aj / V;

	rdfout.precision(5);
	rdfout << "# r\tcurr.\taccu.\t\tdV\tNpair(curr.)\tNpair(accu.)\t\tnorm(curr.)\tnorm(accu.)\n";
	rdfout << "# \n# ctr_i: " << _globalCtr[i] << "\n# ctr_j: " << _globalCtr[j]
	       << "\n# V: " << V << "\n# _universalRDFTimesteps: " << _numberOfRDFTimesteps
	       << "\n# _universalAccumulatedTimesteps: " << _accumulatedNumberOfRDFTimesteps
	       << "\n# rho_i: " << rho_i << " (acc. " << rho_Ai << ")"
	       << "\n# rho_j: " << rho_j << " (acc. " << rho_Aj << ")"
	       << "\n# \n";

	for(unsigned l=0; l < this->_bins; l++) {
		double rmin = l * _intervalLength;
		double rmid = (l+0.5) * _intervalLength;
		double rmax = (l+1.0) * _intervalLength;
		double r3min = rmin*rmin*rmin;
		double r3max = rmax*rmax*rmax;
		double dV = (4.0 / 3.0) * M_PI * (r3max - r3min);

		unsigned long N_pair = _globalDistribution[i][j-i][l] / _numberOfRDFTimesteps;
		unsigned long N_Apair = _globalAccumulatedDistribution[i][j-i][l]
		                                                               / _accumulatedNumberOfRDFTimesteps;
		double N_pair_norm;
		double N_Apair_norm;
		if(i == j)
		{
			N_pair_norm = 0.5*N_i*rho_i*dV;
			N_Apair_norm = 0.5*N_Ai*rho_Ai*dV;
		}
		else
		{
			N_pair_norm = N_i*rho_j*dV;
			N_Apair_norm = N_Ai*rho_Aj*dV;
		}

		rdfout << rmid << "\t" << N_pair/N_pair_norm
				<< "\t" << N_Apair/N_Apair_norm
				<< "\t\t" << dV << "\t" << N_pair << "\t" << N_Apair
				<< "\t\t" << N_pair_norm << "\t" << N_Apair_norm << "\n";
	}
	rdfout.close();
}
