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
#include <map>

using namespace Log;

RDF::RDF(double intervalLength, unsigned int bins, const std::vector<Component>& components) :
	_intervalLength(intervalLength),
	_bins(bins),
	_numberOfComponents(components.size()),
	_components(components),
	_doCollectSiteRDF(false),
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

	this->_localSiteDistribution = new unsigned long****[_numberOfComponents];
	this->_globalSiteDistribution = new unsigned long****[_numberOfComponents];
	this->_globalAccumulatedSiteDistribution = new unsigned long****[_numberOfComponents];
	for(unsigned i = 0; i < _numberOfComponents; i++) {
		this->_globalCtr[i] = 0;
		this->_globalAccumulatedCtr[i] = 0;

		this->_localDistribution[i] = new unsigned long*[_numberOfComponents-i];
		this->_globalDistribution[i] = new unsigned long*[_numberOfComponents-i];
		this->_globalAccumulatedDistribution[i] = new unsigned long*[_numberOfComponents-i];

		unsigned ni = _components[i].numSites();
		this->_localSiteDistribution[i] = new unsigned long***[_numberOfComponents-i];
		this->_globalSiteDistribution[i] = new unsigned long***[_numberOfComponents-i];
		this->_globalAccumulatedSiteDistribution[i] = new unsigned long***[_numberOfComponents-i];

		for(unsigned k=0; i+k < _numberOfComponents; k++) {
			this->_localDistribution[i][k] = new unsigned long[bins];
			this->_globalDistribution[i][k] = new unsigned long[bins];
			this->_globalAccumulatedDistribution[i][k] = new unsigned long[bins];

			for(unsigned l=0; l < bins; l++) {
				this->_localDistribution[i][k][l] = 0;
				this->_globalDistribution[i][k][l] = 0;
				this->_globalAccumulatedDistribution[i][k][l] = 0;
			}

			unsigned nj = _components[i+k].numSites();
			if(ni+nj > 2) {
				this->_doCollectSiteRDF = true;

				this->_localSiteDistribution[i][k] = new unsigned long**[ni];
				this->_globalSiteDistribution[i][k] = new unsigned long**[ni];
				this->_globalAccumulatedSiteDistribution[i][k] = new unsigned long**[ni];

				for(unsigned m=0; m < ni; m++) {
					this->_localSiteDistribution[i][k][m] = new unsigned long*[nj];
					this->_globalSiteDistribution[i][k][m] = new unsigned long*[nj];
					this->_globalAccumulatedSiteDistribution[i][k][m] = new unsigned long*[nj];
					for(unsigned n=0; n < nj; n++) {
						this->_localSiteDistribution[i][k][m][n] = new unsigned long[bins];
						this->_globalSiteDistribution[i][k][m][n] = new unsigned long[bins];
						this->_globalAccumulatedSiteDistribution[i][k][m][n] = new unsigned long[bins];
						for(unsigned l=0; l < bins; l++) {
							this->_localSiteDistribution[i][k][m][n][l] = 0;
							this->_globalSiteDistribution[i][k][m][n][l] = 0;
							this->_globalAccumulatedSiteDistribution[i][k][m][n][l] = 0;
							// cout << "init " << i << "\t" << k << "\t" << m << "\t" << n << "\t" << l << "\n";
						}
					}
				}
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
		unsigned ni = _components[i].numSites();
		for(unsigned k=0; i+k < _numberOfComponents; k++) {
			unsigned nj = _components[i+k].numSites();
			for(unsigned l=0; l < _bins; l++) {
				this->_globalAccumulatedDistribution[i][k][l] += this->_globalDistribution[i][k][l];
				if(ni + nj > 2) {
					for(unsigned m=0; m < ni; m++) {
						for(unsigned n=0; n < nj; n++) {
							this->_globalAccumulatedSiteDistribution[i][k][m][n][l] += this->_globalSiteDistribution[i][k][m][n][l];
						}
					}
				}
			}
		}
	}
}

void RDF::collectRDF(DomainDecompBase* dode) {
	// Communicate component-component RDFs
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

	// communicate site-site RDFs
	for(unsigned i=0; i < _numberOfComponents; i++) {
		unsigned ni = _components[i].numSites();
		for(unsigned k=0; i+k < _numberOfComponents; k++) {
			unsigned nj = _components[i+k].numSites();
			if(ni+nj > 2) {
				dode->collCommInit(_bins * ni * nj);
				for(unsigned l=0; l < _bins; l++) {
					for(unsigned m=0; m < ni; m++) {
						for(unsigned n=0; n < nj; n++) {
							// d1 = this->_localSiteDistribution[i][k][m][n][l];
							dode->collCommAppendUnsLong(_localSiteDistribution[i][k][m][n][l]);
						}
					}
				}

				// domainDecomp->reducevalues(&dZ, &d0, &d1, &d2);
				dode->collCommAllreduceSum();

				for(unsigned l=0; l < _bins; l++) {
					for(unsigned m=0; m < ni; m++) {
						for(unsigned n=0; n < nj; n++) {
							// if(!this->_localRank) this->_globalSiteDistribution[i][k][m][n][l] = d1;
							_globalSiteDistribution[i][k][m][n][l] = dode->collCommGetUnsLong();
						}
					}
				}
				dode->collCommFinalize();

			}
		}
	}
}


void RDF::reset() {
	_numberOfRDFTimesteps = 0;
	for(unsigned i=0; i < _numberOfComponents; i++) {
		_globalCtr[i] = 0;

		unsigned ni = _components[i].numSites();
		for(unsigned k=0; i+k < _numberOfComponents; k++) {
			unsigned nj = _components[i+k].numSites();
			for(unsigned l=0; l < _bins; l++) {
				_localDistribution[i][k][l] = 0;
				_globalDistribution[i][k][l] = 0;

				if(ni+nj > 2) {
					for(unsigned m=0; m < ni; m++) {
						for(unsigned n=0; n < nj; n++) {
							this->_localSiteDistribution[i][k][m][n][l] = 0;
							this->_globalSiteDistribution[i][k][m][n][l] = 0;
						}
					}
				}
			}
		}
	}
}


void RDF::setOutputTimestep(unsigned int timestep) {
	_RDFOutputTimesteps = timestep;
}

void RDF::setOutputPrefix(std::string prefix) {
	_RDFOutputPrefix = prefix;
}


void RDF::doOutput(DomainDecompBase* domainDecomposition, const Domain* domain, unsigned long simStep) {
	if (simStep > 0 && simStep % _RDFOutputTimesteps == 0) {
		collectRDF(domainDecomposition);

		
		if( domainDecomposition->getRank() == 0 ) {
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

	unsigned ni = _components[i].numSites();
	unsigned nj = _components[j].numSites();

	double V = domain->getGlobalVolume();
	double N_i = _globalCtr[i] / (double)_numberOfRDFTimesteps;
	double N_Ai = _globalAccumulatedCtr[i] / (double)_accumulatedNumberOfRDFTimesteps;
	double N_j = _globalCtr[j] / (double)_numberOfRDFTimesteps;
	double N_Aj = _globalAccumulatedCtr[j] / (double)_accumulatedNumberOfRDFTimesteps;
	double rho_i = N_i / V;
	double rho_Ai = N_Ai / V;
	double rho_j = N_j / V;
	double rho_Aj = N_Aj / V;

	rdfout.precision(5);
	rdfout << "# r\tcurr.{loc, int}\taccu.{loc, int}\t\tdV\tNpair(curr.)\tNpair(accu.)\t\tnorm(curr.)\tnorm(accu.)";
	if(ni+nj > 2) {
		for(unsigned m=0; m < ni; m++) {
			rdfout << "\t";
			for(unsigned n=0; n < nj; n++) {
				rdfout << "\t(" << m << ", " << n << ")_curr{loc, int}   (" << m << ", " << n << ")_accu{loc, int}";
			}
		}
	}
	rdfout << "\n";
	rdfout << "# \n# ctr_i: " << _globalCtr[i] << "\n# ctr_j: " << _globalCtr[j]
	       << "\n# V: " << V << "\n# _universalRDFTimesteps: " << _numberOfRDFTimesteps
	       << "\n# _universalAccumulatedTimesteps: " << _accumulatedNumberOfRDFTimesteps
	       << "\n# rho_i: " << rho_i << " (acc. " << rho_Ai << ")"
	       << "\n# rho_j: " << rho_j << " (acc. " << rho_Aj << ")"
	       << "\n# \n";

	double N_pair_int = 0.0;
	double N_Apair_int = 0.0;
	std::map< unsigned, std::map<unsigned, double> > Nsite_pair_int, Nsite_Apair_int;
	for(unsigned l=0; l < this->_bins; l++) {
		double rmin = l * _intervalLength;
		double rmid = (l+0.5) * _intervalLength;
		double rmax = (l+1.0) * _intervalLength;
		double r3min = rmin*rmin*rmin;
		double r3max = rmax*rmax*rmax;
		double dV = (4.0 / 3.0) * M_PI * (r3max - r3min);

		double N_pair = _globalDistribution[i][j-i][l] / (double) _numberOfRDFTimesteps;
		N_pair_int += N_pair;
		double N_Apair = _globalAccumulatedDistribution[i][j-i][l]
		                                                               / (double) _accumulatedNumberOfRDFTimesteps;
		N_Apair_int += N_Apair;
		double N_pair_norm = 0.0;
		double N_Apair_norm = 0.0;
		double N_pair_int_norm = 0.0;
		double N_Apair_int_norm = 0.0;

		if(i == j)
		{
			N_pair_norm = 0.5*N_i*(N_i-1.0) * dV/V;
			N_Apair_norm = 0.5*N_Ai*(N_Ai-1.0) * dV/V;
			N_pair_int_norm = 0.5*N_i*(N_i-1.0) * 4.1887902*r3max/V;
			N_Apair_int_norm = 0.5*N_Ai*(N_Ai-1.0) * 4.1887902*r3max/V;
		}
		else
		{
			N_pair_norm = N_i*N_j * dV/V;
			N_Apair_norm = N_Ai*N_Aj * dV/V;
			N_pair_int_norm = N_i*N_j * 4.1887902*r3max/V;
			N_Apair_int_norm = N_Ai*N_Aj * 4.1887902*r3max/V;
		}

		rdfout << rmid << "\t" << N_pair/N_pair_norm << " " << N_pair_int/N_pair_int_norm
				<< "\t" << N_Apair/N_Apair_norm << " " << N_Apair_int/N_Apair_int_norm
				<< "\t\t" << dV << "\t" << N_pair << "\t" << N_Apair
				<< "\t\t" << N_pair_norm << "\t" << N_Apair_norm;

		if(ni+nj > 2)
		{
			for(unsigned m=0; m < ni; m++)
			{
				rdfout << "\t";
				for(unsigned n=0; n < nj; n++)
				{
					double p = _globalSiteDistribution[i][j-i][m][n][l] / (double)_numberOfRDFTimesteps;
					Nsite_pair_int[m][n] += p;
					double ap = _globalAccumulatedSiteDistribution[i][j-i][m][n][l] / (double)_accumulatedNumberOfRDFTimesteps;
					Nsite_Apair_int[m][n] += ap;
					rdfout << "\t" << p/N_pair_norm << " " << Nsite_pair_int[m][n]/N_pair_int_norm
							<< "   " << ap/N_Apair_norm << " " << Nsite_Apair_int[m][n]/N_Apair_int_norm;
				}
			}
		}
		rdfout << "\n";
	}
	rdfout.close();
}
