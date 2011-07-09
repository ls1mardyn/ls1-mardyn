/*
 * StatisticsWriter.cpp
 *
 * @Date: 01.02.2011
 * @Author: eckhardw
 */

#include "StatisticsWriter.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/LinkedCells.h"
#include "particleContainer/Cell.h"
#include "molecules/Molecule.h"
#include "utils/Logger.h"

using namespace std;
using namespace Log;

StatisticsWriter::StatisticsWriter(unsigned int writeFrequency, const std::string& fileName,
		const LinkedCells& container)
: _outputPrefix(fileName),
  _writeFrequency(writeFrequency),
  _container(container)
{	}


StatisticsWriter::~StatisticsWriter() {
}


void StatisticsWriter::initOutput(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain) {

#ifndef NDEBUG
	if (static_cast<const void*>(&_container) != static_cast<const void*>(particleContainer)) {
		global_log->error() << "VTKGridWriter works only with PlottableLinkCells!" << std::endl;
		exit(1);
	}
#endif

	std::stringstream fileNameStream;
	fileNameStream << _outputPrefix << "_LC" << ".stat";
	ofstream outfile(fileNameStream.str().c_str());

	outfile << "Number of Cells (including halo): [" << _container._cellsPerDimension[0]
	        << "," << _container._cellsPerDimension[1] << "," << _container._cellsPerDimension[2] << "]" << endl;

	outfile << "HaloWidth in num cells: [" << _container._haloWidthInNumCells[0]
		        << "," << _container._haloWidthInNumCells[1] << "," << _container._haloWidthInNumCells[2] << "]" << endl;

	outfile << "Cell Width: [" << _container._cellLength[0] << "," << _container._cellLength[1] << "," << _container._cellLength[2] << "]" << endl;

	outfile << endl;
	outfile << "SizeOf(Molecule) = " << sizeof(Molecule) << endl;

	std::vector<bool> count;

	Molecule* tmpMolecule = particleContainer->begin();
	while (tmpMolecule != particleContainer->end()) {
		if (tmpMolecule->componentid() >= count.size()) {
			count.resize(tmpMolecule->componentid()+1, false);
		}

		if (count[tmpMolecule->componentid()] == false) {
			outfile << "Molecule(cid=" << tmpMolecule->componentid() << ") " << tmpMolecule->totalMemsize() << endl;
			count[tmpMolecule->componentid()] = true;
		}
		tmpMolecule = particleContainer->next();
	}

	outfile.close();

}


void StatisticsWriter::doOutput(
		ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
		Domain* domain, unsigned long simstep,
		std::list<ChemicalPotential>* lmu
) {

	if (simstep % _writeFrequency != 0) {
		return;
	}

	std::vector<int> count;

	if(domainDecomp->getRank()==0){

		for (unsigned int i = 0; i < _container._cells.size(); i++) {
			if (_container._cells[i].isHaloCell()) {
				continue;
			}

			unsigned int numberOfMolecules = _container._cells[i].getMoleculeCount();
			if (numberOfMolecules >= count.size()) {
				count.resize(numberOfMolecules+1, 0);
			}
			count[numberOfMolecules]++;

		}

		std::stringstream fileNameStream;
		fileNameStream << _outputPrefix << "_" << simstep << ".stat";
		ofstream outfile(fileNameStream.str().c_str());

		int sum = 0;
		for (unsigned int i = 0; i < count.size(); i++) {
			outfile << i << " " << count[i] << endl;
			sum += count[i] * i;
		}
		outfile.close();
		global_log->info() << "StatisticsWriter counted " << sum << " molecules" << endl;

	} else {
		global_log->warning() << "Writing Cell Occupancy only for root node!" << endl;
	}

#ifdef ENABLE_MPI
		// do a reduction
#endif
}


void StatisticsWriter::finishOutput(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain) {
}

