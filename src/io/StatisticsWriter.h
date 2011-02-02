/*
 * StatisticsWriter.h
 *
 * @Date: 01.02.2011
 * @Author: eckhardw
 */

#ifndef STATISTICSWRITER_H_
#define STATISTICSWRITER_H_

#include "io/OutputBase.h"
#include "particleContainer/LinkedCells.h"

#include <string>
#include <fstream>
#include <vector>

/**
 * Plot statistics to a file. At the moment, only a histogramm showing the
 * occupancy of the cells is produces.
 *
 * x: number of Molecules
 * y: number of cells with the respective number of molecules
 */
class StatisticsWriter : public OutputBase{

public:

	StatisticsWriter(unsigned int writeFrequency, const std::string& fileName, const LinkedCells& container);

	virtual ~StatisticsWriter();


	virtual void doOutput(
			ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			Domain* domain, unsigned long simstep,
			std::list<ChemicalPotential>* lmu
	);


	//! Write basic properties of the linked cell container to a file
	virtual void initOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain) ;

	//! NOP
	virtual void finishOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);

private:

	//! prefix for the names of the output file
	std::string _outputPrefix;

	std::ofstream _resultStream;

	const unsigned int _writeFrequency;

	const LinkedCells& _container;

};

#endif /* STATISTICSWRITER_H_ */
