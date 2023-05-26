/*
 * VTKMoleculeWriter.h
 *
 * @Date: 24.08.2010
 * @Author: Wolfgang Eckhardt
 */

#ifndef VTKWRITER_H_
#define VTKWRITER_H_

#include <utility>

#include "plugins/PluginBase.h"
#include "io/vtk/VTKMoleculeWriterImplementation.h"

/**
 * This class is an implementation of the OutputBase for the VTK file format.
 *
 * @TODO Think about a way how to handle the setup of particleConainers, Domain, etc...
 *       Maybe some kind of factory?
 */
class VTKMoleculeWriter: public PluginBase {

private:

	/**
	 * Interval how often output shall be written.
	 */
	unsigned int _writeFrequency{50ul};
	/**
	 * When the first output Shall be written.
	 *
	 * This can be used to skip some iterations
	 * or offset the write frequency from the MPI rebalancing to better visualize tuning.
	 */
	unsigned int _writeFrequencyOffset{0ul};
	/**
	 * Whether there should also be an output for the initial state of the system.
	 * If _firstWriteIteration == 0 this has no effect.
	 */
	bool _writeInitialState{true};

	std::string _fileName{};

    /**
     * Counts how many particle containers exist in total.
     * Is initialized by multiple calls to init.
     * Added to support multiple particle containers.
     * */
    unsigned int _particleContainerCount{0ul};

    /**
     *
     * Added to support multiple particle containers
     * */
    unsigned int _endStepCount{0ul};


public:
	VTKMoleculeWriter() = default;

	VTKMoleculeWriter(unsigned int frequency, std::string name):
		_writeFrequency(frequency), _fileName(std::move(name)) {}

	virtual ~VTKMoleculeWriter() {}

    /**
     * With every call to init of this instance _particleContainerCount will rise by one.
     * This was added to support multiple particle containers.
     * */
	virtual void init(ParticleContainer *particleContainer,
                      DomainDecompBase *domainDecomp, Domain *domain);

    /**
     * Writes the final vtk file.
     * Depending on the amount of particle containers, multiple calls to endStep are needed to actually write the file.
     * The file only gets written if all particle containers have been handled.
     * */
	virtual void endStep(
            ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
            Domain *domain, unsigned long simstep
    );

	//! @todo document me!
	virtual void finish(ParticleContainer *particleContainer,
						DomainDecompBase *domainDecomp, Domain *domain);

	std::string getPluginName() {
		return std::string("VTKMoleculeWriter");
	}
	static PluginBase* createInstance() { return new VTKMoleculeWriter(); }

	void readXML(XMLfileUnits& xmlconfig);

private:
	void outputParallelVTKFile(unsigned int numProcs, unsigned long simstep,
			VTKMoleculeWriterImplementation& impl);

};

#endif /* VTKWRITER_H_ */
