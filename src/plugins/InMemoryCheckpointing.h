/*
 * InMemoryCheckpointing.h
 *
 *  Created on: 3 Jul 2018
 *      Author: tchipevn
 */

#ifndef SRC_PLUGINS_INMEMORYCHECKPOINTING_H_
#define SRC_PLUGINS_INMEMORYCHECKPOINTING_H_

#include "PluginBase.h"
#include "molecules/MoleculeForwardDeclaration.h"
#include "molecules/Molecule.h"

#include <vector>

class Snapshot;

class InMemoryCheckpointing: public PluginBase {
public:
	InMemoryCheckpointing() {}
	virtual ~InMemoryCheckpointing() {}

	void init(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain) {}

    void readXML(XMLfileUnits& xmlconfig);

    /**
     * @brief restarting takes place here
     */
	void beforeEventNewTimestep(
			ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			unsigned long simstep
	);

    void beforeForces(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            unsigned long simstep
    ){}

    /**
     * @brief writing takes place here
     */
    void afterForces(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            unsigned long simstep
    ) {}

    void endStep(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            Domain* domain, unsigned long simstep);

    void finish(ParticleContainer* particleContainer,
                              DomainDecompBase* domainDecomp, Domain* domain){}

    std::string getPluginName() {
    	return std::string("InMemoryCheckpointing");
    }

	static PluginBase* createInstance() { return new InMemoryCheckpointing(); }

	class Snapshot {
	public:
		void addMolecule(const Molecule& m) {
			_molecules.push_back(m);
		}

		double getCurrentTime() const {
			return _currentTime;
		}

		void setCurrentTime(double currentTime) {
			_currentTime = currentTime;
		}

		const std::array<double, 3>& getBoxDims() const {
			return _boxDims;
		}

		void setBoxDims(const std::array<double, 3>& boxDims) {
			_boxDims = boxDims;
		}

		unsigned long getGlobalNumberOfMolecules() const {
			return _globalNumberOfMolecules;
		}

		void setGlobalNumberOfMolecules(unsigned long globalNumberOfMolecules) {
			_globalNumberOfMolecules = globalNumberOfMolecules;
		}

		double getTemperature() const {
			return _temperature;
		}

		void setTemperature(double temperature) {
			_temperature = temperature;
		}

		int getRank() const {
			return _rank;
		}

		void setRank(int rank) {
			_rank = rank;
		}

		const std::vector<Molecule>& getMolecules() const {
			return _molecules;
		}

		void clearMolecules() {
			_molecules.clear();
		}

	private:
		std::vector<Molecule> _molecules;
		double _currentTime;
		int _rank; // who do these molecules belong to?

		// the following fields are maybe unnecessary, but leaving them here now for consistency to written headers in file-checkpoints
		unsigned long _globalNumberOfMolecules;
		double _temperature; // maybe not necessary; for consistency to currently written headers
		std::array<double, 3> _boxDims; // maybe not necessary; for consistency to currently written headers

	};

private:
	Snapshot _snapshot; // make an std::vector eventually
	unsigned long _writeFrequency;
	unsigned long _restartAtIteration;
};

#endif /* SRC_PLUGINS_INMEMORYCHECKPOINTING_H_ */
