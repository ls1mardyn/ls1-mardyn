#ifndef SRC_IO_VECTORIZATIONTUNER_H_
#define SRC_IO_VECTORIZATIONTUNER_H_

/// An enum, that describes, whether the molecule count should be increased exponentially or linearly.
enum MoleculeCntIncreaseTypeEnum{
	linear,      //!< linear, the molecule count is increased linearly.
	exponential, //!< exponential, the molecule counts are distributed exponentially.
	both         //!< both, do linear and exponential measurements, linear measurements will stop after 32 molecules.
};


#include <vector>
#include <string>

#include "particleContainer/adapter/CellProcessor.h"
#include "particleContainer/adapter/FlopCounter.h"
#include "plugins/PluginBase.h"
#include "ensemble/EnsembleBase.h"
#include "parallel/LoadCalc.h"

class Component;
#include "particleContainer/ParticleCellForwardDeclaration.h"


/** @brief VectorizationTuner class.
 *
 * This class is used to get detailed information about the performance of the VectorizedCellProcessor.
 * For different scenarios, the performance is evaluated and output.
 * Later this could be used to actually use this class as a tuner, i.e. to use the best possible vectorization method for the actual computation.
 */
class VectorizationTuner: public PluginBase {

public:
	/**
	 * @brief Constructor of VectorizationTuner for the xml input mode.
	 */
    VectorizationTuner() = default;

	/**
	 * Destructor of the VectorizationTuner class.
	 */
	~VectorizationTuner() override = default;

	//documentation in PluginBase
	void init(ParticleContainer *particleContainer,
              DomainDecompBase *domainDecomp, Domain *domain) override;

	/** @brief Read in XML configuration for the VectorizationTuner.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <plugin type="VectorizationTuner">
	     <outputprefix>STRING</outputprefix>
	     <minmoleculecnt>INTEGER</minmoleculecnt>
	     <maxmoleculecnt>INTEGER</maxmoleculecnt>
	     <numRepetitionsMax>INTEGER</numRepetitionsMax>
	     <moleculecntincreasetype>linear OR exponential OR both</moleculecntincreasetype> <!--default: both-->
	     <mode>stdout OR file</mode> <!--print to stdout or file-->
	   </plugin>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig) override;

	static PluginBase* createInstance() { return new VectorizationTuner(); }

	//documentation in PluginBase, does nothing.
	void endStep(ParticleContainer * /*particleContainer*/, DomainDecompBase * /*domainDecomp*/,
                 Domain * /*domain*/, unsigned long /*simstep*/) override {}

	//documentation in PluginBase, does nothing.
	void finish(ParticleContainer * /*particleContainer*/,
				DomainDecompBase * /*domainDecomp*/, Domain * /*domain*/) override;

	//documentation in PluginBase.
	std::string getPluginName() override {
		return std::string("VectorizationTuner");
	}

	/**
	 * used in the KDDecomposition to fill the TunerTimes object
	 *
	 * particleNums stores the number of particles until which the tuner should measure the values
	 *
	 * useExistingFiles denotes whether the tuner is allowed to read the tuner values from existing tuner files (if they exist, otherwise they are still calculated)
	 *
	 * generateNewFiles denotes whether the tuner should write his measured values to files, which also overwrites old files.
	 * If both are true and tuner files are available then the values that were read are written again,
	 */
	void tune(std::vector<Component>& ComponentList, TunerLoad& times, std::vector<int> particleNums, bool generateNewFiles, bool useExistingFiles, bool allowMPIReduce);

private:
	/// The output prefix, that should be prefixed before the output files.
	std::string _outputPrefix{"MarDyn"};

	/// The minimal amount of molecules
	unsigned int _minMoleculeCnt{2};

	/// The maximal amount of molecules
	unsigned int _maxMoleculeCnt{512};

	/// An enum, that describes, whether the molecule count should be increased exponentially or linearly.
	MoleculeCntIncreaseTypeEnum _moleculeCntIncreaseType{both};

	/// The CellProcessor, that should be used to iterate over the cells.
	CellProcessor* _cellProcessor{nullptr};

	/// The cutoff radius
	double _cutoffRadius{};

	/// The cutoff Radius for the LJ potential
	double _LJCutoffRadius{};

	/// Maximal number of repetitions for each measurement.
	/// The actual number of repetitions is _numRepetitionsMax / numMols^2, but always at least 50.
	unsigned int _numRepetitionsMax{4000000};

	/// The cutoff radius
	static constexpr double _cutoffRadiusBig{5.};

	/// The cutoff Radius for the LJ potential
	static constexpr double _LJCutoffRadiusBig{5.};

	/// FlopCounter that utilizes a big cutoff radius
	std::unique_ptr<FlopCounter> _flopCounterBigRc;

	/// FlopCounter that utilizes a normal cutoff radius
	std::unique_ptr<FlopCounter> _flopCounterNormalRc;

	/// FlopCounter for zero cutoff radius
	std::unique_ptr<FlopCounter> _flopCounterZeroRc;

	/*
	 * Writes the given TunerTimes to a file
	 */
	void writeFile(const TunerLoad& vecs);

	/*
	 * Fills the TunerTimes object with values read from a file
	 *
	 * returns true on success
	 */
	bool readFile(TunerLoad& times);

	/**
	 * This function is the main routine of this plugin. Multiple simulations are started from here.
	 * @param ComponentList
	 * @param vcp
	 * @param fc
	 */
	void tune(std::vector<Component>& ComponentList);

	/**
	 * Performs multiple iterations of the selected simulation, that is also set up here.
	 * @param ComponentList
	 * @param vcp
	 * @param fc
	 * @param numMols
	 * @param gflopsOwn
	 * @param gflopsPair
	 */
	void iterate(std::vector<Component>& ComponentList, unsigned int numMols, double& gflopsOwnBig, double& gflopsPairBig, double& gflopsOwnNormal, double& gflopsPairNormalFace,
			double& gflopsPairNormalEdge, double& gflopsPairNormalCorner, double& gflopsOwnZero, double& gflopsPairZero);

	void iteratePair (unsigned int numRepetitions,
			ParticleCell& firstCell, ParticleCell& secondCell, double& time);

	void iterateOwn (unsigned int numRepetitions,
			ParticleCell& cell, double& time);

	/**
	 * @brief Calculation of the molecule interactions within a single cell.
	 * Initializes the Calculation and preprocesses the cell. Once it is preprocessed, multiple (numRepetitions) iterations are performed on that cell. Afterwards it is postprecessed.
	 * @param cp
	 * @param cell1
	 * @param numRepetitions
	 */
	void runOwn(CellProcessor& cp, ParticleCell& cell1, unsigned int numRepetitions);

	/** @brief Calculation of the molecule interactions between two neighboring cells.
	 * Initializes the Calculation and preprocesses both cells. Once they are preprocessed, multiple (numRepetitions) iterations are performed on the cells. Afterwards they are postprecessed.
	 * @param cp
	 * @param cell1
	 * @param cell2
	 * @param numRepetitions The amount of repetitions, that should be performed on the single cell.
	 */
	void runPair(CellProcessor& cp, ParticleCell& cell1, ParticleCell& cell2, unsigned int numRepetitions);

	/**
	 * @brief Initializes the Molecules in an equidistant mesh within the box.
	 *
	 * @param boxMin
	 * @param boxMax
	 * @param comp
	 * @param cell1
	 * @param cell2
	 */
	void initMeshOfMolecules(const double boxMin[3], const double boxMax[3], Component& comp, ParticleCell& cell1, ParticleCell& cell2);

	/**
	 * @brief Initializes the molecules uniformly randomly distributed within the box. The box is set using boxMin and boxMax.
	 * @param boxMin
	 * @param boxMax
	 * @param comp
	 * @param cell
	 * @param numMols
	 */
	void initUniformRandomMolecules(Component& comp, ParticleCell& cell, unsigned int numMols);

	void initUniformRandomMolecules(Component& comp1, Component& comp2, ParticleCell& cell, unsigned int numMolsFirst, unsigned int numMolsSecond);


	/**
	 * @brief Moves all molecules of the cell by the vector specified by direction.
	 * @param direction
	 * @param cell
	 * @param numMols
	 */
	void moveMolecules(double direction[3], ParticleCell& cell);

	/**
	 * @brief Initializes the molecules normally distributed within each cell.
	 * Be careful with this. I have no idea, whether this works correctly and what happens, if initialized particles are outside of the boundary of cells (note: Steffen Seckler)
	 *
	 * @param boxMin
	 * @param boxMax
	 * @param comp
	 * @param cell1
	 * @param cell2
	 * @param numMols
	 */
	void initNormalRandomMolecules(double boxMin[3], double boxMax[3], Component& comp, ParticleCell& cell1, ParticleCell& cell2, unsigned int numMols);

	/**
	 * Removes all molecules from the given cell.
	 *
	 * @param cell
	 */
	void clearMolecules(ParticleCell& cell);

	void initCells(ParticleCell& main, ParticleCell& face, ParticleCell& edge, ParticleCell& corner);

	void iterateOwn(unsigned int numRepetitions,
			ParticleCell& cell,
			double& gflopsPair, FlopCounter& flopCounter);
	void iteratePair(unsigned int numRepetitions,
			ParticleCell& firstCell, ParticleCell& secondCell,
			double& gflopsPair, FlopCounter& flopCounter);

	class VTWriterI {
	public:
		virtual ~VTWriterI() = default;
		virtual void initWrite(const std::string& outputPrefix, double cutoffRadius, double LJCutoffRadius,
				double cutoffRadiusBig, double LJCutoffRadiusBig)=0;
		virtual void writeHeader(const std::string& distributionTypeString)=0;
		virtual void write(unsigned int numMols, double gflopsOwnBig, double gflopsPairBig, double gflopsOwnNormal,
				double gflopsPairNormalFace, double gflopsPairNormalEdge, double gflopsPairNormalCorner,
				double gflopsOwnZero, double gflopsPairZero)=0;
		virtual void close()=0;
	protected:
		int _rank;
	};

	class VTWriter: public VTWriterI {
	public:
		void initWrite(const std::string& outputPrefix, double cutoffRadius, double LJCutoffRadius,
				double cutoffRadiusBig, double LJCutoffRadiusBig) override;
		void writeHeader(const std::string& distributionTypeString) override;
		void write(unsigned int numMols, double gflopsOwnBig, double gflopsPairBig, double gflopsOwnNormal,
				double gflopsPairNormalFace, double gflopsPairNormalEdge, double gflopsPairNormalCorner,
				double gflopsOwnZero, double gflopsPairZero) override;
		void close() override;
	private:
		ofstream _myfile;
	};

	class VTWriterStatistics: public VTWriterI {
	public:
		void initWrite(const std::string& outputPrefix, double cutoffRadius, double LJCutoffRadius,
				double cutoffRadiusBig, double LJCutoffRadiusBig) override;
		void writeHeader(const std::string& distributionTypeString) override;
		void write(unsigned int numMols, double gflopsOwnBig, double gflopsPairBig, double gflopsOwnNormal,
				double gflopsPairNormalFace, double gflopsPairNormalEdge, double gflopsPairNormalCorner,
				double gflopsOwnZero, double gflopsPairZero) override;
		void close() override;
	private:
		void writeStatistics(double input, const std::string& name);
	};

	std::unique_ptr<VTWriterI> vtWriter{new VTWriter()};
};

#endif  // SRC_IO_VECTORIZATIONTUNER_H_
