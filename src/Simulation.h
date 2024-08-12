#ifndef SIMULATION_H_
#define SIMULATION_H_

#include <memory>
#include <any>

#include "io/TimerProfiler.h"
#include "thermostats/VelocityScalingThermostat.h"
#include "utils/FixedSizeQueue.h"
#include "utils/FunctionWrapper.h"
#include "utils/SysMon.h"
#include "utils/Random.h"



// plugins
#include "plugins/PluginFactory.h"

#if !defined (SIMULATION_SRC) or defined (IN_IDE_PARSER)
class Simulation;
/** Global pointer to the simulation object. Do not use directly. Instead use the reference. */
extern Simulation* global_simulation;
#endif

class ParticleInsertion;
class Ensemble;



/** Reference to the global simulation object */
#define _simulation (*global_simulation)

#include <list>
#include <vector>
#include <string>
#include <io/TaskTimingProfiler.h>

#ifdef STEEREO
class SteereoSimSteering;
class SteereoCouplingSim;
#endif

class Domain;
class ParticleContainer;
class ParticlePairsHandler;
class CellProcessor;
class Integrator;
class PluginBase;
class DomainDecompBase;
class InputBase;
class RDF;
class FlopCounter;
class LongRangeCorrection;
class Homogeneous;
class Planar;
class TemperatureControl;
class SphericalTemperatureControl;
class MemoryProfiler;

// by Stefan Becker
const int VELSCALE_THERMOSTAT = 1;

namespace bhfmm {
class FastMultipoleMethod;
} // bhfmm


/** @brief Controls the simulation process
 *  @author Martin Bernreuther <bernreuther@hlrs.de> et al. (2010)
 *
 * Simulation parameters are provided via a xml config file or can be set directly via the corresponding methods.
 */
class Simulation {
private:
	Simulation(Simulation &simulation);
	Simulation& operator=(Simulation &simulation);
	void resetVirials();
	void updateForces();

public:
	/** Instantiate simulation object */
	Simulation();

	/** destruct simulation object */
	~Simulation();

	/** @brief Read in XML configuration for simulation and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <simulation>
	     <integrator type="STRING"><!-- see Integrator class documentation --></integrator>
	     <run>
	       <production>
	         <steps>INTEGER</steps>
	       </production>
	       <equilibration>
	         <steps>INTEGER</steps>
	       </equilibration>
	       <currenttime>DOUBLE</currenttime>
	     </run>
	     <ensemble type=STRING> <!-- see Ensemble class documentation --></ensemble>
	     <algorithm>
	       <cutoffs>
	          <radiusLJ>DOUBLE</radiusLJ>
	       </cutoffs>
	       <electrostatic type='ReactionField'>
	         <epsilon>DOUBLE</epsilon>
	       </electrostatic>
	       <datastructure type=STRING><!-- see ParticleContainer class documentation --></datastructure>
	       <parallelisation type=STRING><!-- see DomainDecompBase class documentation -->
	         <timerForLoad>STRING</timerForLoad><!-- Timer to use as load. Requires valid timer name! -->
	         <timerForLoad_AveragingLength>UINT</timerForLoad_AveragingLength><!-- Defines how many time measurements should be averaged as input for the load balancing.-->
	       </parallelisation>
	       <thermostats>
	         <thermostat type='VelocityScaling' componentId=STRING><!-- componentId can be component id or 'global' -->
	           <temperature>DOUBLE</temperature>
	         </thermostat>
	       </thermostats>
	     </algorithm>
	     <output>
	       <outputplugin name=STRING enabled="yes|no"><!-- see OutputBase class and specific plugin documentation --></outputplugin>
	     </output>
	     <plugin name=STRING enabled="yes|no" (default yes)><!-- see PluginBase class and specific plugin documentation --></plugin>
	   </simulation>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig);

	/** @brief Terminate simulation with given exit code.
	 *
	 * The exit method takes care over the right way to terminate the application in a correct way
	 * for the different parallelization schemes. e.g. terminating other processes in MPI parallel
	 * execution mode.
	 */
	static void exit(int exitcode);

	/** @brief process configuration file
	 *
	 * calls initConfigXML
	 * @param[in]  filename filename of the input file
	 */
	void readConfigFile(std::string filename);

	/** @brief Opens given XML file and reads in parameters for the simulaion.
	 * @param[in]  inputfilename filename of the XML input file
	 */
	void initConfigXML(const std::string& inputfilename);

	/** @brief calculate all values for the starting timepoint
	 *
	 * After the input file has been read in, only the information which is
	 * directly in the inputfile is known at time step zero.
	 * This includes e.g. Molecule position and velocities, but not forces
	 * which act on the molecules or macroscopic values at the potential.
	 * This method has to ensure, that all values are available at
	 * time step zero
	 */
	void prepare_start();

	/** @brief Controls the main loop of the simulation.
	 *
	 * precondition for this method is that initialize() has been called.
	 * The main loop calls all methods which have to be called in each
	 * iteration step, e.g. initializes the molecule Container, calculates
	 * the forces, ... \n
	 * For the integration, a seperate integration object is used.
	 * This method is written in a way that should make it possible to use
	 * different integrators without having to change anything.
	 * Whenever something is done which might make it necessary for an integrator
	 * to do something, the integrator is informed about that using
	 * the integrators corresponding method (see integrator documentation)
	 * It follows a coarse outline of what has to be done:
	 * - Do all output that has to be done in each time step
	 * - Inform the iterator that all information is available (new timestep)
	 * - As the position of the molecules changed, the domain decomposition and
	 *     the datastructure have to be updated
	 * - calculate forces (and afterwards delete halo)
	 * - Inform the iterator that forces for the new positions have been calculated
	 * - The iterator should now have finished everything to be done in this time step,
	 *     so the macroscopic values can be calculated
	 * - velocity and angular momentum have to be scaled
	 */
	void simulate();

	/** @brief call plugins every nth-simstep
	 *
	 * The present method serves as a redirection to the actual plugins.
	 * That includes
	 * a) particular plugin objects included in _plugins,
	 *
	 * @param[in]  simstep timestep of the plugins
	 */
	void pluginEndStepCall(unsigned long simstep);

	/** @brief clean up simulation */
	void finalize();

	/** The following things have to be done here:
	 * - bring all molecules to the corresponding processes (including copies for halo)
	 * - update the caches of the molecules
	 * - update the ParticleContainer
	 */
	void updateParticleContainerAndDecomposition(double lastTraversalTime, bool useTimers);

	/**
	 * Performs both the decomposition and the cell traversal in an overlapping way.
	 * The overlapping is needed to speed up the overall computation. The order of cells
	 * traversed will be different, than for the non-overlapping case, slightly different results are possible.
	 */
	void performOverlappingDecompositionAndCellTraversalStep(double etime);

	/**
	 * Set the private _domainDecomposition variable to a new pointer.
	 * @param domainDecomposition the new va
	 */
	void setDomainDecomposition(DomainDecompBase* domainDecomposition);

	/** Return a reference to the domain decomposition used in the simulation */
	DomainDecompBase& domainDecomposition() { return *_domainDecomposition; }

	/** Get pointer to the domain */
	Domain* getDomain() { return _domain; }

	/** Get pointer to the integrator */
	Integrator* getIntegrator() {return _integrator; }

	/** Get pointer to the molecule container */
	ParticleContainer* getMoleculeContainer() { return _moleculeContainer; }

	/** Set the number of time steps to be performed in the simulation */
	void setNumTimesteps( unsigned long steps ) { _numberOfTimesteps = steps; }
	/** Get the number of time steps to be performed in the simulation */
	unsigned long getNumTimesteps() { return _numberOfTimesteps; }
	/** Get initial number of steps */
	unsigned long getNumInitTimesteps() { return _initSimulation; }
	/** Get the number of the actual time step currently processed in the simulation. */
	unsigned long getSimulationStep() { return _simstep; }
	/** Set Loop Time Limit in seconds */
	void setLoopAbortTime(double time) {
		Log::global_log->info() << "Max loop-abort-time set: " << time << "\n";
		_wallTimeEnabled = true;
		_maxWallTime = time;
	}

	double getcutoffRadius() const { return _cutoffRadius; }
	void setcutoffRadius(double cutoffRadius) { _cutoffRadius = cutoffRadius; }
	double getLJCutoff() const { return _LJCutoffRadius; }
	void setLJCutoff(double LJCutoffRadius) { _LJCutoffRadius = LJCutoffRadius; }
	unsigned long getTotalNumberOfMolecules() const;

	/** @brief Temperature increase factor function during automatic equilibration.
	 * @param[in]  current simulation time step
	 *
	 * The present version of Mardyn provides the option of
	 * equilibrating the system at an increased temperature. This has,
	 * of course, the advantage that the equilibration is accelerated.
	 * In a special case that is of particular interest to a relevant
	 * subset of the developers, namely MD simulation of nucleation,
	 * it also reduces the amount of clusters formed during
	 * equilibration. Then, the actual nucleation process can be
	 * investigated for the equilibrated system.
	 *
	 * Equilibration at an increased temperature is nothing new, in
	 * principle it could be considered a variant of simulated
	 * annealing. For the purpose of accelerating the relaxation,
	 * starting from an unrealistic initial state, it was e.g.
	 * proposed in the thesis of M. Kreitmeir.
	 *
	 * The key parameter to regulating the equilibration is
	 * this->_initCanonical. The system is gradually heated and cooled
	 * again, and from this->_initCanonical on the temperature is
	 * constant (i.e. this method returns 1.0).
	 *
	 * Thus, if temperature increase is to be avoided, the user should
	 * set the "initCanonical" parameter to "0".
	 */
	double Tfactor(unsigned long simstep);

	void initCanonical(unsigned long t) { this->_initCanonical = t; }
	void initGrandCanonical(unsigned long t) { this->_initGrandCanonical = t; }
	void initStatistics(unsigned long t) { this->_initStatistics = t; }
	unsigned long getInitStatistics() const { return this->_initStatistics; }

	void setSimulationTime(double curtime) { _simulationTime = curtime; }
	void advanceSimulationTime(double timestep) { _simulationTime += timestep; }
	double getSimulationTime() { return _simulationTime; }

	void setEnsemble(Ensemble *ensemble) { _ensemble = ensemble; }
	Ensemble* getEnsemble() { return _ensemble; }

	std::shared_ptr<MemoryProfiler> getMemoryProfiler() {
		return _memoryProfiler;
	}

	TimerProfiler* timers() {
		return &_timerProfiler;
	}

	//! get Planck constant
	double getH() {return h;}
	//! set Planck constant
	void setH(double h_extern) {h = h_extern;}

private:


	double _simulationTime; /**< Simulation time t in reduced units */


	/** maximum id of particles */
	/** @todo remove this from the simulation class */
	unsigned long _maxMoleculeId;

	/** maximum distance at which the forces between two molecules still have to be calculated. */
	double _cutoffRadius;

	/** LJ cutoff (may be smaller than the RDF/electrostatics cutoff) */
	double _LJCutoffRadius;

	/** A thermostat can be specified to account for the directed
	 * motion, which means that only the undirected kinetic energy is
	 * maintained constant. In many cases, this is an essential to
	 * simulating the system in the desired way, most notably for
	 * applying MD to nanoscopic fluid dynamics.
	 *
	 * Obviously, that requires calculating the directed velocity.
	 * Since the velocity is often fairly constant over time, it is
	 * not necessary to redetermine it every time the thermostat is
	 * applied. The property this->_collectThermostatDirectedVelocity
	 * regulates how often (in number of time steps) the directed
	 * velocity is evaluated.
	 */
	unsigned _collectThermostatDirectedVelocity;

	//! by Stefan Becker: the Type of the thermostat(velocity scaling or Andersen or...)
	//! appropriate tokens stored as constants at the top of this file
	int _thermostatType;

	unsigned long _numberOfTimesteps;
public:
	unsigned long getNumberOfTimesteps() const;

private:
	/**< Number of discrete time steps to be performed in the simulation */

	unsigned long _simstep;             /**< Actual time step in the simulation. */

	/** initial number of steps */
	unsigned long _initSimulation;
	/** step number for the end of the configurational equilibration */
	unsigned long _initCanonical;
	/** step number for activation of the muVT ensemble */
	unsigned long _initGrandCanonical;
	/** step number for activation of all sorts of statistics */
	unsigned long _initStatistics;

	Ensemble* _ensemble;

	/** Datastructure for finding neighbours efficiently */
	ParticleContainer* _moleculeContainer;

	/** Handler describing what action is to be done for each particle pair */
	ParticlePairsHandler* _particlePairsHandler;

	/** New cellhandler, which will one day replace the particlePairsHandler here completely. */
	CellProcessor* _cellProcessor;

	/** module which handles the domain decomposition */
	DomainDecompBase* _domainDecomposition;

	/** numerical solver for the particles equations of motion */
	Integrator* _integrator;

	/** all macroscopic (local and global) information */
	Domain* _domain;

	/** responsible for reading in the phasespace (header+data) */
	InputBase* _inputReader;

	/** prefix for the names of all output files */
	std::string _outputPrefix;

	//! number of time steps after which the canceling is carried outline
	unsigned _momentumInterval;

	//! random number generator
	Random _rand;

	/** Long Range Correction */
	LongRangeCorrection* _longRangeCorrection;

	/** Temperature Control (Slab Thermostat) */
	TemperatureControl* _temperatureControl;

	/** SphericalTemperatureControl (Spherical Thermostat) */
	SphericalTemperatureControl* _sphericalTemperatureControl;
	
	/** The Fast Multipole Method object */
	bhfmm::FastMultipoleMethod* _FMM;

	/** manager for all timers in the project except the MarDyn main timer */
	TimerProfiler _timerProfiler;

	//! used to get information about the memory consumed by the process and the overall system.
	std::shared_ptr<MemoryProfiler> _memoryProfiler;

#ifdef TASKTIMINGPROFILE
	/** Used to track what thread worked on which task for how long and plot it **/
	TaskTimingProfiler* _taskTimingProfiler;
#endif
public:
#ifdef TASKTIMINGPROFILE
	TaskTimingProfiler* getTaskTimingProfiler(){
		return _taskTimingProfiler;
	}
#endif
	//! computational time for one execution of traverseCell
	double getAndResetOneLoopCompTime() {
		if(_loopCompTimeSteps==0){
			return 1.;
		}
		double t = _loopCompTime/_loopCompTimeSteps;
		_loopCompTime = 0.;
		_loopCompTimeSteps = 0;
		return t;
	}
	void setOutputPrefix( std::string prefix ) { _outputPrefix = prefix; }
	void setOutputPrefix( char *prefix ) { _outputPrefix = std::string( prefix ); }
	std::string getOutputPrefix() { return _outputPrefix; }

	void enableFinalCheckpoint() { _finalCheckpoint = true; }
	void disableFinalCheckpoint() { _finalCheckpoint = false; }

	void useLegacyCellProcessor() { _legacyCellProcessor = true; }

	void enableMemoryProfiler() {
		_memoryProfiler = std::make_shared<MemoryProfiler>();
		_memoryProfiler->registerObject(reinterpret_cast<MemoryProfilable**>(&_moleculeContainer));
		_memoryProfiler->registerObject(reinterpret_cast<MemoryProfilable**>(&_domainDecomposition));
	}

	void setForcedCheckpointTime(double time) { _forced_checkpoint_time = time; }

	/** initialize all member variables with a suitable value */
	void initialize();

	/** @brief get plugin
	 * @return pointer to the plugin if it is active, otherwise nullptr
	 */
	PluginBase* getPlugin(const std::string& name);

	std::list<PluginBase*>* getPluginList(){
		return &_plugins;
	}

	/** Global energy log */
	void initGlobalEnergyLog();
	void writeGlobalEnergyLog(const double& globalUpot, const double& globalT, const double& globalPressure);

	CellProcessor *getCellProcessor() const;

	/** @brief Refresh particle IDs to continuous numbering starting at zero*/
	void refreshParticleIDs();

	/** @brief Checks if Simsteps or MaxWallTime are reached */
	bool keepRunning();

private:

	/**
	 * Parse everything in the XML path </mardyn/simulation/options/option>
	 * and store it in _miscOptions.
	 * @param xmlconfig
	 */
	void parseMiscOptions(XMLfileUnits& xmlconfig);

	/// the timer used for the load calculation.
	Timer* _timerForLoad{nullptr};

	Timer _timeFromStart;
	double _maxWallTime = -1;
	bool _wallTimeEnabled = false;

	/** Enable final checkpoint after simulation run. */
	bool _finalCheckpoint;

	/** use legacyCellProcessor instead of vectorizedCellProcessor */
	bool _legacyCellProcessor = false;

	/**
	 * Specifies whether to use overlapping p2p (peer-to-peer) communication or not.
	 * If false: overlapping is only performed for unpacking and packing of particles.
	 * If true: overlapping of the leaving and halo particle exchange with the actual force calculation is performed.
	 */
	bool _overlappingP2P {false};

	/** List of plugins to use */
	std::list<PluginBase*> _plugins;

	/** Map of all call backs.
	 * The key is the name of the callback.
	 * Each element contains a std::function object.
	 * Please check the specific plugins for the actual signature of the function and use an appropriate any_cast!
	 */
	std::map<std::string, FunctionWrapper> _callbacks;

	VelocityScalingThermostat _velocityScalingThermostat;

	/** This is Planck's constant. (Required for the Metropolis
	 * criterion which is used for the grand canonical ensemble).
	 * Of course, what is actually specified here is not the value
	 * of h in and of itself, since that is a universal constant, but
	 * the value of h AS EXPRESSED IN REDUCED UNITS, i.e. for the
	 * internal use of the program.
	 */
	double h;

	/** Time after which the application should write a checkpoint in seconds. */
	double _forced_checkpoint_time;

	//! computational time for one loop
	double _loopCompTime;

	int _loopCompTimeSteps;

	/** Global energy log */
	unsigned long _nWriteFreqGlobalEnergy;
	std::string _globalEnergyLogFilename;

	/**
	 * Mechanism to easily add arbitrary options via the XML path <mardyn/simulation/options/option>
	 * Any option is identified via the label "name"
	 *
	 * @note Currently, only bool options are supported.
	 *
	 * Example:
	 *   <option name="refreshIDs">true</option>
	 */
	std::map<std::string, bool> _miscOptions;

	FixedSizeQueue<double> _lastTraversalTimeHistory;


public:
	/*** @brief Performs one time step. Is called as many times as timesteps are needed for the full simulation.
	 * 
	 * In library mode, this function is used to externally control the simulation. In normal usage, this function is 
	 * only called in the simulate() function, within the simulation loop.
	 * preSimLoopSteps() needs to be called before this function is called.
	*/
	void simulateOneTimestep();

	/*** @brief Performs necessary setup before the first timestep, which includes initializing all timers. */
	void preSimLoopSteps();

	/*** @brief Performs immediate cleanup after the simulation ends.
	 * 
	 * The final checkpoint is written, the final plugin call is made, and all timers are stopped. 
	 * Relevant timer information is printed.
	 */
	void postSimLoopSteps();

	/*** @brief Used to flag simulation as finished and to proceed to cleanup. Only used when ls1 is compiled and used
	 * as a library, and not in normal usage.
	 * 
	 * Normally, the simulation is controlled entirely by the simulate() function, which determines whether the 
	 * simulation is complete by checking the keepRunning() function.
	 * However in library mode, it is expected that the external code will manually control the simulation by calling
	 * simulateOneTimeStep(). As such, the keepRunning() function will never be called, and the exit conditions will 
	 * never be checked on ls1's side. \n
	 * This function acts as a way for the external code to perform some cleanup steps, before ending the simulation
	 * on the ls1 side. For now, it only sets the simulationDone boolean to true, which is used as a check within
	 * the postSimLoopSteps() function.
	 * */
	void markSimAsDone();

private:
	// stores the timing info for the previous load. This is used for the load calculation and the rebalancing.
	double previousTimeForLoad = 0.;
	/*** @brief Act as safeguards for the preSimLoopSteps(), simulateOneTimestep() and postSimLoopSteps() functions.
	 * 
	 * These three functions are public, since they need to be reachable by external code when ls1 is compiled
	 * as a library. However it is possible to call them out of order, causing unexpected behaviour. As such,
	 * these bool values help keep track of the simulation and can be used to verify the state the simulation is in. 
	 */
	bool preSimLoopStepsDone = false, simulationDone = false, postSimLoopStepsDone = false;

};
#endif /*SIMULATION_H_*/
