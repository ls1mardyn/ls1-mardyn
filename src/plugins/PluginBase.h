#ifndef PLUGINBASE_H_
#define PLUGINBASE_H_

#include <any>
#include <list>
#include <map>
#include <string>
#include <functional>

#include "utils/FunctionWrapper.h"

class ParticleContainer;
class DomainDecompBase;
class Domain;
class XMLfileUnits;


/** @todo Mark all parameters as const: output plugins should not modify the state of the simulation. */
/** @todo get rid of the domain parameter */
/** @todo clean up all classes implementing this interface */


/** @brief The PluginBase class provides the interface for any kind of output/plugin classes - called "(output) plugins".
 *
 * There are a lot of different things that one might want to write out during a simulation,
 * e.g. thermodynamic values, graphical information, time measurements, ...
 * For all cases in which this output happens regularly at the end of each time step
 * the PluginBase class provides a common interface. The interface provides access to
 * the most important data: the particle container and the domain decomposition.
 *
 * Of course, several plugins plugins will be needed in some cases. So the idea is, that
 * all available plugins are registered in the PluginFactory and initialized
 * in the Simulation at runtime as requested by the input file. The plugin will then be
 * called at the respective points in the simulation automatically.
 *
 * Therefore, each plugin has to implement at least the following five methods:
 * - init: will be called once in the beginning
 * - readXML: reads in the plugin configuration from config.xml
 * - endStep: will be called each time step
 * - finish: will be called at the end
 * - getPluginName: returning the output pulugin name
 * - createInstance: returning an instance object as follows
 * \code{.cpp}
 *   static PluginBase* createInstance() { return new MyPlugin(); }   // class name is MyPlugin
 * \endcode
 */
class PluginBase {
public:
    //! @brief Subclasses should use their constructur to pass parameters (e.g. filenames)
    PluginBase(){}

    virtual ~PluginBase(){}

    /** @brief Method init will be called at the begin of the simulation.
     *
     * This method will be called once at the begin of the simulation just
     * right before the main time step loop.
     * It can be used e.g. to open output files or initialize statistics.
     * @param particleContainer  particle container storing the (local) molecules
     * @param domainDecomp       domain decomposition in use
     * @param domain
     */
    virtual void init(ParticleContainer* particleContainer,
                            DomainDecompBase* domainDecomp, Domain* domain) = 0;

    /** @brief Method readXML will be called once for each plugin section in the input file.
     *
     * This method can be used to read in parameters from the corresponding plugin section in
     * the xml config file. The method will be called once after an instance of the plugin
     * is created.
     *
     * @note The same plugins may be specified multiple times in the xml config file.
     *       It is the responsibility of the plugin to handle this case in a propper way.
     *
     * The following xml object structure will be provided to the plugin:
     * \code{.xml}
       <plugin name="plugin name">
         <!-- options for the specific plugin -->
       </plugin>
       \endcode
     *
     * @param xmlconfig  section of the xml file
     */
    virtual void readXML(XMLfileUnits& xmlconfig) = 0;


    /** @brief Method will be called first thing in a new timestep. */
	virtual void beforeEventNewTimestep(
			ParticleContainer* /* particleContainer */, DomainDecompBase* /* domainDecomp */,
			unsigned long /* simstep */
	) {};

    /** @brief Method beforeForces will be called before forcefields have been applied
     * no alterations w.r.t. Forces shall be made here
     *
     */

    virtual void beforeForces(
            ParticleContainer* /* particleContainer */, DomainDecompBase* /* domainDecomp */,
            unsigned long /* simstep */
    ) {};

    /** @brief Method siteWiseForces will be called before forcefields have been applied
     *  alterations to sitewise forces and fullMolecule forces can be made here
     */

    virtual void siteWiseForces(
            ParticleContainer* /* particleContainer */, DomainDecompBase* /* domainDecomp */,
            unsigned long /* simstep */
    ) {};

    /** @brief Method afterForces will be called after forcefields have been applied
     *  no sitewise Forces can be applied here
     */
    virtual void afterForces(
            ParticleContainer* /* particleContainer */, DomainDecompBase* /* domainDecomp */,
            unsigned long /* simstep */
    ) {};


    // make pure virtual?
    /** @brief Method endStep will be called at the end of each time step.
     *
     * This method will be called every time step passing the simstep as an additional parameter.
     * It can be used e.g. to write per time step data to a file or perform additional computations.
     * @param particleContainer  particle container storing the (local) molecules
     * @param domainDecomp       domain decomposition in use
     * @param domain
     */
    virtual void endStep(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            Domain* domain, unsigned long simstep) = 0;

    /** @brief Method finish will be called at the end of the simulation
     *
     * This method will be called once at the end of the simulation.
     * It can be used e.g. to closing output files or writing final statistics.
     * @param particleContainer  particle container storing the (local) molecules
     * @param domainDecomp       domain decomposition in use
     * @param domain
     */
    virtual void finish(ParticleContainer* particleContainer,
                              DomainDecompBase* domainDecomp, Domain* domain) = 0;

    /** @brief return the name of the plugin */
    virtual std::string getPluginName()  = 0;

	/**
	 * Register callbacks to callbackMap.
	 * This allows to make functions of a plugin accessible to other plugins.
	 * New callbacks should be added to callbackMap.
	 * Example syntax:
	 * - register a function that returns a local value:
	 * \code
	 *   callbackMap["getMyLocalValue"] = [this] { return _myLocalValue; };
	 * \endcode
	 * - register a function that calls a local function and returns its return value:
	 * \code
	 *   callbackMap["callMyFunct"] = [this] { return myFunct(); };
	 * \endcode
	 * @param callbackMap Add callbacks to this map.
	 */
	virtual void registerCallbacks(std::map<std::string, FunctionWrapper>& callbackMap) {
		// Empty by default.
	}

	/**
	 * Save callbacks from the callbackMap locally.
	 * This allows a plugin to call functions from other plugins.
	 * Example syntax:
	 * - store a function that returns an unsigned long:
	 * \code
	 *   std::function<unsigned long(void)> myFunction;
	 *   myFunction = callbackMap.at("getSomeLocalValue").get<unsigned long>();
	 * \endcode
	 * - store a function that calls some function of another plugin with an input value (int):
	 * \code
	 *   std::function<void(int)> myFunction;
	 *   myFunction = callbackMap.at("doSth").get<void, int>();
	 * \endcode
	 * @param callbackMap Get callbacks from this map.
	 */
	virtual void accessAllCallbacks(const std::map<std::string, FunctionWrapper>& callbackMap) {
		// Empty by default.
	}
};

#endif /* PLUGINBASE_H */
