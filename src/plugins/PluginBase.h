#ifndef PLUGINBASE_H_
#define PLUGINBASE_H_

#include <list>
#include <map>
#include <string>

#include "../ensemble/GrandCanonical.h"
#include "../ensemble/CavityEnsemble.h"

class ParticleContainer;
class DomainDecompBase;
class Domain;
class XMLfileUnits;


/** @todo Mark all parameters as const: output plugins should not modify the state of the simulation. */
/** @todo get rid of the domain parameter */
/** @todo get rid of lmu and mcav as well, if possible. */
/** @todo clean up all classes implementing this interface */


/** @brief The outputBase class provides the interface for any kind of output classes - called "output plugins".
 *
 * There are a lot of different things that one might want to write out during a simulation,
 * e.g. thermodynamic values, graphical information, time measurements, ...
 * For all cases in which this output happens regularly at the end of each time step
 * the outputBase class provides a common interface. The interface provides access to
 * the most important data: the particle container and the domain decomposition.
 *
 * Of course, several output plugins will be needed in some cases. So the idea is, that
 * all available output plugins are registered in the OutputPluginFactory and initialized
 * in the Simulation at runtime as requested by the input file. The plugin will then be
 * called at the respective points in the simulation automatically.
 *
 * Therefore, each output plugin has to implement at least the following five methods:
 * - initOutput: will be called once in the beginning
 * - doOutput: will be called each time step
 * - finishOutput: will be called at the end
 * - getPluginName: returning the output pulugin name
 * - createInstance: returning an instance object as follows
 * \code{.cpp}
 *   static OutputBase* createInstance() { return new MyPlugin(); }   // class name is MyPlugin
 * \endcode
 *
 * Optionally each plugin can provide its own implemenation of the following method to
 * read in a corresoonding section from the xml config file:
 * - readXML: read in xml section from config file
 */
class PluginBase {
public:
    //! @brief Subclasses should use their constructur to pass parameters (e.g. filenames)
    PluginBase(){}

    virtual ~PluginBase(){}

    /** @brief Method initOutput will be called at the begin of the simulation.
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

    /** @brief Method readXML will be called once for each output plugin section in the input file.
     *
     * This method can be used to read in parameters from the corresponding output plugin section in
     * the xml config file. The method will be called once after an instance of the output plugin
     * is created.
     *
     * @note The same output plugins may be specified multiple times in the xml config file.
     *       It is the responsibility of the output plugin to handle this case in a propper way.
     *
     * The following xml object structure will be provided to the output plugin:
     * \code{.xml}
       <outputplugin name="plugin name">
         <!-- options for the specific plugin -->
       </outputplugin>
       \endcode
     *
     * @param xmlconfig  section of the xml file
     */
    virtual void readXML(XMLfileUnits& xmlconfig) = 0;


    /** @brief Method beforeForces will be called before forcefields have been applied
     *
     * make pure Virtual ?
     */

    virtual void beforeForces(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            unsigned long simstep
    ) {};

    /** @brief Method afterForces will be called after forcefields have been applied
     *
     * make pure Virtual ?
     */
    virtual void afterForces(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            unsigned long simstep
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

    /** @brief Method finalOutput will be called at the end of the simulation
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
};

#endif /* PLUGINBASE_H */
