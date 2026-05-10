/*
 * FlopRateWriter.h
 *
 *  Created on: 7 Feburary 2018
 *      Author: kruegener
 */

#ifndef SRC_UTILS_TESTPLUGIN_H_
#define SRC_UTILS_TESTPLUGIN_H_

#include "PluginBase.h"

class TestPlugin: public PluginBase {
public:
    TestPlugin() {}
    ~TestPlugin() {}

    //! @brief will be called at the beginning of the simulation
    //!
    //! Some OutputPlugins will need some initial things to be done before
    //! the output can start, e.g. opening some files. This method will
    //! be called once at the beginning of the simulation (see Simulation.cpp)
    void init(ParticleContainer* particleContainer,
                    DomainDecompBase* domainDecomp, Domain* domain) {
        Log::global_log->debug()  << "[TESTPLUGIN] TESTPLUGIN INIT" << std::endl;
    }

    void readXML(XMLfileUnits& xmlconfig) {
        Log::global_log -> debug() << "[TESTPLUGIN] READING XML" << std::endl;
    }

    /** @brief Method will be called first thing in a new timestep. */
	void beforeEventNewTimestep(
			ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			unsigned long simstep, bool signalled
	) {
		Log::global_log -> debug() << "[TESTPLUGIN] BEFORE EVENT NEW TIMESTEP" << std::endl;
        if (signalled)
		    Log::global_log -> debug() << "[TESTPLUGIN] SIGUSR1 RECEIVED" << std::endl;
	};

    /** @brief Method beforeForces will be called before forcefields have been applied
            *
            * make pure Virtual ?
    */

    void beforeForces(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            unsigned long simstep
    ) {
        Log::global_log -> debug() << "[TESTPLUGIN] BEFORE FORCES" << std::endl;
    }

    /** @brief Method afterForces will be called after forcefields have been applied
     *
     * make pure Virtual ?
     */
    void afterForces(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            unsigned long simstep
    ) {
        Log::global_log->debug()  << "[TESTPLUGIN] TESTPLUGIN AFTER FORCES" << std::endl;
    }

    /** @brief Method endStep will be called at the end of each time step.
     *
     * This method will be called every time step passing the simstep as an additional parameter.
     * It can be used e.g. to write per time step data to a file or perform additional computations.
     * @param particleContainer  particle container storing the (local) molecules
     * @param domainDecomp       domain decomposition in use
     * @param domain
     */
    void endStep(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            Domain* domain, unsigned long simstep
    ) {
        Log::global_log->debug()  << "[TESTPLUGIN] ENDSTEP" << std::endl;
    }

    /** @brief Method finalOutput will be called at the end of the simulation
     *
     * This method will be called once at the end of the simulation.
     * It can be used e.g. to closing output files or writing final statistics.
     * @param particleContainer  particle container storing the (local) molecules
     * @param domainDecomp       domain decomposition in use
     * @param domain
     */
    void finish(ParticleContainer* particleContainer,
                        DomainDecompBase* domainDecomp, Domain* domain) {
        Log::global_log->debug()  << "[TESTPLUGIN] FINISHING" << std::endl;
    }

    /** @brief return the name of the plugin */
    std::string getPluginName()  {
        Log::global_log->debug()  << "[TESTPLUGIN] GETTING NAME" << std::endl;
        return "TestPlugin";}
    static PluginBase* createInstance() {
        Log::global_log->debug()  << "[TESTPLUGIN] CREATE INSTANCE" << std::endl;
        return new TestPlugin(); }

private:

    std::ofstream _fileStream;
    unsigned long _writeFrequency;
    std::string _outputPrefix;
};

#endif /* SRC_UTILS_TESTPLUGIN_H_ */
