/*
 * FlopRateWriter.h
 *
 *  Created on: 7 Feburary 2018
 *      Author: kruegener
 */

#ifndef SRC_IO_FLOPRATEWRITER_H_
#define SRC_IO_FLOPRATEWRITER_H_

#include "PluginBase.h"

class testPlugin: public PluginBase {
public:
    testPlugin() {}
    ~testPlugin() {}

    //! @brief will be called at the beginning of the simulation
    //!
    //! Some OutputPlugins will need some initial things to be done before
    //! the output can start, e.g. opening some files. This method will
    //! be called once at the beginning of the simulation (see Simulation.cpp)
    void init(ParticleContainer* particleContainer,
                    DomainDecompBase* domainDecomp, Domain* domain) {};

    void readXML(XMLfileUnits& /*xmlconfig*/) {};

    /** @brief Method beforeForces will be called before forcefields have been applied
            *
            * make pure Virtual ?
    */

    virtual void beforeForces(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            unsigned long simstep
    ) {}

    /** @brief Method afterForces will be called after forcefields have been applied
     *
     * make pure Virtual ?
     */
    virtual void afterForces(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            unsigned long simstep
    ) {
        global_log->info()  << "[TESTPLUGIN] TESTPLUGIN AFTER FORCES" << std::endl;
    }

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
            Domain* domain
    ) {};

    /** @brief Method finalOutput will be called at the end of the simulation
     *
     * This method will be called once at the end of the simulation.
     * It can be used e.g. to closing output files or writing final statistics.
     * @param particleContainer  particle container storing the (local) molecules
     * @param domainDecomp       domain decomposition in use
     * @param domain
     */
    virtual void finish(ParticleContainer* particleContainer,
                        DomainDecompBase* domainDecomp, Domain* domain) {};

    /** @brief return the name of the plugin */
    virtual std::string getPluginName()  {return "testPlugin";};
    static PluginBase* createInstance() { return new testPlugin(); }

private:
    void setPrefix(double f_in, double& f_out, char& prefix) const;

    std::ofstream _fileStream;
    unsigned long _writeFrequency;
    std::string _outputPrefix;
};

#endif /* SRC_IO_FLOPRATEWRITER_H_ */
