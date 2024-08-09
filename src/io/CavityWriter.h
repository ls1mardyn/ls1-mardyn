#ifndef CAVITYWRITER_H_
#define CAVITYWRITER_H_

#include <string>

#include "plugins/PluginBase.h"

class CavityEnsemble;

/**@brief IO-plugin: Used to create instances of CavityEnsemble
 * and output the corresponding files containing the detected Cavity coordinates.
 */
class CavityWriter final : public PluginBase {
public:
    /** @brief Initialisation of all necessary field is done in readXML and init.*/
    CavityWriter() = default;

    ~CavityWriter() final = default;

    /** @brief reads in parameters for CavityEnsemble and output
     *<br>
     *	writefrequency: output timestep frequency AFTER initStatistics timesteps<br>
     *	componentid: can be multiple. Starts separate CavityEnsembles for each component.<br>
     *	radius: search radius around pseudo-molecules to look for neighbors with componentid=0/1...<br>
     *	maxNeighbours: number of neighbours allowed in radius for a pseudo-molecule to count as a cavity<br>
     *	Nx,Ny,Nz: grid numbers for grid of pseudo-molecules, sampling cavity properties at their grid position<br>
     *	<br>
     *	for non-overlapping coverage: radius*N_xyz = domainSize_xyz<br>
     *	<br>
     * \code{.xml}
     * <plugin name="CavityWriter" enabled="yes">
                <writefrequency>10</writefrequency>
                <componentid>0</componentid>
                <componentid>1</componentid>
                <radius>2.0</radius>
                <maxNeighbours>1</maxNeighbours>
                <Nx>40</Nx>
                <Ny>40</Ny>
                <Nz>40</Nz>
                <ControlVolume>
                    <x0> LowerBound (1.0) </x0>
                    <x1> UpperBound (2.0) </x1>
                    <y0> LowerBound (1.0) </y0>
                    <y1> UpperBound (2.0) </y1>
                    <z0> LowerBound (1.0) </z0>
                    <z1> UpperBound (2.0) </z1>
                </ControlVolume>
            </plugin>
     *\endcode
     *
     * @param xmlconfig
     */
    void readXML(XMLfileUnits &xmlconfig) final;

    /** @brief
     *
     * Set all necessary variables for each CavityEnsemble. One is used per Component.
     * Ensembles are assigned subDomains depending on the rank of the MPI process.
     * Also, if a smaller controlVolume is specified it gets passed here.
     *
     * @param particleContainer
     * @param domainDecomp
     * @param domain
     */
    void init(ParticleContainer *particleContainer,
              DomainDecompBase *domainDecomp, Domain *domain) final;

    /** @brief Method will be called first thing in a new timestep.
     * If the writeFrequency is reached and this is the first time, this will
     * set the required Quaternions for the pseudo-molecules.*/
    void beforeEventNewTimestep(
            ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
            unsigned long simstep
    ) final;

    /** @brief Method afterForces will be called after forcefields have been applied
     *  no sitewise Forces can be applied here.
     *
     * If the writeFrequency is reached, this will trigger a loop through all pseudo-molecules
     * to determine their status as either cavity or non-cavity and start the communication.
     */
    void afterForces(
            ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
            unsigned long simstep
    ) final;

    /** @brief
     *
     * Output the sampled cavity information from all ensembles to file.
     *
     * @param particleContainer
     * @param domainDecomp
     * @param domain
     * @param simstep
     */
    void endStep(
            ParticleContainer *particleContainer,
            DomainDecompBase *domainDecomp, Domain *domain,
            unsigned long simstep
    ) final;

    /** @brief
     *
     * Nothing done here. Necessary for inheritance from PluginBase.
     *
     * @param particleContainer
     * @param domainDecomp
     * @param domain
     */
    void finish(ParticleContainer *particleContainer,
                DomainDecompBase *domainDecomp, Domain *domain);

    std::string getPluginName() final {
        return std::string("CavityWriter");
    }

    static PluginBase *createInstance() { return new CavityWriter(); }

    std::map<unsigned, CavityEnsemble *> getMcav() { return _mcav; }

private:
    std::string _outputPrefix;
    unsigned long _writeFrequency;
    bool _appendTimestamp;
    bool _incremental;
    int _Nx = 0, _Ny = 0, _Nz = 0;
    int _maxNeighbors = 0;
    float _radius = 0.0f;
    std::map<unsigned, CavityEnsemble *> _mcav;
    double _controlVolume[6];
};

#endif /* CAVITYWRITER_H_ */
