//Calculation of the Fluid Wall interaction by a function

#ifndef MIRROR_H_
#define MIRROR_H_

#include "PluginBase.h"

#include <string>
#include <map>
#include <list>

enum MirrorDirections {
	MD_LEFT_MIRROR = 0,
	MD_RIGHT_MIRROR = 1
};

class ParticleContainer;
class DomainDecompBase;
class Domain;
class ChemicalPotential;
class CavityEnsemble;
class Mirror : public PluginBase
{
public:
	// constructor and destructor
	Mirror();
	~Mirror();

	/** @brief Read in XML configuration for Mirror and all its included objects.
	 *
	 * The following XML object structure is handled by this method:
	 * \code{.xml}
		<plugin name="Mirror">
			<yPos> <float> </yPos>                     <!-- mirror position -->
			<forceConstant> <float> </forceConstant>   <!-- strength of redirection -->
			<direction> <int> </direction>             <!-- 0|1 , i.e. left |<-- or right -->| mirror -->
		</plugin>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig) override;

	void init(ParticleContainer *particleContainer,
			  DomainDecompBase *domainDecomp, Domain *domain) override;

    /** @brief Method afterForces will be called after forcefields have been applied
     *
     * make pure Virtual ?
     */
    void afterForces(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            unsigned long simstep
    ) override;

	void endStep(
			ParticleContainer *particleContainer,
			DomainDecompBase *domainDecomp, Domain *domain,
			unsigned long simstep) override {}

	void finish(ParticleContainer *particleContainer,
				DomainDecompBase *domainDecomp, Domain *domain) override {}

	std::string getPluginName() override {return std::string("Mirror");}
	static PluginBase* createInstance() {return new Mirror();}

private:
		void VelocityChange(ParticleContainer* particleContainer);
		void readNormDistr();

private:
	double _yPos;
	double _forceConstant;
	int _direction;
	bool _bReflect;
	struct NormMB{
		bool enabled;
		struct NormFnames{
			std::string vxz;
			std::string vy;
		} fname;
		std::list<double> vxz;
		std::list<double> vy;
	} _norm;
};

#endif /*MIRROR_H_*/
