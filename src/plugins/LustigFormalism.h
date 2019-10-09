/*
 * LustigFormalism.h
 *
 *  Created on: 18.05.2016
 *      Author: mheinen
 */

#ifndef LUSTIGFORMALISM_H_
#define LUSTIGFORMALISM_H_

#include "PluginBase.h"
#include "molecules/Comp2Param.h"
#include <string>
#include <cstdint>
#include <vector>
#include <array>

class Domain;
class DomainDecompBase;

// restart
struct RestartCtrl
{
	bool isRestart;
	std::string filenameSums;
	uint64_t timestep;
};

class LustigFormalism : public PluginBase
{
public:
	LustigFormalism();
	virtual ~LustigFormalism();

	/** @brief Method init will be called at the begin of the simulation.
	 *
	 * This method will be called once at the begin of the simulation just
	 * right before the main time step loop.
	 * It can be used e.g. to open output files or initialize statistics.
	 * @param particleContainer  particle container storing the (local) molecules
	 * @param domainDecomp       domain decomposition in use
	 * @param domain
	 */
	void init(ParticleContainer* particleContainer,
							DomainDecompBase* domainDecomp, Domain* domain) override;

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
	void readXML(XMLfileUnits& xmlconfig) override;


	/** @brief Method will be called first thing in a new timestep. */
	void beforeEventNewTimestep(
			ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			unsigned long simstep
	) override;

	/** @brief Method beforeForces will be called before forcefields have been applied
	 * no alterations w.r.t. Forces shall be made here
	 *
	 */

	void beforeForces(
			ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			unsigned long simstep
	) override;

	/** @brief Method siteWiseForces will be called before forcefields have been applied
	 *  alterations to sitewise forces and fullMolecule forces can be made here
	 */

	void siteWiseForces(
			ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			unsigned long simstep
	) override;

	/** @brief Method afterForces will be called after forcefields have been applied
	 *  no sitewise Forces can be applied here
	 */
	void afterForces(
			ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			unsigned long simstep
	) override;


	// make pure virtual?
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
			Domain* domain, unsigned long simstep) override;

	/** @brief Method finish will be called at the end of the simulation
	 *
	 * This method will be called once at the end of the simulation.
	 * It can be used e.g. to closing output files or writing final statistics.
	 * @param particleContainer  particle container storing the (local) molecules
	 * @param domainDecomp       domain decomposition in use
	 * @param domain
	 */
	void finish(ParticleContainer* particleContainer,
							  DomainDecompBase* domainDecomp, Domain* domain) override;

	/** @brief return the name of the plugin */
	std::string getPluginName() override {return std::string("LustigFormalism");}
	static PluginBase* createInstance() {return new LustigFormalism();}

	void SetWriteFreq(unsigned long nWriteFreq, unsigned long nStart, unsigned long nStop, unsigned long nWriteFreqSums)
	{_nWriteFreq = nWriteFreq; _nStart = nStart; _nStop = nStop; _nWriteFreqSums = nWriteFreqSums;}
	void InitRestart();
	unsigned long GetStart() {return _nStart;}
	unsigned long GetStop()  {return _nStop;}
	void InitNVT(Domain* domain, unsigned long N, double V, double T, double cutoffRadiusLJ);
	void Init(const double& U6, const double& dUdV, const double& d2UdV2);
	void InitWidom(const double& DU, const double& T);
	void InitMSD(const unsigned long& nID, double const * dr);
	void CalcGlobalValues(DomainDecompBase* domainDecomp);
	void CalcDerivatives();
	void WriteHeader(DomainDecompBase* domainDecomp, Domain* domain);
	void WriteData(DomainDecompBase* domainDecomp, unsigned long simstep);
	void TriggerNewSimstep(unsigned long simstep){_bSimstepTrigger = true; _simstep = simstep;}
	void ResetTriggerNewSimstep(){_bSimstepTrigger = false;}
	bool IsNewSimstep(){return _bSimstepTrigger;}
	void SetNumWidomTests(unsigned int nNumWidomTests){_nNumWidomTestsGlobal = nNumWidomTests;}
	void EventConfigurationSampled(Domain* domain);
	void StoreDomainDecompPointer(DomainDecompBase* domainDecomp) {_domainDecomp = domainDecomp;}
	void StoreTimestepLength(double dTimestepLength) { _Dt_ts = dTimestepLength;}
	void StartSamplingMSD(){_bSampleMSD = true;}
	bool SamplingStartedMSD(){return _bSampleMSD;}

private:
	// reset local values
	void ResetSums();
	void ResetLocalValues();
	void InitDatastructures();
	void InitSums();


private:
	unsigned long _N;
	double _InvN;
	double _V;
	double _InvV;
	double _mInvV;
	double _InvV2;
	double _T;
	unsigned long _tsBufferIndex;
	unsigned long _nWriteFreq;
	unsigned long _nWriteFreqSums;
	unsigned long _nStart;
	unsigned long _nStop;
	unsigned long _nNumConfigs;
	double _rho;
	double _v;
	double _v2;
	double _beta;
	double _beta2;
	double _beta3;

	// local
	std::vector<double> _ULocal;
	std::vector<double> _dUdVLocal;
	std::vector<double> _d2UdV2Local;

	// global
	std::vector<double> _UGlobal;
	std::vector<double> _dUdVGlobal;
	std::vector<double> _d2UdV2Global;

	std::vector<double> _U2Global;
	std::vector<double> _U3Global;
	std::vector<double> _dUdV2Global;
	std::vector<double> _UdUdVGlobal;
	std::vector<double> _U2dUdVGlobal;
	std::vector<double> _UdUdV2Global;
	std::vector<double> _Ud2UdV2Global;

	// sums
	double _UGlobalSum;
	double _U2GlobalSum;
	double _U3GlobalSum;
	double _dUdVGlobalSum;
	double _d2UdV2GlobalSum;
	double _dUdV2GlobalSum;
	double _UdUdVGlobalSum;
	double _U2dUdVGlobalSum;
	double _UdUdV2GlobalSum;
	double _Ud2UdV2GlobalSum;

	// LRC
	double _U_LRC;
	double _dUdV_LRC;
	double _d2UdV2_LRC;
	double _dU_LRC;

	// ideal
	double _A00i;
	double _A10i;
	double _A01i;
	double _A20i;
	double _A11i;
	double _A02i;
	double _A30i;
	double _A21i;
	double _A12i;

	// residual
	double _A00r;
	double _A10r;
	double _A01r;
	double _A20r;
	double _A11r;
	double _A02r;
	double _A30r;
	double _A21r;
	double _A12r;
	// thermodynamic properties
	double _Cv;

	//! parameter streams for each possible pair of molecule-types
	Comp2Param _comp2params;

	// mu
	std::vector<double> _WidomEnergyLocal;
	std::vector<double> _WidomEnergyGlobal;
	double _WidomEnergyGlobalSum;
	unsigned int _nNumWidomTestsLocal;
	unsigned int _nNumWidomTestsGlobal;
	double _mu_res;

	bool _bSimstepTrigger;
	unsigned long _simstep;

	DomainDecompBase* _domainDecomp;

	// restart
	RestartCtrl _restartCtrl;
//	std::string _strRestartFilenameSums;
//	unsigned long _nRestartTimestep;

	// MSD
	//double _D;           // self-diffusion coefficient
	double _MSD;          // Mean Squared Displacement
	double _Dt_ts;       // time step length
	double _dInv6NDt_ts;  // 1 / (6*N*Dt_ts)
	std::vector<std::array<double,3> > _dDisplacementVecLocal;
	std::vector<std::array<double,3> > _dDisplacementVecGlobal;
	double _dDisplacementGlobalSum;
	bool _bSampleMSD;

};  // class LustigFormalism



#endif /* LUSTIGFORMALISM_H_ */
