#include "io/ResultWriter.h"

#include "Domain.h"
#include "Simulation.h"
#include "parallel/DomainDecompBase.h"
#include "utils/Logger.h"

#include <fstream>

using Log::global_log;

void ResultWriter::readXML(XMLfileUnits& xmlconfig) {
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	global_log->info() << "[ResultWriter] Write frequency: " << _writeFrequency << endl;

	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	global_log->info() << "[ResultWriter] Output prefix: " << _outputPrefix << endl;

	xmlconfig.getNodeValue("writeprecision", _writePrecision);
	global_log->info() << "[ResultWriter] Write precision: " << _writePrecision << endl;
}

void ResultWriter::init(ParticleContainer * /*particleContainer*/,
                        DomainDecompBase *domainDecomp, Domain * /*domain*/) {

	// Only main rank writes data to file
	if(domainDecomp->getRank() == 0) {
		const string resultfile(_outputPrefix+".res");
		std::ofstream resultStream;
		resultStream.open(resultfile.c_str(), std::ios::out);
		resultStream << std::setw(10) << "#step"
			<< std::setw(_writePrecision+15) << "time"
			<< std::setw(_writePrecision+15) << "U_pot_avg"
			<< std::setw(_writePrecision+15) << "U_kin_avg"
			<< std::setw(_writePrecision+15) << "U_kinTrans_avg"
			<< std::setw(_writePrecision+15) << "U_kinRot_avg"
			<< std::setw(_writePrecision+15) << "p_avg"
			<< std::setw(_writePrecision+15) << "c_v"
			<< std::setw(_writePrecision+15) << "N"
			<< endl;
		resultStream.close();
	}
}

void ResultWriter::endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
                           unsigned long simstep) {

	const uint64_t globalNumMolecules = domain->getglobalNumMolecules(true, particleContainer, domainDecomp);

	_uPot_acc += domain->getGlobalUpot();
	_p_acc += domain->getGlobalPressure();

	double uKinLocal = 0.0F;
	double uKinGlobal = 0.0F;
	double uKinTransLocal = 0.0F;
	double uKinTransGlobal = 0.0F;
	double uKinRotLocal = 0.0F;
	double uKinRotGlobal = 0.0F;
	#if defined(_OPENMP)
	#pragma omp parallel reduction(+:uKinLocal,uKinTransLocal,uKinRotLocal)
	#endif
	{
		for (auto moleculeIter = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
			 moleculeIter.isValid(); ++moleculeIter) {
			uKinLocal += moleculeIter->U_kin();
			uKinTransLocal += moleculeIter->U_trans();
			uKinRotLocal += moleculeIter->U_rot();
		}
	}
#ifdef ENABLE_MPI
    MPI_Allreduce(&uKinLocal, &uKinGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&uKinTransLocal, &uKinTransGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&uKinRotLocal, &uKinRotGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    uKinGlobal = uKinLocal;
	uKinTransGlobal = uKinTransLocal;
	uKinRotGlobal = uKinRotLocal;
#endif
	_uKin_acc += uKinGlobal;
	_uKinTrans_acc += uKinTransGlobal;
	_uKinRot_acc += uKinRotGlobal;

	_numSamples++;

	if ((simstep % _writeFrequency == 0) and (simstep > 0UL)) {
		// Only main rank writes data to file
		if (domainDecomp->getRank() == 0) {
			const string resultfile(_outputPrefix+".res");
			std::ofstream resultStream;
			resultStream.open(resultfile.c_str(), std::ios::app);
			resultStream << std::setw(10) << simstep << std::setw(_writePrecision+15) << std::scientific << std::setprecision(_writePrecision) << _simulation.getSimulationTime()
				<< std::setw(_writePrecision+15) << std::scientific << std::setprecision(_writePrecision) << _uPot_acc/_numSamples
				<< std::setw(_writePrecision+15) << std::scientific << std::setprecision(_writePrecision) << _uKin_acc/_numSamples
				<< std::setw(_writePrecision+15) << std::scientific << std::setprecision(_writePrecision) << _uKinTrans_acc/_numSamples
				<< std::setw(_writePrecision+15) << std::scientific << std::setprecision(_writePrecision) << _uKinRot_acc/_numSamples
				<< std::setw(_writePrecision+15) << std::scientific << std::setprecision(_writePrecision) << _p_acc/_numSamples
				<< std::setw(_writePrecision+15) << std::scientific << std::setprecision(_writePrecision) << domain->cv()
				<< std::setw(_writePrecision+15) << std::scientific << std::setprecision(_writePrecision) << globalNumMolecules
				<< endl;
			resultStream.close();
		}

		// Reset values
		_numSamples = 0UL;
		_uPot_acc = 0.0F;
		_uKin_acc = 0.0F;
		_uKinTrans_acc = 0.0F;
		_uKinRot_acc = 0.0F;
		_p_acc = 0.0F;
	}
}
