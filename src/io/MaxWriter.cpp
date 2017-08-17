
#include "io/MaxWriter.h"

#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "Simulation.h"
#include "utils/Logger.h"

#include <iomanip>

using Log::global_log;
using namespace std;

MaxWriter::MaxWriter()
	:
	_writeFrequency(1000),
	_outputPrefix("maxvals"),
	_dMaxValuesLocal(nullptr),
	_dMaxValuesGlobal(nullptr),
	_numQuantities(3),
	_numValsPerQuantity(4),
	_numValsPerComponent(3*4),
	_numComponents(2),
	_numVals(3*4*2)
{
	_numComponents = global_simulation->getEnsemble()->getComponents()->size()+1;  // 0: all components
}

MaxWriter::~MaxWriter()
{
	delete[] _dMaxValuesLocal;
	delete[] _dMaxValuesGlobal;
}

void MaxWriter::readXML(XMLfileUnits& xmlconfig)
{
	global_log->info() << "------------------------------------------------------------------------" << std::endl;
	global_log->info() << "MaxWriter" << std::endl;

	_writeFrequency = 1000;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	global_log->info() << "Write frequency: " << _writeFrequency << endl;

	_outputPrefix = "maxvals";
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	global_log->info() << "Output prefix: " << _outputPrefix << endl;

	global_log->info() << "------------------------------------------------------------------------" << std::endl;
}

void MaxWriter::initOutput(ParticleContainer* /*particleContainer*/,
			      DomainDecompBase* domainDecomp, Domain* /*domain*/)
{
	// init data structures
	this->initDataStructures();

	// initialize files (root process only)
	if( 0 != domainDecomp->getRank() )
		return;

	std::stringstream sstrFilename[3];
	sstrFilename[0] << _outputPrefix << "_max_veloc.dat";
	sstrFilename[1] << _outputPrefix << "_max_angmo.dat";
	sstrFilename[2] << _outputPrefix << "_max_force.dat";

	std::stringstream sstrOutput[3];
	for(uint8_t fi=0; fi<3; ++fi)
		sstrOutput[fi] << "                 simstep";

	for(uint32_t cid=0; cid<_numComponents; ++cid)
	{
		// velocity
		sstrOutput[0] << "                 vabs_c" << cid;
		sstrOutput[0] << "                   vx_c" << cid;
		sstrOutput[0] << "                   vy_c" << cid;
		sstrOutput[0] << "                   vz_c" << cid;

		// angular momentum
		sstrOutput[1] << "                 Labs_c" << cid;
		sstrOutput[1] << "                   Lx_c" << cid;
		sstrOutput[1] << "                   Ly_c" << cid;
		sstrOutput[1] << "                   Lz_c" << cid;

		// force
		sstrOutput[2] << "                 Fabs_c" << cid;
		sstrOutput[2] << "                   Fx_c" << cid;
		sstrOutput[2] << "                   Fy_c" << cid;
		sstrOutput[2] << "                   Fz_c" << cid;
	}

	for(uint8_t fi=0; fi<3; ++fi)
		sstrOutput[fi] << endl;

	for(uint8_t fi=0; fi<3; ++fi)
	{
		ofstream ofs(sstrFilename[fi].str().c_str(), ios::out);
		ofs << sstrOutput[fi].str();
		ofs.close();
	}
}

void MaxWriter::doOutput( ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain,
	unsigned long simstep, list<ChemicalPotential>* /*lmu*/, map<unsigned, CavityEnsemble>* mcav )
{
	this->doSampling(particleContainer);

	if(simstep % _writeFrequency != 0)
		return;

	this->calculateGlobalValues();
	this->resetLocalValues();
	this->writeData(domainDecomp);

}

void MaxWriter::finishOutput(ParticleContainer* /*particleContainer*/,
				DomainDecompBase* /*domainDecomp*/, Domain* /*domain*/)
{
}

void MaxWriter::initDataStructures()
{
	_numQuantities = 3;  // velocity, angular momentum, force
	_numValsPerQuantity = 4;  // Quantity (Q): Qabs, Qx, Qy, Qz
	_numValsPerComponent = _numQuantities * _numValsPerQuantity;
	_numVals = _numValsPerComponent * _numComponents;

	_dMaxValuesLocal  = new double[_numVals];
	_dMaxValuesGlobal = new double[_numVals];

	this->resetLocalValues();
}

void MaxWriter::doSampling(ParticleContainer* particleContainer)
{
	for (ParticleIterator pit = particleContainer->iteratorBegin();
			pit != particleContainer->iteratorEnd(); ++pit)
	{
		uint32_t cid = pit->componentid()+1;  // 0: all components
		uint32_t nOffsetComponent = cid*_numValsPerComponent;
		double v2 = pit->v2();
		double L2 = pit->L2();
		double F2 = pit->F2();

		for(uint32_t cid=0; cid<_numComponents; ++cid)
		{
			// velocity
			if(v2 > _dMaxValuesLocal[QO_VELOCITY_OFFSET])
			{
				_dMaxValuesLocal[QO_VELOCITY_OFFSET] = v2;
				for(uint8_t dim=0; dim<3; ++dim)
					_dMaxValuesLocal[QO_VELOCITY_OFFSET+1+dim] = pit->v(dim);
			}
			if(v2 > _dMaxValuesLocal[QO_VELOCITY_OFFSET+nOffsetComponent])
			{
				_dMaxValuesLocal[QO_VELOCITY_OFFSET+nOffsetComponent] = v2;
				for(uint8_t dim=0; dim<3; ++dim)
					_dMaxValuesLocal[QO_VELOCITY_OFFSET+1+dim+nOffsetComponent] = pit->v(dim);
			}

			// angular momentum
			if(L2 > _dMaxValuesLocal[QO_ANGULAR_MOMENTUM_OFFSET])
			{
				_dMaxValuesLocal[QO_ANGULAR_MOMENTUM_OFFSET] = L2;
				for(uint8_t dim=0; dim<3; ++dim)
					_dMaxValuesLocal[QO_ANGULAR_MOMENTUM_OFFSET+1+dim] = pit->v(dim);
			}
			if(L2 > _dMaxValuesLocal[QO_ANGULAR_MOMENTUM_OFFSET+nOffsetComponent])
			{
				_dMaxValuesLocal[QO_ANGULAR_MOMENTUM_OFFSET+nOffsetComponent] = L2;
				for(uint8_t dim=0; dim<3; ++dim)
					_dMaxValuesLocal[QO_ANGULAR_MOMENTUM_OFFSET+1+dim+nOffsetComponent] = pit->v(dim);
			}

			// force
			if(F2 > _dMaxValuesLocal[QO_FORCE_OFFSET])
			{
				_dMaxValuesLocal[QO_FORCE_OFFSET] = F2;
				for(uint8_t dim=0; dim<3; ++dim)
					_dMaxValuesLocal[QO_FORCE_OFFSET+1+dim] = pit->v(dim);
			}
			if(F2 > _dMaxValuesLocal[QO_FORCE_OFFSET+nOffsetComponent])
			{
				_dMaxValuesLocal[QO_FORCE_OFFSET+nOffsetComponent] = F2;
				for(uint8_t dim=0; dim<3; ++dim)
					_dMaxValuesLocal[QO_FORCE_OFFSET+1+dim+nOffsetComponent] = pit->v(dim);
			}
		}
	}
}

void MaxWriter::calculateGlobalValues()
{
#ifdef ENABLE_MPI

	MPI_Reduce( _dMaxValuesLocal, _dMaxValuesGlobal, _numVals, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

#else
	// Scalar quantities
	for(uint32_t vi=0; vi<_numVals; ++vi)
		_dMaxValuesGlobal[vi] = _dMaxValuesLocal[vi];
#endif
}

void MaxWriter::resetLocalValues()
{
	for(uint32_t vi=0; vi<_numVals; ++vi)
		_dMaxValuesLocal[vi] = 0.0;
}

void MaxWriter::writeData(DomainDecompBase* domainDecomp)
{
	// write data
	std::stringstream sstrFilename[3];
	sstrFilename[0] << _outputPrefix << "_max_veloc.dat";
	sstrFilename[1] << _outputPrefix << "_max_angmo.dat";
	sstrFilename[2] << _outputPrefix << "_max_force.dat";

	if( 0 != domainDecomp->getRank() )
		return;

	std::stringstream sstrOutput[3];

	// velocity
	sstrOutput[0] << std::setw(24) << global_simulation->getSimulationStep();
	for(uint32_t cid=0; cid<_numComponents; ++cid)
	{
		uint32_t nOffsetComponent = cid*_numValsPerComponent;
		double vmax = sqrt(_dMaxValuesGlobal[QO_VELOCITY_OFFSET+nOffsetComponent]);
		sstrOutput[0] << FORMAT_SCI_MAX_DIGITS << vmax;
		for(uint8_t vi=1; vi<_numValsPerQuantity; ++vi)
			sstrOutput[0] << FORMAT_SCI_MAX_DIGITS << _dMaxValuesGlobal[QO_VELOCITY_OFFSET+nOffsetComponent+vi];
	}
	sstrOutput[0] << endl;

	// angular momentum
	sstrOutput[1] << std::setw(24) << global_simulation->getSimulationStep();
	for(uint32_t cid=0; cid<_numComponents; ++cid)
	{
		uint32_t nOffsetComponent = cid*_numValsPerComponent;
		double Lmax = sqrt(_dMaxValuesGlobal[QO_ANGULAR_MOMENTUM_OFFSET+nOffsetComponent]);
		sstrOutput[1] << FORMAT_SCI_MAX_DIGITS << Lmax;
		for(uint8_t vi=1; vi<_numValsPerQuantity; ++vi)
			sstrOutput[1] << FORMAT_SCI_MAX_DIGITS << _dMaxValuesGlobal[QO_ANGULAR_MOMENTUM_OFFSET+nOffsetComponent+vi];
	}
	sstrOutput[1] << endl;

	// force
	sstrOutput[2] << std::setw(24) << global_simulation->getSimulationStep();
	for(uint32_t cid=0; cid<_numComponents; ++cid)
	{
		uint32_t nOffsetComponent = cid*_numValsPerComponent;
		double Fmax = sqrt(_dMaxValuesGlobal[QO_FORCE_OFFSET+nOffsetComponent]);
		sstrOutput[2] << FORMAT_SCI_MAX_DIGITS << Fmax;
		for(uint8_t vi=1; vi<_numValsPerQuantity; ++vi)
			sstrOutput[2] << FORMAT_SCI_MAX_DIGITS << _dMaxValuesGlobal[QO_FORCE_OFFSET+nOffsetComponent+vi];
	}
	sstrOutput[2] << endl;

	for(uint8_t fi=0; fi<3; ++fi)
	{
		ofstream ofs(sstrFilename[fi].str().c_str(), ios::app);
		ofs << sstrOutput[fi].str();
		ofs.close();
	}
}






















