
#include "io/MaxWriter.h"

#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "Simulation.h"
#include "utils/Logger.h"

#include <iomanip>
#include <vector>
#include <array>


MaxWriter::MaxWriter()
	:
	_writeFrequency(1000),
	_outputPrefix("maxvals"),
	_dMaxValuesLocal(),
	_dMaxValuesGlobal(),
	_numQuantities(7),
	_numValsPerQuantity(4),
	_numValsPerComponent(7*4),
	_numComponents(2),
	_numVals(7*4*2)
{
	_numComponents = global_simulation->getEnsemble()->getComponents()->size()+1;  // 0: all components
}

MaxWriter::~MaxWriter() {}

void MaxWriter::readXML(XMLfileUnits& xmlconfig)
{
	Log::global_log->info() << "------------------------------------------------------------------------" << std::endl;
	Log::global_log->info() << "MaxWriter" << std::endl;

	_writeFrequency = 1000;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	Log::global_log->info() << "Write frequency: " << _writeFrequency << std::endl;

	_outputPrefix = "maxvals";
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	Log::global_log->info() << "Output prefix: " << _outputPrefix << std::endl;

	Log::global_log->info() << "------------------------------------------------------------------------" << std::endl;
}

void MaxWriter::init(ParticleContainer * /*particleContainer*/,
                     DomainDecompBase *domainDecomp, Domain * /*domain*/)
{
	// init data structures
	this->initDataStructures();

	// initialize files (root process only)
	if( 0 != domainDecomp->getRank() )
		return;

	std::stringstream sstrFilename[4];
	sstrFilename[0] << _outputPrefix << "_max_veloc.dat";
	sstrFilename[1] << _outputPrefix << "_max_angmo.dat";
	sstrFilename[2] << _outputPrefix << "_max_force.dat";
	sstrFilename[3] << _outputPrefix << "_max_torqe.dat";

	std::stringstream sstrOutput[4];
	for(uint32_t qi=0; qi<_numQuantities; ++qi)
		sstrOutput[qi] << "                 simstep";

	for(uint32_t cid=0; cid<_numComponents; ++cid)
	{
		// velocity
		sstrOutput[0] << "                 vabs_c" << cid;
		sstrOutput[0] << "                  v+x_c" << cid;
		sstrOutput[0] << "                  v+y_c" << cid;
		sstrOutput[0] << "                  v+z_c" << cid;
		sstrOutput[0] << "                  v-x_c" << cid;
		sstrOutput[0] << "                  v-y_c" << cid;
		sstrOutput[0] << "                  v-z_c" << cid;

		// angular momentum
		sstrOutput[1] << "                 Labs_c" << cid;
		sstrOutput[1] << "                  L+x_c" << cid;
		sstrOutput[1] << "                  L+y_c" << cid;
		sstrOutput[1] << "                  L+z_c" << cid;
		sstrOutput[1] << "                  L-x_c" << cid;
		sstrOutput[1] << "                  L-y_c" << cid;
		sstrOutput[1] << "                  L-z_c" << cid;

		// force
		sstrOutput[2] << "                 Fabs_c" << cid;
		sstrOutput[2] << "                  F+x_c" << cid;
		sstrOutput[2] << "                  F+y_c" << cid;
		sstrOutput[2] << "                  F+z_c" << cid;
		sstrOutput[2] << "                  F-x_c" << cid;
		sstrOutput[2] << "                  F-y_c" << cid;
		sstrOutput[2] << "                  F-z_c" << cid;

		// torque
		sstrOutput[3] << "                 Mabs_c" << cid;
		sstrOutput[3] << "                  M+x_c" << cid;
		sstrOutput[3] << "                  M+y_c" << cid;
		sstrOutput[3] << "                  M+z_c" << cid;
		sstrOutput[3] << "                  M-x_c" << cid;
		sstrOutput[3] << "                  M-y_c" << cid;
		sstrOutput[3] << "                  M-z_c" << cid;
	}

	for(uint32_t qi=0; qi<_numQuantities; ++qi)
	{
		std::ofstream ofs(sstrFilename[qi].str().c_str(), std::ios::out);
		sstrOutput[qi] << std::endl;
		ofs << sstrOutput[qi].str();
		ofs.close();
	}
}

void MaxWriter::endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
                        unsigned long simstep)
{
	this->doSampling(particleContainer);

	if(simstep % _writeFrequency != 0)
		return;

	this->calculateGlobalValues(domainDecomp);
	this->resetLocalValues();
	this->writeData(domainDecomp);

}

void MaxWriter::finish(ParticleContainer * /*particleContainer*/,
					   DomainDecompBase * /*domainDecomp*/, Domain * /*domain*/)
{
}

void MaxWriter::initDataStructures()
{
	_numQuantities = 4;  // velocity, angular momentum, force, torque
	_numValsPerQuantity = 7;  // Quantity (Q): Qabs, Qx_max, Qy_max, Qz_max, Qx_min, Qy_min, Qz_min
	_numValsPerComponent = _numQuantities * _numValsPerQuantity;
	_numVals = _numValsPerComponent * _numComponents;

	_dMaxValuesLocal.resize(_numVals);
	_dMaxValuesGlobal.resize(_numVals);

	this->resetLocalValues();
}

void MaxWriter::doSampling(ParticleContainer* particleContainer)
{
	for (auto pit = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
			pit.isValid(); ++pit)
	{
		uint32_t cid = pit->componentid()+1;  // 0: all components
		uint32_t nOffsetComponent = cid*_numValsPerComponent;
		std::array<std::array<double,7>,4> arrQuantities;
		// squared absolute values
		arrQuantities.at(0).at(0) = pit->v2();
		arrQuantities.at(1).at(0) = pit->L2();
		arrQuantities.at(2).at(0) = pit->F2();
		arrQuantities.at(3).at(0) = pit->M2();
		// x-component
		arrQuantities.at(0).at(1) = pit->v(0);
		arrQuantities.at(1).at(1) = pit->D(0);
		arrQuantities.at(2).at(1) = pit->F(0);
		arrQuantities.at(3).at(1) = pit->M(0);
		// y-component
		arrQuantities.at(0).at(2) = pit->v(1);
		arrQuantities.at(1).at(2) = pit->D(1);
		arrQuantities.at(2).at(2) = pit->F(1);
		arrQuantities.at(3).at(2) = pit->M(1);
		// z-component
		arrQuantities.at(0).at(3) = pit->v(2);
		arrQuantities.at(1).at(3) = pit->D(2);
		arrQuantities.at(2).at(3) = pit->F(2);
		arrQuantities.at(3).at(3) = pit->M(2);

		for(uint32_t qi=0; qi<_numQuantities; ++qi)
		{
			uint32_t nOffsetQuantity = _numValsPerQuantity*qi;

			// all components
			_dMaxValuesLocal[nOffsetQuantity] = std::max(_dMaxValuesLocal[nOffsetQuantity], arrQuantities.at(qi).at(0));

			for(uint32_t dim=1; dim<4; ++dim) {
				// positive direction (+)
				_dMaxValuesLocal[nOffsetQuantity+dim] = std::max(_dMaxValuesLocal[nOffsetQuantity+dim], arrQuantities.at(qi).at(dim));
				// negative direction (-)
				_dMaxValuesLocal[nOffsetQuantity+dim+3] = std::min(_dMaxValuesLocal[nOffsetQuantity+dim+3], arrQuantities.at(qi).at(dim));
			}

			// specific component
			_dMaxValuesLocal[nOffsetComponent+nOffsetQuantity] = std::max(_dMaxValuesLocal[nOffsetComponent+nOffsetQuantity], arrQuantities.at(qi).at(0));

			for(uint32_t dim=1; dim<4; ++dim) {
				// positive direction (+)
				_dMaxValuesLocal[nOffsetComponent+nOffsetQuantity+dim] = std::max(_dMaxValuesLocal[nOffsetComponent+nOffsetQuantity+dim], arrQuantities.at(qi).at(dim));

				// negative direction (-)
				_dMaxValuesLocal[nOffsetComponent+nOffsetQuantity+dim+3] = std::min(_dMaxValuesLocal[nOffsetComponent+nOffsetQuantity+dim+3], arrQuantities.at(qi).at(dim));
			}
		}
	}
}

void MaxWriter::calculateGlobalValues(DomainDecompBase *domainDecomp)
{
	domainDecomp->collCommInit(_numVals);
	for (auto val : _dMaxValuesLocal) {
		domainDecomp->collCommAppendDouble(val);
	}
	domainDecomp->collCommAllreduceCustom(ReduceType::MAX);
	for (auto & val : _dMaxValuesGlobal) {
		val = domainDecomp->collCommGetDouble();
	}
	domainDecomp->collCommFinalize();
}

void MaxWriter::resetLocalValues()
{
	for(uint32_t vi=0; vi<_numVals; ++vi)
		_dMaxValuesLocal[vi] = 0.0;
}

void MaxWriter::writeData(DomainDecompBase* domainDecomp)
{
	// only root process writes data to files
	if( 0 != domainDecomp->getRank() )
		return;

	std::stringstream sstrFilename[4];
	sstrFilename[0] << _outputPrefix << "_max_veloc.dat";
	sstrFilename[1] << _outputPrefix << "_max_angmo.dat";
	sstrFilename[2] << _outputPrefix << "_max_force.dat";
	sstrFilename[3] << _outputPrefix << "_max_torqe.dat";

	std::stringstream sstrOutput[4];

	// write data to streams
	for(uint32_t qi=0; qi<_numQuantities; ++qi)
	{
		uint32_t nOffsetQuantity = _numValsPerQuantity*qi;

		sstrOutput[qi] << std::setw(24) << global_simulation->getSimulationStep();
		for(uint32_t cid=0; cid<_numComponents; ++cid)
		{
			uint32_t nOffsetComponent = cid*_numValsPerComponent;
			double vmax = sqrt(_dMaxValuesGlobal[nOffsetComponent+nOffsetQuantity]);
			sstrOutput[qi] << FORMAT_SCI_MAX_DIGITS << vmax;
			for(uint32_t vi=1; vi<_numValsPerQuantity; ++vi)
				sstrOutput[qi] << FORMAT_SCI_MAX_DIGITS << _dMaxValuesGlobal[nOffsetComponent+nOffsetQuantity+vi];
		}
		sstrOutput[qi] << std::endl;
	}

	// write streams to files
	for(uint32_t qi=0; qi<_numQuantities; ++qi)
	{
		std::ofstream ofs(sstrFilename[qi].str().c_str(), std::ios::app);
		ofs << sstrOutput[qi].str();
		ofs.close();
	}
}
