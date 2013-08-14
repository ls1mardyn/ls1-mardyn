#include "io/PovWriter.h"

#include <fstream>
#include <sstream>

#include "Common.h"
#include "Domain.h"
#include "ensemble/EnsembleBase.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "Simulation.h"
#include "utils/Logger.h"

using Log::global_log;
using namespace std;


PovWriter::PovWriter(unsigned long writeFrequency, string outputPrefix, bool incremental) {
	_outputPrefix = outputPrefix;
	_writeFrequency = writeFrequency;
	_incremental = incremental;

	if (outputPrefix == "default") {
		_appendTimestamp = true;
	}
	else {
		_appendTimestamp = false;
	}
}

PovWriter::~PovWriter() {}


void PovWriter::readXML(XMLfileUnits& xmlconfig) {
	_writeFrequency = 1;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	global_log->info() << "Write frequency: " << _writeFrequency << endl;

	_outputPrefix = "mardyn";
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	global_log->info() << "Output prefix: " << _outputPrefix << endl;
	
	int incremental = 1;
	xmlconfig.getNodeValue("incremental", incremental);
	_incremental = (incremental != 0);
	global_log->info() << "Incremental numbers: " << _incremental << endl;
	
	int appendTimestamp = 0;
	xmlconfig.getNodeValue("appendTimestamp", appendTimestamp);
	if(appendTimestamp > 0) {
		_appendTimestamp = true;
	}
	global_log->info() << "Append timestamp: " << _appendTimestamp << endl;
}

void PovWriter::initOutput(ParticleContainer* particleContainer,
                           DomainDecompBase* domainDecomp, Domain* domain) {
}

void PovWriter::doOutput(ParticleContainer* particleContainer,
                         DomainDecompBase* domainDecomp, Domain* domain,
                         unsigned long simstep, list<ChemicalPotential>* lmu) {
	if (simstep % _writeFrequency == 0) {
		stringstream filenamestream;
		filenamestream << _outputPrefix;

		if(_incremental) {
			/* align file numbers with preceding '0's in the required range from 0 to _numberOfTimesteps. */
			unsigned long numTimesteps = _simulation.getNumTimesteps();
			int num_digits = (int) ceil( log( double( numTimesteps / _writeFrequency ) ) / log(10.) );
			filenamestream << "-" << aligned_number( simstep / _writeFrequency, num_digits, '0' );
		}
		if(_appendTimestamp) {
			filenamestream << "-" << gettimestring();
		}
		filenamestream << ".pov";

		ofstream ostrm(filenamestream.str().c_str());

		ostrm << "// " << filenamestream.str() << endl;
		ostrm << "// moldy" << endl;
		time_t now;
		now = time(NULL);
		ostrm << "// " << ctime(&now) << endl;

		ostrm << "// bb: [0," << domain->getGlobalLength(0) << "]^3" << endl;
		ostrm << "//*PMRawBegin" << endl;
		ostrm << "background {rgb <1,1,1>}" << endl;
		ostrm << "//*PMRawEnd" << endl;
		vector<Component>* dcomponents = _simulation.getEnsemble()->components();
		for (unsigned int i = 0; i < dcomponents->size(); ++i) {
			ostringstream osstrm;
			osstrm.clear();
			osstrm.str("");
			osstrm << " pigment {color rgb <" << (i + 1) % 2 << "," << (i + 1) / 2 % 2 << "," << (i + 1) / 4 % 2 << ">}";
			osstrm << " finish{ambient 0.5 diffuse 0.4 phong 0.3 phong_size 3}";
			ostrm << "#declare T" << i << " = ";
			dcomponents->at(i).writePOVobjs(ostrm, osstrm.str());
			ostrm << endl;
		}
		ostrm << endl;
		ostrm << "camera { perspective" << endl;
		float xloc = -.1 * domain->getGlobalLength(0);
		float yloc = 1.1 * domain->getGlobalLength(1);
		float zloc = -1.5 * domain->getGlobalLength(2);
		ostrm << " location <" << xloc << ", " << yloc << ", " << zloc << ">" << endl;
		ostrm << " look_at <" << .5 * domain->getGlobalLength(0) << ", " << .5 * domain->getGlobalLength(1) << ", " << .5 * domain->getGlobalLength(2) << ">" << endl;
		ostrm << "}" << endl;
		ostrm << endl;
		ostrm << "light_source { <" << xloc << ", " << yloc << ", " << zloc << ">, color rgb <1,1,1> }" << endl;
		ostrm << "light_source { <0,0,0>, color rgb <1,1,1> }" << endl;
		ostrm << "light_source { <0,0," << domain->getGlobalLength(2) << ">, color rgb <1,1,1> }" << endl;
		ostrm << "light_source { <0," << domain->getGlobalLength(1) << ",0>, color rgb <1,1,1> }" << endl;
		ostrm << "light_source { <0," << domain->getGlobalLength(1) << "," << domain->getGlobalLength(2) << ">, color rgb <1,1,1> }" << endl;
		ostrm << "light_source { <" << domain->getGlobalLength(0) << ",0,0>, color rgb <1,1,1> }" << endl;
		ostrm << "light_source { <" << domain->getGlobalLength(0) << ",0," << domain->getGlobalLength(2) << ">, color rgb <1,1,1> }" << endl;
		ostrm << "light_source { <" << domain->getGlobalLength(0) << "," << domain->getGlobalLength(1) << ",0>, color rgb <1,1,1> }" << endl;
		ostrm << "light_source { <" << domain->getGlobalLength(0) << "," << domain->getGlobalLength(1) << "," << domain->getGlobalLength(2) << ">, color rgb <1,1,1> }" << endl;
		ostrm << endl;
		ostrm << "// " << dcomponents->size() << " objects for the atoms following..." << endl;
		double mrot[3][3];
		for (Molecule* pos = particleContainer->begin(); pos != particleContainer->end(); pos = particleContainer->next()) {
			(pos->q()).getRotinvMatrix(mrot);
			//cout << "object { T0 rotate <0,0,0> translate <0,0,0>}" << endl;
			ostrm << "object { T" << pos->componentid();
			ostrm << " matrix <"
			      << mrot[0][0] << "," << mrot[0][1] << "," << mrot[0][2] << ","
			      << mrot[1][0] << "," << mrot[1][1] << "," << mrot[1][2] << ","
			      << mrot[2][0] << "," << mrot[2][1] << "," << mrot[2][2] << ","
			      << pos->r(0) << "," << pos->r(1) << "," << pos->r(2)
			      << ">";
			ostrm << "}" << endl;
		}
		ostrm.close();
	}
}

void PovWriter::finishOutput(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) {}
