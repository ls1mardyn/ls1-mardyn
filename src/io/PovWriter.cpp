#include "io/PovWriter.h"

#include <fstream>
#include <sstream>
#include <vector>

#include "Common.h"
#include "Domain.h"
#include "ensemble/EnsembleBase.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "Simulation.h"
#include "utils/Logger.h"

using Log::global_log;
using namespace std;



static void writePOVobjs(Component const &component, std::ostream& ostrm, string para) {
	if (component.numLJcenters() <= 0) return;
	if (component.numLJcenters() == 1) {
		LJcenter LJsite = component.ljcenter(0);
		ostrm << "sphere {<" << LJsite.rx() << "," << LJsite.ry() << "," << LJsite.rz() << ">," << .5 * LJsite.sigma() << " " << para << "}";
	}
	else {
		ostrm << "blob { threshold 0.01 ";
		for (auto LJsite : component.ljcenters()) {
			ostrm << "sphere {<" << LJsite.rx() << "," << LJsite.ry() << "," << LJsite.rz() << ">," << .5 * LJsite.sigma() << ", strength 1 } ";
		}
		ostrm << para << "}";
	}
}


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
	} else{
		_appendTimestamp = false;
	}
	global_log->info() << "Append timestamp: " << _appendTimestamp << endl;
}

void PovWriter::init(ParticleContainer * /*particleContainer*/,
                     DomainDecompBase * /*domainDecomp*/, Domain * /*domain*/) {
}

void PovWriter::endStep(ParticleContainer *particleContainer,
                        DomainDecompBase * /*domainDecomp*/, Domain *domain,
                        unsigned long simstep) {
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
		const auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
		tm unused{};
		ostrm << "// " << std::put_time(localtime_r(&now, &unused), "%c") << endl;

		ostrm << "// bb: [0," << domain->getGlobalLength(0) << "]^3" << endl;
		ostrm << "//*PMRawBegin" << endl;
		ostrm << "background {rgb <1,1,1>}" << endl;
		ostrm << "//*PMRawEnd" << endl;
		vector<Component>* dcomponents = _simulation.getEnsemble()->getComponents();
		for (unsigned int i = 0; i < dcomponents->size(); ++i) {
			ostringstream osstrm;
			osstrm.clear();
			osstrm.str("");
			osstrm << " pigment {color rgb <" << (i + 1) % 2 << "," << (i + 1) / 2 % 2 << "," << (i + 1) / 4 % 2 << ">}";
			osstrm << " finish{ambient 0.5 diffuse 0.4 phong 0.3 phong_size 3}";
			ostrm << "#declare T" << i << " = ";
			writePOVobjs(dcomponents->at(i), ostrm, osstrm.str());
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
		for (auto pos = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); pos.isValid(); ++pos) {
			(pos->q()).getRotMatrix(mrot);
			ostrm << "object { T" << pos->componentid();
			ostrm << " matrix <"
			      << mrot[0][0] << "," << mrot[0][1] << "," << mrot[0][2] << ","
			      << mrot[1][0] << "," << mrot[1][1] << "," << mrot[1][2] << ","
			      << mrot[2][0] << "," << mrot[2][1] << "," << mrot[2][2] << ","
			      << pos->r(0) << "," << pos->r(1) << "," << pos->r(2)
			      << ">";
			ostrm << "}" << "\n";
		}
		ostrm.close();
	}
}

void PovWriter::finish(ParticleContainer * /*particleContainer*/, DomainDecompBase * /*domainDecomp*/,
					   Domain * /*domain*/) {}
