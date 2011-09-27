// PovWriter.cpp

#include "io/PovWriter.h"
#include "Common.h"
#include "Domain.h"
#include "particleContainer/ParticleContainer.h"
#include "molecules/Molecule.h"

#include <ctime>
#include <sstream>
#include <fstream>

using namespace std;

PovWriter::PovWriter(unsigned long writeFrequency, string filename, unsigned long numberOfTimesteps, bool incremental) {
	_filename = filename;
	_writeFrequency = writeFrequency;
	_incremental = incremental;
	_numberOfTimesteps = numberOfTimesteps;

	if (filename == "default")
		_filenameisdate = true;
	else
		_filenameisdate = false;
}

PovWriter::~PovWriter() {
}

void PovWriter::initOutput(ParticleContainer* particleContainer,
                           DomainDecompBase* domainDecomp, Domain* domain) {
}

void PovWriter::doOutput(ParticleContainer* particleContainer,
                         DomainDecompBase* domainDecomp, Domain* domain,
                         unsigned long simstep, list<ChemicalPotential>* lmu) {
	if (simstep % _writeFrequency == 0) {

		stringstream filenamestream;
		if (_filenameisdate) {
			filenamestream << "mardyn" << gettimestring();
		}
		else {
			filenamestream << _filename;
		}

		if (_incremental) {
			filenamestream << "-";
			/* align file numbers with preceding '0's in the required range from 0 to _numberOfTimesteps. */
			int num_digits = ceil(log(double(_numberOfTimesteps / _writeFrequency)) / log(10.));
			filenamestream << aligned_number(simstep / _writeFrequency, num_digits, '0');
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
		vector<Component> dcomponents = domain->getComponents();
		for (unsigned int i = 0; i < dcomponents.size(); ++i) {
			ostringstream osstrm;
			osstrm.clear();
			osstrm.str("");
			osstrm << " pigment {color rgb <" << (i + 1) % 2 << "," << (i + 1) / 2 % 2 << "," << (i + 1) / 4 % 2 << ">}";
			osstrm << " finish{ambient 0.5 diffuse 0.4 phong 0.3 phong_size 3}";
			ostrm << "#declare T" << i << " = ";
			//ostrm << "sphere {<0,0,0>,0.5 pigment {color rgb<1,0,0>} finish{ambient 0.5 diffuse 0.4 phong 0.3 phong_size 3} scale 1.}";
			dcomponents.at(i).writePOVobjs(ostrm, osstrm.str());
			ostrm << endl;
		}
		ostrm << endl;
		ostrm << "camera { perspective" << endl;
		float xloc = -.1 * domain->getGlobalLength(0);
		float yloc = 1.1 * domain->getGlobalLength(1);
		float zloc = -1.5 * domain->getGlobalLength(2);
		ostrm << " location <" << xloc << ", " << yloc << ", " << zloc << ">" << endl;
		//ostrm << " direction <0, 0, 1>" << endl;
		//ostrm << " right <1.33333, 0, 0>" << endl;
		//ostrm << " up <0, 1, 0>" << endl;
		//ostrm << " sky <0, 1, 0>" << endl;
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
		ostrm << "// " << dcomponents.size() << " objects for the atoms following..." << endl;
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

		// RK
		/* map cluster ID to color */
		//if ((pos->clusterid() != -1) && (m_clusters.find(pos->clusterid())->second > MINCLUSTERSIZE)) {
		//	ostrm << "  pigment {color rgb<" <<
		//			(0.25*(pos->clusterid()%5)) << "," <<
		//			(0.25*((pos->clusterid()/5)%5)) << "," <<
		//			(0.25*((pos->clusterid()/125)%5)) << "," <<
		//			">}";
		//}
		//else {
		//	ostrm << "  pigment {color rgb<0.9,0.9,0.9>}";
		//}
		// /RK

		ostrm.close();
	}
}

void PovWriter::finishOutput(ParticleContainer* particleContainer,
                             DomainDecompBase* domainDecomp, Domain* domain) {
}
