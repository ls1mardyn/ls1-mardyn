#include "io/GridGenerator.h"

#include "Domain.h"
#include "Simulation.h"
#include "ensemble/EnsembleBase.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "utils/Logger.h"
#include "utils/Random.h"
#include "utils/xmlfileUnits.h"

#include <string>
#include <map>

using Log::global_log;
using namespace std;

void GridGenerator::readXML(XMLfileUnits& xmlconfig) {

	/* Lattice */
	string latticeSystem;
	string latticeCentering;
	xmlconfig.getNodeValue("lattice@type", latticeSystem);
	global_log->info() << "Lattice system type: " << latticeSystem << endl;
	xmlconfig.getNodeValue("lattice@centering", latticeCentering);
	global_log->info() << "Lattice centering: " << latticeCentering << endl;

	double a[3];
	double b[3];
	double c[3];
	long dims[3];
	xmlconfig.changecurrentnode("lattice");
	xmlconfig.getNodeValueReduced("vec[@id='a']/x", a[0]);
	xmlconfig.getNodeValueReduced("vec[@id='a']/y", a[1]);
	xmlconfig.getNodeValueReduced("vec[@id='a']/z", a[2]);
	global_log->info() << "Vec a: " << a[0] << ", " << a[1] << ", " << a[2] << endl;
	xmlconfig.getNodeValueReduced("vec[@id='b']/x", b[0]);
	xmlconfig.getNodeValueReduced("vec[@id='b']/y", b[1]);
	xmlconfig.getNodeValueReduced("vec[@id='b']/z", b[2]);
	global_log->info() << "Vec b: " << b[0] << ", " << b[1] << ", " << b[2] << endl;
	xmlconfig.getNodeValueReduced("vec[@id='c']/x", c[0]);
	xmlconfig.getNodeValueReduced("vec[@id='c']/y", c[1]);
	xmlconfig.getNodeValueReduced("vec[@id='c']/z", c[2]);
	global_log->info() << "Vec c: " << c[0] << ", " << c[1] << ", " << c[2] << endl;
	dims[0] = xmlconfig.getNodeValue_long("dims@a");
	dims[1] = xmlconfig.getNodeValue_long("dims@b");
	dims[2] = xmlconfig.getNodeValue_long("dims@c");
	global_log->info() << "Dims: " << dims[0] << ", " << dims[1] << ", " << dims[2] << endl;
	xmlconfig.changecurrentnode("..");

	std::map<string, LatticeSystem> latticeSystemName2Enum;
	latticeSystemName2Enum.insert(pair<string,LatticeSystem>("triclinic",   triclinic));
	latticeSystemName2Enum.insert(pair<string,LatticeSystem>("monoclinic",  monoclinic));
	latticeSystemName2Enum.insert(pair<string,LatticeSystem>("orthorombic", orthorombic));
	latticeSystemName2Enum.insert(pair<string,LatticeSystem>("tetragonal",  tetragonal));
	latticeSystemName2Enum.insert(pair<string,LatticeSystem>("rhomboedral", rhomboedral));
	latticeSystemName2Enum.insert(pair<string,LatticeSystem>("hexagonal",   hexagonal));
	latticeSystemName2Enum.insert(pair<string,LatticeSystem>("cubic",cubic));

	std::map<string, LatticeCentering> latticeCenteringName2Enum;
	latticeCenteringName2Enum.insert(pair<string,LatticeCentering>("primitive",  primitive));
	latticeCenteringName2Enum.insert(pair<string,LatticeCentering>("body",  body));
	latticeCenteringName2Enum.insert(pair<string,LatticeCentering>("face",  face));
	latticeCenteringName2Enum.insert(pair<string,LatticeCentering>("base A",  base_A));
	latticeCenteringName2Enum.insert(pair<string,LatticeCentering>("base B",  base_B));
	latticeCenteringName2Enum.insert(pair<string,LatticeCentering>("base C",  base_C));

	_lattice.init(latticeSystemName2Enum[latticeSystem], latticeCenteringName2Enum[latticeCentering], a, b, c, dims);

	/* Basis */
	xmlconfig.changecurrentnode("basis");
	XMLfile::Query query = xmlconfig.query("site");
	string oldpath = xmlconfig.getcurrentnodepath();
	XMLfile::Query::const_iterator siteIter;
	for(siteIter = query.begin(); siteIter; siteIter++) {
		molecule_t m;
		xmlconfig.changecurrentnode(siteIter);
		string componentid;
		xmlconfig.getNodeValue("componentid", componentid);
		m.cid = _simulation.getEnsemble()->component(componentid)->ID();
		xmlconfig.getNodeValueReduced("coordinate@x", m.r[0]);
		xmlconfig.getNodeValueReduced("coordinate@y", m.r[1]);
		xmlconfig.getNodeValueReduced("coordinate@z", m.r[2]);
		global_log->info() << "Adding molecule cid=" << componentid << "(" << m.cid << "), (x,y,z)=(" << m.r[0] << "," << m.r[1] << "," << m.r[2] << ")" << endl;
		_basis.addMolecule(m);
	}
	xmlconfig.changecurrentnode(oldpath);
	xmlconfig.changecurrentnode("..");
	/* Generator */
	xmlconfig.changecurrentnode("origin");
	xmlconfig.getNodeValueReduced("x", _origin[0]);
	xmlconfig.getNodeValueReduced("y", _origin[1]);
	xmlconfig.getNodeValueReduced("z", _origin[2]);
	global_log->info() << "Origin: " << _origin[0] << ", " << _origin[1] << ", " << _origin[2] << endl;
	xmlconfig.changecurrentnode("..");
	_generator.init(_lattice, _basis, _origin);
}

long unsigned int GridGenerator::readPhaseSpace(ParticleContainer* particleContainer, list< ChemicalPotential >* lmu, Domain* domain, DomainDecompBase* domainDecomp) {
	unsigned long numMolecules = 0;
	molecule_t m; /* molecule type as provided by the generator */

	Ensemble* ensemble = _simulation.getEnsemble();
	Random rng;

	_simulation.getDomain()->disableComponentwiseThermostat();
	while(_generator.getMolecule(&m) > 0) {
		Component* component = ensemble->component(m.cid);
		Molecule molecule(0, component); /* Molecule type as provided by mardyn */
		double v_abs = sqrt(/*kB=1*/ ensemble->T() / molecule.component()->m());
		double phi, theta;
		phi = rng.rnd();
		theta = rng.rnd();
		double v[3];
		v[0] = v_abs * sin(phi);
		v[1] = v_abs * cos(phi) * sin(theta);
		v[2] = v_abs * cos(phi) * cos(theta);
		molecule.setid(numMolecules);
		for(int d = 0; d < 3; d++) {
			molecule.setr(d, m.r[d]);
			molecule.setv(d, v[d]);
		}
		particleContainer->addParticle(molecule);
		numMolecules++;
	}
	global_log->info() << "Number of inserted molecules: " << numMolecules << endl;
	particleContainer->updateMoleculeCaches();
	return numMolecules;
}
