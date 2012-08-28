/*
 * RDFTest.cpp
 *
 * @Date: 15.02.2011
 * @Author: eckhardw
 */

#include "RDFTest.h"

#include "RDF.h"
#include "Domain.h"
#include "Simulation.h"
#include "parallel/DomainDecompBase.h"
#include "parallel/DomainDecompDummy.h"
#include "particleContainer/LinkedCells.h"
#include "particleContainer/adapter/ParticlePairs2PotForceAdapter.h"

#ifdef ENABLE_MPI
#include "parallel/DomainDecomposition.h"
#endif

#include <sstream>

using namespace std;
TEST_SUITE_REGISTRATION(RDFTest);

RDFTest::RDFTest() {
}

RDFTest::~RDFTest() {
}

void RDFTest::testRDFCountSequential12_LinkedCell() {
	delete _domainDecomposition;
	// will be freed by the tearDown()-method.
	_domainDecomposition = new DomainDecompDummy();

	ParticleContainer* moleculeContainer =
			initializeFromFile(ParticleContainerFactory::LinkedCell,
					"1clj-regular-2x2x3.inp", 1.8);
	testRDFCountSequential12(moleculeContainer);
	delete moleculeContainer;
}

void RDFTest::testRDFCountSequential12_AdaptiveCell() {
	delete _domainDecomposition;
	// will be freed by the tearDown()-method.
	_domainDecomposition = new DomainDecompDummy();

	ParticleContainer* moleculeContainer = initializeFromFile(
			ParticleContainerFactory::AdaptiveSubCell,
			"1clj-regular-2x2x3.inp", 1.8);
	testRDFCountSequential12(moleculeContainer);
	delete moleculeContainer;
}

void RDFTest::testRDFCountSequential12(ParticleContainer* moleculeContainer) {
	ParticlePairs2PotForceAdapter handler(*_domain);
	const vector<Component>& components = _domain->getComponents();
	ASSERT_EQUAL((size_t) 1, components.size());

	moleculeContainer->update();
	moleculeContainer->updateMoleculeCaches();

	/* The number of pairs counted by the RDF also depends on the particles in the halo.
	 * So count first with the halo being empty, and then being populated. */
	RDF rdf(0.018, 100, components);
	handler.setRDF(&rdf);
	rdf.tickRDF();
	moleculeContainer->traversePairs(&handler);
	rdf.collectRDF(_domainDecomposition);

	for (int i = 0; i < 100; i++) {
		if (i == 55) {
			ASSERT_EQUAL(20ul, rdf._globalDistribution[0][0][i]);
		} else if (i == 78) {
			ASSERT_EQUAL(22ul, rdf._globalDistribution[0][0][i]);
		} else if (i == 96) {
			ASSERT_EQUAL(8ul, rdf._globalDistribution[0][0][i]);
		} else {
			ASSERT_EQUAL(0ul, rdf._globalDistribution[0][0][i]);
		}
	}

	rdf.accumulateRDF();
	rdf.reset();
	rdf.tickRDF();

	// now the same with halo particles present.
	_domainDecomposition->exchangeMolecules(moleculeContainer,
			_domain->getComponents(), _domain);
	moleculeContainer->traversePairs(&handler);
	rdf.collectRDF(_domainDecomposition);
	rdf.accumulateRDF();

	for (int i = 0; i < 100; i++) {
		if (i == 55) {
			ASSERT_EQUAL(20ul, rdf._globalDistribution[0][0][i]);
			ASSERT_EQUAL(40ul, rdf._globalAccumulatedDistribution[0][0][i]);
		} else if (i == 78) {
			ASSERT_EQUAL(22ul, rdf._globalDistribution[0][0][i]);
			ASSERT_EQUAL(44ul, rdf._globalAccumulatedDistribution[0][0][i]);
		} else if (i == 83) {
			ASSERT_EQUAL(4ul, rdf._globalDistribution[0][0][i]);
			ASSERT_EQUAL(4ul, rdf._globalAccumulatedDistribution[0][0][i]);
		} else if (i == 96) {
			ASSERT_EQUAL(8ul, rdf._globalDistribution[0][0][i]);
			ASSERT_EQUAL(16ul, rdf._globalAccumulatedDistribution[0][0][i]);
		} else {
			ASSERT_EQUAL(0ul, rdf._globalDistribution[0][0][i]);
			ASSERT_EQUAL(0ul, rdf._globalAccumulatedDistribution[0][0][i]);
		}
	}
}

void RDFTest::testRDFCountLinkedCell() {
	ParticleContainer* moleculeContainer = initializeFromFile(
			ParticleContainerFactory::LinkedCell, "1clj-regular-12x12x12.inp",
			1.8);
	testRDFCount(moleculeContainer);
	delete moleculeContainer;
}

void RDFTest::testRDFCountAdaptiveCell() {
	ParticleContainer* moleculeContainer = initializeFromFile(
			ParticleContainerFactory::AdaptiveSubCell,
			"1clj-regular-12x12x12.inp", 1.8);
	testRDFCount(moleculeContainer);
	delete moleculeContainer;
}

void RDFTest::testRDFCount(ParticleContainer* moleculeContainer) {
	ParticlePairs2PotForceAdapter handler(*_domain);
	const vector<Component>& components = _domain->getComponents();
	ASSERT_EQUAL((size_t) 1, components.size());

	_domainDecomposition->balanceAndExchange(true, moleculeContainer,
			_domain->getComponents(), _domain);
	moleculeContainer->updateMoleculeCaches();

	RDF rdf(0.018, 100, components);
	handler.setRDF(&rdf);
	rdf.tickRDF();
	moleculeContainer->traversePairs(&handler);
	rdf.collectRDF(_domainDecomposition);

	// assert number of pairs counted
	for (int i = 0; i < 100; i++) {
		if (i == 55) {
			ASSERT_EQUAL(4752ul, rdf._globalDistribution[0][0][i]);
		} else if (i == 78) {
			ASSERT_EQUAL(8712ul, rdf._globalDistribution[0][0][i]);
		} else if (i == 83) {
			ASSERT_EQUAL(432ul, rdf._globalDistribution[0][0][i]);
		} else if (i == 96) {
			ASSERT_EQUAL(5324ul, rdf._globalDistribution[0][0][i]);
		} else {
			ASSERT_EQUAL(0ul, rdf._globalDistribution[0][0][i]);
		}
	}

	rdf.accumulateRDF();
	rdf.reset();

	rdf.tickRDF();
	moleculeContainer->traversePairs(&handler);
	rdf.collectRDF(_domainDecomposition);
	rdf.accumulateRDF();

	for (int i = 0; i < 100; i++) {
		stringstream msg;
		msg << "at index " << i;
		if (i == 55) {
			ASSERT_EQUAL(4752ul, rdf._globalDistribution[0][0][i]);
		} else if (i == 78) {
			ASSERT_EQUAL(8712ul, rdf._globalDistribution[0][0][i]);
		} else if (i == 83) {
			ASSERT_EQUAL(432ul, rdf._globalDistribution[0][0][i]);
		} else if (i == 96) {
			ASSERT_EQUAL(5324ul, rdf._globalDistribution[0][0][i]);
		} else {
			ASSERT_EQUAL_MSG(msg.str(), 0ul, rdf._globalDistribution[0][0][i]);
		}

		// the accumulated global distribution must be now twice the global distribution
		ASSERT_EQUAL(rdf._globalAccumulatedDistribution[0][0][i], 2 * rdf._globalDistribution[0][0][i]);
	}
}

void RDFTest::testSiteSiteRDFLinkedCell() {
	ParticleContainer* moleculeContainer = initializeFromFile(
			ParticleContainerFactory::LinkedCell, "2clj-regular.inp", 3.5);
	testSiteSiteRDF(moleculeContainer);
	delete moleculeContainer;
}

void RDFTest::testSiteSiteRDF(ParticleContainer* moleculeContainer) {

	if (_domainDecomposition->getNumProcs() > 8) {
		ASSERT_FAIL("RUN THIS TEST WITH <= 8 PROCESSORS!");
	}

	ParticlePairs2PotForceAdapter handler(*_domain);
	const vector<Component>& components = _domain->getComponents();
	ASSERT_EQUAL((size_t) 1, components.size());

	_domainDecomposition->balanceAndExchange(true, moleculeContainer,
			_domain->getComponents(), _domain);
	moleculeContainer->update();
	moleculeContainer->updateMoleculeCaches();

	RDF rdf(0.05, 101, components);
	handler.setRDF(&rdf);
	rdf.tickRDF();
	moleculeContainer->traversePairs(&handler);
	rdf.collectRDF(_domainDecomposition);

	for (int i = 0; i < 101; i++) {
		//		std::cout << "Bin " << i << ": " << rdf._globalSiteDistribution[0][0][0][0][i] <<
		//					", " << rdf._globalSiteDistribution[0][0][0][1][i] <<
		//					", " << rdf._globalSiteDistribution[0][0][1][0][i] <<
		//					", " << rdf._globalSiteDistribution[0][0][1][1][i] << std::endl;
		if (i == 20) {
			ASSERT_EQUAL(0ul, rdf._globalSiteDistribution[0][0][0][0][i]);
			ASSERT_EQUAL(16ul, rdf._globalSiteDistribution[0][0][0][1][i]);
			ASSERT_EQUAL(0ul, rdf._globalSiteDistribution[0][0][1][0][i]);
			ASSERT_EQUAL(0ul, rdf._globalSiteDistribution[0][0][1][1][i]);
		} else if (i == 60) {
			ASSERT_EQUAL(16ul, rdf._globalSiteDistribution[0][0][0][0][i]);
			ASSERT_EQUAL(0ul, rdf._globalSiteDistribution[0][0][0][1][i]);
			ASSERT_EQUAL(0ul, rdf._globalSiteDistribution[0][0][1][0][i]);
			ASSERT_EQUAL(16ul, rdf._globalSiteDistribution[0][0][1][1][i]);
		} else if (i == 100) {
			ASSERT_EQUAL(0ul, rdf._globalSiteDistribution[0][0][0][0][i]);
			ASSERT_EQUAL(16ul, rdf._globalSiteDistribution[0][0][0][1][i]);
			ASSERT_EQUAL(0ul, rdf._globalSiteDistribution[0][0][1][0][i]);
			ASSERT_EQUAL(0ul, rdf._globalSiteDistribution[0][0][1][1][i]);
		} else {
			ASSERT_EQUAL(0ul, rdf._globalSiteDistribution[0][0][0][0][i]);
			ASSERT_EQUAL(0ul, rdf._globalSiteDistribution[0][0][0][1][i]);
			ASSERT_EQUAL(0ul, rdf._globalSiteDistribution[0][0][1][0][i]);
			ASSERT_EQUAL(0ul, rdf._globalSiteDistribution[0][0][1][1][i]);
		}
	}

	rdf.accumulateRDF();
	rdf.reset();
	rdf.tickRDF();

	// test the accumulation of counts...
	moleculeContainer->traversePairs(&handler);
	rdf.collectRDF(_domainDecomposition);
	rdf.accumulateRDF();

	for (int i = 0; i < 101; i++) {
		if (i == 20) {
			ASSERT_EQUAL(0ul, rdf._globalSiteDistribution[0][0][0][0][i]);
			ASSERT_EQUAL(16ul, rdf._globalSiteDistribution[0][0][0][1][i]);
			ASSERT_EQUAL(0ul, rdf._globalSiteDistribution[0][0][1][0][i]);
			ASSERT_EQUAL(0ul, rdf._globalSiteDistribution[0][0][1][1][i]);
			// accumulated counts
			ASSERT_EQUAL(0ul, rdf._globalAccumulatedSiteDistribution[0][0][0][0][i]);
			ASSERT_EQUAL(32ul, rdf._globalAccumulatedSiteDistribution[0][0][0][1][i]);
			ASSERT_EQUAL(0ul, rdf._globalAccumulatedSiteDistribution[0][0][1][0][i]);
			ASSERT_EQUAL(0ul, rdf._globalAccumulatedSiteDistribution[0][0][1][1][i]);
		} else if (i == 60) {
			ASSERT_EQUAL(16ul, rdf._globalSiteDistribution[0][0][0][0][i]);
			ASSERT_EQUAL(0ul, rdf._globalSiteDistribution[0][0][0][1][i]);
			ASSERT_EQUAL(0ul, rdf._globalSiteDistribution[0][0][1][0][i]);
			ASSERT_EQUAL(16ul, rdf._globalSiteDistribution[0][0][1][1][i]);
			// accumulated counts
			ASSERT_EQUAL(32ul, rdf._globalAccumulatedSiteDistribution[0][0][0][0][i]);
			ASSERT_EQUAL(0ul, rdf._globalAccumulatedSiteDistribution[0][0][0][1][i]);
			ASSERT_EQUAL(0ul, rdf._globalAccumulatedSiteDistribution[0][0][1][0][i]);
			ASSERT_EQUAL(32ul, rdf._globalAccumulatedSiteDistribution[0][0][1][1][i]);
		} else if (i == 100) {
			ASSERT_EQUAL(0ul, rdf._globalSiteDistribution[0][0][0][0][i]);
			ASSERT_EQUAL(16ul, rdf._globalSiteDistribution[0][0][0][1][i]);
			ASSERT_EQUAL(0ul, rdf._globalSiteDistribution[0][0][1][0][i]);
			ASSERT_EQUAL(0ul, rdf._globalSiteDistribution[0][0][1][1][i]);
			// accumulated counts
			ASSERT_EQUAL(0ul, rdf._globalAccumulatedSiteDistribution[0][0][0][0][i]);
			ASSERT_EQUAL(32ul, rdf._globalAccumulatedSiteDistribution[0][0][0][1][i]);
			ASSERT_EQUAL(0ul, rdf._globalAccumulatedSiteDistribution[0][0][1][0][i]);
			ASSERT_EQUAL(0ul, rdf._globalAccumulatedSiteDistribution[0][0][1][1][i]);
		} else {
			ASSERT_EQUAL(0ul, rdf._globalSiteDistribution[0][0][0][0][i]);
			ASSERT_EQUAL(0ul, rdf._globalSiteDistribution[0][0][0][1][i]);
			ASSERT_EQUAL(0ul, rdf._globalSiteDistribution[0][0][1][0][i]);
			ASSERT_EQUAL(0ul, rdf._globalSiteDistribution[0][0][1][1][i]);
			// accumulated counts
			ASSERT_EQUAL(0ul, rdf._globalAccumulatedSiteDistribution[0][0][0][0][i]);
			ASSERT_EQUAL(0ul, rdf._globalAccumulatedSiteDistribution[0][0][0][1][i]);
			ASSERT_EQUAL(0ul, rdf._globalAccumulatedSiteDistribution[0][0][1][0][i]);
			ASSERT_EQUAL(0ul, rdf._globalAccumulatedSiteDistribution[0][0][1][1][i]);
		}
	}

}

void RDFTest::evaluateForcesHalo(double*** forces, std::string inp_file) {

	delete _domainDecomposition;
	delete _domain;
	_domain = new Domain(_rank, NULL);
	_domainDecomposition = new DomainDecompDummy();
	ParticleContainer* moleculeContainer_halo = initializeFromFile(
			ParticleContainerFactory::LinkedCell, inp_file.c_str(), 32.1254);//29.5333);
	int numParticles = moleculeContainer_halo->getNumberOfParticles();
	moleculeContainer_halo->update();
	//_domainDecomposition->exchangeMolecules(_moleculeContainer, _domain->getComponents(), _domain);
	_domainDecomposition->balanceAndExchange(true, moleculeContainer_halo,
			_domain->getComponents(), _domain);
	// The cache of the molecules must be updated/build after the exchange process,
	// as the cache itself isn't transferred
	moleculeContainer_halo->updateMoleculeCaches();
	ParticlePairs2PotForceAdapter handler_halo(*_domain);
	const vector<Component>& components_halo = _domain->getComponents();
	ASSERT_EQUAL((size_t) 1, components_halo.size());
	//RDF rdf_halo(0.321254, 100, components_halo);
	//rdf_halo.collectRDF(_domainDecomposition);
	//moleculeContainer->traversePairs(&handler);
	moleculeContainer_halo->traversePairs(&handler_halo);
	moleculeContainer_halo->deleteOuterParticles();
	Molecule* moleculePtr;

	for (int i = 0; i < numParticles; i++)
		(*forces)[i][0] = 0;

	int num = 0;
	for (moleculePtr = moleculeContainer_halo->begin(); moleculePtr
			!= moleculeContainer_halo->end() && num < numParticles; moleculePtr
			= moleculeContainer_halo->next()) {

		moleculePtr->calcFM();

		double halo[3] = { 0, 0, 0 };

		double rmin[3], rmax[3], halo_L[3], low_limit, high_limit;
		for (int i = 0; i < 3; i++) {
			rmin[i] = moleculeContainer_halo->getBoundingBoxMin(i);
			rmax[i] = moleculeContainer_halo->getBoundingBoxMax(i);
			halo_L[i] = moleculeContainer_halo->get_halo_L(i);
			low_limit = rmin[i] + halo_L[i];
			high_limit = rmax[i] - halo_L[i];
			if (moleculePtr->r(i) < low_limit)
				halo[i] = moleculePtr->r(i);
			else if (moleculePtr->r(i) >= high_limit)
				halo[i] = moleculePtr->r(i);
		}

		if ((*forces)[moleculePtr->id()][0] == 0) {
			(*forces)[moleculePtr->id()][0] = moleculePtr->id();
			for (int i = 0; i < 3; i++)
				(*forces)[moleculePtr->id()][i + 1] = moleculePtr->F(i);
			for (int i = 4; i < 7; i++)
				(*forces)[moleculePtr->id()][i] = halo[i - 4];
			for (int i = 7; i < 10; i++)
				(*forces)[moleculePtr->id()][i] = moleculePtr->getLeftxF()[i - 7];
			num++;
		} else
		cout << "problem " << endl;

	}

}

void RDFTest::evaluateForcesRDF(double*** forces, std::string inp_file,
		std::string rdf_file) {
	delete _domainDecomposition;
	//delete _domain;
	_domain = new Domain(_rank, NULL);
	_domainDecomposition = new RDFDummyDecomposition();
	ParticleContainer* moleculeContainer = initializeFromFile(
			ParticleContainerFactory::LinkedCell, inp_file.c_str(), 32.1254);//29.5333);
	moleculeContainer->update();
	//_domainDecomposition->exchangeMolecules(_moleculeContainer, _domain->getComponents(), _domain);
	//_domainDecomposition->balanceAndExchange(true, moleculeContainer,
	//		_domain->getComponents(), _domain);
	// The cache of the molecules must be updated/build after the exchange process,
	// as the cache itself isn't transferred
	moleculeContainer->updateMoleculeCaches();
	ParticlePairs2PotForceAdapter handler(*_domain);
	const vector<Component>& components = _domain->getComponents();
	ASSERT_EQUAL((size_t) 1, components.size());
	//RDF rdf(0.321254, 100, components);
	//rdf.collectRDF(_domainDecomposition);
	std::vector<std::string> file_names;
	file_names.push_back(rdf_file);
	moleculeContainer->traversePairs(&handler, file_names);
	//moleculeContainer->traversePairs(&handler);
	Molecule* moleculePtr;

	int num = 0;
	for (moleculePtr = moleculeContainer->begin(); moleculePtr
			!= moleculeContainer->end(); moleculePtr
			= moleculeContainer->next()) {
		moleculePtr->calcFM();
		(*forces)[moleculePtr->id()][0] = moleculePtr->id();
		for (int i = 0; i < 3; i++)
			(*forces)[moleculePtr->id()][i + 1] = moleculePtr->F(i);
		for (int i = 4; i < 7; i++)
			(*forces)[moleculePtr->id()][i] = moleculePtr->getLeftxRdfF()[i - 4];
		num++;

	}

}

void RDFTest::compareFiles() {
	FILE* file_rdf = fopen("forces_10000molecules_2steps_rdf.txt", "r");
	FILE* file_periodic = fopen("forces_10000molecules_2steps_periodic.txt", "r");
	FILE* file_comp = fopen("file_comparison_rdf_periodic.txt", "w");
	int read1 = 1, read2 = 1;
	double total_abs_diff = 0;
	int ok = 0;
	while (read1 != EOF && read2 != EOF) {
		double rdf0, rdf1, rdf2, rdf3, rdf4, rdf5, periodic0, periodic1, periodic2, periodic3, periodic4, periodic5;

		int num1, num2;
		read1 = fscanf(file_rdf, "%d %lg %lg %lg %lg %lg %lg \n", &num1, &rdf0, &rdf1, &rdf2, &rdf3, &rdf4, &rdf5);
		read2 = fscanf(file_periodic, "%d %lg %lg %lg %lg %lg %lg \n", &num2, &periodic0, &periodic1, &periodic2, &periodic3, &periodic4, &periodic5);
		if (num1 == num2 && rdf0 != -1) {
			double diff[3];
			ok++;
			diff[0] = periodic0 - rdf0;
			diff[1] = periodic1 - rdf1;
			diff[2] = periodic2 - rdf2;
			double abs_diff = std::sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);
			total_abs_diff += abs_diff;
			fprintf(file_comp, "%d		%g %g %g		%g %g %g	%g %g %g	%g %g %g 	%g\n", num1, periodic0, periodic1, periodic2, rdf0, rdf1, rdf2, periodic3, periodic4, periodic5, rdf3, rdf4, rdf5, abs_diff);

		} else cout<<"problem, not the same id!"<<endl;
	}
	fclose(file_rdf);
	fclose(file_periodic);
	fclose(file_comp);

	cout<<"comparison: total abs diff " << total_abs_diff/ok << endl;
}

void RDFTest::testRDFPressureCondition() {
	std::string inp_file = "../test_input/Ethan_10k_equilibrated.inp";//"../test_input/Ethane_eqv.inp";
	std::string rdf_file = "../test_input/rdf_Ethan_10k_0-0.000052000.rdf";//"../test_input/rdf_Ethan_16k_0-0.000100100.rdf";
	double allowed_tolerance = 0.001;
	ParticleContainer* moleculeContainer = initializeFromFile(
			ParticleContainerFactory::LinkedCell, inp_file.c_str(), 32.1254);//29.5333);
	cout<<"number of particles: "<<moleculeContainer->getNumberOfParticles()<<endl;
	double** forces_rdf = new double*[moleculeContainer->getNumberOfParticles() + 1];
			//new double*[moleculeContainer->getNumberOfParticles()];
	for (unsigned int i = 0; i < 159015; i++)
		forces_rdf[i] = (double*) malloc(7 * sizeof(double));//new double[4];

	double** forces_halo =
			new double*[moleculeContainer->getNumberOfParticles()];
	for (unsigned int i = 0; i < 159015; i++)
		forces_halo[i] = new double[10];
	cout<<"here"<<endl;
	evaluateForcesRDF(&forces_rdf, inp_file, rdf_file);
	evaluateForcesHalo(&forces_halo, inp_file);

	FILE* file = fopen("rdf_pressure_condition.txt", "w");


	double sum_abs_diff = 0;
	int halos = 0;
	for (unsigned int num = 1; num <= moleculeContainer->getNumberOfParticles(); num++) {
		double diff[3];
		if (forces_halo[num][0] == forces_rdf[num][0]) {
			for (int i = 1; i < 4; i++) {
				diff[i - 1] = forces_halo[num][i] - forces_rdf[num][i];
			}
			double local_diff = 0;
			if ((forces_halo[num][4]) != 0 || (forces_halo[num][5]) != 0 ||
					(forces_halo[num][6]) != 0) {
				local_diff = std::sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);
				sum_abs_diff += local_diff;
				halos++;
			}
			fprintf(file, "%7d %11g %11g %11g		%11g %11g %11g		%11g %11g %11g		%11g %11g %11g	%7g %7g %7g %11g\n",
					(int) forces_halo[num][0], forces_halo[num][1],
					forces_halo[num][2], forces_halo[num][3],
					forces_rdf[num][1], forces_rdf[num][2], forces_rdf[num][3],
					forces_halo[num][7], forces_halo[num][8], forces_halo[num][9],
					forces_rdf[num][4], forces_rdf[num][5], forces_rdf[num][6],
					forces_halo[num][4], forces_halo[num][5],
					forces_halo[num][6], local_diff);
		}
		else
			fprintf(file, "diff\n");
	}
	fclose(file);
	printf("average abs difference: %g\n", sum_abs_diff/halos);
}

void RDFTest::readOutputFile() {
	FILE* file = fopen("comparison_ethan160k_rc4.5_1000steps_leftBoundary_other.txt", "r");
	int read = 1;
	double avg = 0;
	double num = 0;
	double max = 0;
	int d, less90=0;
	int maxd;
	double f1, f2, f3, f4, f5, f6, f7, f8, f9, val, val1, val2;
	read = fscanf(file, "total_abs: %lg less than 90 percent for: %d molecules \n", &f1, &d);
	cout<<" f1 "<<f1<<endl;
	read = fscanf(file, "total periodic: %lg 	total_rdf: %lg \n", &f1, &f2);
	while (read != EOF) {
		read = fscanf(file, "%d %lg %lg %lg %lg %lg %lg %lg %lg %lg 		%lg %lg %lg\n", &d, &f1, &f2, &f3, &f4, &f5, &f6, &f7, &f8, &f9, &val, &val1, &val2);
		if (val != -1 && abs(val) < 100000) {
			avg += val;
			num ++;
			if (val > max)  {
				max = val;
				maxd = d;
			}
			if (val < 90) less90++;
		} else {
			if (abs(val) >= 100000)
				cout<<"val, d"<<val<<" "<<d<<endl;
		}
	}
	cout<<"actual average is "<<avg/num<<endl;
	cout<<"max is "<<max<<" for molecule id "<<maxd<<endl;
	cout<<"less90 "<<less90<<endl;
}

void RDFTest::testReadRDFFile() {

	delete _domainDecomposition;
	// will be freed by the tearDown()-method.
	_domainDecomposition = new DomainDecompDummy();

	// to get proper components
	initializeFromFile(ParticleContainerFactory::LinkedCell,
			"Ethane_eqv.inp", 36.916625);
	ParticlePairs2PotForceAdapter handler(*_domain);
	const vector<Component>& components = _domain->getComponents();
	ASSERT_EQUAL((size_t) 1, components.size());

	RDF rdf(0.321254, 100, components);
	vector<double> rmids;
	vector<double> globalDist;
	vector<double> globalADist;
	vector<vector<double> > globalSiteDist;
	vector<vector<double> > globalSiteADist;
	//../test_input/rdf_Ethan_16k_0-0.000100100.rdf
	rdf.readRDFInputFile("rdf_Ethan_20k_0a_rc4_0-0.000090000.rdf", 0, 0, 2, 2,
			&rmids, &globalDist, &globalADist, &globalSiteDist,
			&globalSiteADist);

	FILE* rmidsFile = fopen("rmids", "w");
	FILE* globalDistFile = fopen("globalDistRDF.txt", "w");
	FILE* globalADistFile = fopen("globalADistRDF.txt", "w");
	FILE* globalSiteADistFile = fopen("globalSiteADistRDF.txt", "w");
	FILE* globalSiteDistFile = fopen("globalSiteDistRDF.txt", "w");
	for (unsigned int i = 0; i < rmids.size(); i++) {
		fprintf(rmidsFile, "%g \n", rmids[i]);
		fprintf(globalDistFile, "%g \n", globalDist[i]);
		fprintf(globalADistFile, "%g \n", globalADist[i]);
		for (unsigned int m = 0; m < components[0].numSites(); m++) {
			for (unsigned int n = 0; n < components[0].numSites(); n++) {
				fprintf(globalSiteDistFile, "%g ", globalSiteDist[m
						* components[0].numSites() + n][i]);
				fprintf(globalSiteADistFile, "%g ", globalSiteADist[m
						* components[0].numSites() + n][i]);
			}
		}
		fprintf(globalSiteDistFile, "\n");
		fprintf(globalSiteADistFile, "\n");
	}
	fclose(rmidsFile);
	fclose(globalDistFile);
	fclose(globalADistFile);
	fclose(globalSiteDistFile);
	fclose(globalSiteADistFile);

	FILE* expected = fopen("globalSiteDistRDF_expected.txt", "r");
	for (unsigned int i = 0; i < globalSiteDist.size(); i++) {
		double curr1, curr2, curr3, curr4;
		fscanf(expected, "%lg %lg %lg %lg\n", &curr1, &curr2, &curr3, &curr4);
		ASSERT_EQUAL(curr1, globalSiteDist[0][i]);
		ASSERT_EQUAL(curr2, globalSiteDist[1][i]);
		ASSERT_EQUAL(curr3, globalSiteDist[2][i]);
		ASSERT_EQUAL(curr4, globalSiteDist[3][i]);
	}
	fclose(expected);

}
