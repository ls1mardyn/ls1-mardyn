/*
 * WignerRotationTest.cpp
 *
 *  Created on: Jun 16, 2015
 *      Author: uwe
 */

#include "WignerRotationTest.h"
#include <cmath>
#include <iostream>
#include "bhfmm/pseudoParticles/SHMultipoleParticle.h"
#include "bhfmm/pseudoParticles/SHLocalParticle.h"
#include "bhfmm/expansions/SolidHarmonicsExpansion.h"
#include "bhfmm/utils/WignerMatrix.h"
#include "bhfmm/utils/RotationParameterLookUp.h"
#include "bhfmm/utils/Vector3.h"

TEST_SUITE_REGISTRATION(WignerRotationTest);

WignerRotationTest::WignerRotationTest() {
}

WignerRotationTest::~WignerRotationTest() {
}

void WignerRotationTest::testM2MWignerRotation() {

	/* order of expansion */
	const int order = 5;

	/* setup a SHMultipoleParticle */
	bhfmm::SHMultipoleParticle MP(order); // center is implicitly (0.0, 0.0, 0.0)
	MP.setRadius(4.0);

	bhfmm::Vector3<double> particle_position; // source particle
	double charge;


	particle_position[0] = 1.2;
	particle_position[1] = -1.3;
	particle_position[2] = 1.4;
	charge = 1.0;
	MP.addSource(particle_position, charge);

	particle_position[0] = -1.5;
	particle_position[1] = -0.6;
	particle_position[2] = -1.2;
	charge = 1.0;
	MP.addSource(particle_position, charge);

	/* setup translation vector */
	bhfmm::Vector3<double> translation_vector;
	translation_vector[0] = -8.0;
	translation_vector[1] = -8.0;
	translation_vector[2] = -8.0;

	/* setup rotation angles and z-translation vector */
	double projXY_len = sqrt(
			translation_vector[0] * translation_vector[0]
					+ translation_vector[1] * translation_vector[1]);

	double len = translation_vector.L2Norm();
	bhfmm::Vector3<double> translation_vector_Z;
	translation_vector_Z[0] = 0.0;
	translation_vector_Z[1] = 0.0;
	translation_vector_Z[2] = len;

	double phi = atan2(translation_vector[1], translation_vector[0]);
	double theta = atan2(projXY_len, translation_vector[2]);

	/* setup Sin(phi)/Cos(phi) lookup */
	double* CosSin = new double[(order+1)*2];
	for (int m = 0; m <= order; ++m) {
		CosSin[2*m] 		= cos(m*phi);
		CosSin[2*m + 1] 	= sin(m*phi);
	}

	/* translate without Wigner-rotation */
	bhfmm::SolidHarmonicsExpansion L_trafo = convoluteLL(MP.getExpansion(), evaluateLOfR(order, translation_vector));

	/* translate using Wigner-rotation */
	bhfmm::WignerMatrix W_pos(order, true);
	bhfmm::WignerMatrix W_neg(order, true);

	W_pos.evaluate(theta);
	W_neg.evaluate(-theta);

	/* initialize rotation parameter table */
	bhfmm::RotationParameterLookUp::tab = new bhfmm::RotationParameterLookUp(order);
	bhfmm::RotationParameterLookUp::tab->initFromDirEval();

	// pre-multiply prefactors
	for (int l = 0; l <= order; ++l) {
		for (int m = 0; m <= l; ++m) {
			for (int k = -l; k<=l; ++k) {
				const double factor = bhfmm::RotationParameterLookUp::tab->acc_c(l,m,k);
				W_pos.acc(l,m,k) *= factor;
				W_neg.acc(l,m,k) *= factor;
			}
		}
	}

	/* clean up */
	delete bhfmm::RotationParameterLookUp::tab;


	bhfmm::SolidHarmonicsExpansion L_trafo_Wigner = rotatePhi(MP.getExpansion(), CosSin, 1);

	L_trafo_Wigner = rotateThetaL(L_trafo_Wigner, W_pos);
	L_trafo_Wigner = convoluteLL_Z(L_trafo_Wigner, evaluateLOfR(order, translation_vector_Z));
	L_trafo_Wigner = rotateThetaL(L_trafo_Wigner, W_neg);
	L_trafo_Wigner = rotatePhi(L_trafo_Wigner, CosSin, -1);

	/* calculate difference */
	double abs_diff = 0.0;
	const unsigned num_entries = L_trafo.getNumEntries();

	double* C_L_trafo = &L_trafo.getC(0, 0);
	double* S_L_trafo = &L_trafo.getS(0, 0);

	double* C_L_trafo_Wigner = &L_trafo_Wigner.getC(0, 0);
	double* S_L_trafo_Wigner = &L_trafo_Wigner.getS(0, 0);

	// compare C parts
	for(unsigned i=0; i<num_entries/2; ++i) {
		abs_diff += std::abs(C_L_trafo[i] - C_L_trafo_Wigner[i]);
	}
	// compare S parts
	for(unsigned i=0; i<num_entries/2; ++i) {
		abs_diff += std::abs(S_L_trafo[i] - S_L_trafo_Wigner[i]);
	}

//	MP.getExpansion().print();
//	L_trafo_Wigner.print();
//	L_trafo.print();

	/* assertion */
	ASSERT_DOUBLES_EQUAL(abs_diff, 0.0, 1e-10);

}

void WignerRotationTest::testL2LWignerRotation() {

	/* order of expansion */
	const int order = 5;

	/* setup a SHMultipoleParticle */
	bhfmm::SHLocalParticle MP(order); // center is implicitly (0.0, 0.0, 0.0)
	MP.setRadius(1.0);

	bhfmm::Vector3<double> particle_position; // source particle
	double charge;


	particle_position[0] = 1.2;
	particle_position[1] = -1.3;
	particle_position[2] = 1.4;
	charge = 1.0;
	MP.addSource(particle_position, charge);

	particle_position[0] = -1.5;
	particle_position[1] = -0.6;
	particle_position[2] = -1.2;
	charge = 1.0;
	MP.addSource(particle_position, charge);

	/* setup translation vector */
	bhfmm::Vector3<double> translation_vector;
	translation_vector[0] = -1.0;
	translation_vector[1] = -1.0;
	translation_vector[2] = -1.0;

	/* setup rotation angles and z-translation vector */
	double projXY_len = sqrt(
			translation_vector[0] * translation_vector[0]
					+ translation_vector[1] * translation_vector[1]);

	double len = translation_vector.L2Norm();
	bhfmm::Vector3<double> translation_vector_Z;
	translation_vector_Z[0] = 0.0;
	translation_vector_Z[1] = 0.0;
	translation_vector_Z[2] = len;

	double phi = atan2(translation_vector[1], translation_vector[0]);
	double theta = atan2(projXY_len, translation_vector[2]);

	/* setup Sin(phi)/Cos(phi) lookup */
	double* CosSin = new double[(order+1)*2];
	for (int m = 0; m <= order; ++m) {
		CosSin[2*m] 		= cos(m*phi);
		CosSin[2*m + 1] 	= sin(m*phi);
	}

	/* translate without Wigner-rotation */
	bhfmm::SolidHarmonicsExpansion M_trafo = convoluteLM(evaluateLOfR(order, translation_vector), MP.getExpansion());

	/* translate using Wigner-rotation */
	bhfmm::WignerMatrix W_pos(order, true);
	bhfmm::WignerMatrix W_neg(order, true);

	W_pos.evaluate(theta);
	W_neg.evaluate(-theta);

	/* initialize rotation parameter table */
	bhfmm::RotationParameterLookUp::tab = new bhfmm::RotationParameterLookUp(order);
	bhfmm::RotationParameterLookUp::tab->initFromDirEval();

	// pre-multiply prefactors
	for (int l = 0; l <= order; ++l) {
		for (int m = 0; m <= l; ++m) {
			for (int k = -l; k<=l; ++k) {
				const double factor = bhfmm::RotationParameterLookUp::tab->acc_c(l,k,m);
				W_pos.acc(l,m,k) *= factor;
				W_neg.acc(l,m,k) *= factor;
			}
		}
	}

	/* clean up */
	delete bhfmm::RotationParameterLookUp::tab;

	bhfmm::SolidHarmonicsExpansion M_trafo_Wigner = rotatePhi(MP.getExpansion(), CosSin, +1);

	M_trafo_Wigner = rotateThetaM(M_trafo_Wigner, W_pos);
	M_trafo_Wigner = convoluteL_ZM(evaluateLOfR(order, translation_vector_Z), M_trafo_Wigner);
	M_trafo_Wigner = rotateThetaM(M_trafo_Wigner, W_neg);
	M_trafo_Wigner = rotatePhi(M_trafo_Wigner, CosSin, -1);

	/* calculate difference */
	double abs_diff = 0.0;
	const unsigned num_entries = M_trafo.getNumEntries();

	double* C_M_trafo = &M_trafo.getC(0, 0);
	double* S_M_trafo = &M_trafo.getS(0, 0);

	double* C_M_trafo_Wigner = &M_trafo_Wigner.getC(0, 0);
	double* S_M_trafo_Wigner = &M_trafo_Wigner.getS(0, 0);

	// compare C parts
	for(unsigned i=0; i<num_entries/2; ++i) {
		abs_diff += std::abs(C_M_trafo[i] - C_M_trafo_Wigner[i]);
	}
	// compare S parts
	for(unsigned i=0; i<num_entries/2; ++i) {
		abs_diff += std::abs(S_M_trafo[i] - S_M_trafo_Wigner[i]);
	}

//	MP.getExpansion().print();
//	M_trafo_Wigner.print();
//	M_trafo.print();

	/* assertion */
	ASSERT_DOUBLES_EQUAL(abs_diff, 0.0, 1e-10);

}

void WignerRotationTest::testM2LWignerRotation() {
	/* order of expansion */
	const int order = 5;

	/* setup a SHMultipoleParticle */
	bhfmm::SHMultipoleParticle MP(order); // center is implicitly (0.0, 0.0, 0.0)
	MP.setRadius(4.0);

	bhfmm::Vector3<double> particle_position; // source particle
	double charge;


	particle_position[0] = 1.2;
	particle_position[1] = -1.3;
	particle_position[2] = 1.4;
	charge = 1.0;
	MP.addSource(particle_position, charge);

//	particle_position[0] = -1.5;
//	particle_position[1] = -0.6;
//	particle_position[2] = -1.2;
//	charge = 1.0;
//	MP.addSource(particle_position, charge);

	/* setup translation vector */
	bhfmm::Vector3<double> translation_vector;
	translation_vector[0] = 0.0;
	translation_vector[1] = 60.4713;
	translation_vector[2] = 0.0;

	/* setup rotation angles and z-translation vector */
	double projXY_len = sqrt(
			translation_vector[0] * translation_vector[0]
					+ translation_vector[1] * translation_vector[1]);

	double len = translation_vector.L2Norm();
	bhfmm::Vector3<double> translation_vector_Z;
	translation_vector_Z[0] = 0.0;
	translation_vector_Z[1] = 0.0;
	translation_vector_Z[2] = len;

	double phi = atan2(translation_vector[1], translation_vector[0]);
	double theta = atan2(projXY_len, translation_vector[2]);

//	std::cout << "phi: " << phi*180.0/M_PI << std::endl;
//	std::cout << "theta: " << theta*180.0/M_PI << std::endl;

	/* setup Sin(phi)/Cos(phi) lookup */
	double* CosSin = new double[(order+1)*2];
	for (int m = 0; m <= order; ++m) {
		CosSin[2*m] 		= cos(m*phi);
		CosSin[2*m + 1] 	= sin(m*phi);
	}

	/* translate without Wigner-rotation */
	bhfmm::SolidHarmonicsExpansion LM_trafo = convoluteLM(setAtMinusR(MP.getConstExpansion()), evaluateMOfR(order, translation_vector));

	/* translate using Wigner-rotation */
	bhfmm::WignerMatrix W_pos(order, true);
	bhfmm::WignerMatrix W_neg(order, true);

	W_pos.evaluate(theta);
	W_neg.evaluate(-theta);

	/* initialize rotation parameter table */
	bhfmm::RotationParameterLookUp::tab = new bhfmm::RotationParameterLookUp(order);
	bhfmm::RotationParameterLookUp::tab->initFromDirEval();

	// pre-multiply prefactors
	for (int l = 0; l <= order; ++l) {
		for (int m = 0; m <= l; ++m) {
			for (int k = -l; k<=l; ++k) {
				W_pos.acc(l,m,k) *= bhfmm::RotationParameterLookUp::tab->acc_c(l,m,k);
				W_neg.acc(l,m,k) *= bhfmm::RotationParameterLookUp::tab->acc_c(l,k,m);
			}
		}
	}

	/* clean up */
	delete bhfmm::RotationParameterLookUp::tab;

	bhfmm::SolidHarmonicsExpansion LM_trafo_Wigner = rotatePhi(setAtMinusR(MP.getConstExpansion()), CosSin, +1);
	LM_trafo_Wigner = rotateThetaL(LM_trafo_Wigner, W_pos);
	LM_trafo_Wigner = convoluteLM_Z(LM_trafo_Wigner, evaluateMOfR(order, translation_vector_Z));
	LM_trafo_Wigner = rotateThetaM(LM_trafo_Wigner, W_neg);
	LM_trafo_Wigner = rotatePhi(LM_trafo_Wigner, CosSin, -1);

	/* calculate difference */
	double abs_diff = 0.0;
	const unsigned num_entries = LM_trafo.getNumEntries();

	double* C_M_trafo = &LM_trafo.getC(0, 0);
	double* S_M_trafo = &LM_trafo.getS(0, 0);

	double* C_M_trafo_Wigner = &LM_trafo_Wigner.getC(0, 0);
	double* S_M_trafo_Wigner = &LM_trafo_Wigner.getS(0, 0);

	// compare C parts
	for(unsigned i=0; i<num_entries/2; ++i) {
		abs_diff += std::abs(C_M_trafo[i] - C_M_trafo_Wigner[i]);
	}
	// compare S parts
	for(unsigned i=0; i<num_entries/2; ++i) {
		abs_diff += std::abs(S_M_trafo[i] - S_M_trafo_Wigner[i]);
	}

//	MP.getExpansion().print();

//	LM_trafo.print();
//	LM_trafo_Wigner.print();

	/* assertion */
	ASSERT_DOUBLES_EQUAL(abs_diff, 0.0, 1e-10);
}
