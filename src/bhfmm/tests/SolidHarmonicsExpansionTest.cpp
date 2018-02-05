/*
 * SolidHarmonicsExpansionTest.cpp
 *
 *  Created on: 2 Feb 2018
 *      Author: tchipevn
 */

#include "SolidHarmonicsExpansionTest.h"
#include "bhfmm/pseudoParticles/SHMultipoleParticle.h"
#include "bhfmm/pseudoParticles/SHLocalParticle.h"

using namespace bhfmm;

TEST_SUITE_REGISTRATION(SolidHarmonicsExpansionTest);

SolidHarmonicsExpansionTest::SolidHarmonicsExpansionTest() {
}

SolidHarmonicsExpansionTest::~SolidHarmonicsExpansionTest() {
}

void localpotforce(Vector3<double> source_pos, Vector3<double> target_pos, double q1, double q2, double & u, Vector3<double>& f) {
	Vector3<double> dr = target_pos - source_pos;
	double dr2 = dr.L2NormSquare();

	double invdr2 = 1.0 / dr2;
	double invdr = sqrt(invdr2);

	double q1q2 = q1 * q2;
	u = q1q2 * invdr;
	const double fac = u * invdr2;
	for (unsigned short d = 0; d < 3; d++)
		f[d] = fac * dr[d];
}

double calc_err(double exact_u, Vector3<double> exact_f, double u, Vector3<double> f) {

	double abs_err_u = std::abs(exact_u - u);
	double rel_err_u = std::abs(abs_err_u / exact_u);
//	test_log->info() << "Abs Err Pot SH:   " << abs_err_u << std::endl;
//	test_log->info() << "Rel Err Pot SH:   " << rel_err_u << std::endl;

	double abs_err_f = (exact_f - f).L2Norm();
	double rel_err_f = abs_err_f / f.L2Norm();
//	test_log->info() << "Abs Err Force SH: " << abs_err_u << std::endl;
//	test_log->info() << "Rel Err Force SH: " << rel_err_f << std::endl;

	return abs_err_u;
}

double pot_err_bound(Vector3<double> source_pos, double source_q, Vector3<double> target_pos, int order) {
	// nach Kabadshow's diss, Section 2.7.1.
	double a = source_pos.L2Norm();
	double r = target_pos.L2Norm();
	double temp = std::pow((a/r), order+1);
	double err_bound = std::abs(source_q) / (r - a) * temp;
	return err_bound;
}

void SolidHarmonicsExpansionTest::test_P2M_M2P() {
	//Source Particle
	double source_q = 1.0;
	double source_p[3] = {0.3, -0.2, 0.1};
	Vector3<double> source_p_vec(source_p);

	//Target Particle
	double target_q = -1.0;
	double target_p[3] = {2.0, -2.0, 1.3};
	Vector3<double> target_p_vec(target_p);

	double m_p[3] = {0.4, -0.1, 0.2};
	Vector3<double> m_p_vec(m_p);

//	test_log->info() << "Testing test_P2M_M2P." << std::endl;

	double ex_pot = 0.0;
	Vector3<double> ex_force(0.0);
	localpotforce(source_p_vec, target_p_vec, source_q, target_q, ex_pot, ex_force);

//	test_log->info() << "exact pot & force: " << ex_pot << " " << ex_force << std::endl;

	// create a multipole expansion
	for (int order = 0; order <= 10; ++order) {
		Vector3<double> target_p_force(0.0);
		double target_u = 0.0;

		SHMultipoleParticle multipole(order);
		multipole.setCenter(m_p_vec);

		multipole.addSource(source_p_vec, source_q);

		// compute potential
		multipole.actOnTarget(target_p_vec, target_q, target_u, target_p_force);

//		test_log->info() << "      pot & force: " << target_u << " " << target_p_force << std::endl;
		double err_bound = pot_err_bound(source_p_vec, source_q, target_p_vec, order);
//		test_log->info() << "abs error bound for potential:" << err_bound << std::endl;
		double abs_err = calc_err(ex_pot, ex_force, target_u, target_p_force);
//		test_log->info() << "abs error       for potential:" << abs_err << std::endl;
		ASSERT_TRUE(abs_err <= err_bound);
	}
}

void SolidHarmonicsExpansionTest::test_P2L_L2P() {
	//Source Particle
	double source_q = 1.0;
	double source_p[3] = {0.3, -0.2, 0.1};
	Vector3<double> source_p_vec(source_p);

	//Target Particle
	double target_q = -1.0;
	double target_p[3] = {2.0, -2.0, 1.3};
	Vector3<double> target_p_vec(target_p);

	double l_p[3] = {2.1, -1.9, 1.5};
	Vector3<double> l_p_vec(l_p);

//	test_log->info() << "Testing test_P2L_L2P." << std::endl;

	double ex_pot = 0.0;
	Vector3<double> ex_force(0.0);
	localpotforce(source_p_vec, target_p_vec, source_q, target_q, ex_pot, ex_force);

//	test_log->info() << "exact pot & force: " << ex_pot << " " << ex_force << std::endl;

	// create a multipole expansion
	for (int order = 0; order <= 10; ++order) {
		Vector3<double> target_p_force(0.0);
		double target_u = 0.0;

		SHLocalParticle local(order);
		local.setCenter(l_p_vec);

		local.addSource(source_p_vec, source_q);

		// compute potential
		local.actOnTarget(target_p_vec, target_q, target_u, target_p_force);

//		test_log->info() << "      pot & force: " << target_u << " " << target_p_force << std::endl;
		double err_bound = pot_err_bound(source_p_vec, source_q, target_p_vec, order);
//		test_log->info() << "abs error bound for potential:" << err_bound << std::endl;
		double abs_err = calc_err(ex_pot, ex_force, target_u, target_p_force);
//		test_log->info() << "abs error       for potential:" << abs_err << std::endl;
		ASSERT_TRUE(abs_err <= err_bound);
	}
}
