/*
 * PrincipalAxisTransform.cpp
 *
 * @Date: 17.09.2011
 * @Author: eckhardw
 */

#include "PrincipalAxisTransform.h"
#include "molecules/Component.h"
#include "eig3.h"


/**
 * @param[in] site
 * @param[out] totalMass
 * @param[out] centerOfMass
 */
static void sumProperties(const Site& site, double& totalMass, double centerOfMass[3]);

/**
 * @param [out] matrix 3x3 matrix for moment of inertia
 */
static void shiftCenterAndBuildMatrix(Site& site, const double centerOfMass[3], double matrix[3][3]);

static void transformCoordinatesAndBuildMatrix(Site& s, const double eigenVectors[3][3], double A[3][3]);

void principalAxisTransform(Component& component) {

	double totalMass = 0.;
	double centerOfMass[] = {0., 0., 0.};
	for (size_t i = 0; i < component.numLJcenters(); i++) {
		sumProperties(component.ljcenter(i), totalMass, centerOfMass);
	}
	for (size_t i = 0; i < component.numCharges(); i++) {
		sumProperties(component.charge(i), totalMass, centerOfMass);
	}

	for (int i = 0; i < 3; i++) {
		centerOfMass[i] = centerOfMass[i] / totalMass;
	}

	std::cout << "Moving center of mass [" << centerOfMass[0] << "," << centerOfMass[1] << "," << centerOfMass[2] << "]" << std::endl;

	// adjust coordinates wrt center of mass and assemble mass matrix
	double momentMatrix[3][3] = { {0., 0., 0.} , {0., 0., 0.} , {0., 0., 0.} };
	for (size_t i = 0; i < component.numLJcenters(); i++) {
		shiftCenterAndBuildMatrix(component.ljcenter(i), centerOfMass, momentMatrix);
	}
	for (size_t i = 0; i < component.numCharges(); i++) {
		shiftCenterAndBuildMatrix(component.charge(i), centerOfMass, momentMatrix);
	}
	for (size_t i = 0; i < component.numDipoles(); i++) {
		shiftCenterAndBuildMatrix(component.dipole(i), centerOfMass, momentMatrix);
	}
	for (size_t i = 0; i < component.numQuadrupoles(); i++) {
		shiftCenterAndBuildMatrix(component.quadrupole(i), centerOfMass, momentMatrix);
	}

	std::cout << "momentMatrix:" << std::endl;
	for (int i = 0; i < 3; i++) {
		std::cout << "[ " << std::endl;
		for (int j = 0; j < 3; j++) {
			std::cout << " " << momentMatrix[i][j];
		}
		std::cout << "]" << std::endl;
	}

	double eigenVectors[3][3];
	double eigenValues[3];
	eigen_decomposition(momentMatrix, eigenVectors, eigenValues);

	// the code I use has a different sign in the second column
	// this is original: at least second and third vector are exchanged
	// eigenVectors[0][1] = - eigenVectors[0][1];
	// eigenVectors[1][1] = - eigenVectors[1][1];
	// eigenVectors[2][1] = - eigenVectors[2][1];
	for (int i = 0; i < 3; i++) {
		double tmp = eigenVectors[1][i];
		eigenVectors[1][i] = eigenVectors[2][i];
		eigenVectors[2][i] = tmp;
	}

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			momentMatrix[i][j] = 0.;
		}
	}

	std::cout << "Eigenvectors:" << std::endl;
	for (int i = 0; i < 3; i++) {
		std::cout << "[ ";
		for (int j = 0; j < 3; j++) {
			std::cout << " " << eigenVectors[i][j];
		}
		std::cout << "]" << std::endl;
	}

	// adjust coordinates wrt. principal axis
	for (size_t i = 0; i < component.numLJcenters(); i++) {
		transformCoordinatesAndBuildMatrix(component.ljcenter(i), eigenVectors, momentMatrix);
	}
	for (size_t i = 0; i < component.numCharges(); i++) {
		transformCoordinatesAndBuildMatrix(component.charge(i), eigenVectors, momentMatrix);
	}
	for (size_t i = 0; i < component.numDipoles(); i++) {
		transformCoordinatesAndBuildMatrix(component.dipole(i), eigenVectors, momentMatrix);
	}
	for (size_t i = 0; i < component.numQuadrupoles(); i++) {
		transformCoordinatesAndBuildMatrix(component.quadrupole(i), eigenVectors, momentMatrix);
	}
	component.setI11(momentMatrix[0][0]);
	component.setI22(momentMatrix[1][1]);
	component.setI33(momentMatrix[2][2]);


	std::cout << "Diagonal: " << momentMatrix[0][0] << "," << momentMatrix[1][1] << "," << momentMatrix[2][2] <<std::endl;
}

static void sumProperties(const Site& site, double& totalMass, double centerOfMass[3]) {
	totalMass += site.m();
	centerOfMass[0] += site.m() * site.rx();
	centerOfMass[1] += site.m() * site.ry();
	centerOfMass[2] += site.m() * site.rz();
}

static void shiftCenterAndBuildMatrix(Site& s, const double centerOfMass[3], double A[3][3]) {
	const std::array<double, 3> r = s.r();
	for (int j = 0; j < 3; j++) {
		s.setR(j, r[j] - centerOfMass[j]);
	}
	A[0][0]+= (s.ry() * s.ry() + s.rz() * s.rz()) * s.m(); // (s.y*s.y+s.z*s.z)*s.mass
	A[0][1]-= s.rx() * s.ry() * s.m(); // s.x*s.y*s.mass
	A[0][2]-= s.rx() * s.rz() * s.m();// s.x*s.z*s.mass
	A[1][1]+= (s.rx() * s.rx() + s.rz() * s.rz()) * s.m();// (s.x*s.x+s.z*s.z)*s.mass
	A[1][2]-= s.ry() * s.rz() * s.m();// s.y*s.z*s.mass
	A[2][2]+= (s.rx() * s.rx() + s.ry() * s.ry()) * s.m();// (s.x*s.x+s.y*s.y)*s.mass
	A[1][0] = A[0][1];
	A[2][0] = A[0][2];
	A[2][1] = A[1][2];
}


static void transformCoordinatesAndBuildMatrix(Site& s, const double eigvec[3][3], double A[3][3]) {
	double x = s.rx();
	double y = s.ry();
	double z = s.rz();
	s.setR(0, eigvec[0][0] * x + eigvec[0][1] * y + eigvec[0][2] * z);
	s.setR(1, eigvec[1][0] * x + eigvec[1][1] * y + eigvec[1][2] * z);
	s.setR(2, eigvec[2][0] * x + eigvec[2][1] * y + eigvec[2][2] * z);
	A[0][0] += (s.ry() * s.ry() + s.rz() * s.rz()) * s.m();
	A[0][1] -= s.rx() * s.ry() * s.m();
	A[0][2] -= s.rx() * s.rz() * s.m();
	A[1][1] += (s.rx() * s.rx() + s.rz() * s.rz()) * s.m();
	A[1][2] -= s.ry() * s.rz() * s.m();
	A[2][2] += (s.rx() * s.rx() + s.ry() * s.ry()) * s.m();
	A[1][0] = A[0][1];
	A[2][0] = A[0][2];
	A[2][1] = A[1][2];
}
