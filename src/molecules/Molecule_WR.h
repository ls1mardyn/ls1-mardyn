/*
 * Molecule_WR.h
 *
 *  Created on: 21 Jan 2017
 *      Author: tchipevn
 */

#ifndef SRC_MOLECULES_MOLECULE_WR_H_
#define SRC_MOLECULES_MOLECULE_WR_H_

#include "MoleculeInterface.h"
#include "particleContainer/adapter/vectorization/SIMD_TYPES.h"

class CellDataSoA_WR;

class Molecule_WR : public MoleculeInterface {
	enum StorageState {
		SOA = 0,
		AOS = 1
	};

public:
	Molecule_WR(unsigned long id = 0, Component *component = nullptr,
        double rx = 0., double ry = 0., double rz = 0.,
        double vx = 0., double vy = 0., double vz = 0.,
        double q0 = 0., double q1 = 0., double q2 = 0., double q3 = 0.,
        double Dx = 0., double Dy = 0., double Dz = 0.
	) {
		_state = AOS;
		_r[0] = rx;
		_r[1] = ry;
		_r[2] = rz;
		_v[0] = vx;
		_v[1] = vy;
		_v[2] = vz;
		_id = id;
		_soa = nullptr;

		if(not _initCalled) {
			initStaticVars();
		}
	}

	Molecule_WR(CellDataSoA_WR * soa, size_t index) {
		_state = SOA;
		_soa = soa;
		_soa_index = index;

		if (not _initCalled) {
			initStaticVars();
		}
	}

	~Molecule_WR() {}

	unsigned long id() const;
	void setid(unsigned long id);
	void setr(unsigned short d, double r);
	void setv(unsigned short d, double v);
	double r(unsigned short d) const;
	double v(unsigned short d) const;


	void setComponent(Component *component) {
		_component = component;
	}

	Component* component() const {
		return _component;
	}

	double F(unsigned short d) const {
		assert(false);
		return 0.0;
	}

	const Quaternion& q() const {
		assert(false);
		return _quaternion;
	}

	void setq(Quaternion q) {
		_quaternion = q;
	}

	double D(unsigned short d) const {
		return 0.0;
	}
	double M(unsigned short d) const {
		return 0.0;
	}
	double Vi(unsigned short d) const {
		return 0.0;
	}

	void setD(unsigned short d, double D) {}

	inline void move(int d, double dr) {
		setr(d, r(d) + dr);
	}

	double gMass() {
		return mass();
	}
	double getI(unsigned short d) const {
		assert(false);
		// TODO: check values for single-centered molecules
		return 0.0;
	}


	double U_rot() {
		return 0.0;
	}

	void setupSoACache(CellDataSoABase * const s, unsigned iLJ, unsigned iC, unsigned iD, unsigned iQ) {
		assert(false);
		// should this ever be like called?
		setSoA(s);
		_soa_index = iLJ;
	}

	void setSoA(CellDataSoABase * const s);

	void setStartIndexSoA_LJ(unsigned i) {
		assert(false);
	}
	void setStartIndexSoA_C(unsigned i) {
		assert(false);
	}
	void setStartIndexSoA_D(unsigned i) {
		assert(false);
	}
	void setStartIndexSoA_Q(unsigned i) {
		assert(false);
	}

	unsigned int numSites() const {
		return 1;
	}
	unsigned int numOrientedSites() const {
		return 0;
	}
	unsigned int numLJcenters() const {
		return 1;
	}
	unsigned int numCharges() const {
		return 0;
	}
	unsigned int numDipoles() const {
		return 0;
	}
	unsigned int numQuadrupoles() const {
		return 0;
	}

	std::array<double, 3> site_d(unsigned int i) const {
		assert(i == 0);
		return std::array<double, 3>( { r(0), r(1), r(2) });
	}

	std::array<double, 3> ljcenter_d(unsigned int i) const {
		assert(i == 0);
		return std::array<double, 3>( { r(0), r(1), r(2) });
	}
	std::array<double, 3> charge_d(unsigned int i) const {
		assert(false);
		return std::array<double, 3>({0.0, 0.0, 0.0});
	}
	std::array<double, 3> dipole_d(unsigned int i) const {
		assert(false);
		return std::array<double, 3>({0.0, 0.0, 0.0});
	}
	std::array<double, 3> quadrupole_d(unsigned int i) const {
		assert(false);
		return std::array<double, 3>({0.0, 0.0, 0.0});
	}

	std::array<double, 3> site_d_abs(unsigned int i) const {
		assert(i == 0);
		return std::array<double, 3>( { r(0), r(1), r(2) });
	}
	std::array<double, 3> ljcenter_d_abs(unsigned int i) const {
		assert(i == 0);
		return std::array<double, 3>( { r(0), r(1), r(2) });
	}
	std::array<double, 3> charge_d_abs(unsigned int i) const {
		assert(false);
		return std::array<double, 3>({0.0, 0.0, 0.0});
	}
	std::array<double, 3> dipole_d_abs(unsigned int i) const {
		assert(false);
		return std::array<double, 3>({0.0, 0.0, 0.0});
	}
	std::array<double, 3> quadrupole_d_abs(unsigned int i) const {
		assert(false);
		return std::array<double, 3>({0.0, 0.0, 0.0});
	}

	std::array<double, 3> dipole_e(unsigned int i) const {
		assert(false);
		return std::array<double, 3>({0.0, 0.0, 0.0});
	}
	std::array<double, 3> quadrupole_e(unsigned int i) const {
		assert(false);
		return std::array<double, 3>({0.0, 0.0, 0.0});
	}

	std::array<double, 3> site_F(unsigned int i) const {
		assert(false);
		return std::array<double, 3>({0.0, 0.0, 0.0});
	}
	std::array<double, 3> ljcenter_F(unsigned int i) const {
		assert(false);
		return std::array<double, 3>({0.0, 0.0, 0.0});
	}
	std::array<double, 3> charge_F(unsigned int i) const {
		assert(false);
		return std::array<double, 3>({0.0, 0.0, 0.0});
	}
	std::array<double, 3> dipole_F(unsigned int i) const {
		assert(false);
		return std::array<double, 3>({0.0, 0.0, 0.0});
	}
	std::array<double, 3> quadrupole_F(unsigned int i) const {
		assert(false);
		return std::array<double, 3>({0.0, 0.0, 0.0});
	}

	void normalizeQuaternion() {}
	void computeLJcenter_d(unsigned int i, double result[3]) const {}
	void computeCharge_d(unsigned int i, double result[3]) const {}
	void computeDipole_d(unsigned int i, double result[3]) const {}
	void computeQuadrupole_d(unsigned int i, double result[3]) const {}
	void computeDipole_e(unsigned int i, double result[3]) const {}
	void computeQuadrupole_e(unsigned int i, double result[3]) const {}

	unsigned long totalMemsize() const {
		//todo: check
		assert(false);
		return sizeof(*this);
	}

	void setF(double F[3]) {}
	void setM(double M[3]) {}
	void setVi(double Vi[3]) {}
	void scale_v(double s) {}
	void scale_F(double s) {}
	void scale_D(double s) {}
	void scale_M(double s) {}
	void Fadd(const double a[]) {}
	void Madd(const double a[]) {}
	void Viadd(const double a[]) {}
	void vadd(const double ax, const double ay, const double az) {}
	void vsub(const double ax, const double ay, const double az) {}
	void Fljcenteradd(unsigned int i, double a[]) {}
	void Fljcentersub(unsigned int i, double a[]) {}
	void Fchargeadd(unsigned int i, double a[]) {}
	void Fchargesub(unsigned int i, double a[]) {}
	void Fdipoleadd(unsigned int i, double a[]) {}
	void Fdipolesub(unsigned int i, double a[]) {}
	void Fquadrupoleadd(unsigned int i, double a[]) {}
	void Fquadrupolesub(unsigned int i, double a[]) {}
	void upd_preF(double dt) {}
	void upd_postF(double dt_halve, double& summv2, double& sumIw2) {}
	void calculate_mv2_Iw2(double& summv2, double& sumIw2) {}
	void calculate_mv2_Iw2(double& summv2, double& sumIw2, double offx, double offy, double offz) {}
	static std::string getWriteFormat(); // TODO
	void write(std::ostream& ostrm) const {}
	void clearFM() {}
	void calcFM() {}
	void check(unsigned long id) {}


private:
	static void initStaticVars();

    static Component *_component;  /**< IDentification number of its component type */
    static Quaternion _quaternion;
    static bool _initCalled;

    StorageState _state;

    // if the state is AOS, the following values are read:
    vcp_real_calc _r[3];  /**< position coordinates */
	vcp_real_calc _v[3];  /**< velocity */
	unsigned long _id;

	// if the state is SOA, the values are read from the SoA:
	CellDataSoA_WR * _soa;
	size_t _soa_index;
};

#endif /* SRC_MOLECULES_MOLECULE_WR_H_ */
