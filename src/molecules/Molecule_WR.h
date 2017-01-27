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
#include "particleContainer/adapter/CellDataSoA_WR.h"

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
	}

	Molecule_WR(CellDataSoA_WR * soa, size_t index) {
		_state = SOA;
		_soa = soa;
		_soa_index = index;
	}

	~Molecule_WR() {}

	double r(unsigned short d) const {
		assert(_state == SOA or _state == AOS);

		if (_state == AOS) {
			return _r[d];
		} else {
			size_t linOffset = _soa->_mol_r.dimensionToOffset(d);
			return _soa->_mol_r.linearCrossAccess(linOffset + _soa_index);
		}
	}

	double v(unsigned short d) const {
		assert(_state == SOA or _state == AOS);

		if (_state == AOS) {
			return _v[d];
		} else {
			size_t linOffset = _soa->_mol_v.dimensionToOffset(d);
			return _soa->_mol_v.linearCrossAccess(linOffset + _soa_index);
		}
	}

	unsigned long id() const {
		assert(_state == SOA or _state == AOS);

		if (_state == AOS) {
			return _id;
		} else {
			return _soa->_mol_uid[_soa_index];
		}
	}

	void setr(unsigned short d, double r) {
		assert(_state == SOA or _state == AOS);

		if (_state == AOS) {
			_r[d] = r;
		} else {
			size_t linOffset = _soa->_mol_r.dimensionToOffset(d);
			_soa->_mol_r.linearCrossAccess(linOffset + _soa_index) = r;
		}
	}

	void setv(unsigned short d, double v) {
		assert(_state == SOA or _state == AOS);

		if (_state == AOS) {
			_v[d] = v;
		} else {
			size_t linOffset = _soa->_mol_r.dimensionToOffset(d);
			_soa->_mol_v.linearCrossAccess(linOffset + _soa_index) = v;
		}
	}

	void setid(unsigned long id) {
		assert(_state == SOA or _state == AOS);

		if (_state == AOS) {
			_id = id;
		} else {
			_soa->_mol_uid[_soa_index] = id;
		}
	}

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

	void setSoA(CellDataSoABase * const s) {
		assert(_state == AOS);
		CellDataSoA_WR * derived;
#ifndef NDEBUG
		derived = nullptr;
		derived = dynamic_cast<CellDataSoA_WR *>(s);
		if(derived == nullptr) {
			global_log->error() << "expected CellDataSoA_WR pointer for m" << _id << std::endl;
			assert(false);
		}
#else
		derived = static_cast<CellDataSoA *>(s);
#endif
		_soa = derived;
	}

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

	unsigned long totalMemsize() const = 0;

	void setF(double F[3]) = 0;
	void setM(double M[3]) = 0;
	void setVi(double Vi[3]) = 0;
	void scale_v(double s) = 0;
	void scale_v(double s, double offx, double offy, double offz);
	void scale_F(double s) = 0;
	void scale_D(double s) = 0;
	void scale_M(double s) = 0;
	void Fadd(const double a[]) = 0;
	void Madd(const double a[]) = 0;
	void Viadd(const double a[]) = 0;
	void vadd(const double ax, const double ay, const double az) = 0;
	void vsub(const double ax, const double ay, const double az) = 0;
	void Fljcenteradd(unsigned int i, double a[]) = 0;
	void Fljcentersub(unsigned int i, double a[]) = 0;
	void Fchargeadd(unsigned int i, double a[]) = 0;
	void Fchargesub(unsigned int i, double a[]) = 0;
	void Fdipoleadd(unsigned int i, double a[]) = 0;
	void Fdipolesub(unsigned int i, double a[]) = 0;
	void Fquadrupoleadd(unsigned int i, double a[]) = 0;
	void Fquadrupolesub(unsigned int i, double a[]) = 0;
	void upd_preF(double dt) = 0;
	void upd_postF(double dt_halve, double& summv2, double& sumIw2) = 0;
	void calculate_mv2_Iw2(double& summv2, double& sumIw2) = 0;
	void calculate_mv2_Iw2(double& summv2, double& sumIw2, double offx, double offy, double offz) = 0;
	static std::string getWriteFormat(); // TODO
	void write(std::ostream& ostrm) const = 0;
	void clearFM() = 0;
	void calcFM() = 0;
	void check(unsigned long id) = 0;


private:
    static Component *_component;  /**< IDentification number of its component type */
    static Quaternion _quaternion;

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
