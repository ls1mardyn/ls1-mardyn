#ifndef SRC_MOLECULES_MOLECULERMM_H_
#define SRC_MOLECULES_MOLECULERMM_H_

#include "MoleculeInterface.h"
#include "particleContainer/adapter/vectorization/SIMD_TYPES.h"

#ifdef UNIT_TESTS
#define DEBUG_FUNCTIONALITY_HACKS
#endif

class CellDataSoARMM;

class MoleculeRMM : public MoleculeInterface {
public:
	enum StorageState {
		STORAGE_SOA = 0,
		STORAGE_AOS = 1
	};

public:
	MoleculeRMM(unsigned long id = 0, Component *component = nullptr,
        double rx = 0., double ry = 0., double rz = 0.,
        double vx = 0., double vy = 0., double vz = 0.,
        double = 1., double = 0., double = 0., double = 0., /*q0, q1, q2, q3*/
        double = 0., double = 0., double = 0. /*Dx, Dy, Dz*/
	) {
		_state = STORAGE_AOS;
		_r[0] = rx;
		_r[1] = ry;
		_r[2] = rz;
		_v[0] = vx;
		_v[1] = vy;
		_v[2] = vz;
		_id = id;
		_soa = nullptr;
		_soa_index = 0;

		if(component != nullptr) {
			_component = component;
			_quaternion = Quaternion(1.0, 0.0, 0.0, 0.0);
			_initCalled = true;
		} else if(not _initCalled) {
			initStaticVars();
		}
	}

	// copy constructor should always create an AOS molecule?
	MoleculeRMM(const MoleculeRMM& other) {
		_state = STORAGE_AOS;
		for (int d = 0; d < 3; ++d) {
			setr(d, other.r(d));
			setv(d, other.v(d));
		}
		_id = other.getID();
		_soa = nullptr;
		_soa_index = 0;
	}

	MoleculeRMM(CellDataSoARMM * soa, size_t index) {
		_state = STORAGE_SOA;
		_soa = soa;
		_soa_index = index;

		if (not _initCalled) {
			initStaticVars();
		}
	}

	~MoleculeRMM() override = default;

	unsigned long getID() const override;
	void setid(unsigned long id) override;
	void setr(unsigned short d, double r) override;
	void setv(unsigned short d, double v) override;
	void setF(unsigned short, double F) override;
	double r(unsigned short d) const override;
	double v(unsigned short d) const override;


	void setComponent(Component *component) override {
		_component = component;
		this->updateMassInertia();
	}

	Component* component() const override {
		return _component;
	}

#ifndef DEBUG_FUNCTIONALITY_HACKS
	double F(unsigned short /*d*/) const override {
		mardyn_assert(false);
		return 0.0;
	}
#else
	double F(unsigned short d) const override { return v(d);}
#endif

	const Quaternion& q() const override {
		return _quaternion;
	}

	void setq(Quaternion q) override {
		_quaternion = q;
	}

	double D(unsigned short /*d*/) const override {
		return 0.0;
	}
	double M(unsigned short /*d*/) const override {
		return 0.0;
	}
	double Vi(unsigned short /*d*/) const override {
		return 0.0;
	}

	void setD(unsigned short /*d*/, double /*D*/) override {}

	inline void move(int d, double dr) override {
		setr(d, r(d) + dr);
	}

	double getI(unsigned short /*d*/) const override {
		mardyn_assert(false);
		// TODO: check values for single-centered molecules
		return 0.0;
	}

	void updateMassInertia() override {}


	double U_rot() override {
		return 0.;
	}

	double U_rot_2() override {
		return 0.;
	}

	void setupSoACache(CellDataSoABase * const s, unsigned iLJ, unsigned /*iC*/, unsigned /*iD*/, unsigned /*iQ*/) override {
		mardyn_assert(false);
		// should this ever be like called?
		setSoA(s);
		_soa_index = iLJ;
	}

	void setSoA(CellDataSoABase * s) override;

	void setStartIndexSoA_LJ(unsigned i) override {
		_soa_index = i;
	}
	void setStartIndexSoA_C(unsigned /*i*/) override {
		mardyn_assert(false);
	}
	void setStartIndexSoA_D(unsigned /*i*/) override {
		mardyn_assert(false);
	}
	void setStartIndexSoA_Q(unsigned /*i*/) override {
		mardyn_assert(false);
	}

	unsigned int numSites() const override {
		return 1;
	}
	unsigned int numOrientedSites() const override {
		return 0;
	}
	unsigned int numLJcenters() const override {
		return 1;
	}
	unsigned int numCharges() const override {
		return 0;
	}
	unsigned int numDipoles() const override {
		return 0;
	}
	unsigned int numQuadrupoles() const override {
		return 0;
	}

	std::array<double, 3> site_d(unsigned int /*i*/) const override { return emptyArray3(); }

	std::array<double, 3> ljcenter_d(unsigned int i) const override { return emptyArray3(); }
	std::array<double, 3> charge_d(unsigned int /*i*/) const override { return emptyArray3(); }
	std::array<double, 3> dipole_d(unsigned int /*i*/) const override { return emptyArray3(); }
	std::array<double, 3> quadrupole_d(unsigned int /*i*/) const override { return emptyArray3(); }

	std::array<double, 3> site_d_abs(unsigned int i) const override { return ljcenter_d_abs(i); }
	std::array<double, 3> ljcenter_d_abs(unsigned int i) const override {
		mardyn_assert(i == 0);
		return r_arr();
	}
	std::array<double, 3> charge_d_abs(unsigned int /*i*/) const override { return emptyArray3(); }
	std::array<double, 3> dipole_d_abs(unsigned int /*i*/) const override { return emptyArray3(); }
	std::array<double, 3> quadrupole_d_abs(unsigned int /*i*/) const override { return emptyArray3(); }

	std::array<double, 3> dipole_e(unsigned int /*i*/) const override { return emptyArray3(); }
	std::array<double, 3> quadrupole_e(unsigned int /*i*/) const override { return emptyArray3(); }

	std::array<double, 3> site_F(unsigned int /*i*/) const override { return emptyArray3(); }
	std::array<double, 3> ljcenter_F(unsigned int /*i*/) const override { return emptyArray3(); }
	std::array<double, 3> charge_F(unsigned int /*i*/) const override { return emptyArray3(); }
	std::array<double, 3> dipole_F(unsigned int /*i*/) const override { return emptyArray3(); }
	std::array<double, 3> quadrupole_F(unsigned int /*i*/) const override { return emptyArray3(); }

	void normalizeQuaternion() override {}
	std::array<double, 3> computeLJcenter_d(unsigned int /*i*/) const override { return emptyArray3(); }
	std::array<double, 3> computeCharge_d(unsigned int /*i*/) const override { return emptyArray3(); }
	std::array<double, 3> computeDipole_d(unsigned int /*i*/) const override { return emptyArray3(); }
	std::array<double, 3> computeQuadrupole_d(unsigned int /*i*/) const override { return emptyArray3(); }
	std::array<double, 3> computeDipole_e(unsigned int /*i*/) const override { return emptyArray3(); }
	std::array<double, 3> computeQuadrupole_e(unsigned int /*i*/) const override { return emptyArray3(); }

	unsigned long totalMemsize() const override {
		//todo: check
		mardyn_assert(false);
		return sizeof(*this);
	}
        
	void setF(double /*F*/ [3]) override {}
	void setM(double /*M*/[3]) override {}
	void setVi(double /*Vi*/[9]) override {}
	void Fadd(const double /*a*/[]) override {}
	void Madd(const double /*a*/[]) override {}
	void Viadd(const double /*a*/[]) override {}
	void vadd(const double ax, const double ay, const double az) override {
		setv(0, v(0) + ax); setv(1, v(1) + ay); setv(2, v(2) + az);
	}
	void vsub(const double ax, const double ay, const double az) override {
		setv(0, v(0) - ax); setv(1, v(1) - ay); setv(2, v(2) - az);
	}

#ifndef DEBUG_FUNCTIONALITY_HACKS
	void Fljcenteradd(unsigned int /*i*/, double /*a*/[]) override {}
	void Fljcentersub(unsigned int /*i*/, double /*a*/[]) override {}
#else
	void Fljcenteradd(unsigned int i, double a[]) override {
		mardyn_assert(i == 0);
		for(int d = 0; d < 3; ++d)
			setv(d, v(d) + a[d]);
	}
	void Fljcentersub(unsigned int i, double a[]) override {
		mardyn_assert(i == 0);
		for(int d = 0; d < 3; ++d)
			setv(d, v(d) - a[d]);
	}
#endif

	void Fchargeadd(unsigned int /*i*/, double /*a*/[]) override {}
	void Fchargesub(unsigned int /*i*/, double /*a*/[]) override {}
	void Fdipoleadd(unsigned int /*i*/, double /*a*/[]) override {}
	void Fdipolesub(unsigned int /*i*/, double /*a*/[]) override {}
	void Fquadrupoleadd(unsigned int /*i*/, double /*a*/[]) override {}
	void Fquadrupolesub(unsigned int /*i*/, double /*a*/[]) override {}
	void upd_preF(double /*dt*/) override {
		mardyn_assert(false);
	}
	void upd_postF(double /*dt_halve*/, double& /*summv2*/, double& /*sumIw2*/) override {
		mardyn_assert(false);
	}
	void calculate_mv2_Iw2(double& summv2, double& /*sumIw2*/) override {
		summv2 += _component->m() * v2();
	}
	void calculate_mv2_Iw2(double& summv2, double& /*sumIw2*/, double offx, double offy, double offz) override {
		double vcx = _v[0] - offx;
		double vcy = _v[1] - offy;
		double vcz = _v[2] - offz;
		summv2 += _component->m() * (vcx*vcx + vcy*vcy + vcz*vcz);
	}
	static std::string getWriteFormat();
	void write(std::ostream& /*ostrm*/) const override;
	void writeBinary(std::ostream& /*ostrm*/) const override {}
	void clearFM() override {}
	void calcFM() override {}
	void check(unsigned long /*id*/) override {}

	static Component * getStaticRMMComponent() {
		return _component;
	}

	void setStorageState(StorageState s) {
		_state = s;
	}

	StorageState getStorageState() const {
		return _state;
	}
        

	void buildOwnSoA() override {
		mardyn_assert(_state == STORAGE_AOS);
	}
	void releaseOwnSoA() override {
		mardyn_assert(_state == STORAGE_AOS);
	}

private:
	static std::array<double, 3> emptyArray3() {
		//mardyn_assert(false);
		std::array<double, 3> ret{0., 0., 0.};
		return ret;
	}

	static void initStaticVars();

    static Component *_component;  /**< IDentification number of its component type */
    static Quaternion _quaternion;
    static bool _initCalled;

    StorageState _state;

    // if the state is AOS, the following values are read:
    vcp_real_calc _r[3];  /**< position coordinates */
    vcp_real_accum _v[3];  /**< velocity */
    unsigned long _id;

	// if the state is SOA, the values are read from the SoA:
	CellDataSoARMM * _soa;
	size_t _soa_index;
};

std::ostream& operator<<( std::ostream& os, const MoleculeRMM& m );

#endif /* SRC_MOLECULES_MOLECULERMM_H_ */
