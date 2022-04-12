/**
 * @file AutoPasFullMolecule.h
 * @author seckler
 * @date 20.09.18
 */

#pragma once

#include <autopas/molecularDynamics/MoleculeLJ.h>
#include <autopas/utils/SoAType.h>
#include <autopas/utils/inBox.h>

#include "FullMolecule.h"

/**
 * class that implements additional functions to make the molecule compatible with autopas
 */
class AutoPasSimpleMolecule final : public MoleculeInterface, public autopas::ParticleFP64 {
public:
	explicit AutoPasSimpleMolecule(unsigned long id = 0, Component* component = nullptr, double rx = 0., double ry = 0.,
								   double rz = 0., double vx = 0., double vy = 0., double vz = 0., double q0 = 1.,
								   double q1 = 1., double q2 = 0., double q3 = 0., double Dx = 0., double Dy = 0.,
								   double Dz = 0.);

	AutoPasSimpleMolecule(const AutoPasSimpleMolecule& m) = default;

	~AutoPasSimpleMolecule() override = default;

    /**
     * Enums used as ids for accessing and creating a dynamically sized SoA.
     */
    enum AttributeNames : int { ptr, id, posX, posY, posZ, forceX, forceY, forceZ, typeId, ownershipState };

    /**
     * The type for the SoA storage.
     *
     * @note The attribute owned is of type float but treated as a bool.
     * This means it shall always only take values 0.0 (=false) or 1.0 (=true).
     * The reason for this is the easier use of the value in calculations (See LJFunctor "energyFactor")
     */
    using SoAArraysType =
        typename autopas::utils::SoAType<AutoPasSimpleMolecule *, size_t /*id*/, double /*x*/, double /*y*/,
        double /*z*/, double /*fx*/, double /*fy*/, double /*fz*/,
        size_t /*typeid*/, autopas::OwnershipState /*ownershipState*/>::Type;

    /**
     * Getter, which allows access to an attribute using the corresponding attribute name (defined in AttributeNames).
     * @tparam attribute Attribute name.
     * @return Value of the requested attribute.
     * @note The value of owned is return as floating point number (true = 1.0, false = 0.0).
     */
    template <AttributeNames attribute>
    constexpr typename std::tuple_element<static_cast<size_t>(attribute), SoAArraysType>::type::value_type get() {
      if constexpr (attribute == AttributeNames::ptr) {
        return this;
      } else if constexpr (attribute == AttributeNames::id) {
        return getID();
      } else if constexpr (attribute == AttributeNames::posX) {
        return getR()[0];
      } else if constexpr (attribute == AttributeNames::posY) {
        return getR()[1];
      } else if constexpr (attribute == AttributeNames::posZ) {
        return getR()[2];
      } else if constexpr (attribute == AttributeNames::forceX) {
        return getF()[0];
      } else if constexpr (attribute == AttributeNames::forceY) {
        return getF()[1];
      } else if constexpr (attribute == AttributeNames::forceZ) {
        return getF()[2];
      } else if constexpr (attribute == AttributeNames::typeId) {
        return getTypeId();
      } else if constexpr (attribute == AttributeNames::ownershipState) {
        return this->_ownershipState;
      } else {
        throw std::runtime_error("AutoPasSimpleMolecule::get() unknown attribute");
      }
    }

    /**
     * Setter, which allows set an attribute using the corresponding attribute name (defined in AttributeNames).
     * @tparam attribute Attribute name.
     * @param value New value of the requested attribute.
     * @note The value of owned is extracted from a floating point number (true = 1.0, false = 0.0).
     */
    template <AttributeNames attribute>
    constexpr void set(
        typename std::tuple_element<static_cast<size_t>(attribute), SoAArraysType>::type::value_type value) {
      if constexpr (attribute == AttributeNames::id) {
        setID(value);
      } else if constexpr (attribute == AttributeNames::posX) {
        _r[0] = value;
      } else if constexpr (attribute == AttributeNames::posY) {
        _r[1] = value;
      } else if constexpr (attribute == AttributeNames::posZ) {
        _r[2] = value;
      } else if constexpr (attribute == AttributeNames::forceX) {
        _f[0] = value;
      } else if constexpr (attribute == AttributeNames::forceY) {
        _f[1] = value;
      } else if constexpr (attribute == AttributeNames::forceZ) {
        _f[2] = value;
      } else if constexpr (attribute == AttributeNames::typeId) {
        setTypeId(value);
      } else if constexpr (attribute == AttributeNames::ownershipState) {
        this->_ownershipState = value;
      } else {
        throw std::runtime_error("AutoPasSimpleMolecule::set() unknown attribute");
      }
    }

	unsigned long getID() const override { return _id; }

	void setid(unsigned long id) override { setID(id); }

	void setComponent(Component* component) override { _component = component; }

	void setr(unsigned short d, double r) override { _r[d] = r; }

	void setv(unsigned short d, double v) override { _v[d] = v; }

	void setF(unsigned short d, double F) override { _f[d] = F; }

	Component* component() const override { return _component; }

	double r(unsigned short d) const override { return getR()[d]; }

	double v(unsigned short d) const override { return getV()[d]; }

	double F(unsigned short d) const override { return getF()[d]; }

	const Quaternion& q() const override { return _quaternion; }

	void setq(Quaternion q) override { _quaternion = q; }

	double D(unsigned short d) const override { return 0.; }

	double M(unsigned short d) const override { return 0.; }

	double Vi(unsigned short d) const override { return 0.; }

	void setD(unsigned short d, double D) override {}

	inline void move(int d, double dr) override { setr(d, r(d) + dr); }

	// by Stefan Becker
	double getI(unsigned short d) const override { return 1.; }

	double U_rot() override { return 0.; }

	double U_rot_2() override { return 0.; }

	double U_pot() override { return 0.; }

	void updateMassInertia() override {}

	void setupSoACache(CellDataSoABase* const s, unsigned iLJ, unsigned iC, unsigned iD, unsigned iQ) override {}

	void setSoA(CellDataSoABase* const s) override{};

	void setStartIndexSoA_LJ(unsigned i) override{};

	void setStartIndexSoA_C(unsigned i) override{};

	void setStartIndexSoA_D(unsigned i) override{};

	void setStartIndexSoA_Q(unsigned i) override{};

	unsigned int numSites() const override { return 1; };

	unsigned int numOrientedSites() const override { return 0; }

	unsigned int numLJcenters() const override { return 1; }

	unsigned int numCharges() const override { return 0; }

	unsigned int numDipoles() const override { return 0; }

	unsigned int numQuadrupoles() const override { return 0; }

	std::array<double, 3> site_d(unsigned int i) const override { return emptyArray3(); }

	std::array<double, 3> ljcenter_d(unsigned int i) const override { return emptyArray3(); }

	std::array<double, 3> charge_d(unsigned int i) const override { return emptyArray3(); }

	std::array<double, 3> dipole_d(unsigned int i) const override { return emptyArray3(); }

	std::array<double, 3> quadrupole_d(unsigned int i) const override { return emptyArray3(); }

	std::array<double, 3> site_d_abs(unsigned int i) const override { return getR(); }

	std::array<double, 3> ljcenter_d_abs(unsigned int i) const override { return getR(); }

	std::array<double, 3> charge_d_abs(unsigned int i) const override { return emptyArray3(); }

	std::array<double, 3> dipole_d_abs(unsigned int i) const override { return emptyArray3(); }

	std::array<double, 3> quadrupole_d_abs(unsigned int i) const override { return emptyArray3(); }

	std::array<double, 3> dipole_e(unsigned int i) const override { return emptyArray3(); }

	std::array<double, 3> quadrupole_e(unsigned int i) const override { return emptyArray3(); }

	std::array<double, 3> site_F(unsigned int i) const override { return getF(); }

	std::array<double, 3> ljcenter_F(unsigned int i) const override { return emptyArray3(); }

	std::array<double, 3> charge_F(unsigned int i) const override { return emptyArray3(); }

	std::array<double, 3> dipole_F(unsigned int i) const override { return emptyArray3(); }

	std::array<double, 3> quadrupole_F(unsigned int i) const override { return emptyArray3(); }

	void normalizeQuaternion() override {}

	std::array<double, 3> computeLJcenter_d(unsigned int i) const override { return emptyArray3(); }

	std::array<double, 3> computeCharge_d(unsigned int i) const override { return emptyArray3(); }

	std::array<double, 3> computeDipole_d(unsigned int i) const override { return emptyArray3(); }

	std::array<double, 3> computeQuadrupole_d(unsigned int i) const override { return emptyArray3(); }

	std::array<double, 3> computeDipole_e(unsigned int i) const override { return emptyArray3(); }

	std::array<double, 3> computeQuadrupole_e(unsigned int i) const override { return emptyArray3(); }

	unsigned long totalMemsize() const override { return sizeof(*this); }

	void setF(double F[3]) override {
		for (unsigned short i = 0; i < 3; i++) _f[i] = F[i];
	}

	void setF(const std::array<double, 3>& f) { autopas::Particle::setF(f); }

	void setM(double M[3]) override {}

	void setVi(double Vi[9]) override {}

	void setU(const double upot) override {}

	void Fadd(const double F[]) override {
		for (unsigned short i = 0; i < 3; i++) _f[i] += F[i];
	}

	void Madd(const double a[]) override {}

	void Viadd(const double a[]) override {}

	void vadd(const double ax, const double ay, const double az) override {
		std::array<double, 3> addV_arr{ax, ay, az};
		addV(addV_arr);
	}

	void vsub(const double ax, const double ay, const double az) override {
		std::array<double, 3> addV_arr{-ax, -ay, -az};
		addV(addV_arr);
	}

	void Uadd(const double upot) override {}

	void Fljcenteradd(unsigned int i, double a[]) override { vadd(a[0], a[1], a[2]); }

	void Fljcentersub(unsigned int i, double a[]) override { vsub(a[0], a[1], a[2]); }

	void Fchargeadd(unsigned int i, double a[]) override {}

	void Fchargesub(unsigned int i, double a[]) override {}

	void Fdipoleadd(unsigned int i, double a[]) override {}

	void Fdipolesub(unsigned int i, double a[]) override {}

	void Fquadrupoleadd(unsigned int i, double a[]) override {}

	void Fquadrupolesub(unsigned int i, double a[]) override {}

	// Leapfrog integration:
	void upd_preF(double dt) override;

	void upd_postF(double dt_halve, double& summv2, double& sumIw2) override;

	void calculate_mv2_Iw2(double& summv2, double& sumIw2) override { summv2 += _component->m() * v2(); }

	void calculate_mv2_Iw2(double& summv2, double& sumIw2, double offx, double offy, double offz) override {
		double vcx = _v[0] - offx;
		double vcy = _v[1] - offy;
		double vcz = _v[2] - offz;
		summv2 += _component->m() * (vcx * vcx + vcy * vcy + vcz * vcz);
	}

	static std::string getWriteFormat() { return "IRV"; }

	void write(std::ostream& ostrm) const override {
		ostrm << getID() << "\t" << r(0) << " " << r(1) << " " << r(2) << "\t" << v(0) << " " << v(1) << " " << v(2)
			  << "\t" << std::endl;
	}

	void writeBinary(std::ostream& ostrm) const override {
		ostrm.write(reinterpret_cast<const char*>(&_id), 8);
		ostrm.write(reinterpret_cast<const char*>(&(_r[0])), 8);
		ostrm.write(reinterpret_cast<const char*>(&(_r[1])), 8);
		ostrm.write(reinterpret_cast<const char*>(&(_r[2])), 8);
		ostrm.write(reinterpret_cast<const char*>(&(_v[0])), 8);
		ostrm.write(reinterpret_cast<const char*>(&(_v[1])), 8);
		ostrm.write(reinterpret_cast<const char*>(&(_v[2])), 8);
	}

	void clearFM() override {
		for (unsigned short d = 0; d < 3; ++d) {
			setF(d, 0.);
		}
	}
	void calcFM() override {}
	void check(unsigned long id) override {}

	void buildOwnSoA() override {}

	void releaseOwnSoA() override {}

	bool inBox(const double rmin[3], const double rmax[3]) const override {
		return MoleculeInterface::inBox(rmin, rmax);
	}

	bool inBox(const std::array<double, 3>& rmin, const std::array<double, 3>& rmax) const {
		return autopas::utils::inBox(this->getR(), rmin, rmax);
	}

	size_t getTypeId() const {
		mardyn_assert(_component != nullptr);
		return _component->ID();
	}

private:
	static std::array<double, 3> emptyArray3() {
		// mardyn_assert(false);
		std::array<double, 3> ret{0., 0., 0.};
		return ret;
	}

	static Component* _component;
	static Quaternion _quaternion;
};

std::ostream& operator<<( std::ostream& os, const AutoPasSimpleMolecule& m );