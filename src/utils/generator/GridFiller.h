/*
 * Copyright (c) 2013-2017 Christoph Niethammer <christoph.niethammer@gmail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#ifndef GRIDFILLER_H_
#define GRIDFILLER_H_

#include <array>
#include <random>

#include "Basis.h"
#include "Lattice.h"
#include "utils/generator/ObjectFillerBase.h"


/** The GridFiller returns molecules within an object placed on a lattice using a specified lattice basis. */
class GridFiller : public ObjectFillerBase {
public:
	   GridFiller() : _lattice(), _basis(), _origin{{0.0, 0.0, 0.0}}, _object(nullptr), _latticeOccupancy(1.0), _dis(0.0, 1.0), _gen(0) {}
	   ~GridFiller(){}

	/** @brief Read in XML configuration for GridFiller and all its included objects.
	 *
	 * If a density is provided a face-centered cubic (or orthorhombic) lattice will be used. If in this case also a lattice occupancy factor (0-1] is provided,
	 * a finer grid will be used and only the specified fraction of points will be used. The lattice vectors will be scaled to
	 * achieve the desired density taking the occupancy factor into account. By default the occupancy factor is 1 (use all lattice points).
	 * Note that the desired density is (probably) not met exactly due to the lattice structure (discrete number of molecules)
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <filler type="GridFiller">
	     <lattice><!-- ... --></lattice>
	     <basis><!-- ... --></basis>
	     <latticeOrigin> <x>DOUBLE</x> <y>DOUBLE</y> <z>DOUBLE</z> </latticeOrigin>
	     <density>DOUBLE</density>
	     <latticeOccupancy>DOUBLE</latticeOccupancy>
	   </filler>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig);

	/** Initialize the generator
		* @param[in]  lattice  The underlying point lattice to be used
		* @param[in]  basis    The molecular basis to be put on each lattice point
		* @param[in]  origin   The origin for the lattice
		* @param[in]  object   Volume object to be filled
		*/
	void init(Lattice& lattice, Basis& basis, double origin[3]);

	/** Initialize the generator with current internal state */
	void init();

	/* Set object to fill */
	void setObject(std::shared_ptr<Object> object) { _object = object; }

	/** Get a single molecule
	 * By subsequent calls all molecules will be returned, one by one.
	 * @param[out] molecule  Pointer to molecule data structure where to store the molecule data (coordinate and component id)
	 * @return     0 if no more molecules can be returned
	 */
	int getMolecule(Molecule *molecule);

	std::string getPluginName() { return std::string("GridFiller"); }
	static ObjectFillerBase* createInstance() { return new GridFiller(); }

private:
	Lattice _lattice;
	Basis _basis;
	std::array<double, 3> _origin;
	std::shared_ptr<Object> _object;
	double _latticeOccupancy;

	std::uniform_real_distribution<> _dis;
	std::mt19937 _gen;

	/* Internal values/counters used during the creation by getMolecule */
	long _baseCount;
	double _lattice_point[3];

	bool useDensity = false;
};

#endif  // GRIDFILLER_H_
