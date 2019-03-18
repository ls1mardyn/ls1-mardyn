#ifndef OBJECTFILLERBASE_H_
#define OBJECTFILLERBASE_H_

#include <string>

#include "molecules/Molecule.h"
#include "utils/generator/Objects.h"
#include "utils/xmlfileUnits.h"

/** FillerBase interface */
class ObjectFillerBase {
public:
	ObjectFillerBase() = default;
	virtual ~ObjectFillerBase() = default;

	/** @brief Read in XML configuration for Filler and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <filler type="STRING">
	     <!-- see derived filler classes -->
	   </filler>
	   \endcode
	 */
	virtual void readXML(XMLfileUnits& xmlconfig) = 0;

	/** Initialize the generator with current internal state */
	virtual void init() = 0;

	/* Set object to fill */
	virtual void setObject(std::shared_ptr<Object> object) = 0;

	/** Get a single molecule
	 * By subsequent calls all molecules will be returned, one by one.
	 * @param[out] molecule  Pointer to molecule data structure where to store the molecule data (coordinate and component id)
	 * @return     0 if no more molecules can be returned
	 */
	virtual int getMolecule(Molecule *molecule) = 0;

	virtual std::string getPluginName() = 0;
};

#endif  // OBJECTFILLERBASE_H_

