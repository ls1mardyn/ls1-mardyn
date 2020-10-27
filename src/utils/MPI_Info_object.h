#ifndef MPI_INFO_OBJECT_H
#define MPI_INFO_OBJECT_H

#ifdef ENABLE_MPI

#include <string>

#include "mpi.h"

#include "utils/xmlfile.h"

/** @brief MPI Info object implements functionalities to handle MPI Info objects, e.g. initialize from a XML file.
 */
class MPI_Info_object{
public:
	MPI_Info_object()
		: _mpi_info(MPI_INFO_NULL)
		{}
	~MPI_Info_object();

	/** @brief read key value pairs to be set in the MPI info object
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <pair> <key>STRING</key> <value>STRING</value> </pair>
	   <pair> ...
	   \endcode
	 */
	void readXML(XMLfile& xmlconfig);

	/** @brief reset MPI info object */
	void reset();
	/** @brief Add a hint key value pair */
	void add_hint(std::string& key, std::string& value);
	/** @brief Obtains the MPI info object as a duplicate of the internal one. */
	void get_dup_MPI_Info(MPI_Info* mpi_info);
	operator MPI_Info() const { return _mpi_info; }

private:
	MPI_Info _mpi_info;
};

#endif  // ENABLE_MPI

#endif  // MPI_INFO_OBJECT_H

