#include "utils/MPI_Info_object.h"

#include "Simulation.h"
#include "utils/Logger.h"

#if ENABLE_MPI

MPI_Info_object::MPI_Info_object() {
	MPI_CHECK( MPI_Info_create(&_mpi_info) );
}

MPI_Info_object::~MPI_Info_object() {
	MPI_CHECK( MPI_Info_free(&_mpi_info) );
}

void MPI_Info_object::readXML(XMLfile& xmlconfig) {
	XMLfile::Query query = xmlconfig.query("pair");
	global_log->debug() << "Number of key value pairs: " << query.card() << endl;

	string oldpath = xmlconfig.getcurrentnodepath();
	for(auto pairIter = query.begin(); pairIter; ++pairIter) {
		xmlconfig.changecurrentnode(pairIter);
		std::string key, value;
		xmlconfig.getNodeValue("key", key);
		xmlconfig.getNodeValue("value", value);
		global_log->debug() << "Found MPI Info key '" << key << "': " << value << std::endl;
		add_pair(key, value);
	}
	xmlconfig.changecurrentnode(oldpath);
}

void MPI_Info_object::add_pair(std::string& key, std::string& value) {
	if(key.size() > MPI_MAX_INFO_KEY) {
		global_log->error() << "MPI Info key name longer than allowed." << std::endl;
		return;
	}
	if(value.size() > MPI_MAX_INFO_VAL) {
		global_log->error() << "MPI Info value longer than allowed." << std::endl;
		return;
	}
	MPI_CHECK( MPI_Info_set(_mpi_info, key.c_str(), value.c_str()) );
}

void MPI_Info_object::get_MPI_Info(MPI_Info* mpi_info) {
	MPI_CHECK( MPI_Info_dup(_mpi_info, mpi_info) );
}

#endif  // ENABLE_MPI
