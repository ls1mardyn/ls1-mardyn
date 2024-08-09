#include "utils/MPI_Info_object.h"

#include "utils/Logger.h"

#ifdef ENABLE_MPI

MPI_Info_object::~MPI_Info_object() {
	reset();	// will also free _mpi_info object if necessary
}

void MPI_Info_object::reset() {
	if(_mpi_info != MPI_INFO_NULL) {
		MPI_CHECK( MPI_Info_free(&_mpi_info) );
	}
}

void MPI_Info_object::readXML(XMLfile& xmlconfig) {
	if(_mpi_info == MPI_INFO_NULL) {
		MPI_CHECK( MPI_Info_create(&_mpi_info) );
	}

	XMLfile::Query query = xmlconfig.query("hint");
	Log::global_log->debug() << "[MPI_Info_object]\tNumber of hint key value pairs: " << query.card() << std::endl;

	std::string oldpath = xmlconfig.getcurrentnodepath();
	for(auto pairIter = query.begin(); pairIter; ++pairIter) {
		xmlconfig.changecurrentnode(pairIter);
		std::string key, value;
		xmlconfig.getNodeValue("key", key);
		xmlconfig.getNodeValue("value", value);
		Log::global_log->debug() << "[MPI_Info_object]\treadXML: found MPI Info hint '" << key << "': " << value << std::endl;
		add_hint(key, value);
	}
	xmlconfig.changecurrentnode(oldpath);
}

void MPI_Info_object::add_hint(std::string& key, std::string& value) {
	if(key.size() > MPI_MAX_INFO_KEY) {
		Log::global_log->error() << "MPI Info key name longer than allowed." << std::endl;
		return;
	}
	if(value.size() > MPI_MAX_INFO_VAL) {
		Log::global_log->error() << "MPI Info value longer than allowed." << std::endl;
		return;
	}
	Log::global_log->info() << "[MPI_Info_object]\tsetting MPI Info hint " << key << "=" << value << std::endl;
	MPI_CHECK( MPI_Info_set(_mpi_info, const_cast<char*>(key.c_str()), const_cast<char*>(value.c_str())) );
}

void MPI_Info_object::get_dup_MPI_Info(MPI_Info* mpi_info) {
	MPI_CHECK( MPI_Info_dup(_mpi_info, mpi_info) );
}

#endif  // ENABLE_MPI
