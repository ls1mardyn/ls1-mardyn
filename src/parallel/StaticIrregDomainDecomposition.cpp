/**
 * @file StaticIrregDomainDecompostion.cpp
 * @author amartyads
 * @date 21.02.24
 */

#include "StaticIrregDomainDecomposition.h"
#include "Domain.h"
#include "utils/Logger.h"
#include <fstream>
#include <sstream>

StaticIrregDomainDecomposition::StaticIrregDomainDecomposition(Domain* domain) : 
	StaticIrregDomainDecomposition(domain, MPI_COMM_WORLD, {std::vector<unsigned int>{},{},{}}) {}

StaticIrregDomainDecomposition::StaticIrregDomainDecomposition(
	Domain* domain, MPI_Comm comm, const std::array<std::vector<unsigned int>, DIMgeom>& subdomainWeights) : 
		DomainDecomposition(comm, {(int)subdomainWeights[0].size(), (int)subdomainWeights[1].size(), (int)subdomainWeights[2].size()}), 
		_subdomainWeights(subdomainWeights),
		_boxMin{0},
		_boxMax{0},
		_domainLength{domain->getGlobalLength(0), domain->getGlobalLength(1), domain->getGlobalLength(2)} {
			if(_subdomainWeights[0].size() == 0 && _subdomainWeights[1].size() == 0 && _subdomainWeights[2].size() == 0) { //default behaviour, regular grid
				for (int i = 0; i < 3; i++) {
					for (int j = 0; j < _gridSize[i]; j++) {
						_subdomainWeights[i].push_back(1);  //equal subdomains, so weight = 1
					}
				}
			}
			updateSubdomainDimensions();
		}

void StaticIrregDomainDecomposition::readXML(XMLfileUnits& xmlconfig) {
	DomainDecompMPIBase::readXML(xmlconfig);
	//bypass DomainDecomposition readXML to avoid reading MPIGridDims

	if(xmlconfig.changecurrentnode("subdomainSizeCSV")) {
		std::string filename = xmlconfig.getNodeValue_string("filename");
		updateSubdomainWeightsFromFile(filename);
		for (int i = 0; I < _subdomainWeights.size(); ++i) {
			_gridSize[i] = static_cast<int>(_subdomainWeights[0].size()); //_gridSize still contains the number of ranks per dimension, just not the actual "size" of subdomains
		}
		xmlconfig.changecurrentnode("..");
		initMPIGridDims();							//to recalculate _coords
		updateSubdomainDimensions();				//recalculate sizes from _coords
	}
}

double StaticIrregDomainDecomposition::getBoundingBoxMin(int dimension, Domain* domain) {
	return _boxMin[dimension];
}

double StaticIrregDomainDecomposition::getBoundingBoxMax(int dimension, Domain* domain) {
	return _boxMax[dimension];
}

void StaticIrregDomainDecomposition::updateSubdomainDimensions() {
	std::array<int, 3> totWeight{0}; //stores the total weight across the dimensions
	for(int i = 0; i < 3; i++) {
		for(unsigned int weight: _subdomainWeights[i]) {
			totWeight[i] += weight;
		}
	}
	Log::global_log->debug() << "totweight: " << totWeight[0] << " " << totWeight[1] << " " << totWeight[2]  << std::endl;

	int backWeight; //stores the cumulative weights of all subdomains previous to current subdomain
	for(int i = 0; i < 3; i++) {
		backWeight = 0;
		for(int j = 0; j < _coords[i]; j++) {
			backWeight += _subdomainWeights[i][j];
		}
		Log::global_log->debug() << "backweight: " << i << " " << backWeight << " coords: " << _coords[0] << ", " <<_coords[1] << ", " << _coords[2] << std::endl;

		//calculate box bounds from cumulative weights of previous ranks, and the weight of the current rank
		_boxMin[i] = static_cast<double>(backWeight) * _domainLength[i] / totWeight[i];
		_boxMax[i] = _boxMin[i] + (static_cast<double>(_subdomainWeights[i][_coords[i]]) * _domainLength[i] / totWeight[i]);
	}
}

void StaticIrregDomainDecomposition::updateSubdomainWeightsFromFile(std::string filename) {
	if(filename.empty()) {
		Log::global_log->fatal() << "CSV filename to read domain decomposition from is empty! Please check config file!";
		Simulation::exit(5000);
	}
	std::ifstream file(filename.c_str());
	if(!file.good()) {
		Log::global_log->fatal() << "CSV file to read domain decomposition from does not exist! Please check config file!";
		Simulation::exit(5001);
	}

	std::string line;
	for(int i = 0; i < 3; i++) { //only reads the first 3 lines, theoretically the rest of the file can contain whatever
		getline(file, line);
		if(line.empty()) {
			Log::global_log->fatal() << "CSV has less than 3 lines! Please check CSV file!";
			Simulation::exit(5002);
		}
		_subdomainWeights[i].clear();
		std::stringstream ss(line);
		while(ss.good()) {
			int temp;
			ss >> temp;
			if (temp <=0) {
				Log::global_log->fatal() << "CSV has non-natural number! Only weights > 0 allowed, please check CSV file!";
				Simulation::exit(5003);
			}
			_subdomainWeights[i].push_back(temp);
			if (ss.peek() == ',' || ss.peek() == ' ') //skip commas and spaces
				ss.ignore();
		}
	}
	Log::global_log->info() << "Weights for subdomains for StaticIrregDomainDecomposition have been read" << std::endl;
	for (int i = 0; i < 3; i++) {
		std::stringstream ss;
		for (int w: _subdomainWeights[i]) {
			ss << w << " ";
		}
		Log::global_log->info() << "Weights for axis " << i << ": " << ss.str() << std::endl;
	}
}