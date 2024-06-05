/**
 * @file StaticIrregDomainDecompostion.cpp
 * @author amartyads
 * @date 21.02.24
 */

#include "StaticIrregDomainDecomposition.h"
#include "Domain.h"
#include "utils/Logger.h"
#include <fstream>
#include <numeric>
#include <sstream>

StaticIrregDomainDecomposition::StaticIrregDomainDecomposition(Domain *domain)
    : StaticIrregDomainDecomposition(domain, MPI_COMM_WORLD,
                                     {std::vector<unsigned int>{}, {}, {}}) {}

StaticIrregDomainDecomposition::StaticIrregDomainDecomposition(
    Domain *domain, MPI_Comm comm,
    const std::array<std::vector<unsigned int>, DIMgeom> &subdomainWeights)
    : DomainDecomposition(comm, {(int)subdomainWeights[0].size(),
                                 (int)subdomainWeights[1].size(),
                                 (int)subdomainWeights[2].size()}),
      _subdomainWeights(subdomainWeights), _domainLength{
                                               domain->getGlobalLength(0),
                                               domain->getGlobalLength(1),
                                               domain->getGlobalLength(2)} {
  if (_subdomainWeights[0].size() == 0 && _subdomainWeights[1].size() == 0 &&
      _subdomainWeights[2].size() == 0) { // default behaviour, regular grid
    /* If we have empty vectors for subdomain weights, we successfully
    initialize a regular grid by leveraging DomainDecomposition's default
    behaviour. However, we now need to use the updated _coords and _gridSize to
    write these weights back to _subdomainWeights, since otherwise the box
    dimension calculation will be broken. So for each axis, we insert a number
    of ones equal to the number of subdomains, to denote that we have equal
    subdomains on this axis, amounting to _gridSize.
    */
    for (int i = 0; i < 3; i++) {
      _subdomainWeights[i].reserve(_gridSize[i]);
      for (int j = 0; j < _gridSize[i]; j++) {
        _subdomainWeights[i].push_back(1); // equal subdomains, so weight = 1
      }
    }
  }
  updateSubdomainDimensions();
}

void StaticIrregDomainDecomposition::readXML(XMLfileUnits &xmlconfig) {
  DomainDecompMPIBase::readXML(xmlconfig);
  // bypass DomainDecomposition readXML to avoid reading MPIGridDims

  std::string filename = xmlconfig.getNodeValue_string("subdomainWeightsCSV");
  if (!filename.empty()) {
    updateSubdomainWeightsFromFile(filename);
    for (int i = 0; i < _subdomainWeights.size(); ++i) {
      _gridSize[i] = static_cast<int>(
          _subdomainWeights[i]
              .size()); //_gridSize still contains the number of ranks per
                        // dimension, just not the actual "size" of subdomains
    }
    initMPIGridDims();           // to recalculate _coords
    updateSubdomainDimensions(); // recalculate sizes from _coords
  }
}

double StaticIrregDomainDecomposition::getBoundingBoxMin(int dimension,
                                                         Domain *domain) {
  return _boxMin[dimension];
}

double StaticIrregDomainDecomposition::getBoundingBoxMax(int dimension,
                                                         Domain *domain) {
  return _boxMax[dimension];
}

void StaticIrregDomainDecomposition::updateSubdomainDimensions() {
  for (int i = 0; i < 3; i++) {
    const auto backWeight =
        std::reduce(_subdomainWeights[i].begin(),
                    _subdomainWeights[i].begin() + _coords[i], 0u);
    const auto totalWeight =
        std::reduce(_subdomainWeights[i].begin() + _coords[i],
                    _subdomainWeights[i].end(), backWeight);
    Log::global_log->debug()
        << "Dim: " << i << " totalWeight: " << totalWeight
        << " backWeight: " << backWeight << " coords: " << _coords[0] << ", "
        << _coords[1] << ", " << _coords[2] << std::endl;

    // calculate box bounds from cumulative weights of previous ranks, and the
    // weight of the current rank
    _boxMin[i] =
        static_cast<double>(backWeight) * _domainLength[i] / totalWeight;
    _boxMax[i] =
        _boxMin[i] + (static_cast<double>(_subdomainWeights[i][_coords[i]]) *
                      _domainLength[i] / totalWeight);
  }
}

void StaticIrregDomainDecomposition::updateSubdomainWeightsFromFile(
    const std::string &filename) {
  std::ifstream file(filename.c_str());
  if (!file.good()) {
    Log::global_log->fatal() << "CSV file to read domain decomposition from "
                                "does not exist! Please check config file!";
    Simulation::exit(5001);
  }

  std::string line;
  for (int i = 0; i < 3; i++) { // only reads the first 3 lines, theoretically
                                // the rest of the file can contain whatever
    getline(file, line);
    if (line.empty()) {
      Log::global_log->fatal()
          << "CSV has less than 3 lines! Please check CSV file!";
      Simulation::exit(5002);
    }
    _subdomainWeights[i].clear();
    std::stringstream ss(line);
    if (!ss.good()) {
      Log::global_log->fatal() << "EOF or I/O error occured on line " << i
                               << " of CSV. Please check CSV file!";
      Simulation::exit(5004);
    }
    while (ss.good()) {
      int temp;
      ss >> temp;
      if (temp <= 0) {
        Log::global_log->fatal() << "CSV has non-natural number! Only weights "
                                    "> 0 allowed, please check CSV file!";
        Simulation::exit(5003);
      }
      _subdomainWeights[i].push_back(temp);
      if (ss.peek() == ',' || ss.peek() == ' ') // skip commas and spaces
        ss.ignore();
    }
    if (_subdomainWeights[i].empty()) {
      Log::global_log->fatal()
          << "Weights empty, failed reading operation, please check CSV file!";
      Simulation::exit(5005);
    }
  }
  Log::global_log->info() << "Weights for subdomains for "
                             "StaticIrregDomainDecomposition have been read"
                          << std::endl;
  for (int i = 0; i < 3; i++) {
    std::stringstream ss;
    for (auto w : _subdomainWeights[i]) {
      ss << w << " ";
    }
    Log::global_log->info()
        << "Weights for axis " << i << ": " << ss.str() << std::endl;
  }
}
