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

  if (xmlconfig.changecurrentnode("subdomainWeights")) {
    _gridSize[0] = _gridSize[1] = _gridSize[2] =
        1; // ensure that, if the xml for a direction is empty, that axis is not
           // decomposed
    for (short i = 0; i < 3; i++) {
      _subdomainWeights[i].clear();
    }

    std::string weightsx = xmlconfig.getNodeValue_string("x");
    if (!weightsx.empty()) {
      std::stringstream ss(weightsx);
      while (ss.good()) {
        int temp;
        ss >> temp;
        if (temp <= 0) {
          Log::global_log->fatal()
              << "Weights in x axis have a non-natural number! Only weights > "
                 "0 allowed, please check XML file!";
          Simulation::exit(5003);
        }
        _subdomainWeights[0].push_back(temp);
        if (ss.peek() == ',' || ss.peek() == ' ') // skip commas and spaces
          ss.ignore();
      }
    } else {
      _subdomainWeights[0].push_back(
          1); // no decomposition, whole length spanned
    }

    std::string weightsy = xmlconfig.getNodeValue_string("y");
    if (!weightsy.empty()) {
      std::stringstream ss(weightsy);
      while (ss.good()) {
        int temp;
        ss >> temp;
        if (temp <= 0) {
          Log::global_log->fatal()
              << "Weights in y axis have a non-natural number! Only weights > "
                 "0 allowed, please check XML file!";
          Simulation::exit(5003);
        }
        _subdomainWeights[1].push_back(temp);
        if (ss.peek() == ',' || ss.peek() == ' ') // skip commas and spaces
          ss.ignore();
      }
    } else {
      _subdomainWeights[1].push_back(
          1); // no decomposition, whole length spanned
    }

    std::string weightsz = xmlconfig.getNodeValue_string("z");
    if (!weightsz.empty()) {
      std::stringstream ss(weightsz);
      while (ss.good()) {
        int temp;
        ss >> temp;
        if (temp <= 0) {
          Log::global_log->fatal()
              << "Weights in z axis have a non-natural number! Only weights > "
                 "0 allowed, please check XML file!";
          Simulation::exit(5003);
        }
        _subdomainWeights[2].push_back(temp);
        if (ss.peek() == ',' || ss.peek() == ' ') // skip commas and spaces
          ss.ignore();
      }
    } else {
      _subdomainWeights[2].push_back(
          1); // no decomposition, whole length spanned
    }

    Log::global_log->info() << "Weights for subdomains for "
                               "StaticIrregDomainDecomposition have been read"
                            << std::endl;
    for (short i = 0; i < _subdomainWeights.size(); i++) {
      std::stringstream ss;
      for (auto w : _subdomainWeights[i]) {
        ss << w << " ";
      }
      Log::global_log->info()
          << "Weights for axis " << i << ": " << ss.str() << std::endl;
      _gridSize[i] = static_cast<int>(
          _subdomainWeights[i]
              .size()); //_gridSize still contains the number of ranks per
                        // dimension, just not the actual "size" of subdomains
    }
    initMPIGridDims();           // to recalculate _coords
    updateSubdomainDimensions(); // recalculate sizes from _coords
    xmlconfig.changecurrentnode("..");
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
