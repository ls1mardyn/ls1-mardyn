/**
 * @file StaticIrregDomainDecompostion.h
 * @author amartyads
 * @date 21.02.24
 */

#pragma once

#include "DomainDecompMPIBase.h"
#include "DomainDecomposition.h"

/** @brief Extends DomainDecomposition to implement static irregular grids
 *
 * */
class StaticIrregDomainDecomposition : public DomainDecomposition {
public:
  /** @brief Default constructor. Passes default values to the main constructor.
   *
   * The _subdomainWeights passed here is empty, so that the size() of
   * individual vectors is zero, to trigger the default MPI coord setup from
   * DomainDecomposition.
   *
   */
  StaticIrregDomainDecomposition(Domain *domain);

  /** @brief Main constructor. Passes values to DomainDecomposition, and inits
   * class members.
   *
   * In case _subdomainWeights is blank (to trigger initial coord breakdown from
   * DomainDecomposition), the generated default _gridSize is taken and used to
   * populate the weights. The behaviour becomes identical to
   * DomainDecomposition, with a regular equally-spaced grid.
   *
   */
  StaticIrregDomainDecomposition(
      Domain *domain, MPI_Comm comm,
      const std::array<std::vector<unsigned int>, DIMgeom> &subdomainWeights);

  /** @brief Reads in XML configuration for StaticIrregDomainDecomposition.
   *
   * The only configuration allowed right now is a CSV file, which contains the
   actual domain breakdown.
   * Even though this class subclasses DomainDecomposition, it bypasses the
   readXML() mehod of DomainDecomposition
   * because MPIGridDims is supposed to be calculated from the CSV
   configuration.
   *
   * The following xml object structure is handled by this method:
   * \code{.xml}
     <parallelisation type="StaticIrregDomainDecomposition">
       <!-- structure handled by DomainDecompMPIBase -->
       <subdomainSizeCSV> <filename>STRING.csv</filename> </subdomainSizeCSV>
     </parallelisation>
     \endcode
   */
  void readXML(XMLfileUnits &xmlconfig) override;

  double getBoundingBoxMin(int dimension, Domain *domain) override;

  double getBoundingBoxMax(int dimension, Domain *domain) override;

  /** @brief Assuming _subdomainWeights is up-to-date, calculates bounds of
   * current subdomain and updates _boxMin and _boxMax.
   *
   * For an explanation on what the weights signify, please see the
   * documentation for the member _subdomainWeights.
   *
   */
  void updateSubdomainDimensions();

  /** @brief Reads in the CSV file given by the XML config, and updates
   * _subdomainWeights.
   *
   * The CSV file is expected to have 3 lines of comma-separated integers, with
   * the integer signifying the "weight" (relative width) of the subdomain. The
   * lines, in order, are expected to be the weights for the x, y and z
   * dimension. Consequently, the number of integers for a dimension signify the
   * number of ranks in that dimension, and is calculated as such.
   *
   */
  void updateSubdomainWeightsFromFile(const std::string &filename);

private:
  /** @brief Stores the weights from the given CSV file.
   *
   * The weights denote the relative width of that subdomain relative to the
   * others in the same dimension. Ex: if the weights in the x dimension are
   * 1,3,1, and the size of the domain in x dimension is 100, the total weight
   * is 1+3+1=5 and the divisions are 1/5, 3/5 and 1/5. Hence, all subdomains
   * with coords (0,y,z) are 20 units in x direction, (1,y,z) are 60 units and
   * (2,y,z) are 20 units. Weights are only relative for the dimension, and
   * weights for different dimensions are independent.
   *
   */
  std::array<std::vector<unsigned int>, 3> _subdomainWeights{{{}, {}, {}}};

  /** @brief Stores the start of the subdomain. Calculated by
   * updateSubdomainDimensions().
   *
   */
  std::array<double, 3> _boxMin{0, 0, 0};

  /** @brief Stores the end of the subdomain. Calculated by
   * updateSubdomainDimensions().
   *
   */
  std::array<double, 3> _boxMax{0, 0, 0};

  /** @brief Stores the domain lengths. Set in the constructor.
   *
   */
  std::array<double, 3> _domainLength;
};