/**
 * @file StaticIrregDomainDecompostion.h
 * @author amartyads
 * @date 21.02.24
 */

#pragma once

#include "DomainDecompMPIBase.h"
#include "DomainDecomposition.h"

/** Extends DomainDecomposition to implement static irregular grids
 *
 * */
class StaticIrregDomainDecomposition : public DomainDecomposition {
public:
  /** Default constructor. Passes default values to the main constructor.
   *
   * The constructor passes the default value of _subDomainWeights to the main
   * constructor, so that the size() of individual vectors is zero, to trigger
   * the default MPI coord setup from DomainDecomposition.
   *
   * @param *domain The domain object defined in Simulation, needed to extract
   * the global simulation bounds.
   *
   */
  StaticIrregDomainDecomposition(Domain *domain);

  /** Main constructor. Passes values to DomainDecomposition, and inits
   * class members.
   *
   * In case _subdomainWeights is blank (to trigger initial coord breakdown from
   * DomainDecomposition), the generated default _gridSize is taken and used to
   * populate the weights. The behaviour becomes identical to
   * DomainDecomposition, with a regular equally-spaced grid.
   *
   * @note This constructor is directly called when MaMiCo coupling is used. In
   * this case the communicator and the weights are provided by MaMiCo.
   *
   * @param *domain The domain object defined in Simulation, needed to extract
   * the global simulation bounds.
   * @param comm The local communicator for the simulation
   * @param _subdomainWeights An array containing 3 vectors, each containing the
   * numeric weights of each subdomain, ordered by axes.
   *
   */
  StaticIrregDomainDecomposition(
      Domain *domain, MPI_Comm comm,
      const std::array<std::vector<unsigned int>, DIMgeom> &subdomainWeights);

  /** Reads in XML configuration for StaticIrregDomainDecomposition.
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
       <subdomainWeightsCSV> STRING.csv </subdomainWeightsCSV>
     </parallelisation>
     \endcode
    *
    * A file with its first line being 1,2,1 defines an x axis with subdomain
   lengths in the 1:2:1 ratio
    * If file not given, default behaviour is an equally spaced grid, same as
   DomainDecomposition
    *
    * @param &xmlconfig The xml node from which to read the CSV filename with
   weights
    *
   */
  void readXML(XMLfileUnits &xmlconfig) override;

  double getBoundingBoxMin(int dimension, Domain *domain) override;

  double getBoundingBoxMax(int dimension, Domain *domain) override;

  /** Assuming _subdomainWeights is up-to-date, calculates bounds of
   * current subdomain and updates _boxMin and _boxMax.
   *
   * For an explanation on what the weights signify, please see the
   * documentation for the member _subdomainWeights.
   *
   */
  void updateSubdomainDimensions();

  /** Reads in the CSV file given by the XML config, and updates
   * _subdomainWeights.
   *
   * The CSV file is expected to have 3 lines of comma-separated integers, with
   * the integer signifying the "weight" (relative width) of the subdomain. The
   * lines, in order, are expected to be the weights for the x, y and z
   * dimension. Consequently, the number of integers for a dimension signify the
   * number of ranks in that dimension, and is calculated as such.
   *
   * @param &filename The CSV file from which to read weights. Obtained from the
   * XML config.
   *
   */
  void updateSubdomainWeightsFromFile(const std::string &filename);

private:
  /** Stores the weights from the given CSV file.
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

  /** Stores the start of the subdomain. Calculated by
   * updateSubdomainDimensions().
   *
   */
  std::array<double, 3> _boxMin{0, 0, 0};

  /** Stores the end of the subdomain. Calculated by
   * updateSubdomainDimensions().
   *
   */
  std::array<double, 3> _boxMax{0, 0, 0};

  /** Stores the domain lengths. Set in the constructor.
   *
   */
  std::array<double, 3> _domainLength;
};