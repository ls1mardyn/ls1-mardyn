//
// Created by Kruegener on 8/19/2018.
//

#ifndef MARDYN_TRUNK_PROFILEBASE_H
#define MARDYN_TRUNK_PROFILEBASE_H

#include "../../Domain.h"
#include "../../parallel/DomainDecompBase.h"

class SpatialProfile;

struct SamplingInformation {
	double universalInvProfileUnit[3]; // Inv. Bin Sizes
	double universalProfileUnit[3]; // Bin Sizes
	double globalLength[3]; // Size of Domain
	double segmentVolume; // Size of one Sampling grid bin
	double universalCentre[3]; // Centre coords for cylinder system
	unsigned long globalNumMolecules; // number of molecules in system
	unsigned long numMolFixRegion; // number of molecules in Fix Region
	bool cylinder; // Cartesian or Cylinder output
};

/** @brief Base class for all Profile outputs used by KartesianProfile.
 *
 * The major steps for all profiles are <b>recording</b> the profile data, <b>communication</b>, writing the <b>output file</b>
 * and <b>resetting</b> everything for the next recording period. Each of these steps has a function associated with it that
 * needs to be implemented by all profiles inheriting this class to be able to work with KartesianProfile.
 * KartesianProfile also needs to know the number of items added to the communication per
 * Molecule. So a DensityProfile should return 1 via Profile::comms(). A Velocity3dProfile should return 3. <br>
 * A very simple <b>example</b> of how to use this class is the DensityProfile.
 *
 */
class ProfileBase {

public:

	virtual ~ProfileBase() {};

	/** @brief Init function is given a pointer to the KartesianProfile object handling this profile. Same for all profiles.
	 *
	 * @param kartProf Pointer to KartesianProfile. Grants access to necessary global params.
	 */
	virtual void init(SamplingInformation& samplingInformation) { _samplInfo = samplingInformation; };

	/** @brief The recording step defines what kind of data needs to be recorded for a single molecule with a corresponding uID.
	 *
	 * @param mol Reference to Molecule, needed to extract info such as velocity or Virial.
	 * @param uID uID of molecule in sampling grid, needed to put data in right spot in the profile arrays.
	 */
	virtual void record(Molecule& mol, unsigned long uID) = 0;

	/** @brief Append all necessary communication per bin to the DomainDecomposition. Append from e.g. _localProfile.
	 *
	 * @param domainDecomp DomainDecomposition handling the communication.
	 * @param uID uID of molecule in sampling grid, needed to put data in right spot in the profile arrays.
	 */
	virtual void collectAppend(DomainDecompBase* domainDecomp, unsigned long uID) = 0;

	/** @brief Get global values after AllReduceSum per bin. Write to e.g. _globalProfile.
	 *
	 * @param domainDecomp DomainDecomposition handling the communication.
	 * @param uID uID of molecule in sampling grid, needed to put data in right spot in the profile arrays.
	 */
	virtual void collectRetrieve(DomainDecompBase* domainDecomp, unsigned long uID) = 0;

	/** @brief Whatever is necessary to output for this profile.
	 *
	 * This function varies wildly between profiles. The Profile should output to its desired format here and handle all
	 * file IO for one profile writing step.
	 * @param prefix File prefix including the global _outputPrefix for all profiles and the current timestep. Should be
	 * appended by some specific file ending for this specific profile.
	 */
	virtual void output(std::string prefix, long unsigned accumulatedDatasets) = 0;

	/** @brief Used to reset all array contents for a specific uID in order to start the next recording timeframe.
	 *
	 * @param uID uID of molecule in sampling grid, needed to put data in right spot in the profile arrays.
	 */
	virtual void reset(unsigned long uID) = 0;

	/** @brief 1D profiles like a number density profile should return 1 here. 3D profiles that have 3 entries per bin
	 * that need to be communicated would need to return 3. Adjust as needed. Same number as commAppends in collectAppend.
	 *
	 * @return Number of nedded communications per bin so the communicator can be setup correctly.
	 */
	virtual int comms() = 0;

protected:
	// output file prefix
	std::string _profilePrefix;

	SamplingInformation _samplInfo;
	long _accumulatedDatasets = -1; // Number of Datasets between output writes / profile resets // -1 if not set properly

	/** @brief Write Single Data Entry for Matrix with given uID to outfile
	 *
	 * @param uID unique ID of bin
	 * @param outfile outfile to write to
	 */
	virtual void writeDataEntry(unsigned long uID, std::ofstream& outfile) const = 0;

	/**@brief Matrix writing routine to avoid code duplication
	 *
	 * @param outfile opened filestream from Profile
	 */
	void writeMatrix(std::ofstream& outfile);

	void writeKartMatrix(std::ofstream& outfile);

	/**@brief STUB for simple Matrix output without headers
	 *
	 * @param outfile opened filestream from Profile
	 */
	// TODO: implement if needed
	void writeSimpleMatrix(std::ofstream& outfile);

	/**@brief cylinder Matrix output
	 *
	 * @param outfile opened filestream from Profile
	 */
	void writeCylMatrix(std::ofstream& outfile);
};

#endif //MARDYN_TRUNK_PROFILEBASE_H
