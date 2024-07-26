#ifndef ASCIIREADER_H_
#define ASCIIREADER_H_

#include "io/InputBase.h"
#include <string>
#include <fstream>

/** @brief The ASCIIReader reads in phasespace information using ls1-MarDyn's (old) ASCII input file format.
 */
class ASCIIReader : public InputBase {
public:
	ASCIIReader();

	~ASCIIReader() {}

	void readXML(XMLfileUnits& xmlconfig);

	/** @brief Set the phase space file to be read.
	 *
	 * @param filename full path to the phase space file
	 */
	void setPhaseSpaceFile(std::string filename);

	//! @brief For this class, header and data are in the same file, so there is no separate header file
	void setPhaseSpaceHeaderFile(std::string filename);

	//! DEPRECATED!
	//! The information stored in the header must be provided via the xml config file
	//! @brief reads in header of the input file (including component description)
	//!
	//! The Header in the input file consists of several elements. An element starts
	//! in a new line with the element name followed by some whitespace (e.g. "Temperature ").
	//! After that, there can be any number of tokens belonging to that element.
	//! A description of the different elements follows below. But first
	//! some rules dealing with the order of the elements:
	//! \li The first header element must start with the string "MDProject"
	//!     followed by a timestamp.
	//! \li The last header element must start with the string NumberOfMolecules.
	//!     After that, the header is over and the molecule data follows
	//! \li The order of the remaining header lines is not important.
	//! \li Header elements beginning with "#" are ignored.
	//!
	//! Now the description of the different elements:
	//! \li MDProject: One token follows with a version number (see description of _inpversion);
	//! \li currenttime: One token follows with the start time
	//! \li Temperature: One token follows with the temperature
	//! \li Length: Three tokens follow with the length of the simulation box
	//!     in the three dimensions (x, y and z)
	//! \li NumberOfComponents: Here follow several tokens, the first one is the
	//!     actual Number of Components.
	//!     Then the values describing the components have to follow, seperated by
	//!     whitespace. For each component, the following values have to be provided:
	//!     - Number of Lennard-Jones-Centers, Number of Dipoles, Number of Quadrupoles (all int)
	//!     - For each LJ-Center: x-coord., y-coord., z-coord., mass, epsilon, sigma, tcutoff, do_shift (all double)
	//!     - For each Dipole: x-coord., y-coord., z-coord., eMyx, eMyy, eMyz, absMy (all double)
	//!     - For each Quadrupole: x-coord., y-coord., z-coord., eQx, eQy, eQz, absQ (all double)
	//!     - moments of inertia for principal axes: I11, I22, I33 (all double)
	//!     - For each pair of different components: xi, eta (both double)
	//!     - epsilonRF (double)
	//! \li NumberOfMolecules: One token follows with the number of molecules
	void readPhaseSpaceHeader(Domain* domain, double timestep) {};

	//! @brief reads in the data of all molecules
	//!
	//! The Molecule Data starts in a new line with the string "MoleculeFormat"
	//! followed by whitespace and a string representing the format.
	//! In the standard case (format ICRVQD), the following values are provided for each molecule:
	//! \li id of the molecule (int)
	//! \li id of the component of the molecule (int)
	//! \li Coordinates: x, y, z (all double)
	//! \li velocities: vx, vy, vz (all double)
	//! \li Orientation (quaternion): q0, q1, q2, q3 (all double)
	//! \li Angular Momentum: Dx, Dy, Dz (all double)
	//!
	//! @param particleContainer Here the Molecules from the input file are stored
	//! @return Number of molecules read in from the input phase space file
	unsigned long readPhaseSpace(ParticleContainer* particleContainer, Domain* domain, DomainDecompBase* domainDecomp);

private:

	std::string _phaseSpaceFile;
	std::string _phaseSpaceHeaderFile;
	std::fstream _phaseSpaceFileStream;
	std::fstream _phaseSpaceHeaderFileStream;

};

#endif  // ASCIIREADER_H_
