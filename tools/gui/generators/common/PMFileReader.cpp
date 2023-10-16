/*
 * PMFileReader.cpp
 *
 * @Date: 12.07.2011
 * @Author: eckhardw
 */

#include "PMFileReader.h"
#include "Parameters/ParameterWithIntValue.h"
#include "Conversions.h"
#include <fstream>
#include <cstdlib>


const std::string nSiteTypesTag = "NSiteTypes";
const std::string nSitesTag = "NSites";

const std::string siteTypeTag = "SiteType";
const std::string siteTypeLJ126Tag = "LJ126";
const std::string siteTypeChargeTag = "Charge";
const std::string siteTypeDipoleTag = "D";
const std::string siteTypeQuadrupoleTag = "Q";

const std::string LJ126Tags[] = {"x", "y", "z", "sigma", "epsilon", "mass"};
const std::string ChargeTags[] = {"x", "y", "z", "charge", "mass", "shielding"};
const std::string DipoleTags[] = {"x", "y", "z", "theta", "phi", "dipole", "mass", "shielding"};
const std::string QuadrupoleTags[] = {"x", "y", "z", "theta", "phi", "quadrupole", "mass", "shielding"};



template <typename T>
void readTag(const std::string tagToRead, T& value, std::ifstream& input) {
	std::string line;
	do {
		getline(input, line);
		std::cout << "Read Line: " << line << std::endl;
	} while (line.size() == 0 || line[0] == '#');

	std::istringstream linestream(line);
	std::string tag;
	linestream >> tag;

	if (tag != tagToRead) {
		std::cout << "Error in PMFile: expected tag " << tagToRead <<
				" but got " << tag << std::endl;
		// todo: we'll replace this by throwing an exception in the case we
		//       build it with gui, otherwise we simply exit!
		exit(-1);
	}

	linestream >> tag;
	if (tag != "=") {
		std::cout << "Error in PMFile: expected \"=\" " << " but got " << tag << std::endl;
		if (tag.size() > 0 && tag[0] == '=') {
			std::cout << "\t YOU PROBABLY FORGOT A WHITESPACE?" << std::endl;
		}
		exit(-1);
	}

	linestream >> value;
}

void readLJ126(ComponentParameters* parameters, Generator* generator, std::ifstream& input, const int nSites) {
	for (int i = 0; i < nSites; i++) {
		double LJ126Value = 0;
		std::string baseName = parameters->getNameId() + ".LJCenter" + convertToString(i);

		for (int k = 0; k < 6; k++) {
			readTag(LJ126Tags[k], LJ126Value, input);
			ParameterWithDoubleValue* LJ126Parameter =
					dynamic_cast<ParameterWithDoubleValue*> (ParameterCollection::findParameter(parameters, baseName + "." + LJ126Tags[k]));
			LJ126Parameter->setValue(LJ126Value);
			generator->setParameter(LJ126Parameter);
		}
	}
}

void readCharge(ComponentParameters* parameters, Generator* generator, std::ifstream& input, const int nSites) {
	for (int i = 0; i < nSites; i++) {
		double value = 0;
		std::string baseName = parameters->getNameId() + ".Charge" + convertToString(i);

		for (int k = 0; k < 6; k++) {
			readTag(ChargeTags[k], value, input);

			if (ChargeTags[k] == "shielding") {
				std::cout << "PMFileReader: Ignoring value " << "\"shielding\"!" << std::endl;
				continue;
			}

			ParameterWithDoubleValue* chargeParameter =
					dynamic_cast<ParameterWithDoubleValue*> (ParameterCollection::findParameter(parameters, baseName + "." + ChargeTags[k]));
			chargeParameter->setValue(value);
			generator->setParameter(chargeParameter);
		}
	}
}


void readDipole(ComponentParameters* parameters, Generator* generator, std::ifstream& input, const int nSites) {
	for (int i = 0; i < nSites; i++) {
		double value = 0;
		std::string baseName = parameters->getNameId() + ".Dipole" + convertToString(i);

		for (int k = 0; k < 8; k++) {
			readTag(DipoleTags[k], value, input);

			if (DipoleTags[k] == "shielding" || DipoleTags[k] == "phi" || DipoleTags[k] == "theta" || DipoleTags[k] == "mass") {
				std::cout << "PMFileReader: Ignoring value " << DipoleTags[k] << std::endl;
				continue;
			}

			ParameterWithDoubleValue* dipoleParameter =
					dynamic_cast<ParameterWithDoubleValue*> (ParameterCollection::findParameter(parameters, baseName + "." + DipoleTags[k]));
			dipoleParameter->setValue(value);
			generator->setParameter(dipoleParameter);
		}

		ParameterWithDoubleValue* dipoleParameter =
				dynamic_cast<ParameterWithDoubleValue*> (ParameterCollection::findParameter(parameters, baseName + ".eMyz"));
		dipoleParameter->setValue(1.0);
		generator->setParameter(dipoleParameter);
	}

}


void readQuadrupole(ComponentParameters* parameters, Generator* generator, std::ifstream& input, const int nSites) {
	for (int i = 0; i < nSites; i++) {
		double value = 0;
		std::string baseName = parameters->getNameId() + ".Quadrupole" + convertToString(i);

		for (int k = 0; k < 8; k++) {
			readTag(QuadrupoleTags[k], value, input);

			if (QuadrupoleTags[k] == "shielding" || QuadrupoleTags[k] == "phi" || QuadrupoleTags[k] == "theta" || QuadrupoleTags[k] == "mass") {
				std::cout << "PMFileReader: Ignoring value " << QuadrupoleTags[k] << std::endl;
				continue;
			}

			ParameterWithDoubleValue* quadrupoleParameter =
					dynamic_cast<ParameterWithDoubleValue*> (ParameterCollection::findParameter(parameters, baseName + "." + QuadrupoleTags[k]));
			quadrupoleParameter->setValue(value);
			generator->setParameter(quadrupoleParameter);
		}

		ParameterWithDoubleValue* quadrupoleParameter =
				dynamic_cast<ParameterWithDoubleValue*> (ParameterCollection::findParameter(parameters, baseName + ".eQz"));
		quadrupoleParameter->setValue(1.0);
		generator->setParameter(quadrupoleParameter);
	}
}

void reset(Generator* generator, ComponentParameters* parameters) {
	std::string numberSitesNames[] = { ".NumberOfLJCenters", ".NumberOfCharges",
								  ".NumberOfDipoles", ".NumberOfQuadrupoles" };

	for (int i = 0; i < 4; i++) {
		ParameterWithIntValue* numLJParameters =
				dynamic_cast<ParameterWithIntValue*> (ParameterCollection::findParameter(parameters, parameters->getNameId() + numberSitesNames[i]));
		numLJParameters->setValue(0);
		generator->setParameter(numLJParameters);
	}
}

void PMFileReader::readPMFile(const std::string& filename, Generator* generator, ComponentParameters* parameters) {
	reset(generator, parameters);

	std::ifstream input(filename.c_str());

	int nSiteTypes = 0;
	readTag(nSiteTypesTag, nSiteTypes, input);

	std::string siteType;
	int nSites = 0;
	std::string numberSitesName;

	for (int i = 0; i < nSiteTypes; i++) {
		readTag(siteTypeTag, siteType, input);
		readTag(nSitesTag, nSites, input);

		if (siteType == siteTypeLJ126Tag) {
			numberSitesName = ".NumberOfLJCenters";
		} else if (siteType == siteTypeChargeTag) {
			numberSitesName = ".NumberOfCharges";
		} else if (siteType == siteTypeDipoleTag) {
			numberSitesName = ".NumberOfDipoles";
		} else if (siteType == siteTypeQuadrupoleTag) {
			numberSitesName = ".NumberOfQuadrupoles";
		} else {
			std::cout << "Error in PMFile: unexpected siteType " << siteType << std::endl;
			exit(-1);
		}

		// set the number of XXX centers parameter and get the adapted parameter model
		ParameterWithIntValue* numLJParameters =
				dynamic_cast<ParameterWithIntValue*> (ParameterCollection::findParameter(parameters, parameters->getNameId() + numberSitesName));
		numLJParameters->setValue(nSites);
		generator->setParameter(numLJParameters);

		std::vector<ParameterCollection*> newGeneratorParameters = generator->getParameters();
		ComponentParameters* newParameters = NULL;
		for (size_t i = 0; i < newGeneratorParameters.size(); i++) {
			newParameters =
					dynamic_cast<ComponentParameters*> (ParameterCollection::findParameter(newGeneratorParameters[i], parameters->getNameId()));
			if (newParameters != NULL) {
				break;
			}
		}

		if (siteType == siteTypeLJ126Tag) {
			readLJ126(newParameters, generator, input, nSites);
		} else if (siteType == siteTypeChargeTag) {
			readCharge(newParameters, generator, input, nSites);
		} else if (siteType == siteTypeDipoleTag) {
			readDipole(newParameters, generator, input, nSites);
		} else if (siteType == siteTypeQuadrupoleTag) {
			readQuadrupole(newParameters, generator, input, nSites);
		} else {
			std::cout << "Error in PMFile: unexpected siteType " << siteType << std::endl;
			exit(-1);
		}


		for (size_t i = 0; i < newGeneratorParameters.size(); i++) {
			ParameterCollection::deleteParameter(newGeneratorParameters[i]);
		}

	}
}
