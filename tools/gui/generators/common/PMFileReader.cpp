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

using namespace std;

const string nSiteTypesTag = "NSiteTypes";
const string nSitesTag = "NSites";

const string siteTypeTag = "SiteType";
const string siteTypeLJ126Tag = "LJ126";
const string siteTypeChargeTag = "Charge";
const string siteTypeDipoleTag = "D";
const string siteTypeQuadrupoleTag = "Q";

const string LJ126Tags[] = {"x", "y", "z", "sigma", "epsilon", "mass"};
const string ChargeTags[] = {"x", "y", "z", "charge", "mass", "shielding"};
const string DipoleTags[] = {"x", "y", "z", "theta", "phi", "dipole", "mass", "shielding"};
const string QuadrupoleTags[] = {"x", "y", "z", "theta", "phi", "quadrupole", "mass", "shielding"};



template <typename T>
void readTag(const std::string tagToRead, T& value, ifstream& input) {
	string line;
	do {
		getline(input, line);
		std::cout << "Read Line: " << line << endl;
	} while (line.size() == 0 || line[0] == '#');

	istringstream linestream(line);
	string tag;
	linestream >> tag;

	if (tag != tagToRead) {
		std::cout << "Error in PMFile: expected tag " << tagToRead <<
				" but got " << tag << endl;
		// todo: we'll replace this by throwing an exception in the case we
		//       build it with gui, otherwise we simply exit!
		exit(-1);
	}

	linestream >> tag;
	if (tag != "=") {
		std::cout << "Error in PMFile: expected \"=\" " << " but got " << tag << endl;
		if (tag.size() > 0 && tag[0] == '=') {
			std::cout << "\t YOU PROBABLY FORGOT A WHITESPACE?" << endl;
		}
		exit(-1);
	}

	linestream >> value;
}

void readLJ126(ComponentParameters* parameters, Generator* generator, ifstream& input, const int nSites) {
	for (int i = 0; i < nSites; i++) {
		double LJ126Value = 0;
		string baseName = parameters->getNameId() + ".LJCenter" + convertToString(i);

		for (int k = 0; k < 6; k++) {
			readTag(LJ126Tags[k], LJ126Value, input);
			ParameterWithDoubleValue* LJ126Parameter =
					dynamic_cast<ParameterWithDoubleValue*> (ParameterCollection::findParameter(parameters, baseName + "." + LJ126Tags[k]));
			LJ126Parameter->setValue(LJ126Value);
			generator->setParameter(LJ126Parameter);
		}
	}
}

void readCharge(ComponentParameters* parameters, Generator* generator, ifstream& input, const int nSites) {
	for (int i = 0; i < nSites; i++) {
		double value = 0;
		string baseName = parameters->getNameId() + ".Charge" + convertToString(i);

		for (int k = 0; k < 6; k++) {
			readTag(ChargeTags[k], value, input);

			if (ChargeTags[k] == "shielding") {
				std::cout << "PMFileReader: Ignoring value " << "\"shielding\"!" << endl;
				continue;
			}

			ParameterWithDoubleValue* chargeParameter =
					dynamic_cast<ParameterWithDoubleValue*> (ParameterCollection::findParameter(parameters, baseName + "." + ChargeTags[k]));
			chargeParameter->setValue(value);
			generator->setParameter(chargeParameter);
		}
	}
}


void readDipole(ComponentParameters* parameters, Generator* generator, ifstream& input, const int nSites) {
	for (int i = 0; i < nSites; i++) {
		double value = 0;
		string baseName = parameters->getNameId() + ".Dipole" + convertToString(i);

		for (int k = 0; k < 6; k++) {
			readTag(DipoleTags[k], value, input);

			if (DipoleTags[k] == "shielding" || DipoleTags[k] == "phi" || DipoleTags[k] == "theta") {
				std::cout << "PMFileReader: Ignoring value " << DipoleTags[k] << endl;
				continue;
			}

			ParameterWithDoubleValue* dipoleParameter =
					dynamic_cast<ParameterWithDoubleValue*> (ParameterCollection::findParameter(parameters, baseName + "." + DipoleTags[k]));
			dipoleParameter->setValue(value);
			generator->setParameter(dipoleParameter);
		}
	}
}


void readQuadrupole(ComponentParameters* parameters, Generator* generator, ifstream& input, const int nSites) {
	for (int i = 0; i < nSites; i++) {
		double value = 0;
		string baseName = parameters->getNameId() + ".Quadrupole" + convertToString(i);

		for (int k = 0; k < 6; k++) {
			readTag(QuadrupoleTags[k], value, input);

			if (QuadrupoleTags[k] == "shielding" || QuadrupoleTags[k] == "phi" || QuadrupoleTags[k] == "theta") {
				std::cout << "PMFileReader: Ignoring value " << QuadrupoleTags[k] << endl;
				continue;
			}

			ParameterWithDoubleValue* quadrupoleParameter =
					dynamic_cast<ParameterWithDoubleValue*> (ParameterCollection::findParameter(parameters, baseName + "." + QuadrupoleTags[k]));
			quadrupoleParameter->setValue(value);
			generator->setParameter(quadrupoleParameter);
		}
	}
}


void PMFileReader::readPMFile(const std::string& filename, Generator* generator, ComponentParameters* parameters) {
	ifstream input(filename.c_str());

	int nSiteTypes = 0;
	readTag(nSiteTypesTag, nSiteTypes, input);

	string siteType;
	int nSites = 0;
	string numberSitesName;

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
			std::cout << "Error in PMFile: unexpected siteType " << siteType << endl;
			exit(-1);
		}

		// set the number of XXX centers parameter and get the adapted parameter model
		ParameterWithIntValue* numLJParameters =
				dynamic_cast<ParameterWithIntValue*> (ParameterCollection::findParameter(parameters, parameters->getNameId() + numberSitesName));
		numLJParameters->setValue(nSites);
		generator->setParameter(numLJParameters);

		vector<ParameterCollection*> newGeneratorParameters = generator->getParameters();
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
			std::cout << "Error in PMFile: unexpected siteType " << siteType << endl;
			exit(-1);
		}


		for (size_t i = 0; i < newGeneratorParameters.size(); i++) {
			ParameterCollection::deleteParameter(newGeneratorParameters[i]);
		}

	}
}
