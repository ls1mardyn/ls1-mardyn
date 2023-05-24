
#include "MS2RestartReader.h"

#include <cmath>
#include <sstream>
#include <string>
#include <iomanip>


int MS2RestartReader::readMolecules(ifstream& inputfilestream, int cid, int startid, MoleculeData* molecules, bool rotationalDOF) {
	char dummy;
	std::string line;
	int numMolecules = 0;

	getline(inputfilestream, line);
	std::stringstream lineStream(line);
	lineStream >> numMolecules;
	std::cout << "Reading numIons= " << line << std::endl;
	std::cout << "numMolecules= " << numMolecules << std::endl;

	int id = startid;
    for (int i = 0; i < numMolecules; i++) {
        molecules[id].id = id;
        molecules[id].cid = cid;
        getline(inputfilestream, line);
        //cout << "read line=" << line << std::endl;
        std::stringstream lineStream(line);
        lineStream >> molecules[id].x[0] >> dummy >> molecules[id].x[1] >> dummy >> molecules[id].x[2];
        id++;
    }

    id = startid;
	for (int i = 0; i < numMolecules; i++) {
	    getline(inputfilestream, line);
	    std::stringstream lineStream(line);
		lineStream >> molecules[id].v[0] >> dummy >> molecules[id].v[1] >> dummy >> molecules[id].v[2];
		id++;
	}

	/***  ATOMS OF SECOND COMPONENT  ***/
	for (int i = 0; i < 4*numMolecules; i++) {
	    getline(inputfilestream, line);
	}

	if (rotationalDOF) {
		id = startid;
		for (int i = 0; i < numMolecules; i++) {
		    getline(inputfilestream, line);
		    std::stringstream lineStream(line);
		    lineStream >> molecules[id].q[0] >> dummy >> molecules[id].q[1] >> dummy >> molecules[id].q[2] >> dummy >> molecules[id].q[3];
		    id++;
		}

		for (int i = 0; i < 4*numMolecules; i++) {
			getline(inputfilestream, line);
		}

		id = startid;
		for (int i = 0; i < numMolecules; i++) {
		    getline(inputfilestream, line);
		    std::stringstream lineStream(line);
			lineStream >> molecules[id].d[0] >> dummy >> molecules[id].d[1] >> dummy >> molecules[id].d[2];
			id++;
		}
	}
	return numMolecules;
}

MS2RestartReader::MoleculeData* MS2RestartReader::readMS2RestartFile(const std::string& file, const int numberOfComponents,
		const int numberOfMolecules, std::vector<bool> hasRotationalDOF) {

    std::ifstream inputfilestream(file.c_str());
    if (!inputfilestream.is_open()) {
        std::cout << "ERROR: Could not open file " << file << std::endl;
        return NULL;
    }

    MoleculeData* molecules = new MoleculeData[numberOfMolecules];
    std::cout << std::setprecision(9);
    std::string line;
    for (int i = 0; i < 4; i++) {
        getline(inputfilestream, line);
        std::cout << " skipped line: " << line << std::endl;
    }

    double boxVolume = 0;
    getline(inputfilestream, line);
	std::stringstream lineStream(line);
	lineStream >> boxVolume;
	std::cout << "BoxVolume = " << boxVolume << std::endl;

    for (int i = 0; i < 6; i++) {
        getline(inputfilestream, line);
        std::cout << " skipped line: " << line << std::endl;
    }

    int start_idx = 0;
    for (int i = 0; i < numberOfComponents; i++) {
    	start_idx = readMolecules(inputfilestream, i, start_idx, molecules, hasRotationalDOF[i]);
    }
    
//    for (int i = 0; i < numberOfMolecules; i++) {
//        molecules[i].print(cout);
//    }

    // shift and scale molecule coordinates (ms2 centers box around origin)
    double boxLength = pow(boxVolume, 1.0 / 3.0);
    for (int i = 0; i < numberOfMolecules; i++) {
    	for (int j = 0; j < 3; j++) {
    		molecules[i].x[j] = ( molecules[i].x[j] + 0.5 ) * boxLength;
    	}
    }
    
    return molecules;
}

