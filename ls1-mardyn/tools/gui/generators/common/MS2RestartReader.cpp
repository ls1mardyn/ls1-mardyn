
#include "MS2RestartReader.h"

#include <cmath>
#include <sstream>
#include <string>
#include <iomanip>

using namespace std;


int MS2RestartReader::readMolecules(ifstream& inputfilestream, int cid, int startid, MoleculeData* molecules, bool rotationalDOF) {
	char dummy;
	string line;
	int numMolecules = 0;

	getline(inputfilestream, line);
	stringstream lineStream(line);
	lineStream >> numMolecules;
	cout << "Reading numIons= " << line << endl;
	cout << "numMolecules= " << numMolecules << endl;

	int id = startid;
    for (int i = 0; i < numMolecules; i++) {
        molecules[id].id = id;
        molecules[id].cid = cid;
        getline(inputfilestream, line);
        //cout << "read line=" << line << endl;
        stringstream lineStream(line);
        lineStream >> molecules[id].x[0] >> dummy >> molecules[id].x[1] >> dummy >> molecules[id].x[2];
        id++;
    }

    id = startid;
	for (int i = 0; i < numMolecules; i++) {
	    getline(inputfilestream, line);
	    stringstream lineStream(line);
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
		    stringstream lineStream(line);
		    lineStream >> molecules[id].q[0] >> dummy >> molecules[id].q[1] >> dummy >> molecules[id].q[2] >> dummy >> molecules[id].q[3];
		    id++;
		}

		for (int i = 0; i < 4*numMolecules; i++) {
			getline(inputfilestream, line);
		}

		id = startid;
		for (int i = 0; i < numMolecules; i++) {
		    getline(inputfilestream, line);
		    stringstream lineStream(line);
			lineStream >> molecules[id].d[0] >> dummy >> molecules[id].d[1] >> dummy >> molecules[id].d[2];
			id++;
		}
	}
	return numMolecules;
}

MS2RestartReader::MoleculeData* MS2RestartReader::readMS2RestartFile(const string& file, const int numberOfComponents,
		const int numberOfMolecules, std::vector<bool> hasRotationalDOF) {

    ifstream inputfilestream(file.c_str());
    if (!inputfilestream.is_open()) {
        cout << "ERROR: Could not open file " << file << endl;
        return NULL;
    }

    MoleculeData* molecules = new MoleculeData[numberOfMolecules];
    cout << setprecision(9);
    string line;
    for (int i = 0; i < 4; i++) {
        getline(inputfilestream, line);
        std::cout << " skipped line: " << line << endl;
    }

    double boxVolume = 0;
    getline(inputfilestream, line);
	stringstream lineStream(line);
	lineStream >> boxVolume;
	std::cout << "BoxVolume = " << boxVolume << endl;

    for (int i = 0; i < 6; i++) {
        getline(inputfilestream, line);
        std::cout << " skipped line: " << line << endl;
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

