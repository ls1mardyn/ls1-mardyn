/*
 * MS2RestartReader.h
 *
 *  Created on: 26.07.2013
 *      Author: eckhardw
 */

#ifndef MS2RESTARTREADER_H_
#define MS2RESTARTREADER_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

class MS2RestartReader {

public:

struct MoleculeData {
    int id;
    int cid;
    double x[3];
    double v[3];
    double q[4];
    double d[3];

    MoleculeData() {
        id = -1;
        cid = -1;
        x[0] = x[1] = x[2] = 0.0;
        v[0] = v[1] = v[2] = 0.0;
        q[0] = 1.0; q[1] = q[2] = q[3] = 0.0;
        d[0] = d[1] = d[2] = 0.0;
     }

    void print(std::ostream& str) {
	str << " Mol: " << id << " " << cid << " "
            << x[0] << " " << x[1] << " " << x[2] << "\t" << v[0] << " " << v[1] << " " << v[2] << "\t"
            << q[0] << " " << q[1] << " " << q[2] << " " <<q[3] << "\t" << d[0] << " " << d[1] << " " << d[2] << std::endl;
    }

};

static MoleculeData* readMS2RestartFile(const std::string& file, const int numberOfComponents,
		const int numberOfMolecules, std::vector<bool> hasRotationalDOF);

private:

static int readMolecules(std::ifstream& inputfilestream, int cid, int startid, MoleculeData* molecules, bool rotationalDOF);

};

#endif /* MS2RESTARTREADER_H_ */
