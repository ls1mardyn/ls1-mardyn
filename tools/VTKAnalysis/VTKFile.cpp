//
// Created by alex on 22.05.23.
//

#include <unordered_map>
#include <sstream>
#include "utils/xmlfileUnits.h"
#include "utils/xmlfile.h"
#include "VTKFile.h"

void VTKFile::readXML(const std::string &file) {
    try{
        XMLfileUnits inp(file);
        if(inp.changecurrentnode("/VTKFile") <= 0 ||
           inp.changecurrentnode("UnstructuredGrid") <= 0 ||
           inp.changecurrentnode("Piece") <= 0) {
            throw std::runtime_error{"VTK file has wrong format!"};
        }

        long numPoints = 0;
        inp.getNodeValue("@NumberOfPoints", numPoints);
        if(numPoints == 0) throw std::runtime_error{"Contains no molecules!"};

        //need to map each index to the correct molecule, as one molecule has multiple points and these point are not necessarily sorted
        std::unordered_map<uint64_t, uint64_t> index_to_mol;
        //here we count for each molecule how many particles it contains to later compute the center position
        std::unordered_map<uint64_t, uint64_t> mol_to_nParticles;
        // handle point data
        {
            XMLfile::Query query = inp.query("PointData/DataArray");
            long numArrays = query.card();
            if(numArrays < 2) throw std::runtime_error{"VTK file has wrong format!"};
            XMLfile::Query::const_iterator iterator;
            std::string oldPath = inp.getcurrentnodepath();
            for(iterator = query.begin(); iterator; iterator++) {
                inp.changecurrentnode(iterator);
                // handle IDs
                if(inp.getNodeValue_string("@Name") == "id")
                {
                    std::string data = inp.getNodeValue_string("./");
                    std::stringstream stream {data};
                    uint64_t index = 0;
                    while(!stream.eof()) {
                        uint64_t value;
                        stream >> value;
                        index_to_mol[index++] = value;
                        if(value >= _molecules.size()) {
                            _molecules.resize(value+1);
                            _molecules[value].id = value;
                            _molecules[value].r = {0,0,0};
                        }
                        mol_to_nParticles[value] += 1;
                    }
                }
                // handle CIDs
                if(inp.getNodeValue_string("@Name") == "component-id")
                {
                    std::string data = inp.getNodeValue_string("./");
                    std::stringstream stream {data};
                    uint64_t index = 0;
                    while(!stream.eof()) {
                        int value;
                        stream >> value;
                        _molecules[index_to_mol[index++]].cid = value;
                    }
                }

            }
            inp.changecurrentnode(oldPath);
        }

        // handle points
        {
            std::string data = inp.getNodeValue_string("Points/DataArray");
            std::stringstream stream {data};
            uint64_t index = 0;
            while(!stream.eof()) {
                float x,y,z;
                stream >> x;
                stream >> y;
                stream >> z;
                _molecules[index_to_mol[index]].r[0] += x;
                _molecules[index_to_mol[index]].r[1] += y;
                _molecules[index_to_mol[index++]].r[2] += z;
            }

            for(auto& mol : _molecules) {
                mol.r[0] /= static_cast<float>(mol_to_nParticles[mol.id]);
                mol.r[1] /= static_cast<float>(mol_to_nParticles[mol.id]);
                mol.r[2] /= static_cast<float>(mol_to_nParticles[mol.id]);
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Error in XML config. Please check your input file!" << std::endl;
        std::cerr << "Exception: " << e.what() << std::endl;
        exit(7);
    }
}

const std::vector<Molecule> &VTKFile::getMolecules() {
    return _molecules;
}
