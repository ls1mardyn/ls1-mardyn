//
// Created by alex on 03.06.23.
//

#ifndef MARDYN_INTERACTIONLOGGER_H
#define MARDYN_INTERACTIONLOGGER_H


#include <string>
#include <utility>
#include <vector>
#include "molecules/Molecule.h"

/**
 * Logger for interactions between molecules.
 * Stores the interaction pair and writes it at a later point into a file.
 * */
class InteractionLogger {
public:
    //! @brief interaction representation
    using interaction_t = std::pair<Molecule, Molecule>;
    //! @brief constructor
    explicit InteractionLogger(std::string  filename = "InteractionLog") : _outputFile(std::move(filename)), _interactions(), fileID(0) { }

    //! @brief log this pair of molecules as one interaction
    void log(Molecule& mi, Molecule& mj) {
        _interactions.emplace_back(mi, mj);
    }

    //! @brief write output file of all interactions
    void writeLogFile() {
        std::ofstream file;
        file.open(_outputFile + std::to_string(fileID) + ".log");
        if(!file.is_open()) exit(250);

        for(auto& [mi, mj] : _interactions) {
            writeMolecule(mi, file);
            file << "#\n";
            writeMolecule(mj, file);
            file << "##\n";
        }
        file.flush();
        file.close();
        fileID++;
    }

    //! @brief clears the recorded interactions
    void clearLog() {
        _interactions.clear();
    }
private:
    //! @brief output file name
    std::string _outputFile;
    //! @brief storage for all interactions
    std::vector<interaction_t> _interactions;
    //! @brief running file ID
    int fileID;

    //! @brief write one molecule to file
    static void writeMolecule(Molecule& m, std::ofstream& file) {
        file << m.getID() << "|";
        file << m.r(0) << "|";
        file << m.r(1) << "|";
        file << m.r(2) << "|";
        file << m.v(0) << "|";
        file << m.v(1) << "|";
        file << m.v(2) << "\n";
    }
};


#endif //MARDYN_INTERACTIONLOGGER_H
