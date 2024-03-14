//
// Created by alex on 22.05.23.
//

#ifndef MARDYN_VTKFILE_H
#define MARDYN_VTKFILE_H


#include <string>
#include "Molecule.h"

/**
 * Representation of input vtk file.
 * Can read in file and store molecules.
 * */
class VTKFile {
public:
    /**
     * @brief reads in standard X_nodeX_XXX.vtu file and loads in all molecules.
     * */
    void readXML(const std::string& file);

    /**
     * @brief returns reference to buffer containing all molecules. is only populated after xml file has been read in.
     * */
    const std::vector<Molecule>& getMolecules();

private:
    //! @brief buffer in which all molecules are stored in
    std::vector<Molecule> _molecules;
};


#endif //MARDYN_VTKFILE_H
