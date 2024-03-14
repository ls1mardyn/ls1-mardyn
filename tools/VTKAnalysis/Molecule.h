//
// Created by alex on 22.05.23.
//

#ifndef MARDYN_MOLECULE_H
#define MARDYN_MOLECULE_H


#include <array>
#include <cstdint>

/**
 * Minimal representation of molecules based on vtk input.
 * */
struct Molecule {
    //! @brief id
    uint64_t id;
    //! @brief component id
    int32_t cid;
    //! @brief position
    std::array<float, 3> r;
};


#endif //MARDYN_MOLECULE_H
