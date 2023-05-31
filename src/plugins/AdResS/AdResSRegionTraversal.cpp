//
// Created by alex on 16.05.23.
//

#include "AdResSRegionTraversal.h"
#include "utils/threeDimensionalMapping.h"
#include "AdResSForceAdapter.h"

AdResSRegionTraversal::AdResSRegionTraversal(std::array<double, 3> checkLow, std::array<double, 3> checkHigh,
                                             ParticleContainer* particleContainer,
                                             std::unordered_map<unsigned long, Resolution>& compMap) :
                                             _cellPairOffsets8Pack(), _particleContainer(particleContainer),
                                             _cells(), _end(), _dims(), _comp_to_res(compMap), _cutoff(0), _cutoff2(0),
                                             _ljcutoff2(0), _cellSize(0), _globLen({0,0,0}){
    using threeDimensionalMapping::threeToOneD;
    using std::make_pair;

    // compute the dimensions of the region that has to be handled (size in cells)
    _cutoff = _simulation.getcutoffRadius();
    _cutoff2 = _cutoff * _cutoff;
    _ljcutoff2 = _simulation.getLJCutoff();
    _ljcutoff2 *= _ljcutoff2;
    _cellSize = _cutoff;

    for(int d = 0; d < 3; d++) {
        _globLen[d] = _simulation.getDomain()->getGlobalLength(d);
    }

    for (int d = 0; d < 3; ++d) {
        _dims[d] = static_cast<long>((checkHigh[d] - checkLow[d]) / _cellSize);
        if(_dims[d] == 0) _dims[d] = 1;
        _dims[d] += 1; // need 1 dummy cell in every direction to get remaining interactions
        _end[d] = _dims[d] - 1;
    }

    // generate offsets in our created cells structure
    {
        long int o   = threeToOneD(0l, 0l, 0l, _dims); // origin
        long int x   = threeToOneD(1l, 0l, 0l, _dims); // displacement to the right
        long int y   = threeToOneD(0l, 1l, 0l, _dims); // displacement ...
        long int z   = threeToOneD(0l, 0l, 1l, _dims);
        long int xy  = threeToOneD(1l, 1l, 0l, _dims);
        long int yz  = threeToOneD(0l, 1l, 1l, _dims);
        long int xz  = threeToOneD(1l, 0l, 1l, _dims);
        long int xyz = threeToOneD(1l, 1l, 1l, _dims);

        int i = 0;
        // if incrementing along X, the following order will be more cache-efficient:
        _cellPairOffsets8Pack[i++] = make_pair(o, o  );
        _cellPairOffsets8Pack[i++] = make_pair(o, y  );
        _cellPairOffsets8Pack[i++] = make_pair(y, z  );
        _cellPairOffsets8Pack[i++] = make_pair(o, z  );
        _cellPairOffsets8Pack[i++] = make_pair(o, yz );

        _cellPairOffsets8Pack[i++] = make_pair(x, yz );
        _cellPairOffsets8Pack[i++] = make_pair(x, y  );
        _cellPairOffsets8Pack[i++] = make_pair(x, z  );
        _cellPairOffsets8Pack[i++] = make_pair(o, x  );
        _cellPairOffsets8Pack[i++] = make_pair(o, xy );
        _cellPairOffsets8Pack[i++] = make_pair(xy, z );
        _cellPairOffsets8Pack[i++] = make_pair(y, xz );
        _cellPairOffsets8Pack[i++] = make_pair(o, xz );
        _cellPairOffsets8Pack[i++] = make_pair(o, xyz);
    }

    //need to create iterators in single thread context, as otherwise iterators would only work partially
    _cells.resize(_dims[0]*_dims[1]*_dims[2]);
    for(long indX = 0; indX < _dims[0]; indX++) {
        for(long indY = 0; indY < _dims[1]; indY++) {
            for(long indZ = 0; indZ < _dims[2]; indZ++) {
                std::array<double,3> low_offset = {_cellSize * indX, _cellSize * indY, _cellSize * indZ};
                std::array<double,3> high_offset = {_cellSize * (indX + 1), _cellSize * (indY + 1), _cellSize * (indZ + 1)};
                std::array<double,3> low = {checkLow[0] + low_offset[0], checkLow[1] + low_offset[1], checkLow[2] + low_offset[2]};
                std::array<double,3> high = {checkLow[0] + high_offset[0], checkLow[1] + high_offset[1], checkLow[2] + high_offset[2]};
                if(indX == _dims[0] - 2) high[0] = checkHigh[0]; // - 1 for last index and - 1 for one dummy cell
                if(indY == _dims[1] - 2) high[1] = checkHigh[1]; // - 1 for last index and - 1 for one dummy cell
                if(indZ == _dims[2] - 2) high[2] = checkHigh[2]; // - 1 for last index and - 1 for one dummy cell
                if(indX == _dims[0] - 1 || indY == _dims[1] - 1 || indZ == _dims[2] - 1) high = low; // make empty dummy cells
                _cells[threeToOneD(indX, indY, indZ, _dims)] = _particleContainer->regionIterator(low.data(), high.data(), ParticleIterator::ALL_CELLS);
            }
        }
    }
}

void AdResSRegionTraversal::traverse(AdResSForceAdapter& forceAdapter, FPRegion& region, bool invert) {
    #if defined(_OPENMP)
    #pragma omp parallel
    #endif
    {
        for (long col = 0; col < 8; ++col) {
            std::array<long, 3> begin = threeDimensionalMapping::oneToThreeD(col, {2,2,2});
            #if defined(_OPENMP)
            #pragma omp for collapse(3) nowait
            #endif
            for (long cell_z = begin[2]; cell_z < _end[2]; cell_z += 2) {
                for (long cell_y = begin[1]; cell_y < _end[1]; cell_y += 2) {
                    for (long cell_x = begin[0]; cell_x < _end[0]; cell_x += 2) {
                        long baseIndex = threeDimensionalMapping::threeToOneD(cell_x, cell_y, cell_z, _dims);
                        const int num_pairs = _cellPairOffsets8Pack.size();
                        for(int j = 0; j < num_pairs; ++j) {
                            std::pair<long, long> current_pair = _cellPairOffsets8Pack[j];

                            unsigned offset1 = current_pair.first;
                            unsigned cellIndex1 = baseIndex + offset1;

                            unsigned offset2 = current_pair.second;
                            unsigned cellIndex2 = baseIndex + offset2;
                            RegionParticleIterator cell1 = RegionParticleIterator(this->_cells[cellIndex1]) ;
                            RegionParticleIterator cell2 = RegionParticleIterator(this->_cells[cellIndex2]);

                            if(cellIndex1 == cellIndex2) processCell(cell1, forceAdapter, region, invert);
                            else processCellPair(cell1, cell2, forceAdapter, region, invert);

                        }
                    }
                }
            }
            #if defined(_OPENMP)
            #pragma omp barrier
            #endif
            // this barrier is needed, since we have a nowait
        } // end for color
    } // end omp parallel
}

void AdResSRegionTraversal::processCell(RegionParticleIterator &cellIt, AdResSForceAdapter& forceAdapter, FPRegion& region, bool invert) {
    std::array<double, 3> dist = {0,0,0};
    // let every molecule in the box created by checkLow and checkHigh interact with the hybrid molecules in this region
    for(; cellIt.isValid(); ++cellIt) {
        Molecule& m1 = *cellIt; // this can be of any type

        auto itInner = cellIt;
        ++itInner;
        for(; itInner.isValid(); ++itInner) {
            Molecule& m2 = *itInner; // must be hybrid
            mardyn_assert(&m1 != &m2);

            //check if inner is FP or CG -> skip
            if(_comp_to_res[m2.componentid()] != Hybrid) continue;

            //check distance
            double dd = m1.dist2(m2, dist.data());
            if(dd < _cutoff2) {
                PairType type = MOLECULE_MOLECULE;
                if(!FPRegion::isInnerPoint(m1.r_arr(), {0,0,0}, _globLen) ||
                   !FPRegion::isInnerPoint(m2.r_arr(), {0,0,0}, _globLen)) {
                    type = MOLECULE_HALOMOLECULE;
                }

                //recompute force and invert it -> last bool param is true
                forceAdapter.processPair(m1, m2, dist, type, dd, (dd < _ljcutoff2), invert,
                                          _comp_to_res, region);
            }
        }
    }
}

void AdResSRegionTraversal::processCellPair(RegionParticleIterator &cell1, RegionParticleIterator& cell2, AdResSForceAdapter& forceAdapter, FPRegion& region, bool invert) {
    std::array<double, 3> dist = {0,0,0};
    // let every molecule in the box created by checkLow and checkHigh interact with the hybrid molecules in this region
    for(auto cell1It = RegionParticleIterator(cell1); cell1It.isValid(); ++cell1It) {
        Molecule& m1 = *cell1It; // this can be of any type

        for(auto cell2It = RegionParticleIterator(cell2); cell2It.isValid(); ++cell2It) {
            Molecule& m2 = *cell2It; // must be hybrid
            if(_comp_to_res[m2.componentid()] != Hybrid) continue;
            //check distance
            double dd = m1.dist2(m2, dist.data());
            if(dd < _cutoff2) {
                PairType type = MOLECULE_MOLECULE;
                if(!FPRegion::isInnerPoint(m1.r_arr(), {0,0,0}, _globLen) ||
                   !FPRegion::isInnerPoint(m2.r_arr(), {0,0,0}, _globLen)) {
                    type = MOLECULE_HALOMOLECULE;
                }
                //recompute force and invert it -> last bool param is true
                forceAdapter.processPair(m1, m2, dist, type, dd, (dd < _ljcutoff2), invert,
                                         _comp_to_res, region);
            }
        }
    }

    for(auto cell1It = RegionParticleIterator(cell1); cell1It.isValid(); ++cell1It) {
        Molecule& m1 = *cell1It; // this can be of any type
        if(_comp_to_res[m1.componentid()] != Hybrid) continue;

        for(auto cell2It = RegionParticleIterator(cell2); cell2It.isValid(); ++cell2It) {
            Molecule& m2 = *cell2It; // must be hybrid
            if(_comp_to_res[m2.componentid()] == Hybrid) continue;
            //check distance
            double dd = m1.dist2(m2, dist.data());
            if(dd < _cutoff2) {
                PairType type = MOLECULE_MOLECULE;
                if(!FPRegion::isInnerPoint(m1.r_arr(), {0,0,0}, _globLen) ||
                   !FPRegion::isInnerPoint(m2.r_arr(), {0,0,0}, _globLen)) {
                    type = MOLECULE_HALOMOLECULE;
                }
                //recompute force and invert it -> last bool param is true
                forceAdapter.processPair(m1, m2, dist, type, dd, (dd < _ljcutoff2), invert,
                                         _comp_to_res, region);
            }
        }
    }
}
