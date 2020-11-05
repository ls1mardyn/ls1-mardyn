/**
 * @file ODFCellProcessor.h
 * @author F. Gratl
 * @date 10/11/19
 */

#pragma once

class ODF;

#include <array>
#include "CellProcessor.h"

/**
 * Class to calculate an orientation distribution function.
 * It calculates the distribution of relative orientations of particles using the ODF class.
 */
class ODFCellProcessor : public CellProcessor {

 public:
  ODFCellProcessor(double cutoffRadius,
                   ODF *odf,
                   const std::array<double, 3> &simBoxSize);

  void initTraversal() override;
  void preprocessCell(ParticleCell &cell) override;
  void processCellPair(ParticleCell &cell1, ParticleCell &cell2, bool sumAll) override;
  void processCell(ParticleCell &cell) override;
  double processSingleMolecule(Molecule *m1, ParticleCell &cell2) override;
  void postprocessCell(ParticleCell &cell) override;
  void endTraversal() override;

 private:

  std::array<double, 3> calcOrientationVector(const Molecule & molecule);

  ODF * const _odf;
  const std::array<double, 3> _simBoxSize;
};



