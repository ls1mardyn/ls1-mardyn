/**
 * @file ODFCellProcessor.h
 * @author F. Gratl
 * @date 10/11/19
 */

#pragma once

#include <io/ODF.h>
#include "CellProcessor.h"

class ODFCellProcessor : public CellProcessor {

 public:

  ODFCellProcessor(const double cutoffRadius,
                   ODF *const odf,
                   const array<double, 3> &simBoxSize);

  void initTraversal() override;
  void preprocessCell(ParticleCell &cell) override;
  void processCellPair(ParticleCell &cell1, ParticleCell &cell2, bool sumAll) override;
  void processCell(ParticleCell &cell) override;
  double processSingleMolecule(Molecule *m1, ParticleCell &cell2) override;
  void postprocessCell(ParticleCell &cell) override;
  void endTraversal() override;

 private:

  std::array<double, 3> calcUpVec1(const Molecule & molecule1);

  ODF * const _odf;
  const array<double, 3> _simBoxSize;
};



