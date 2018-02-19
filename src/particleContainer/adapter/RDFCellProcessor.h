/*
 * RDFCellProcessor.h
 *
 * @Date: 19.02.18
 * @Author: tchipevn
 */

#ifndef RDFCELLPROCESSOR_H_
#define RDFCELLPROCESSOR_H_

#include "particleContainer/adapter/CellProcessor.h"
#include "particleContainer/ParticleCellForwardDeclaration.h"

class RDF;

class RDFCellProcessor : public CellProcessor {

private:
	RDF* const _rdf;

public:
	RDFCellProcessor& operator=(const RDFCellProcessor&) = delete;

	RDFCellProcessor(const double cutoffRadius, RDF* rdf) :
			CellProcessor(cutoffRadius, cutoffRadius), _rdf(rdf) {
	}

	~RDFCellProcessor() {}

	void initTraversal() {}

	void preprocessCell(ParticleCell& /*cell*/) {}

	void processCellPair(ParticleCell& cell1, ParticleCell& cell2, bool sumAll = false);

	double processSingleMolecule(Molecule* /* m1 */, ParticleCell& /*cell2*/) { return 0.0;}

	void processCell(ParticleCell& cell);

	void postprocessCell(ParticleCell& /*cell*/) {}

	void endTraversal() {}
};

#endif /* RDFCELLPROCESSOR_H_ */
