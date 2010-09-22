/*
 * VTKGridWriter.h
 *
 * @Date: 24.08.2010
 * @Author: eckhardw
 */

#ifndef VTKGRIDWRITER_H_
#define VTKGRIDWRITER_H_

#include "io/OutputBase.h"
#include "io/vtk/VTKGridCell.h"
#include "io/vtk/VTKGridVertex.h"

class LinkedCells;
class VTKGridWriterImplementation;

class VTKGridWriter : public OutputBase {

private:

	const unsigned int _writeFrequency;

	const std::string _fileName;

	const LinkedCells& _container;

	VTKGridCell* _cells;

	VTKGridVertex* _vertices;

	//! length of the cells array
	int _numCells;

	//! length of the vertices array
	int _numVertices;

public:

	VTKGridWriter(unsigned int writeFrequency, const std::string& fileName, const LinkedCells& container);
	virtual ~VTKGridWriter();

	virtual void initOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);

	virtual void doOutput(
			ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			Domain* domain, unsigned long simstep,
			std::list<ChemicalPotential>* lmu
	);

	virtual void finishOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);

private:

	void setupVTKGrid();

	void releaseVTKGrid();

	void getCellData(VTKGridCell& cell);

	void outputParallelVTKFile(unsigned int numProcs, unsigned long simstep,
			VTKGridWriterImplementation& impl);

};

#endif /* VTKGRIDWRITER_H_ */
