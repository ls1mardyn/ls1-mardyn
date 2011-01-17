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

/**
 * This class acts as adapter to the VTKGridWriterImplementation, which handles
 * the actual xml writing. It is a friend class of LinkedCells, but reads only
 * its internal data to generate the vtk output.
 */
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

	/**
	 * @param container the LinkedCells particle container. It has to be the same container
	 *                  which is handed in to the methods initOutput() / doOutput() / finishOutput()!
	 */
	VTKGridWriter(unsigned int writeFrequency, const std::string& fileName, const LinkedCells& container);
	virtual ~VTKGridWriter();

	virtual void initOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);

	/**
	 * creates the VTKGrid and sets all the data, which is then written out.
	 *
	 * @note As the gridwriter isn't notified when the particleContainer is changed
	 * (i.e. rebuild()), it is not sufficient to build up the vtkgrid once and only
	 * renew the information every iteration.
	 * Thus every time output is done, the methods setupVTKGrid() and releaseVTKGrid()
	 * are called.
	 */
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
