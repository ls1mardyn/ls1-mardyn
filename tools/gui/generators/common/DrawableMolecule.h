/*
 * DrawableMolecule.h
 *
 *  Created on: Mar 25, 2011
 *      Author: kovacevt
 */
#include "Objects/Object.h"
#include <molecules/Molecule.h>
#include <vector>

#ifndef DRAWABLEMOLECULE_H_
#define DRAWABLEMOLECULE_H_

const int MAXINT = 32767;
/**
 * Class DrawableMolecule depicts a molecule used in the simulation
 */
class DrawableMolecule : public Object {
private:

	Position _x;
	Vector _v;
	unsigned long _id;
	int _numComponents;
	int _cid; // component-id

public:
	/**
	 * Constructors
	 */
	DrawableMolecule(const Molecule& molecule, int numComponents);

	DrawableMolecule();

	/**
	 * Destructor
	 */
	virtual ~DrawableMolecule();

	/**
	 * Draws the molecule
	 * @param numMols - number of molecules to be drawn, used for being able to specify the ID
	 */

	std::vector<std::string> getDrawableValues() const;

	vtkSmartPointer<vtkActor> draw(std::string valueName);

};

#endif /* DrawableMolecule_H_ */
