#include "tests/datastructures/CellTest.h"
#include "datastructures/Cell.h"
#include "molecules/SimpleMolecule.h"

#include <cmath>
#include <iostream>
using namespace std;

utils::Log CellTest::_log("CellTest");

CellTest::CellTest():
  TestCase( "CellTest" ) {
}

CellTest::~CellTest(){
}


void CellTest::allTests(){
  
  datastructures::Cell<SimpleMolecule> cell = datastructures::Cell<SimpleMolecule>();
  // list of molecule pointers
  list<SimpleMolecule*>& moleculePointers = cell.getParticlePointers();;
  std::list<SimpleMolecule*>::iterator MoleculeIterator;
    
  // check whether the list is empty
  for(MoleculeIterator=moleculePointers.begin();MoleculeIterator!=moleculePointers.end();++MoleculeIterator){
    _log.error("allTests", "newly created list must be empty");
  }
  // check whether the list is still empty after removing all molecules 
  cell.removeAllParticles();
  moleculePointers = cell.getParticlePointers();
  for(MoleculeIterator=moleculePointers.begin();MoleculeIterator!=moleculePointers.end();++MoleculeIterator){
    _log.error("allTests", "newly created list must still be empty after removing molecules");
  }
  // insert some molecules into the cell
  SimpleMolecule molecule1(8 ,3, 3.01, 4.42,-1.6, 7.11,-1.35, -0.63);
  SimpleMolecule molecule2(2 ,5, 2.00, 2.51, 0.1, 1.61, 3.12,  3.15);
  SimpleMolecule molecule3(15,1,-5.01, 8.42,-1.6, 0.38,-0.42,  2.13);
  SimpleMolecule molecule4(23,5, 1.24,12.51, 0.1,-5.02,-5.83,  7.45);
  
  cell.addParticle(&molecule1);
  cell.addParticle(&molecule2);
  cell.addParticle(&molecule3);
  cell.addParticle(&molecule4);
  // retrieve the list of molecule pointers from the cell,
  // count the molecules and check whether the correct molecules are in the cell
  moleculePointers = cell.getParticlePointers();
  bool moleculeExistent[4]; // will be set to true when the corresponding molecule is in the cell
  int moleculeCount = 0;
  for(int i=0; i<4; i++){
    moleculeExistent[i] = false; 
  }
  for(MoleculeIterator=moleculePointers.begin();MoleculeIterator!=moleculePointers.end();++MoleculeIterator){
    moleculeCount++;
    if ((*MoleculeIterator)->getId() == 8){
      moleculeExistent[0] = true;
    }
    if ((*MoleculeIterator)->getId() == 2){
      moleculeExistent[1] = true;
    }
    if ((*MoleculeIterator)->getId() == 15){
      moleculeExistent[2] = true;
    }
    if ((*MoleculeIterator)->getId() == 23){
      moleculeExistent[3] = true;
    }
  }
  // error if one of the needed molecules doesn't exist in the cell
  bool correctMoleculesExistent = true;
  for(int i=0; i<4; i++){
    correctMoleculesExistent = (correctMoleculesExistent && moleculeExistent[i]); 
  }  
  if(moleculeCount != 4){
    _log.error("allTests","Wrong number of Molecules in Cell"); 
  }
  else if(!correctMoleculesExistent){
    _log.error("allTests","Wrong Molecules in Cell"); 
  }
  
  // check whether the list is empty after removing all molecules 
  cell.removeAllParticles();
  moleculePointers = cell.getParticlePointers();
  for(MoleculeIterator=moleculePointers.begin();MoleculeIterator!=moleculePointers.end();++MoleculeIterator){
    _log.error("allTests", "newly created list must still be empty after removing molecules");
  } 
}

void CellTest::run(){
  allTests();
}
