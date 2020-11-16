#include "MPIKDNode.h"

#include "utils/mardyn_assert.h"

#include "utils/Logger.h"

#define assertion(x) mardyn_assert(x)
#define assertion1(x,y) { if (!(x)) std::cerr << (y) << std::endl; mardyn_assert(x);}

using Log::global_log;

MPIKDNodePacked::MPIKDNodePacked(const std::bitset<3>& coversWholeDomain, const int& numProcs,
		const int* lowCorner, const int* highCorner, const int& nodeID, const int& owningProc,
		const int& firstChildID, const int& secondChildID, const int& nextSendingProcess,
		const double& load, const double& OptimalLoadPerProcess, const double& expectedDeviation,
		const double& deviation, const int& level):
_nodeID(nodeID),
_owningProc(owningProc),
_firstChildID(firstChildID),
_secondChildID(secondChildID),
_nextSendingProcess(nextSendingProcess),
_load(load),
_OptimalLoadPerProcess(OptimalLoadPerProcess),
_deviationLowerBound(expectedDeviation),
_deviation(deviation),
_level(level) {
  setCoversWholeDomain(coversWholeDomain);
  setNumProcs(numProcs);
  for (int i = 0; i < 3; i++) {
    setLowCorner(i, lowCorner[i]);
    setHighCorner(i, highCorner[i]);
  }
  assertion((31 < (8 * sizeof(int))));
  
}

MPIKDNodePacked::~MPIKDNodePacked() { }


std::bitset<3> MPIKDNodePacked::getCoversWholeDomain() const {
  int mask = (int) (1 << (3)) - 1 ;
  mask = static_cast<int>(mask << (0));
  int tmp = static_cast<int>(_packedRecords0 & mask);
  tmp = static_cast<int>(tmp >> (0));
  std::bitset<3> result = tmp;
  return result;
}


void MPIKDNodePacked::setCoversWholeDomain(const std::bitset<3>& coversWholeDomain) {
  int mask = (int) (1 << (3)) - 1 ;
  mask = static_cast<int>(mask << (0));
  _packedRecords0 = static_cast<int>(_packedRecords0 & ~mask);
  _packedRecords0 = static_cast<int>(_packedRecords0 | coversWholeDomain.to_ulong() << (0));
}



bool MPIKDNodePacked::getCoversWholeDomain(int elementIndex) const {
  assertion(elementIndex>=0);
  assertion(elementIndex<3);
  int mask = 1 << (0);
  mask = mask << elementIndex;
  return (_packedRecords0& mask);
}



void MPIKDNodePacked::setCoversWholeDomain(int elementIndex, const bool& coversWholeDomain) {
  assertion(elementIndex>=0);
  assertion(elementIndex<3);
  //assertion(!coversWholeDomain || coversWholeDomain==1);
  int shift        = 0 + elementIndex; 
  int mask         = 1     << (shift);
  int shiftedValue = coversWholeDomain << (shift);
  _packedRecords0 = _packedRecords0 & ~mask;
  _packedRecords0 = _packedRecords0 |  shiftedValue;
}



void MPIKDNodePacked::flipCoversWholeDomain(int elementIndex) {
  assertion(elementIndex>=0);
  assertion(elementIndex<3);
  int mask = 1 << (0);
  mask = mask << elementIndex;
  _packedRecords0^= mask;
}



int MPIKDNodePacked::getNumProcs() const {
  int mask =  (1 << (28)) - 1;
  mask = static_cast<int>(mask << (3));
  int tmp = static_cast<int>(_packedRecords0 & mask);
  tmp = static_cast<int>(tmp >> (3));
  assertion(( tmp >= 0 &&  tmp <= 268000000));
  return (int) tmp;
}



void MPIKDNodePacked::setNumProcs(const int& numProcs) {
  assertion((numProcs >= 0 && numProcs <= 268000000));
  int mask =  (1 << (28)) - 1;
  mask = static_cast<int>(mask << (3));
  _packedRecords0 = static_cast<int>(_packedRecords0 & ~mask);
  _packedRecords0 = static_cast<int>(_packedRecords0 | numProcs << (3));
}


const int* MPIKDNodePacked::getLowCorner() const {
  return _lowCorner;
}


int MPIKDNodePacked::getLowCorner(int elementIndex) const {
  assertion(elementIndex>=0);
  assertion(elementIndex<3);
  return _lowCorner[elementIndex];
  
}


void MPIKDNodePacked::setLowCorner(int elementIndex, const int& lowCorner) {
  assertion(elementIndex>=0);
  assertion(elementIndex<3);
  _lowCorner[elementIndex]= lowCorner;
  
}


const int* MPIKDNodePacked::getHighCorner() const {
  return _highCorner;
}


int MPIKDNodePacked::getHighCorner(int elementIndex) const {
  assertion(elementIndex>=0);
  assertion(elementIndex<3);
  return _highCorner[elementIndex];
  
}



void MPIKDNodePacked::setHighCorner(int elementIndex, const int& highCorner) {
  assertion(elementIndex>=0);
  assertion(elementIndex<3);
  _highCorner[elementIndex]= highCorner;
  
}



int MPIKDNodePacked::getNodeID() const {
  return _nodeID;
}



void MPIKDNodePacked::setNodeID(const int& nodeID) {
  _nodeID = nodeID;
}



int MPIKDNodePacked::getOwningProc() const {
  return _owningProc;
}



void MPIKDNodePacked::setOwningProc(const int& owningProc) {
  _owningProc = owningProc;
}



int MPIKDNodePacked::getFirstChildID() const {
  return _firstChildID;
}



void MPIKDNodePacked::setFirstChildID(const int& firstChildID) {
  _firstChildID = firstChildID;
}



int MPIKDNodePacked::getSecondChildID() const {
  return _secondChildID;
}



void MPIKDNodePacked::setSecondChildID(const int& secondChildID) {
  _secondChildID = secondChildID;
}



int MPIKDNodePacked::getNextSendingProcess() const {
  return _nextSendingProcess;
}



void MPIKDNodePacked::setNextSendingProcess(const int& nextSendingProcess) {
  _nextSendingProcess = nextSendingProcess;
}



double MPIKDNodePacked::getLoad() const {
  return _load;
}



void MPIKDNodePacked::setLoad(const double& load) {
  _load = load;
}



double MPIKDNodePacked::getOptimalLoadPerProcess() const {
  return _OptimalLoadPerProcess;
}



void MPIKDNodePacked::setOptimalLoadPerProcess(const double& OptimalLoadPerProcess) {
  _OptimalLoadPerProcess = OptimalLoadPerProcess;
}


double MPIKDNodePacked::getDeviationLowerBound() const {
	return _deviationLowerBound;
}

double MPIKDNodePacked::getDeviation() const {
	return _deviation;
}

int MPIKDNodePacked::getLevel() const {
	return _level;
}



std::string MPIKDNodePacked::toString() const {
  std::ostringstream stringstr;
  toString(stringstr);
  return stringstr.str();
}

void MPIKDNodePacked::toString (std::ostream& out) const {
  out << "("; 
  out << "coversWholeDomain:[";
  for (int i = 0; i < 3-1; i++) {
    out << getCoversWholeDomain(i) << ",";
  }
  out << getCoversWholeDomain(3-1) << "]";
  out << ",";
  out << "numProcs:" << getNumProcs();
  out << ",";
  out << "lowCorner:[";
  for (int i = 0; i < 3-1; i++) {
    out << getLowCorner(i) << ",";
  }
  out << getLowCorner(3-1) << "]";
  out << ",";
  out << "highCorner:[";
  for (int i = 0; i < 3-1; i++) {
    out << getHighCorner(i) << ",";
  }
  out << getHighCorner(3-1) << "]";
  out << ",";
  out << "nodeID:" << getNodeID();
  out << ",";
  out << "owningProc:" << getOwningProc();
  out << ",";
  out << "firstChildID:" << getFirstChildID();
  out << ",";
  out << "secondChildID:" << getSecondChildID();
  out << ",";
  out << "nextSendingProcess:" << getNextSendingProcess();
  out << ",";
  out << "load:" << getLoad();
  out << ",";
  out << "OptimalLoadPerProcess:" << getOptimalLoadPerProcess();
  out <<  ")";
}



MPI_Datatype MPIKDNodePacked::Datatype = 0;

/**
 * TODO: incorporate changes of rev. 1000
 * However, at the moment I'm not quite sure how that works...
 */
void MPIKDNodePacked::initDatatype() {
	MPIKDNodePacked dummyMPIKDNodePacked;

	const int Attributes = 13;
	MPI_Datatype subtypes[Attributes] = {
			MPI_INT,		 //lowCorner
			MPI_INT,		 //highCorner
			MPI_INT,		 //nodeID
			MPI_INT,		 //owningProc
			MPI_INT,		 //firstChildID
			MPI_INT,		 //secondChildID
			MPI_INT,		 //nextSendingProcess
			MPI_DOUBLE,		 //load
			MPI_DOUBLE,		 //OptimalLoadPerProcess
			MPI_DOUBLE,		 // expectedDeviation;
			MPI_DOUBLE,		 // deviation;
			MPI_INT,		 // level;
			MPI_INT,		 //_packedRecords0
	};

	int blocklen[Attributes] = {
			3,		 //lowCorner
			3,		 //highCorner
			1,		 //nodeID
			1,		 //owningProc
			1,		 //firstChildID
			1,		 //secondChildID
			1,		 //nextSendingProcess
			1,		 //load
			1,		 //OptimalLoadPerProcess
			1,		 // expectedDeviation;
			1,		 // deviation;
			1,		 // level;
			1,		 //_packedRecords0
	};

	MPI_Aint     disp[Attributes];

	MPI_Aint base;
	MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyMPIKDNodePacked))), &base);
	MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyMPIKDNodePacked._lowCorner[0]))), 		&disp[0] );
	MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyMPIKDNodePacked._highCorner[0]))), 		&disp[1] );
	MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyMPIKDNodePacked._nodeID))), 		&disp[2] );
	MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyMPIKDNodePacked._owningProc))), 		&disp[3] );
	MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyMPIKDNodePacked._firstChildID))), 		&disp[4] );
	MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyMPIKDNodePacked._secondChildID))), 		&disp[5] );
	MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyMPIKDNodePacked._nextSendingProcess))), 		&disp[6] );
	MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyMPIKDNodePacked._load))), 		&disp[7] );
	MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyMPIKDNodePacked._OptimalLoadPerProcess))), 		&disp[8] );
	MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyMPIKDNodePacked._deviationLowerBound))), 		&disp[9] );
	MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyMPIKDNodePacked._deviation))), 		&disp[10] );
	MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyMPIKDNodePacked._level))), 		&disp[11] );
	MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyMPIKDNodePacked._packedRecords0))), 		&disp[12] );

	for (int i=1; i<Attributes; i++) {
		if (!(disp[i] > disp[i-1])) {
			global_log->debug() << "i=" << i << std::endl;
		}
		assertion1( disp[i] > disp[i-1], i );
	}
	for (int i=0; i<Attributes; i++) {
		disp[i] -= base;
	}
	MPI_Type_create_struct( Attributes, blocklen, disp, subtypes, &MPIKDNodePacked::Datatype );
	MPI_Type_commit( &MPIKDNodePacked::Datatype );
}


void MPIKDNodePacked::shutdownDatatype() {
	MPI_Type_free( &MPIKDNodePacked::Datatype );
}

