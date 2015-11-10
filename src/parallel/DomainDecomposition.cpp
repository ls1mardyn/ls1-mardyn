#include "DomainDecomposition.h"

#include "Domain.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "parallel/ParticleData.h"
#include "utils/xmlfileUnits.h"
#include "utils/Logger.h"

using Log::global_log;
using namespace std;

DomainDecomposition::DomainDecomposition() {

	int period[DIM]; // 1(true) when using periodic boundary conditions in the corresponding dimension
	int reorder; // 1(true) if the ranking may be reordered by MPI_Cart_create


	// We create a torus topology, so all boundary conditions are periodic
	for (int d = 0; d < DIM; d++)
		period[d] = 1;
	// Allow reordering of process ranks
	reorder = 1;
	// Find out appropriate grid dimensions
	MPI_CHECK( MPI_Comm_size(MPI_COMM_WORLD, &_numProcs) );

	for (int i = 0; i < DIM; i++) {
		_gridSize[i] = 0;
	}
	MPI_CHECK( MPI_Dims_create( _numProcs, DIM, (int *) &_gridSize ) );

	// Create the communicator
	MPI_CHECK( MPI_Cart_create(MPI_COMM_WORLD, DIM, _gridSize, period, reorder, &_comm) );
	global_log->info() << "MPI grid dimensions: " << _gridSize[0]<<", "<<_gridSize[1]<<", "<<_gridSize[2] << endl;

	// introduce coordinates
	MPI_CHECK( MPI_Comm_rank(_comm, &_rank) );
	MPI_CHECK( MPI_Cart_coords(_comm, _rank, DIM, _coords) );
	global_log->info() << "MPI coordinate of current process: " << _coords[0]<<", "<<_coords[1]<<", "<<_coords[2] << endl;


	for (int d = 0; d < DIM; ++d) {
		if(_gridSize[d] == 1) {
			_coversWholeDomain[d] = true;
		} else {
			_coversWholeDomain[d] = false;
		}
	}

	// Initialize MPI Dataype for the particle exchange once at the beginning.
	ParticleData::setMPIType(_mpi_Particle_data);

//	initCommunicationPartners();
}

DomainDecomposition::~DomainDecomposition() {
	MPI_Type_free(&_mpi_Particle_data);

	if(_comm != MPI_COMM_WORLD) {
		MPI_Comm_free(&_comm);
	}
}

void DomainDecomposition::initCommunicationPartners(ParticleContainer * moleculeContainer, Domain * domain) {

	// corners of the process-specific domain
	double rmin[DIM]; // lower corner
	double rmax[DIM]; // higher corner
	double halo_L[DIM]; // width of the halo strip

	for (int d = 0; d < DIM; d++) {
		rmin[d] = moleculeContainer->getBoundingBoxMin(d);
		rmax[d] = moleculeContainer->getBoundingBoxMax(d);
		halo_L[d] = moleculeContainer->get_halo_L(d);
	}

	int direction;

	for (unsigned short d = 0; d < DIM; d++) {
		// set the ranks
		MPI_CHECK( MPI_Cart_shift(_comm, d, 1, &_partners[d][LOWER]._rank, &_partners[d][HIGHER]._rank) );

		// when moving a particle across a periodic boundary, the molecule position has to change
		// these offset specify for each dimension (x, y and z) and each direction ("left"/lower
		// neighbour and "right"/higher neighbour, how the paritcle coordinates have to be changed.
		// e.g. for dimension x (d=0) and a process on the left boundary of the domain, particles
		// moving to the left get the length of the whole domain added to their x-value
		double offsetLower[DIM];
		double offsetHigher[DIM];
		offsetLower[d] = 0.0;
		offsetHigher[d] = 0.0;

		// process on the left boundary
		if (_coords[d] == 0)
			offsetLower[d] = domain->getGlobalLength(d);
		// process on the right boundary
		if (_coords[d] == _gridSize[d] - 1)
			offsetHigher[d] = -domain->getGlobalLength(d);

		double regToSendLow[DIM]; // Region that belongs to a neighbouring process
		double regToSendHigh[DIM]; // -> regToSendLow

		for (direction = LOWER; direction <= HIGHER; direction++) {
			_CommunicationPartner& partner = _partners[d][direction];

			// set the regions
			for (int i = 0; i < DIM; i++) {
				regToSendLow[i] = rmin[i] - halo_L[i];
				regToSendHigh[i] = rmax[i] + halo_L[i];
			}
			switch (direction) {
			case LOWER:
				regToSendHigh[d] = rmin[d] + halo_L[d];
				break;
			case HIGHER:
				regToSendLow[d] = rmax[d] - halo_L[d];
				break;
			}

			for(int i = 0; i < 3; ++i) {
				partner._regionLow[i] = regToSendLow[i];
				partner._regionHigh[i] = regToSendHigh[i];
			}

			// set the shift
			double shift = 0.0;
			if (direction == LOWER)
				shift = offsetLower[d];
			if (direction == HIGHER)
				shift = offsetHigher[d];

			partner._shift = shift;

		}
	}
}

void DomainDecomposition::readXML(XMLfileUnits& xmlconfig) {
	/* no parameters */
	/* TODO: Maybe add decomposition dimensions, default auto. */
}


void DomainDecomposition::exchangeMolecules(ParticleContainer* moleculeContainer, Domain* domain) {

	// TODO
	static bool initCalled = false;
	if(!initCalled)  {
		initCommunicationPartners(moleculeContainer, domain);
		initCalled = true;
	}

	for (unsigned short d = 0; d < DIM; d++) {
		if (_coversWholeDomain[d]) {
			// use the sequential version
			this->handleDomainLeavingParticles(d, moleculeContainer);
			this->populateHaloLayerWithCopies(d, moleculeContainer);
			continue;
		}

		for (int direction = LOWER; direction <= HIGHER; direction++) {
			_CommunicationPartner& buddy = _partners[d][direction];

			vector<Molecule*> particles;
			moleculeContainer->getRegion(buddy._regionLow, buddy._regionHigh, particles);

			const int n = particles.size();

			// initialize send buffer
			buddy._sendBuf.resize(n);

			for (int i = 0; i < n; ++i) {
				ParticleData::MoleculeToParticleData(buddy._sendBuf[i], *(particles[i]));
				// add offsets for particles transfered over the periodic boundary
				buddy._sendBuf[i].r[d] += buddy._shift;
			}

			MPI_CHECK( MPI_Isend(&(buddy._sendBuf[0]), (int) buddy._sendBuf.size(), _mpi_Particle_data, buddy._rank, 99, _comm, &buddy._sendRequest));
		}

#if 1
		bool allDone = false;
		bool countArrived[2] = {false, false};
		bool messageReceived[2] = {false, false};
		bool messageSent[2] = {false, false};

		while (not allDone) {
			allDone = true;

			// kickstart the send requests
			for(int i = 0; i < 2; ++i) {
				if(not messageSent[i]) {
					_CommunicationPartner& buddy = _partners[d][i];
					int flag = 0;
					MPI_CHECK( MPI_Test(&(buddy._sendRequest), &flag, &buddy._sendStatus) );
					if(flag == 1) {
						messageSent[i] = true;
						buddy._sendBuf.clear();
					}
					allDone &= messageSent[i];
				}
			} // send

			// get the counts
			for(int i = 0; i < 2; ++i) {
				if(not countArrived[i]) {
					_CommunicationPartner& buddy = _partners[d][i];
					int flag = 0;
					MPI_CHECK( MPI_Iprobe(buddy._rank, 99, _comm, &flag, &buddy._recvStatus) );
					if(flag == 1) {
						countArrived[i] = true;
						int numrecv;
						MPI_CHECK( MPI_Get_count(&buddy._recvStatus, _mpi_Particle_data, &numrecv) );
						buddy._recvBuf.resize(numrecv);
						MPI_CHECK( MPI_Irecv(&(buddy._recvBuf[0]), numrecv, _mpi_Particle_data, buddy._rank, 99, _comm, &buddy._recvRequest));
					}
				}
			} // probe recv

			// unpack molecules
			for(int i = 0; i < 2; ++i) {
				if(countArrived[i] and not messageReceived[i]) {
					_CommunicationPartner& buddy = _partners[d][i];
					int flag = 0;
					MPI_CHECK( MPI_Test(&(buddy._recvRequest), &flag, &buddy._recvStatus) );
					if(flag == 1) {
						messageReceived[i] = true;
						int numrecv = buddy._recvBuf.size();

						for (int i = 0; i < numrecv; i++) {
							Molecule *m;
							ParticleData::ParticleDataToMolecule(buddy._recvBuf[i], &m);
							moleculeContainer->addParticlePointer(m);
						}
						buddy._recvBuf.clear();

					}
					allDone &= messageReceived[i];
				}
				allDone &= messageReceived[i];
			} // unpack recv

			// TODO: we can catch deadlocks here due to dying MPI processes (or bugs)
		} // close while loop
#else


		for (int direction = LOWER; direction <= HIGHER; direction++) {
			_CommunicationPartner& buddy = _partners[d][direction];

			int numrecv;
			MPI_CHECK( MPI_Probe(buddy._rank, 99, _comm, &buddy._recvStatus) );
			MPI_CHECK( MPI_Get_count(&buddy._recvStatus, _mpi_Particle_data, &numrecv) );

			buddy._recvBuf.resize(numrecv);
			MPI_CHECK( MPI_Irecv(&(buddy._recvBuf[0]), numrecv, _mpi_Particle_data, buddy._rank, 99, _comm, &buddy._recvRequest));
		}

		// Insert molecules into domain
		for (int direction = LOWER; direction <= HIGHER; direction++) {
			_CommunicationPartner& next_buddy = _partners[d][direction];
			_CommunicationPartner& prev_buddy = _partners[d][1 - direction];

			MPI_CHECK( MPI_Wait(&next_buddy._sendRequest, &next_buddy._sendStatus) );
			MPI_CHECK( MPI_Wait(&prev_buddy._recvRequest, &prev_buddy._recvStatus) );
			// insert received molecules into list of molecules
			int numrecv = prev_buddy._recvBuf.size();

			for (int i = 0; i < numrecv; i++) {
				Molecule *m;
				ParticleData::ParticleDataToMolecule(prev_buddy._recvBuf[i], &m);
				moleculeContainer->addParticlePointer(m);
			}
			// clear memory
			next_buddy._sendBuf.clear();
			prev_buddy._recvBuf.clear();
		}
#endif
	} // close dimension
}

void DomainDecomposition::balanceAndExchange(bool balance, ParticleContainer* moleculeContainer, Domain* domain) {
	exchangeMolecules(moleculeContainer, domain);
}

bool DomainDecomposition::procOwnsPos(double x, double y, double z, Domain* domain) {
	if (x < getBoundingBoxMin(0, domain) || x >= getBoundingBoxMax(0, domain))
		return false;
	else if (y < getBoundingBoxMin(1, domain) || y >= getBoundingBoxMax(1, domain))
		return false;
	else if (z < getBoundingBoxMin(2, domain) || z >= getBoundingBoxMax(2, domain))
		return false;
	else
		return true;
}



double DomainDecomposition::getBoundingBoxMin(int dimension, Domain* domain) {
	return _coords[dimension] * domain->getGlobalLength(dimension) / _gridSize[dimension];
}

double DomainDecomposition::getBoundingBoxMax(int dimension, Domain* domain) {
	return (_coords[dimension] + 1) * domain->getGlobalLength(dimension) / _gridSize[dimension];
}

void DomainDecomposition::printDecomp(string filename, Domain* domain) {

	if (_rank == 0) {
		ofstream povcfgstrm(filename.c_str());
		povcfgstrm << "size " << domain->getGlobalLength(0) << " " << domain->getGlobalLength(1) << " " << domain->getGlobalLength(2) << endl;
		povcfgstrm << "cells " << _gridSize[0] << " " << _gridSize[1] << " " << _gridSize[2] << endl;
		povcfgstrm << "procs " << _numProcs << endl;
		povcfgstrm << "data DomainDecomp" << endl;
		povcfgstrm.close();
	}

	for (int process = 0; process < _numProcs; process++) {
		if (_rank == process) {
			ofstream povcfgstrm(filename.c_str(), ios::app);
			povcfgstrm << _coords[2] * _gridSize[0] * _gridSize[1] + _coords[1] * _gridSize[0] + _coords[0] << " " << _rank << endl;
			povcfgstrm.close();
		}
		barrier();
	}
}

unsigned DomainDecomposition::Ndistribution(unsigned localN, float* minrnd, float* maxrnd) {
	unsigned* moldistribution = new unsigned[_numProcs];
	MPI_CHECK( MPI_Allgather(&localN, 1, MPI_UNSIGNED, moldistribution, 1, MPI_UNSIGNED, _comm) );
	unsigned globalN = 0;
	for (int r = 0; r < _rank; r++)
		globalN += moldistribution[r];
	unsigned localNbottom = globalN;
	globalN += moldistribution[_rank];
	unsigned localNtop = globalN;
	for (int r = _rank + 1; r < _numProcs; r++)
		globalN += moldistribution[r];
	delete[] moldistribution;
	*minrnd = (float) localNbottom / globalN;
	*maxrnd = (float) localNtop / globalN;
	return globalN;
}

void DomainDecomposition::assertIntIdentity(int IX) {
	if (_rank)
		MPI_CHECK( MPI_Send(&IX, 1, MPI_INT, 0, 2 * _rank + 17, _comm) );
	else {
		int recv;
		MPI_Status s;
		for (int i = 1; i < _numProcs; i++) {
			MPI_CHECK( MPI_Recv(&recv, 1, MPI_INT, i, 2 * i + 17, _comm, &s) );
			if (recv != IX) {
				global_log->error() << "IX is " << IX << " for rank 0, but " << recv << " for rank " << i << ".\n";
				MPI_Abort(MPI_COMM_WORLD, 911);
			}
		}
		global_log->info() << "IX = " << recv << " for all " << _numProcs << " ranks.\n";
	}
}

void DomainDecomposition::assertDisjunctivity(TMoleculeContainer* mm) {
	Molecule* m;

	if (_rank) {
		int num_molecules = mm->getNumberOfParticles();
		unsigned long *tids;
		tids = new unsigned long[num_molecules];

		int i = 0;
		for (m = mm->begin(); m != mm->end(); m = mm->next()) {
			tids[i] = m->id();
			i++;
		}
		MPI_CHECK( MPI_Send(tids, num_molecules, MPI_UNSIGNED_LONG, 0, 2674 + _rank, _comm) );
		delete[] tids;
		global_log->info() << "Data consistency checked: for results see rank 0." << endl;
	}
	else {
		map<unsigned long, int> check;

		for (m = mm->begin(); m != mm->end(); m = mm->next())
			check[m->id()] = 0;

		MPI_Status status;
		for (int i = 1; i < _numProcs; i++) {
			int num_recv = 0;
			unsigned long *recv;
			MPI_CHECK( MPI_Probe(i, 2674 + i, _comm, &status) );
			MPI_CHECK( MPI_Get_count(&status, MPI_UNSIGNED_LONG, &num_recv) );
			recv = new unsigned long[num_recv];

			MPI_CHECK( MPI_Recv(recv, num_recv, MPI_UNSIGNED_LONG, i, 2674 + i, _comm, &status) );
			for (int j = 0; j < num_recv; j++) {
				if (check.find(recv[j]) != check.end()) {
					global_log->error() << "Ranks " << check[recv[j]] << " and " << i << " both propagate ID " << recv[j] << endl;
					MPI_Abort(MPI_COMM_WORLD, 1);
				}
				else
					check[recv[j]] = i;
			}
			delete[] recv;
		}
		global_log->info() << "Data consistency checked: No duplicate IDs detected among " << check.size() << " entries." << endl;
	}
}
