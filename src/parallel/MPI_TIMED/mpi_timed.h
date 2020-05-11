//
// Created by jeremyharisch on 11.05.20.
//

#ifndef MARDYN_MPI_TIMED_H
#define MARDYN_MPI_TIMED_H

#include <mpi.h>

#include "parallel/MPI_TIMED/ProcessTimer.h"


public:

	extern "C"
	int MPI_Send(void *buf, int count, MPI_Datatype dtype, int dest, int tag, MPI_Comm comm) {
	  // TODO: Update call counters and start the timer

	  // Call the original MPI function
	  int result = PMPI_Send(buf, count, dtype, dest, tag, comm);

	  // TODO: Stop the timer

	  return result;
	}


private:
	ProcessTimer _processTimer;


#endif // MARDYN_MPI_TIMED_H
