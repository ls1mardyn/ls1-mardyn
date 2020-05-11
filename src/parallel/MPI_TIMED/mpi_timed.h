//
// Created by jeremyharisch on 11.05.20.
//

#ifndef MARDYN_MPI_TIMED_H
#define MARDYN_MPI_TIMED_H

#include <mpi.h>

#include "parallel/MPI_TIMED/ProcessTimer.h"



extern ProcessTimer _processTimer;

extern "C"
int MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest,
						int tag, MPI_Comm comm) {
	std::cout <<"##################send########################" << std::endl;
	// Start the timer
	_processTimer.startTimer();
	// Call the original MPI function
	int result = PMPI_Send(buf, count, datatype, dest, tag, comm);
	// Stop the timer
	_processTimer.stopTimer();
	return result;
}

extern "C"
inline int MPI_Recv (void *buf,int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status) {
	std::cout <<"##################recv########################" << std::endl;
	// Start the timer
	_processTimer.startTimer();
	// Call the original MPI function
	int result = PMPI_Recv(buf, count, datatype, source, tag, comm, status);
	// Stop the timer
	_processTimer.stopTimer();
	return result;
}

extern "C"
inline int MPI_File_open(MPI_Comm comm, const char *filename, int amode, MPI_Info info, MPI_File *fh) {
	std::cout <<"##################open########################" << std::endl;
	// Start the timer
	_processTimer.startTimer();
	// Call the original MPI function
	int result = PMPI_File_open(comm, filename, amode, info, fh);
	// Stop the timer
	_processTimer.stopTimer();
	return result;
}

extern "C"
inline int MPI_File_close(MPI_File * fh){
	std::cout <<"##################close########################" << std::endl;
	// Start the timer
	_processTimer.startTimer();
	// Call the original MPI function
	int result = PMPI_File_close(fh);
	// Stop the timer
	_processTimer.stopTimer();
	return result;
}

extern "C"
inline int MPI_File_write(MPI_File mpi_fh, void *buf, int count, MPI_Datatype datatype, MPI_Status *status) {
	std::cout <<"##################write########################" << std::endl;
	// Start the timer
	_processTimer.startTimer();
	// Call the original MPI function
	int result = PMPI_File_write(mpi_fh, buf, count, datatype, status);
	// Stop the timer
	_processTimer.stopTimer();
	return result;
}

extern "C"
int MPI_Allreduce(void* send_data, void* recv_data, int count, MPI_Datatype datatype, MPI_Op op,
				  MPI_Comm communicator) {
	std::cout <<"##################reduce########################" << std::endl;
	// Start the timer
	_processTimer.startTimer();
	// Call the original MPI function
	int result = PMPI_Allreduce(send_data, recv_data, count, datatype, op, communicator);
	// Stop the timer
	_processTimer.stopTimer();
	return result;
}

extern "C"
inline int MPI_Iprobe(int source, int tag, MPI_Comm comm, int *flag, MPI_Status * status) {
	// Start the timer
	_processTimer.startTimer();
	// Call the original MPI function
	int result = PMPI_Iprobe(source, tag, comm, flag, status);
	// Stop the timer
	_processTimer.stopTimer();
	return result;
}

extern "C"
inline int MPI_Test(MPI_Request * request, int *flag, MPI_Status * status) {
	std::cout <<"##################test########################" << std::endl;
	// Start the timer
	_processTimer.startTimer();
	// Call the original MPI function
	int result = PMPI_Test(request, flag, status);
	// Stop the timer
	_processTimer.stopTimer();
	return result;
}


#endif // MARDYN_MPI_TIMED_H
