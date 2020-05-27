//
// Created by jeremyharisch on 11.05.20.
//

#ifndef MARDYN_MPI_TIMED_H
#define MARDYN_MPI_TIMED_H

#include <mpi.h>
#include "parallel/MPI_TIMED/ProcessTimer.h"

extern ProcessTimer _processTimer;

#ifdef __cplusplus
extern "C"{
#endif

int MPI_Pcontrol(const int level, ...) {
	return _processTimer.switchProfiling(level);
}

int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest,
						int tag, MPI_Comm comm) {
	_processTimer.startTimer();
	int result = PMPI_Send(buf, count, datatype, dest, tag, comm);
	_processTimer.stopTimer();
	return result;
}

int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest,
					 int tag, MPI_Comm comm, MPI_Request *request) {
	_processTimer.startTimer();
	int result = PMPI_Isend(buf, count, datatype, dest, tag, comm, request);
	_processTimer.stopTimer();
	return result;
}

int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
					MPI_Comm comm, MPI_Status *status) {
	_processTimer.startTimer();
	int result = PMPI_Recv(buf, count, datatype, source, tag, comm, status);
	_processTimer.stopTimer();
	return result;
}

int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
					 MPI_Comm comm, MPI_Request *request) {
	_processTimer.startTimer();
	int result = PMPI_Irecv(buf, count, datatype, source, tag, comm, request);
	_processTimer.stopTimer();
	return result;
}

int MPI_File_open(MPI_Comm comm, const char *filename, int amode, MPI_Info info, MPI_File *fh) {
	_processTimer.startTimer();
	int result = PMPI_File_open(comm, filename, amode, info, fh);
	_processTimer.stopTimer();
	return result;
}


int MPI_File_close(MPI_File * fh){
	_processTimer.startTimer();
	int result = PMPI_File_close(fh);
	_processTimer.stopTimer();
	return result;
}


int MPI_File_write(MPI_File mpi_fh, const void *buf, int count, MPI_Datatype datatype, MPI_Status *status) {
	_processTimer.startTimer();
	int result = PMPI_File_write(mpi_fh, buf, count, datatype, status);
	_processTimer.stopTimer();
	return result;
}

int MPI_File_write_at(MPI_File fh, MPI_Offset offset, const void * buf, int count,
							 MPI_Datatype datatype, MPI_Status * status) {
	_processTimer.startTimer();
	int result = PMPI_File_write_at(fh, offset, buf, count, datatype, status);
	_processTimer.stopTimer();
	return result;
}

int MPI_File_get_position(MPI_File fh, MPI_Offset * offset) {
	_processTimer.startTimer();
	int result = PMPI_File_get_position(fh, offset);
	_processTimer.stopTimer();
	return result;
}

int MPI_Exscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
	_processTimer.startTimer();
	int result = PMPI_Exscan(sendbuf, recvbuf, count, datatype, op, comm);
	_processTimer.stopTimer();
	return result;
}

int MPI_Allreduce(const void* send_data, void* recv_data, int count, MPI_Datatype datatype, MPI_Op op,
				  MPI_Comm communicator) {
	_processTimer.startTimer();
	int result = PMPI_Allreduce(send_data, recv_data, count, datatype, op, communicator);
	_processTimer.stopTimer();
	return result;
}

int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype,
					 int root, MPI_Comm comm) {
	_processTimer.startTimer();
	int result = PMPI_Bcast(buffer, count, datatype, root, comm);
	_processTimer.stopTimer();
	return result;
}

int MPI_Probe(int source, int tag, MPI_Comm comm, MPI_Status * status) {
	_processTimer.startTimer();
	int result = PMPI_Probe(source, tag, comm, status);
	_processTimer.stopTimer();
	return result;
}

int MPI_Iprobe(int source, int tag, MPI_Comm comm, int *flag, MPI_Status * status) {
	_processTimer.startTimer();
	int result = PMPI_Iprobe(source, tag, comm, flag, status);
	_processTimer.stopTimer();
	return result;
}


int MPI_Test(MPI_Request * request, int *flag, MPI_Status * status) {
	_processTimer.startTimer();
	int result = PMPI_Test(request, flag, status);
	_processTimer.stopTimer();
	return result;
}

int MPI_Type_create_struct(int count, const int array_of_blocklengths[], const MPI_Aint array_of_displacements[],
								  const MPI_Datatype array_of_types[], MPI_Datatype * newtype) {
	_processTimer.startTimer();
	int result = PMPI_Type_create_struct(count, array_of_blocklengths, array_of_displacements, array_of_types, newtype);
	_processTimer.stopTimer();
	return result;
}

int MPI_Type_commit(MPI_Datatype * datatype) {
	_processTimer.startTimer();
	int result = PMPI_Type_commit(datatype);
	_processTimer.stopTimer();
	return result;
}

int MPI_Type_free(MPI_Datatype *datatype) {
	_processTimer.startTimer();
	int result = PMPI_Type_free(datatype);
	_processTimer.stopTimer();
	return result;
}

int MPI_Cart_create(MPI_Comm comm_old, int ndims, const int dims[], const int periods[], int reorder,
						   MPI_Comm * comm_cart){
	_processTimer.startTimer();
	int result = PMPI_Cart_create(comm_old, ndims, dims, periods, reorder, comm_cart);
	_processTimer.stopTimer();
	return result;
}

int MPI_Cart_coords(MPI_Comm comm, int rank, int maxdims, int coords[]) {
	_processTimer.startTimer();
	int result = PMPI_Cart_coords(comm, rank, maxdims, coords);
	_processTimer.stopTimer();
	return result;
}

int MPI_Cart_rank(MPI_Comm comm, const int coords[], int *rank) {
	_processTimer.startTimer();
	int result = PMPI_Cart_rank(comm, coords, rank);
	_processTimer.stopTimer();
	return result;
}

int MPI_Op_free(MPI_Op * op) {
	_processTimer.startTimer();
	int result = PMPI_Op_free(op);
	_processTimer.stopTimer();
	return result;
}

int MPI_Op_create(MPI_User_function *user_fn, int commute, MPI_Op *op) {
	_processTimer.startTimer();
	int result = PMPI_Op_create(user_fn, commute, op);
	_processTimer.stopTimer();
	return result;
}

int MPI_Get_count(const MPI_Status * status, MPI_Datatype datatype, int *count) {
	_processTimer.startTimer();
	int result = PMPI_Get_count(status, datatype, count);
	_processTimer.stopTimer();
	return result;
}

int MPI_Comm_free(MPI_Comm * comm) {
	_processTimer.startTimer();
	int result = PMPI_Comm_free(comm);
	_processTimer.stopTimer();
	return result;
}

/*
// FIXME: Produces Seg-Fault
int MPI_Comm_rank(MPI_Comm comm, int *rank) {
	_processTimer.startTimer();
	int result = PMPI_Comm_rank(comm, rank);
	_processTimer.stopTimer();
	return result;
}
*/

int MPI_Get_address(const void *location, MPI_Aint * address) {
	_processTimer.startTimer();
	int result = PMPI_Get_address(location, address);
	_processTimer.stopTimer();
	return result;
}

int MPI_Dims_create(int nnodes, int ndims, int dims[]) {
	_processTimer.startTimer();
	int result = PMPI_Dims_create(nnodes, ndims, dims);
	_processTimer.stopTimer();
	return result;
}

int MPI_Finalize(void) {
	int result = PMPI_Finalize();
	// _processTimer.writeProcessTimeLog();
	return result;
}


#ifdef __cplusplus
}
#endif


#endif // MARDYN_MPI_TIMED_H
