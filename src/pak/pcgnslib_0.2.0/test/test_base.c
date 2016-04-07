//! @file test_base.c
//! @author Kyle Horne <horne.kyle@gmail.com>
//! @version 0.2
//!
//! @section LICENSE
//! BSD style license
//!
//! @section DESCRIPTION
//! Test program for pcgns library

#include "pcgnslib.h"

#include "stdio.h"
#include "stdlib.h"
#include "mpi.h"

int main(int argc, char* argv[]) {
	int err;
	int comm_size;
	int comm_rank;
	MPI_Info info;
	int fn;
	int B;
	char basename[100];
	int cell_dim = 3;
	int phys_dim = 3;

	err = MPI_Init(&argc,&argv);
	if(err!=MPI_SUCCESS) cgp_doError;
	err = MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	if(err!=MPI_SUCCESS) cgp_doError;
	err = MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
	if(err!=MPI_SUCCESS) cgp_doError;
	err = MPI_Info_create(&(info));
	if(err!=MPI_SUCCESS) cgp_doError;

	err = cgp_open("test_base.cgns", 0, MPI_COMM_WORLD, &info, &fn);
	if(err!=0) cgp_doError;
	err = cgp_base_write(fn, "Base 1", cell_dim, phys_dim, &B);
	if(err!=0) cgp_doError;
	err = cgp_base_read(fn, B, basename, &cell_dim, &phys_dim);
	if(err!=0) cgp_doError;
	err = cgp_close(fn);
	if(err!=0) cgp_doError;

	err = MPI_Finalize();
	if(err!=MPI_SUCCESS) cgp_doError;
	return err;
	}
