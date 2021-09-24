/* 
   Copyright (C) 2020
   Jesper Larsson Träff <traff@par.tuwien.ac.at>
   Sascha Hunold <hunold@par.tuwien.ac.at>
   This file is part of the MPI Multi-Lane Collectives Library.
 */

#ifndef _LANECOLLECTIVES
#define _LANECOLLECTIVES

#include <mpi.h>

int Get_Lane_comms(MPI_Comm comm, MPI_Comm *nodecomm, MPI_Comm *lanecomm);

int Bcast_lane(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm);


int Gather_lane(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
    MPI_Datatype recvtype, int root, MPI_Comm comm);

int Scatter_lane(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
    MPI_Datatype recvtype, int root, MPI_Comm comm);

int Allgather_lane(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
    MPI_Datatype recvtype, MPI_Comm comm);

int Allgather_lane_zerocopy(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
                   MPI_Datatype recvtype, MPI_Comm comm);

int Reduce_lane(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);

int Allreduce_lane(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);

int Reduce_scatter_block_lane(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);

int Scan_lane(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);

int Exscan_lane(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);

int Alltoall_lane(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
    MPI_Datatype recvtype, MPI_Comm comm);

int Bcast_hier(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm);

int Gather_hier(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
    MPI_Datatype recvtype, int root, MPI_Comm comm);

int Scatter_hier(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
    MPI_Datatype recvtype, int root, MPI_Comm comm);

int Allgather_hier(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
    MPI_Datatype recvtype, MPI_Comm comm);

int Reduce_hier(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);

int Allreduce_hier(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);

int Reduce_scatter_block_hier(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);

int Scan_hier(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);

int Exscan_hier(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);

#endif
