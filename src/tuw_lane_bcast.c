
#include "tuw_lanecoll.h"

int Bcast_lane(void *buffer, int count, MPI_Datatype datatype, int root,
               MPI_Comm comm)
{
  static MPI_Comm decomm   = MPI_COMM_NULL;
  static MPI_Comm nodecomm = MPI_COMM_NULL;
  static MPI_Comm lanecomm = MPI_COMM_NULL;

  int size;
  int lanerank, noderank, nodesize;
  int rootnode, noderoot;
  int i;

  MPI_Aint lb, extent;

  PMPI_Comm_size(comm,&size);
  if (count==0||size==1) return MPI_SUCCESS;

  PMPI_Type_get_extent(datatype,&lb,&extent);

  if (comm!=decomm) {
    Get_Lane_comms(comm,&nodecomm,&lanecomm);
    decomm = comm;
  }

  PMPI_Comm_rank(nodecomm,&noderank);
  PMPI_Comm_size(nodecomm,&nodesize);
  PMPI_Comm_rank(lanecomm,&lanerank);

  rootnode = root/nodesize;
  noderoot = root%nodesize;

  int block, blockcount;
  int counts[nodesize];
  int displs[nodesize];

  block = count/nodesize;

  for (i=0; i<nodesize; i++) counts[i] = block;
  counts[nodesize-1] += count%nodesize;
  displs[0] = 0;
  for (i=1; i<nodesize; i++) displs[i] = displs[i-1]+counts[i-1];
  blockcount = counts[noderank];

  if (lanerank==rootnode) {
    void *recbuf = (noderank==noderoot) ?
                   MPI_IN_PLACE : (char*)buffer+noderank*block*extent;
    if (USEREGCOLL&&count%nodesize==0) {
      PMPI_Scatter(buffer,blockcount,datatype,
                   recbuf,blockcount,datatype,noderoot,nodecomm);
    } else {
      PMPI_Scatterv(buffer,counts,displs,datatype,
                    recbuf,blockcount,datatype,noderoot,nodecomm);
    }
  }
//  PMPI_Bcast((char*)buffer+noderank*block*extent,blockcount,datatype,
//  	    rootnode,lanecomm);

  Bcast_lane_on_lane((char*)buffer+noderank*block*extent,blockcount,datatype,
                     rootnode,lanecomm);

  if (USEREGCOLL&&count%nodesize==0) {
    PMPI_Allgather(MPI_IN_PLACE,0,datatype,
                   buffer,blockcount,datatype,nodecomm);
  } else {
    PMPI_Allgatherv(MPI_IN_PLACE,0,datatype,
                    buffer,counts,displs,datatype,nodecomm);
  }

  return MPI_SUCCESS;
}

int Bcast_lane_on_lane(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)
{
  return PMPI_Bcast(buffer, count, datatype, root, comm);
}
