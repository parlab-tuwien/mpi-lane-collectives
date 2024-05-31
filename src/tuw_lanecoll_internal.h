

#ifndef MPI_LANE_COLLECTIVES_SRC_TUW_LANECOLL_INTERNAL_H_H
#define MPI_LANE_COLLECTIVES_SRC_TUW_LANECOLL_INTERNAL_H_H

#include "tuw_lanecoll.h"

extern const int USEREGCOLL;

int Get_Lane_comms(MPI_Comm comm, MPI_Comm *nodecomm, MPI_Comm *lanecomm);

#define BLOCK(_c,_p) (((_c)%(_p))==0 ? (_c)/(_p) : (_c)/(_p)+1)
#define TYPE MPI_INT
typedef int base_t;
#define ROOT 0
#define OP MPI_SUM

#define LANE_TAG 1024
#define PIPE 50

#endif //MPI_LANE_COLLECTIVES_SRC_TUW_LANECOLL_INTERNAL_H_H
