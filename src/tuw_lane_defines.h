/* 
   Copyright (C) 2020
   Jesper Larsson Träff <traff@par.tuwien.ac.at>
   Sascha Hunold <hunold@par.tuwien.ac.at>
   This file is part of the MPI Multi-Lane Collectives Library.
 */

#ifndef SRC_TUW_LANE_DEFINES_H_
#define SRC_TUW_LANE_DEFINES_H_

#define BLOCK(_c,_p) (((_c)%(_p))==0 ? (_c)/(_p) : (_c)/(_p)+1)
#define TYPE MPI_INT
typedef int base_t;
#define ROOT 0
#define OP MPI_SUM


#define LANE_TAG 1024
#define PIPE 50

#endif /* SRC_TUW_LANE_DEFINES_H_ */
