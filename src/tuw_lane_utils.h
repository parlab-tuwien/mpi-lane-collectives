/* 
   Copyright (C) 2020
   Jesper Larsson Träff <traff@par.tuwien.ac.at>
   Sascha Hunold <hunold@par.tuwien.ac.at>
   This file is part of the MPI Multi-Lane Collectives Library.
 */

#ifndef SRC_TUW_LANE_UTILS_H_
#define SRC_TUW_LANE_UTILS_H_


#include "tuw_lane_defines.h"

void setup_send_recv_buffers(base_t **sendbuf, int num_stype, base_t **recvbuf, int num_rtype, int alloc_both);

int get_max_count(const int *count_ar, const int count_nb);

void create_argv_copy(const int argc, char *argv[], char ***argv_copy);

void free_argv_copy(const int argc, char ***argv_copy);

#endif /* SRC_TUW_LANE_UTILS_H_ */
