/* 
   Copyright (C) 2020
   Jesper Larsson Träff <traff@par.tuwien.ac.at>
   Sascha Hunold <hunold@par.tuwien.ac.at>
   This file is part of the MPI Multi-Lane Collectives Library.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "tuw_lane_utils.h"
#include "tuw_lanecoll.h"

void setup_send_recv_buffers(base_t **sendbuf, int num_stype, base_t **recvbuf, int num_rtype, int alloc_both) {

  *sendbuf = (base_t*) malloc(num_stype * sizeof(base_t));
  assert(*sendbuf!=NULL);
  if (recvbuf != NULL && alloc_both == 1) {
    *recvbuf = (base_t*) malloc(num_rtype * sizeof(base_t));
    assert(*recvbuf!=NULL);
  }

}

/*
 * find maximum count number in count_ar
 */
int get_max_count(const int *count_ar, const int count_nb) {
  int max_count = 0;
  int i;

  if( count_nb > 0 ) {
    max_count = count_ar[0];
  }
  for(i=1; i<count_nb; i++) {
    if( count_ar[i] > max_count ) {
      max_count = count_ar[i];
    }
  }

  return max_count;
}


void create_argv_copy(const int argc, char *argv[], char ***argv_copy) {
  int i;
  *argv_copy = (char **)malloc(argc * sizeof(char*));
  for(i=0; i<argc; i++) {
    (*argv_copy)[i] = strdup(argv[i]);
  }
}

void free_argv_copy(const int argc, char ***argv_copy) {
  int i;
  for(i=0; i<argc; i++) {
    free((*argv_copy)[i]);
  }
  free(*argv_copy);
}
