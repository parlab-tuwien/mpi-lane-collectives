/* 
   Copyright (C) 2020
   Jesper Larsson Träff <traff@par.tuwien.ac.at>
   Sascha Hunold <hunold@par.tuwien.ac.at>
   This file is part of the MPI Multi-Lane Collectives Library.
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <getopt.h>
#include <string.h>

#include <mpi.h>

#include "tuw_lane_defines.h"
#include "tuw_lanecoll.h"
#include "tuw_lanecoll_argparse.h"
#include "tuw_lane_utils.h"


typedef enum mpi_func {
  UNKNOWN = 0,
  BCAST,
  ALLGATHER,
  ALLTOALL,
  ALLREDUCE,
  EXSCAN,
  GATHER,
  SCAN,
  SCATTER,
  REDUCE,
  REDUCE_SCATTER_BLOCK,
  FUNCTIONS_END
} mpi_func_t;


typedef struct fname2ftype {
  const char *func_name;
  const mpi_func_t func_id;
} fname2ftype_t;


fname2ftype_t function_map[] = {
    { "bcast", BCAST },
    { "allgather", ALLGATHER },
    { "allreduce", ALLREDUCE },
    { "scan", SCAN }
};

const int NUM_FUNCS = sizeof(function_map)/sizeof(fname2ftype_t);

enum coll_type {
  DEFAULT,
  LANE,
  HIERARCHICAL
};
typedef enum coll_type coll_type_t;

typedef int (*bcast_func)(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm);
typedef int (*allgather_func)(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
    MPI_Datatype recvtype, MPI_Comm comm);
typedef int (*allreduce_func)(const void *sendbuf, void *recvbuf, int count,
    MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
typedef int (*scan_func)(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
    MPI_Op op, MPI_Comm comm);

typedef struct args {
  mpi_func_t func_id;
  coll_type_t type;
} bench_arg_t;


int parse_bench_args(int argc, char *argv[], bench_arg_t *bench_params) {
  int ret = 0;
  int c;
  extern int opterr, optind;


  int function_found = 0;
  int type_found = 0;

  bench_params->func_id = UNKNOWN;
  bench_params->type = DEFAULT;

  // getopt_long should not complain about unrecognized options
  // we know that this will happen
  opterr = 0;
  optind = 1;

  while (1) {
    static struct option long_options[] =
      {
        {"mpifunc",required_argument, 0, 'm'},
        {"type", required_argument, 0, 't'},
        {0, 0, 0, 0}
      };

    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long(argc, argv, "m:t:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c) {
    case 'm':
    {
      int i;
      for(i=0; i<NUM_FUNCS; i++) {
        if( strcmp(optarg, function_map[i].func_name) == 0 ) {
          function_found = 1;
          bench_params->func_id = function_map[i].func_id;
          break;
        }
      }
    }
    break;

    case 't':
    {
      type_found = 1;
      if( strcmp(optarg, "default") == 0 ) {
        bench_params->type = DEFAULT;
      } else if(strcmp(optarg, "hier") == 0) {
        bench_params->type = HIERARCHICAL;
      } else if(strcmp(optarg, "lane") == 0) {
        bench_params->type = LANE;
      } else {
        fprintf(stderr, "unknown type (supported: default,hier,lane)\n");
        exit(1);
      }
    }
    break;

    case '?':
      /* getopt_long already printed an error message. */
      break;

    default:
      abort();
    }
  }

  if( function_found == 0 ) {
    int i;
    fprintf(stderr, "function unknown; available are: ");
    for(i=0; i<NUM_FUNCS; i++) {
      fprintf(stderr, "%s", function_map[i].func_name);
      if(i < NUM_FUNCS-1) {
        fprintf(stderr, ", ");
      }
    }
    fprintf(stderr, "\n");
    return -1;
  }

  if( type_found == 0 ) {
    fprintf(stderr, "no type (-t) given, use default\n");
    bench_params->type = DEFAULT;
  }

  return ret;
}


int bench_bcast(int argc, char *argv[], const coll_type_t coll_type, const int *count_ar, const int count_nb) {
  int size, i;
  base_t *sendbuf;
  int max_count = 0;
  bcast_func func;
  

  MPI_Comm_size(MPI_COMM_WORLD, &size);

  max_count = get_max_count(count_ar, count_nb);
  setup_send_recv_buffers(&sendbuf, max_count, NULL, 0, 0);

  
  
  

  for (i=0; i<count_nb; i++) {
    int cur_count = count_ar[i];

    switch(coll_type) {
    case LANE:
      
      func = &Bcast_lane;
      break;
    case HIERARCHICAL:
      
      func = &Bcast_hier;
      break;
    case DEFAULT:
      
      func = &MPI_Bcast;
      break;
    default:
      fprintf(stderr, "unknown type... should never happen\n");
      exit(1);
    }

    
    

    

    func(sendbuf, cur_count, TYPE, ROOT, MPI_COMM_WORLD);

    

    

    

    // print_runtime_array name=runtime_coll collective=callname count=cur_count end_time=t2 start_time=t1 type=all
    

  }

  
  

  free(sendbuf);

  return 0;
}

int bench_allgather(int argc, char *argv[], const coll_type_t coll_type, const int *count_ar, const int count_nb) {
  int size, i;
  base_t *sendbuf, *recvbuf;
  int max_count = 0;
  int block_size;
  allgather_func func;
  

  MPI_Comm_size(MPI_COMM_WORLD, &size);

  max_count = get_max_count(count_ar, count_nb);
  block_size = max_count * size;
  setup_send_recv_buffers(&sendbuf, max_count, &recvbuf, block_size, 1);

  
  
  

  for (i=0; i<count_nb; i++) {
    int cur_count = count_ar[i];

    //block_size = BLOCK(cur_count, size);

    switch(coll_type) {
    case LANE:
      
      func = &Allgather_lane;
      break;
    case HIERARCHICAL:
      
      func = &Allgather_hier;
      break;
    case DEFAULT:
      
      func = &MPI_Allgather;
      break;
    default:
      fprintf(stderr, "unknown type... should never happen\n");
      exit(1);
    }

    
    

    

    func(sendbuf, cur_count, TYPE, recvbuf, cur_count, TYPE, MPI_COMM_WORLD);

    

    

    

    // print_runtime_array name=runtime_coll collective=callname count=cur_count end_time=t2 start_time=t1 type=all
    

  }

  
  

  free(sendbuf);
  free(recvbuf);

  return 0;
}

int bench_allreduce(int argc, char *argv[], const coll_type_t coll_type, const int *count_ar, const int count_nb) {
  int size, rank, i;
  base_t *sendbuf, *recvbuf;
  int max_count = 0;
  allreduce_func func;
  

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  max_count = get_max_count(count_ar, count_nb);
  setup_send_recv_buffers(&sendbuf, max_count, &recvbuf, max_count, 1);

  
  
  

  for (i=0; i<count_nb; i++) {
    int cur_count = count_ar[i];
    int j;

    for (j=0; j<cur_count; j++) sendbuf[i] = rank+i;
    for (j=0; j<cur_count; j++) recvbuf[i] = -1;

    switch(coll_type) {
    case LANE:
      
      func = &Allreduce_lane;
      break;
    case HIERARCHICAL:
      
      func = &Allreduce_hier;
      break;
    case DEFAULT:
      
      func = &MPI_Allreduce;
      break;
    default:
      fprintf(stderr, "unknown type... should never happen\n");
      exit(1);
    }

    
    

    

    func(sendbuf, recvbuf, cur_count, TYPE, OP, MPI_COMM_WORLD);

    

    

    

    

  }

  
  

  free(sendbuf);
  free(recvbuf);

  return 0;
}

int bench_scan(int argc, char *argv[], const coll_type_t coll_type, const int *count_ar, const int count_nb) {
  int size, rank, i;
  base_t *sendbuf, *recvbuf;
  int max_count = 0;
  scan_func func;
  

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  max_count = get_max_count(count_ar, count_nb);
  setup_send_recv_buffers(&sendbuf, max_count, &recvbuf, max_count, 1);

  
  
  

  for (i=0; i<count_nb; i++) {
    int cur_count = count_ar[i];
    int j;

    for (j=0; j<cur_count; j++) sendbuf[i] = rank+i;
    for (j=0; j<cur_count; j++) recvbuf[i] = -1;

    switch(coll_type) {
    case LANE:
      
      func = &Scan_lane;
      break;
    case HIERARCHICAL:
      
      func = &Scan_hier;
      break;
    case DEFAULT:
      
      func = &MPI_Scan;
      break;
    default:
      fprintf(stderr, "unknown type... should never happen\n");
      exit(1);
    }

    
    

    

    func(sendbuf, recvbuf, cur_count, TYPE, OP, MPI_COMM_WORLD);

    

    

    

    // print_runtime_array name=runtime_coll collective=callname count=cur_count end_time=t2 start_time=t1 type=all
    

  }

  
  

  free(sendbuf);
  free(recvbuf);

  return 0;
}

int main(int argc, char *argv[]) {

  int *count_ar;
  int count_nb;
  bench_arg_t bench_args;

  char **argv_copy;

  MPI_Init(&argc, &argv);

  create_argv_copy(argc, argv, &argv_copy);
  parse_bench_args(argc, argv_copy, &bench_args);
  free_argv_copy(argc, &argv_copy);

  create_argv_copy(argc, argv, &argv_copy);
  parse_counts(argc, argv_copy, &count_ar, &count_nb);
  free_argv_copy(argc, &argv_copy);

//  printf("count_nb=%d\n", count_nb);
//  for(int i=0; i<count_nb; i++) {
//    printf("count[%d]=%d\n", i, count_ar[i]);
//  }

  switch(bench_args.func_id) {
  case BCAST:
    bench_bcast(argc, argv, bench_args.type, count_ar, count_nb);
    break;
  case ALLGATHER:
    bench_allgather(argc, argv, bench_args.type, count_ar, count_nb);
    break;
  case ALLREDUCE:
    bench_allreduce(argc, argv, bench_args.type, count_ar, count_nb);
    break;
  case SCAN:
    bench_scan(argc, argv, bench_args.type, count_ar, count_nb);
    break;
  default:
    fprintf(stderr, "func_id not supported...aborting\n");
  }

  free(count_ar);

  MPI_Finalize();


  return 0;
}
