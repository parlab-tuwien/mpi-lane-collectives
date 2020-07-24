/* 
   Copyright (C) 2020
   Jesper Larsson Träff <traff@par.tuwien.ac.at>
   Sascha Hunold <hunold@par.tuwien.ac.at>
   This file is part of the MPI Multi-Lane Collectives Library.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#define STRING_SIZE 256
#define COUNT_AR_SIZE 30

static void parse_msize_list(char* msizes, int **count_ar, int *count_nb) {
    char* msizes_tok;
    char* save_str;

    /* Parse the list of message sizes */
    if (msizes != NULL) {
        *count_nb = 0;

        *count_ar = (int*)calloc(COUNT_AR_SIZE, sizeof(int));
        msizes_tok = strtok_r(msizes, ",", &save_str);
        while (msizes_tok != NULL) {
            long msize;
            msize = atol(msizes_tok);
            (*count_ar)[*count_nb] = msize;
            (*count_nb)++;
            msizes_tok = strtok_r(NULL, ",", &save_str);
        }
    }
}


int parse_counts(int argc, char *argv[], int **count_ar, int *count_nb) {

  int ret = 0;
  int c;
  extern int opterr, optind;

  // getopt_long should not complain about unrecognized options
  // we know that this will happen
  opterr = 0;
  optind = 1;

  while (1) {
    static struct option long_options[] =
      {
        {"counts",required_argument, 0, 'c'},
        {0, 0, 0, 0}
      };

    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long(argc, argv, "c:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c) {
    case 'c': {
      parse_msize_list(optarg, count_ar, count_nb);
      break;
    }}
  }

  return ret;
}
