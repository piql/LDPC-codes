/* PRINT-CODE.C - Print a Low Density Parity Check code's matrices. */

/* Copyright (c) 1995-2012 by Radford M. Neal.
 *
 * Permission is granted for anyone to copy, use, modify, and distribute
 * these programs and accompanying documents for any purpose, provided
 * this copyright notice is retained and prominently displayed, and note
 * is made of any changes made to these programs.  These programs and
 * documents are distributed without any warranty, express or implied.
 * As the programs were written for research purposes only, they have not
 * been tested to the degree that would be advisable in any important
 * application.  All use of these programs is entirely at the user's own
 * risk.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "alloc.h"
#include "mod2sparse.h"
#include "mod2dense.h"
#include "mod2convert.h"
#include "rcode.h"

void usage(void);


/* MAIN PROGRAM. */

int main
( int argc,
  char **argv
)
{
  Arena arena;
  char *gen_file;
  int dprint;
  int i, j;

  gen_matrix gm;

  dprint = 0;
  if (argc>1 && strcmp(argv[1],"-d")==0)
  { dprint = 1;
    argc -= 1;
    argv += 1;
  }

  if (!(gen_file = argv[1]) || argv[2])
  { usage();
  }

  arena.size = 16 * 1024 * 1024;
  arena.base = malloc(arena.size);
  arena.used = 0;

  read_gen(&arena, gen_file,0,1, &gm);

  switch (gm.type)
  {
    case 's': 
    { 
      printf("\nGenerator matrix in %s (sparse representation):\n\n",gen_file);

      printf("Column order (message bits at end):\n");
      for (j = 0; j<gm.dim.N; j++) 
      { if (j%20==0) printf("\n");
        printf(" %3d",gm.cols[j]);
      }
      printf("\n\n");

      printf("Row order:\n");
      for (i = 0; i<gm.dim.M; i++) 
      { if (i%20==0) printf("\n");
        printf(" %3d",gm.data.sparse.rows[i]); 
      }
      printf("\n\n");

      if (dprint)
      { mod2dense *Ld, *Ud;
        Ld = mod2dense_allocate(&arena, gm.dim.M,gm.dim.M);
        Ud = mod2dense_allocate(&arena, gm.dim.M,gm.dim.N);
        mod2sparse_to_dense(gm.data.sparse.L,Ld);
        mod2sparse_to_dense(gm.data.sparse.U,Ud);
        printf("L:\n\n");
        mod2dense_print(stdout,Ld);
        printf("\n");
        printf("U:\n\n");
        mod2dense_print(stdout,Ud);
        printf("\n");
      }
      else
      { printf("L:\n\n");
        mod2sparse_print(stdout,gm.data.sparse.L);
        printf("\n");
        printf("U:\n\n");
        mod2sparse_print(stdout,gm.data.sparse.U);
        printf("\n");
      }
     
      break;
    }

    case 'd': case 'm':
    {
      if (gm.type=='d')
      { printf("\nGenerator matrix in %s (dense representation):\n\n",gen_file);
      }
      if (gm.type=='m')
      { printf("\nGenerator matrix in %s (mixed representation):\n\n",gen_file);
      }

      printf("Column order (message bits at end):\n");
      for (j = 0; j<gm.dim.N; j++) 
      { if (j%20==0) printf("\n");
        printf(" %3d",gm.cols[j]);
      }
      printf("\n\n");

      printf (gm.type=='d' ? "Inv(A) X B:\n\n" : "Inv(A):\n\n");
      mod2dense_print(stdout,gm.data.G);
      printf("\n");

      break;
    }

    default: 
    { fprintf(stderr,"Unknown type of generator matrix file\n");
      exit(1);
    }
  }

  return 0;
}


/* PRINT USAGE MESSAGE AND EXIT. */

void usage(void)
{ fprintf(stderr,"Usage: print-gen [ -d ] gen-file\n");
  exit(1);
}
