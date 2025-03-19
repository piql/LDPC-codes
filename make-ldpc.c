/* MAKE-LDPC.C - Make a Low Density Parity Check code's parity check matrix. */

/* Copyright (c) 1995-2012 by Radford M. Neal and Peter Junteng Liu.
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
#include <math.h>

#include "rand.h"
#include "alloc.h"
#include "intio.h"
#include "open.h"
#include "mod2sparse.h"
#include "mod2dense.h"
#include "rcode.h"
#include "distrib.h"


/* METHODS FOR CONSTRUCTING CODES. */

typedef enum 
{ Evencol, 	/* Uniform number of bits per column, with number specified */
  Evenboth 	/* Uniform (as possible) over both columns and rows */
} make_method; 


void make_ldpc (Arena *, int, make_method, distrib *, int, mod2sparse **H, int M, int N);
int *column_partition (Arena *, distrib *, int);
void usage (void);


/* MAIN PROGRAM. */

int main
( int argc,
  char **argv
)
{
  Arena arena;
  make_method method;
  char *file, **meth;
  int seed, no4cycle;
  distrib *d;
  char junk;
  FILE *f;

  mod2sparse *H;
  int M;
  int N;
  
  /* Look at initial arguments. */

  (void)argc;
  if (!(file = argv[1])
   || !argv[2] || sscanf(argv[2],"%d%c",&M,&junk)!=1 || M<=0
   || !argv[3] || sscanf(argv[3],"%d%c",&N,&junk)!=1 || N<=0
   || !argv[4] || sscanf(argv[4],"%d%c",&seed,&junk)!=1)
  { usage();
  }

  /* Look at the arguments specifying the method for producing the code. */

  meth = argv+5;

  if (!meth[0]) usage();

  no4cycle = 0;

  arena.size = 16 * 1024 * 1024;
  arena.base = malloc(arena.size);
  arena.used = 0;

  if (strcmp(meth[0],"evencol")==0 || strcmp(meth[0],"evenboth")==0)
  { method = strcmp(meth[0],"evencol")==0 ? Evencol : Evenboth;
    if (!meth[1])
    { free(arena.base);
      usage();
    }
    d = distrib_create(&arena, meth[1]);
    if (d==0)
    { free(arena.base);
      usage();
    }
    if (meth[2])
    { if (strcmp(meth[2],"no4cycle")==0)
      { no4cycle = 1;
        if (meth[3])
        { free(arena.base);
          usage();
        }
      }
      else
      { free(arena.base);
        usage();
      }
    }
  }
  else
  { free(arena.base);
    usage();
  }

  /* Check for some problems. */

  if (distrib_max(d)>M)
  { fprintf(stderr,
      "At least one checks per bit (%d) is greater than total checks (%d)\n",
      distrib_max(d), M);
    free(arena.base);
    exit(1);
  }

  if (distrib_max(d)==M && N>1 && no4cycle)
  { fprintf(stderr,
      "Can't eliminate cycles of length four with this many checks per bit\n");
    free(arena.base);
    exit(1);
  } 

  /* Make the parity check matrix. */

  make_ldpc(&arena, seed,method,d,no4cycle,&H,M,N);

  /* Write out the parity check matrix. */

  f = open_file_std(file,"wb");
  if (f==NULL) 
  { fprintf(stderr,"Can't create parity check file: %s\n",file);
    free(arena.base);
    exit(1);
  }

  intio_write(f,('P'<<8)+0x80);
  
  if (ferror(f) || !mod2sparse_write(f,H) || fclose(f)!=0)
  { fprintf(stderr,"Error writing to parity check file %s\n",file);
    free(arena.base);
    exit(1);
  }

  free(arena.base);

  return 0;
}


/* PRINT USAGE MESSAGE AND EXIT. */

void usage(void)
{ fprintf(stderr,"Usage:  make-ldpc pchk-file n-checks n-bits seed method\n");
  fprintf(stderr,"Method: evencol  checks-per-col [ \"no4cycle\" ]\n");
  fprintf(stderr,"    or: evencol  checks-distribution [ \"no4cycle\" ]\n");
  fprintf(stderr,"    or: evenboth checks-per-col [ \"no4cycle\" ]\n");
  fprintf(stderr,"    or: evenboth checks-distribution [ \"no4cycle\" ]\n");
  exit(1);
}


/* CREATE A SPARSE PARITY-CHECK MATRIX.  Of size M by N, stored in H. */

void make_ldpc
( Arena *arena,
  int seed,		/* Random number seed */
  make_method method,	/* How to make it */
  distrib *d,		/* Distribution list specified */
  int no4cycle,		/* Eliminate cycles of length four? */
  mod2sparse **H, /* (output parameter) Parity check matrix */
  int M,		/* Number of rows in parity check matrix */
  int N		/* Number of columns in parity check matrix */
)
{
  mod2entry *e, *f, *g, *h;
  int added, uneven, elim4, all_even, n_full, left;
  int i, j, k, t, z, cb_N;
  int *part, *u;

  rand_seed(10*seed+1);

  *H = mod2sparse_allocate(arena, M,N);
  part = column_partition(arena, d,N);

  /* Create the initial version of the parity check matrix. */

  switch (method)
  { 
    case Evencol:
    { 
      z = 0;
      left = part[z];

      for (j = 0; j<N; j++)
      { while (left==0)
        { z += 1;
          if (z>distrib_size(d))
          { abort();
          }
          left = part[z];
        }
        for (k = 0; k<distrib_num(d,z); k++)
        { do
          { i = rand_int(M);
          } while (mod2sparse_find(*H,i,j));
          mod2sparse_insert(arena, *H,i,j);
        }
        left -= 1;
      }

      break;
    }

    case Evenboth:
    {
      cb_N = 0;
      for (z = 0; z<distrib_size(d); z++)
      { cb_N += distrib_num(d,z) * part[z];
      }
      
      u = chk_alloc (arena, cb_N, sizeof *u);

      for (k = cb_N-1; k>=0; k--)
      { u[k] = k%M;
      }
  
      uneven = 0;
      t = 0;
      z = 0;
      left = part[z];

      for (j = 0; j<N; j++)
      { 
        while (left==0)
        { z += 1;
          if (z>distrib_size(d))
          { abort();
          }
          left = part[z];
        }

	for (k = 0; k<distrib_num(d,z); k++)
        { 
          for (i = t; i<cb_N && mod2sparse_find(*H,u[i],j); i++) ;

          if (i==cb_N)
          { uneven += 1;
            do
            { i = rand_int(M);
            } while (mod2sparse_find(*H,i,j));
            mod2sparse_insert(arena, *H,i,j);
          }
          else
          { do
            { i = t + rand_int(cb_N-t);
            } while (mod2sparse_find(*H,u[i],j));
            mod2sparse_insert(arena, *H,u[i],j);
            u[i] = u[t];
            t += 1;
          }
        }

        left -= 1;
      }

      if (uneven>0)
      { fprintf(stderr,"Had to place %d checks in rows unevenly\n",uneven);
      }

      break;
    }

    default: abort();
  }

  /* Add extra bits to avoid rows with less than two checks. */

  added = 0;

  for (i = 0; i<M; i++)
  { e = mod2sparse_first_in_row(*H,i);
    if (mod2sparse_at_end(e))
    { j = rand_int(N);
      e = mod2sparse_insert(arena, *H,i,j);
      added += 1;
    }
    e = mod2sparse_first_in_row(*H,i);
    if (mod2sparse_at_end(mod2sparse_next_in_row(e)) && N>1)
    { do 
      { j = rand_int(N); 
      } while (j==mod2sparse_col(e));
      mod2sparse_insert(arena, *H,i,j);
      added += 1;
    }
  }

  if (added>0)
  { fprintf(stderr,
           "Added %d extra bit-checks to make row counts at least two\n",
           added);
  }

  /* Add extra bits to try to avoid problems with even column counts. */

  n_full = 0;
  all_even = 1;
  for (z = 0; z<distrib_size(d); z++)
  { if (distrib_num(d,z)==M) 
    { n_full += part[z];
    }
    if (distrib_num(d,z)%2==1)
    { all_even = 0;
    }
  }

  if (all_even && N-n_full>1 && added<2)
  { int a;
    for (a = 0; added+a<2; a++)
    { do
      { i = rand_int(M);
        j = rand_int(N);
      } while (mod2sparse_find(*H,i,j));
      mod2sparse_insert(arena, *H,i,j);
    }
    fprintf(stderr,
 "Added %d extra bit-checks to try to avoid problems from even column counts\n",
      a);
  }

  /* Eliminate cycles of length four, if asked, and if possible. */

  if (no4cycle)
  { 
    elim4 = 0;

    for (t = 0; t<10; t++) 
    { k = 0;
      for (j = 0; j<N; j++)
      { for (e = mod2sparse_first_in_col(*H,j);
             !mod2sparse_at_end(e);
             e = mod2sparse_next_in_col(e))
        { for (f = mod2sparse_first_in_row(*H,mod2sparse_row(e));
               !mod2sparse_at_end(f);
               f = mod2sparse_next_in_row(f))
          { if (f==e) continue;
            for (g = mod2sparse_first_in_col(*H,mod2sparse_col(f));
                 !mod2sparse_at_end(g);
                 g = mod2sparse_next_in_col(g))
            { if (g==f) continue;
              for (h = mod2sparse_first_in_row(*H,mod2sparse_row(g));
                   !mod2sparse_at_end(h);
                   h = mod2sparse_next_in_row(h))
              { if (mod2sparse_col(h)==j)
                { do
                  { i = rand_int(M);
                  } while (mod2sparse_find(*H,i,j));
                  mod2sparse_delete(*H,e);
                  mod2sparse_insert(arena, *H,i,j);
                  elim4 += 1;
                  k += 1;
                  goto nextj;
                }
              }
            }
          }
        }
      nextj: ;
      }
      if (k==0) break;
    }

    if (elim4>0)
    { fprintf(stderr,
        "Eliminated %d cycles of length four by moving checks within column\n",
         elim4);
    }

    if (t==10) 
    { fprintf(stderr,
        "Couldn't eliminate all cycles of length four in 10 passes\n");
    }
  }
}


/* PARTITION THE COLUMNS ACCORDING TO THE SPECIFIED PROPORTIONS.  It
   may not be possible to do this exactly.  Returns a pointer to an
   array of integers containing the numbers of columns corresponding 
   to the entries in the distribution passed. */

int *column_partition
( Arena *arena,
  distrib *d,		/* List of proportions and number of check-bits */
  int n			/* Total number of columns to partition */
)
{
  double *trunc;
  int *part;
  int cur, used;
  int i, j;

  part = chk_alloc (arena, distrib_size(d), sizeof(int));
  size_t snapshot = arena->used;
  trunc = chk_alloc (arena, distrib_size(d), sizeof(double));

  used = 0;
  for (i = 0; i<distrib_size(d); i++)
  { cur = floor(distrib_prop(d,i)*n);
    part[i] = cur; 
    trunc[i] = distrib_prop(d,i)*n - cur; 
    used += cur; 
  }

  if (used>n) 
  { abort();
  }
  
  while (used<n)
  { cur = 0;
    for (j = 1; j<distrib_size(d); j++) 
    { if (trunc[j]>trunc[cur])
      { cur = j;
      }
    }
    part[cur] += 1;
    used += 1;
    trunc[cur] = -1;
  }

  arena->used = snapshot; // free trunc
  return part;
}
