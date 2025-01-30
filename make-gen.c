/* MAKE-GEN.C - Make generator matrix from parity-check matrix. */

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
#include "intio.h"
#include "open.h"
#include "mod2sparse.h"
#include "mod2dense.h"
#include "mod2convert.h"
#include "rcode.h"

typedef enum { Sparse, Dense, Mixed } make_method;      /* Ways of making it */

void make_dense_mixed (FILE *, make_method, char *, mod2sparse *H, gen_matrix *gm);     /* Procs to make it */
void make_sparse (FILE *, mod2sparse_strategy, int, int, mod2sparse *H, gen_matrix *gm);
void usage(void);


/* MAIN PROGRAM. */

int main
( int argc,
  char **argv
)
{
  char *pchk_file, *gen_file, *other_gen_file;
  mod2sparse_strategy strategy;
  int abandon_when, abandon_number;
  make_method method;
  char *meth;
  char junk;
  FILE *f;

  mod2sparse *H;
  gen_matrix gm;

  strategy = Mod2sparse_first;
  other_gen_file = NULL;

  /* Look at arguments. */

  (void)argc;
  if (!(pchk_file = argv[1])
   || !(gen_file = argv[2])
   || !(meth = argv[3]))
  { usage();
  }
  
  if (strcmp(meth,"sparse")==0)     
  { method = Sparse;
    strategy = Mod2sparse_minprod;
    abandon_number = 0;
    if (argv[4])
    { if (strcmp(argv[4],"first")==0)        strategy = Mod2sparse_first;
      else if (strcmp(argv[4],"mincol")==0)  strategy = Mod2sparse_mincol;
      else if (strcmp(argv[4],"minprod")==0) strategy = Mod2sparse_minprod;
      else 
      { usage();
      }
      if (argv[5])
      { if (sscanf(argv[5],"%d%c",&abandon_number,&junk)!=1 || abandon_number<=0
         || !argv[6] 
         || sscanf(argv[6],"%d%c",&abandon_when,&junk)!=1 || abandon_when<=0
         || argv[7])
        { usage();
        }
      }
    }
  }
  else if (strcmp(meth,"dense")==0) 
  { method = Dense;
    other_gen_file = argv[4];
    if (other_gen_file && argv[5])
    { usage();
    }
  }
  else if (strcmp(meth,"mixed")==0) 
  { method = Mixed;
    other_gen_file = argv[4];
    if (other_gen_file && argv[5])
    { usage();
    }
  }
  else 
  { usage();
  }

  /* Read parity check matrix. */

  H = read_pchk(pchk_file, &gm.dim);

  if (gm.dim.N<=gm.dim.M)
  { fprintf(stderr,
     "Can't encode if number of bits (%d) isn't greater than number of checks (%d)\n",gm.dim.N,gm.dim.M);
    exit(1);
  }

  /* Create generator matrix file. */

  f = open_file_std(gen_file,"wb");
  if (f==NULL)
  { fprintf(stderr,"Can't create generator matrix file: %s\n",gen_file);
    exit(1);
  }

  /* Allocate space for row and column permutations. */

  gm.cols = chk_alloc (gm.dim.N, sizeof *gm.cols);
  gm.data.sparse.rows = chk_alloc (gm.dim.M, sizeof gm.data.sparse.rows);

  /* Create generator matrix with specified method. */

  switch (method)
  { case Sparse: 
    { make_sparse(f,strategy,abandon_number,abandon_when, H, &gm); 
      break;
    }
    case Dense: case Mixed:
    { make_dense_mixed(f,method,other_gen_file, H, &gm);
      break;
    }
    default: abort();
  }

  /* Check for error writing file. */

  if (ferror(f) || fclose(f)!=0)
  { fprintf(stderr,"Error writing to generator matrix file\n");
    exit(1);
  }

  return 0;
}


/* MAKE DENSE OR MIXED REPRESENTATION OF GENERATOR MATRIX. */

void make_dense_mixed
( FILE *f,
  make_method method,
  char *other_gen_file,
  mod2sparse *H, /* Parity check matrix */
  gen_matrix *gm /* Generator matrix */
)
{ 
  mod2dense *DH, *A, *A2, *AI, *B;
  int i, j, c, c2, n;
  int *rows_inv;

  DH = mod2dense_allocate(gm->dim.M,gm->dim.N);
  AI = mod2dense_allocate(gm->dim.M,gm->dim.M);
  B  = mod2dense_allocate(gm->dim.M,gm->dim.N-gm->dim.M);
  gm->data.G  = mod2dense_allocate(gm->dim.M,gm->dim.N-gm->dim.M);

  mod2sparse_to_dense(H,DH);

  /* If another generator matrix was specified, invert using the set of
     columns it specifies. */

  if (other_gen_file)
  { 
    read_gen(other_gen_file,1,0, gm);

    A = mod2dense_allocate(gm->dim.M,gm->dim.M);
    mod2dense_copycols(DH,A,gm->cols);

    if (!mod2dense_invert(A,AI))
    { fprintf(stderr,
       "Couldn't invert sub-matrix with column order given in other file\n");
      exit(1);
    }

    mod2dense_copycols(DH,B,gm->cols+gm->dim.M);
  }

  /* If no other generator matrix was specified, invert using whatever 
     selection of rows/columns is needed to get a non-singular sub-matrix. */

  if (!other_gen_file)
  {
    A  = mod2dense_allocate(gm->dim.M,gm->dim.N);
    A2 = mod2dense_allocate(gm->dim.M,gm->dim.N);

    n = mod2dense_invert_selected(DH,A2,gm->data.sparse.rows,gm->cols);
    mod2sparse_to_dense(H,DH);  /* DH was destroyed by invert_selected */

    if (n>0)
    { fprintf(stderr,"Note: Parity check matrix has %d redundant checks\n",n);
    }

    rows_inv = chk_alloc (gm->dim.M, sizeof *rows_inv);

    for (i = 0; i<gm->dim.M; i++)
    { rows_inv[gm->data.sparse.rows[i]] = i;
    }

    mod2dense_copyrows(A2,A,gm->data.sparse.rows);
    mod2dense_copycols(A,A2,gm->cols);
    mod2dense_copycols(A2,AI,rows_inv);

    mod2dense_copycols(DH,B,gm->cols+gm->dim.M);
  }

  /* Form final generator matrix. */

  if (method==Dense) 
  { mod2dense_multiply(AI,B,gm->data.G);
  }
  else if (method==Mixed)
  { gm->data.G = AI;
  }
  else
  { abort();
  }

  /* Compute and print number of 1s. */

  if (method==Dense)  
  { c = 0;
    for (i = 0; i<gm->dim.M; i++)
    { for (j = 0; j<gm->dim.N-gm->dim.M; j++)
      { c += mod2dense_get(gm->data.G,i,j);
      }
    }
    fprintf(stderr,
      "Number of 1s per check in Inv(A) X B is %.1f\n", (double)c/gm->dim.M);
  }

  if (method==Mixed)
  { c = 0;
    for (i = 0; i<gm->dim.M; i++)
    { for (j = 0; j<gm->dim.M; j++)
      { c += mod2dense_get(gm->data.G,i,j);
      }
    }
    c2 = 0;
    for (i = gm->dim.M; i<gm->dim.N; i++) 
    { c2 += mod2sparse_count_col(H,gm->cols[i]);
    }
    fprintf(stderr,
     "Number of 1s per check in Inv(A) is %.1f, in B is %.1f, total is %.1f\n",
     (double)c/gm->dim.M, (double)c2/gm->dim.M, (double)(c+c2)/gm->dim.M);
  }

  /* Write the represention of the generator matrix to the file. */

  intio_write(f,('G'<<8)+0x80);

  if (method==Dense)      
  { fwrite ("d", 1, 1, f);
  }
  if (method==Mixed) 
  { fwrite ("m", 1, 1, f);
  }

  intio_write(f,gm->dim.M);
  intio_write(f,gm->dim.N);

  for (i = 0; i<gm->dim.N; i++) 
  { intio_write(f,gm->cols[i]);
  }

  mod2dense_write (f, gm->data.G);
}


/* MAKE SPARSE REPRESENTATION OF GENERATOR MATRIX. */

void make_sparse
( FILE *f,
  mod2sparse_strategy strategy,
  int abandon_number,
  int abandon_when,
  mod2sparse *H, /* Parity check matrix */
  gen_matrix *gm /* Generator matrix */
)
{
  int n, cL, cU, cB;
  int i;

  /* Find LU decomposition. */

  gm->data.sparse.L = mod2sparse_allocate(gm->dim.M,gm->dim.M);
  gm->data.sparse.U = mod2sparse_allocate(gm->dim.M,gm->dim.N);

  n = mod2sparse_decomp(H,gm->dim.M,gm->data.sparse.L,gm->data.sparse.U,gm->data.sparse.rows,gm->cols,strategy,abandon_number,abandon_when);

  if (n!=0 && abandon_number==0)
  { fprintf(stderr,"Note: Parity check matrix has %d redundant checks\n",n);
  }
  if (n!=0 && abandon_number>0)
  { fprintf(stderr,
  "Note: Have %d dependent columns, but this could be due to abandonment.\n",n);
    fprintf(stderr,
  "      Try again with lower abandonment number.\n");
    exit(1);
  }

  /* Compute and print number of 1s. */

  cL = cU = cB = 0;

  for (i = 0; i<gm->dim.M; i++) cL += mod2sparse_count_row(gm->data.sparse.L,i);
  for (i = 0; i<gm->dim.M; i++) cU += mod2sparse_count_row(gm->data.sparse.U,i);
  for (i = gm->dim.M; i<gm->dim.N; i++) cB += mod2sparse_count_col(H,gm->cols[i]);

  fprintf(stderr,
   "Number of 1s per check in L is %.1f, U is %.1f, B is %.1f, total is %.1f\n",
    (double)cU/gm->dim.M, (double)cL/gm->dim.M, (double)cB/gm->dim.M, (double)(cL+cU+cB)/gm->dim.M);

  /* Write it all to the generator matrix file. */

  intio_write(f,('G'<<8)+0x80);

  fwrite ("s", 1, 1, f);

  intio_write(f,gm->dim.M);
  intio_write(f,gm->dim.N);

  for (i = 0; i<gm->dim.N; i++) 
  { intio_write(f,gm->cols[i]);
  }

  for (i = 0; i<gm->dim.M; i++) 
  { intio_write(f,gm->data.sparse.rows[i]);
  }

  mod2sparse_write (f, gm->data.sparse.L);
  mod2sparse_write (f, gm->data.sparse.U);
}


/* PRINT USAGE MESSAGE AND EXIT. */

void usage(void)
{ fprintf (stderr, 
   "Usage:  make-gen pchk-file gen-file method\n");
  fprintf (stderr, 
   "Method: sparse [ \"first\" | \"mincol\" | \"minprod\" ] [ abandon_num abandon_when ]\n");
  fprintf (stderr, 
   "    or: dense [ other-gen-file ]\n");
  fprintf (stderr, 
   "    or: mixed [ other-gen-file ]\n");
  exit(1);
}
