/* RCODE.C - Procedures to read parity check and generator matrices. */

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

#include "alloc.h"
#include "intio.h"
#include "open.h"
#include "mod2sparse.h"
#include "mod2dense.h"
#include "rcode.h"


/* READ PARITY CHECK MATRIX. If an
   error is encountered, a message is displayed on standard error, and the
   program is terminated. */

mod2sparse *read_pchk
( Arena *arena,
  char *pchk_file,
  pchk_dimensions *dim /* (out) Parity check matrix dimensions (rows, columns) */
)
{
  FILE *f;

  f = open_file_std(pchk_file,"rb");
  if (f==NULL)
  { fprintf(stderr,"Can't open parity check file: %s\n",pchk_file);
    exit(1);
  }

  if (intio_read(f)!=('P'<<8)+0x80)
  { fprintf(stderr,"File %s doesn't contain a parity check matrix\n",pchk_file);
    exit(1);
  }

  mod2sparse *H = mod2sparse_read(arena, f);

  if (H==0)
  { fprintf(stderr,"Error reading parity check matrix from %s\n",pchk_file);
    exit(1);
  }

  dim->M = mod2sparse_rows(H);
  dim->N = mod2sparse_cols(H);

  fclose(f);
  return H;
}


/* READ GENERATOR MATRIX. The generator matrix must be 
   compatible with the parity check matrix dimensions.  If the 
   second argument is 1, only the column ordering (the last N-M of which are 
   the indexes of the message bits) is read, into 'cols' .
   Otherwise, everything is read, into gm.  'type' is set to a letter
   indicating which represention is used.

   If an error is encountered, a message is displayed on standard error,
   and the program is terminated. */

void read_gen
( Arena *arena,
  char *gen_file,	/* Name of generator matrix file */
  int cols_only,	/* Read only column ordering? */
  int no_pchk_file,	/* No parity check file used? */
  gen_matrix *gm /* Generator matrix. Dimensions must be set unles no_pchk_file is 1 */
)
{
  int M2, N2;
  FILE *f;
  int i;

  f = open_file_std(gen_file,"rb");
  if (f==NULL)
  { fprintf(stderr,"Can't open generator matrix file: %s\n",gen_file);
    exit(1);
  }

  if (intio_read(f)!=('G'<<8)+0x80)
  { fprintf(stderr,"File %s doesn't contain a generator matrix\n",gen_file);
    exit(1);
  }

  if (fread (&gm->type, 1, 1, f) != 1) goto error;

  M2 = intio_read(f);
  N2 = intio_read(f);

  if (feof(f) || ferror(f)) goto error;

  if (no_pchk_file)
  { gm->dim.M = M2;
    gm->dim.N = N2;
  }
  else 
  { if (M2!=gm->dim.M || N2!=gm->dim.N)
    { fprintf(stderr,
              "Generator matrix and parity-check matrix are incompatible\n");
      exit(1);
    }
  }

  gm->cols = chk_alloc (arena, gm->dim.N, sizeof *gm->cols);
  gm->data.sparse.rows = chk_alloc (arena, gm->dim.M, sizeof *(gm->data.sparse.rows));

  for (i = 0; i<gm->dim.N; i++)
  { gm->cols[i] = intio_read(f);
    if (feof(f) || ferror(f)) goto error;
  }

  if (!cols_only)
  {
    switch (gm->type)
    {
      case 's':
      { 
        for (i = 0; i<gm->dim.M; i++)
        { gm->data.sparse.rows[i] = intio_read(f);
          if (feof(f) || ferror(f)) goto error;
        }

        if ((gm->data.sparse.L = mod2sparse_read(arena, f)) == 0) goto error;
        if ((gm->data.sparse.U = mod2sparse_read(arena, f)) == 0) goto error;
  
        if (mod2sparse_rows(gm->data.sparse.L)!=gm->dim.M || mod2sparse_cols(gm->data.sparse.L)!=gm->dim.M) goto garbled;
        if (mod2sparse_rows(gm->data.sparse.U)!=gm->dim.M || mod2sparse_cols(gm->data.sparse.U)<gm->dim.M) goto garbled;
       
        break;
      }
  
      case 'd':
      {
        if ((gm->data.G = mod2dense_read(arena, f)) == 0) goto error;
  
        if (mod2dense_rows(gm->data.G)!=gm->dim.M || mod2dense_cols(gm->data.G)!=gm->dim.N-gm->dim.M) goto garbled;
  
        break;
      }
  
      case 'm':
      {
        if ((gm->data.G = mod2dense_read(arena, f)) == 0) goto error;
  
        if (mod2dense_rows(gm->data.G)!=gm->dim.M || mod2dense_cols(gm->data.G)!=gm->dim.M) goto garbled;
  
        break;
      }
  
      default: 
      { fprintf(stderr,
         "Unknown type of generator matrix in file %s\n",gen_file);
        exit(1);
      }
    }
  }
  
  fclose(f);

  return;

error:
  fprintf(stderr,"Error reading generator matrix from file %s\n",gen_file);
  exit(1);

garbled:
  fprintf(stderr,"Garbled generator matrix in file %s\n",gen_file);
  exit(1);
}
