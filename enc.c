/* ENC.C - Encoding procedures. */

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
#include "mod2sparse.h"
#include "mod2dense.h"
#include "rcode.h"
#include "enc.h"


/* ENCODE A BLOCK USING A SPARSE REPRESENTATION OF THE GENERATOR MATRIX. */

void sparse_encode
( void *mem,
  size_t mem_size,
  char *sblk,
  char *cblk,
  mod2sparse *H,		/* Parity check matrix */
  gen_matrix *gm		/* Generator matrix */
)
{
  Arena arena;
  int j;

  mod2entry *e;
  char *x, *y;

  arena.base = mem;
  arena.size = mem_size;
  arena.used = 0;

  x = chk_alloc (&arena, gm->dim.M, sizeof *x);
  y = chk_alloc (&arena, gm->dim.M, sizeof *y);

  /* Multiply the vector of source bits by the systematic columns of the 
     parity check matrix, giving x.  Also copy these bits to the coded block. */

  for (j = gm->dim.M; j<gm->dim.N; j++)
  { 
    cblk[gm->cols[j]] = sblk[j-gm->dim.M];

    if (sblk[j-gm->dim.M]==1)
    { for (e = mod2sparse_first_in_col(H,gm->cols[j]);
           !mod2sparse_at_end(e);
           e = mod2sparse_next_in_col(e))
      { x[mod2sparse_row(e)] ^= 1;
      }
    }
  }
 
  /* Solve Ly=x for y by forward substitution, then U(cblk)=y by backward
     substitution. */

  if (!mod2sparse_forward_sub(gm->data.sparse.L,gm->data.sparse.rows,x,y)
   || !mod2sparse_backward_sub(gm->data.sparse.U,gm->cols,y,cblk))
  { 
    abort(); /* Shouldn't occur, even if the parity check matrix has 
                redundant rows */
  }
}


/* ENCODE A BLOCK USING DENSE REPRESENTATION OF GENERATOR MATRIX. */

void dense_encode
( char *sblk,
  char *cblk,
  mod2dense *u,
  mod2dense *v,
  gen_matrix *gm /* Generator matrix */
)
{
  int j;

  /* Copy source bits to the systematic part of the coded block. */

  for (j = gm->dim.M; j<gm->dim.N; j++) 
  { cblk[gm->cols[j]] = sblk[j-gm->dim.M];
  }

  /* Multiply by Inv(A) X B to produce check bits. */

  for (j = gm->dim.M; j<gm->dim.N; j++)
  { mod2dense_set(u,j-gm->dim.M,0,sblk[j-gm->dim.M]); 
  }
  
  mod2dense_multiply(gm->data.G,u,v);

  /* Copy check bits to the right places in the coded block. */

  for (j = 0; j<gm->dim.M; j++)
  { cblk[gm->cols[j]] = mod2dense_get(v,j,0);
  }
}


/* ENCODE A BLOCK USING MIXED REPRESENTATION OF GENERATOR MATRIX. */

void mixed_encode
( char *sblk,
  char *cblk,
  mod2dense *u,
  mod2dense *v,
  mod2sparse *H, /* Parity check matrix */
  gen_matrix *gm /* Generator matrix */
)
{
  mod2entry *e;
  int j;

  /* Multiply the vector of source bits by the message bit columns of the 
     parity check matrix.  Also copy these bits to the coded block.  Take
     account of how columns have been reordered. */

  mod2dense_clear(u);

  for (j = gm->dim.M; j<gm->dim.N; j++)
  { 
    cblk[gm->cols[j]] = sblk[j-gm->dim.M];

    if (sblk[j-gm->dim.M]==1)
    { for (e = mod2sparse_first_in_col(H,gm->cols[j]);
           !mod2sparse_at_end(e);
           e = mod2sparse_next_in_col(e))
      { (void) mod2dense_flip(u,mod2sparse_row(e),0);
      }
    }
  }

  /* Multiply by Inv(A) to produce check bits. */

  mod2dense_multiply(gm->data.G,u,v);

  /* Copy check bits to the right places in the coded block. */

  for (j = 0; j<gm->dim.M; j++)
  { cblk[gm->cols[j]] = mod2dense_get(v,j,0);
  }
}
