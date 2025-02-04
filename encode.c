/* ENCODE.C - Encode message blocks. */

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
#include "blockio.h"
#include "open.h"
#include "mod2sparse.h"
#include "mod2dense.h"
#include "rcode.h"
#include "enc.h"

void usage(void);


/* MAIN PROGRAM. */

int main
( int argc,
  char **argv
)
{
  char *source_file, *encoded_file;
  char *pchk_file, *gen_file;
  mod2dense *u, *v;

  FILE *srcf, *encf;
  char *sblk, *cblk, *chks;
  int i, n;

  mod2sparse *H;
  gen_matrix gm;

  /* Look at initial flag arguments. */

  blockio_flush = 0;

  while (argc>1)
  {
    if (strcmp(argv[1],"-f")==0)
    { if (blockio_flush!=0) usage();
      blockio_flush = 1;
    }
    else 
    { break;
    }

    argc -= 1;
    argv += 1;
  }

  /* Look at remaining arguments. */

  if (!(pchk_file = argv[1])
   || !(gen_file = argv[2])
   || !(source_file = argv[3])
   || !(encoded_file = argv[4])
   || argv[5])
  { usage();
  }

  if ((strcmp(pchk_file,"-")==0) 
    + (strcmp(gen_file,"-")==0) 
    + (strcmp(source_file,"-")==0) > 1)
  { fprintf(stderr,"Can't read more than one stream from standard input\n");
    exit(1);
  }

  /* Read parity check file */

  H = read_pchk(pchk_file, &gm.dim);

  if (gm.dim.N<=gm.dim.M)
  { fprintf(stderr,
 "Can't encode if number of bits (%d) not greater than number of checks (%d)\n",
      gm.dim.N,gm.dim.M);
    exit(1);
  }

  /* Read generator matrix file. */

  read_gen(gen_file,0,0, &gm);

  /* Allocate needed space. */

  if (gm.type=='d')
  { u = mod2dense_allocate(gm.dim.N-gm.dim.M,1);
    v = mod2dense_allocate(gm.dim.M,1);
  }

  else if (gm.type=='m')
  { u = mod2dense_allocate(gm.dim.M,1);
    v = mod2dense_allocate(gm.dim.M,1);
  }

  else
  { u = NULL;
    v = NULL;
  }

  /* Open source file. */

  srcf = open_file_std(source_file,"r");
  if (srcf==NULL)
  { fprintf(stderr,"Can't open source file: %s\n",source_file);
    exit(1);
  }

  /* Create encoded output file. */

  encf = open_file_std(encoded_file,"w");
  if (encf==NULL)
  { fprintf(stderr,"Can't create file for encoded data: %s\n",encoded_file);
    exit(1);
  }

  sblk = chk_alloc (gm.dim.N-gm.dim.M, sizeof *sblk);
  cblk = chk_alloc (gm.dim.N, sizeof *cblk);
  chks = chk_alloc (gm.dim.M, sizeof *chks);

  /* Encode successive blocks. */

  for (n = 0; ; n++)
  { 
    /* Read block from source file. */

    if (blockio_read(srcf,sblk,gm.dim.N-gm.dim.M)==EOF) 
    { break;
    }

    /* Compute encoded block. */

    switch (gm.type)
    { case 's':
      { sparse_encode (sblk, cblk, H, &gm);
        break;
      }
      case 'd':
      { dense_encode (sblk, cblk, u, v, &gm);
        break;
      }
      case 'm':
      { mixed_encode (sblk, cblk, u, v, H, &gm);
        break;
      }
    }

    /* Check that encoded block is a code word. */

    mod2sparse_mulvec (H, cblk, chks);

    for (i = 0; i<gm.dim.M; i++) 
    { if (chks[i]==1)
      { fprintf(stderr,"Output block %d is not a code word!  (Fails check %d)\n",n,i);
        abort(); 
      }
    }

    /* Write encoded block to encoded output file. */

    blockio_write(encf,cblk,gm.dim.N);
    if (ferror(encf))
    { break;
    }
  }

  fprintf(stderr,
    "Encoded %d blocks, source block size %d, encoded block size %d\n",n,gm.dim.N-gm.dim.M,gm.dim.N);

  if (ferror(encf) || fclose(encf)!=0)
  { fprintf(stderr,"Error writing encoded blocks to %s\n",encoded_file);
    exit(1);
  }

  return 0;
}


/* PRINT USAGE MESSAGE AND EXIT. */

void usage(void)
{ fprintf(stderr,
   "Usage:  encode [ -f ] pchk-file gen-file source-file encoded-file\n");
  exit(1);
}
