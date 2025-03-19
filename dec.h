/* DEC.H - Interface to decoding procedures. */

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


/* DECODING METHOD. */

typedef enum 
{ Enum_block, Enum_bit, Prprp
} decoding_method;


/* PROCEDURES RELATING TO DECODING METHODS. */

void enum_decode_setup (Arena *, gen_matrix *gm, int table, char *gen_file);
unsigned enum_decode (void *mem, size_t mem_size, double *, char *, double *, int, mod2sparse *H, gen_matrix *gm, int table, int block_no);

void prprp_decode_setup (int table);
unsigned prprp_decode 
(mod2sparse *, double *, char *, char *, double *, int table, int block_no, int max_iter);

void initprp (mod2sparse *, double *, char *, double *);
void iterprp (mod2sparse *, double *, char *, double *);
