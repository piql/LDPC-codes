/* ALLOC.C - Routine to allocate memory and complain if it doesn't work. */

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

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>

#include "alloc.h"

inline void *chk_alloc_old
( unsigned n,		/* Number of elements */
  unsigned size		/* Size of each element */
)
{ 
  void *p;

  p = calloc(n,size);

  if (p==0)
  { fprintf(stderr,"Ran out of memory (while trying to allocate %d bytes)\n",
      n*size);
    exit(1);
  }

  return p;
}

/* ALLOCATE SPACE AND CHECK FOR ERROR. Displays an error message and exits if the space couldn't be found. */

void *chk_alloc
( Arena *arena,		/* Allocator */
  unsigned n,		/* Number of elements */
  unsigned size		/* Size of each element */
)
{ 
  uint32_t align = size;
  align--;
  align |= align >> 1;
  align |= align >> 2;
  align |= align >> 4;
  align++;
  align = 16 < align ? 16 : align;
  uintptr_t start;
  uintptr_t p;
  size_t allocation_size;
  size_t offset;
  start = (uintptr_t)arena->base + (uintptr_t)arena->used;
  p = start;
  allocation_size = (size_t)n * (size_t)size;
  size_t mask = align - 1;
  p = (p + mask) & ~mask;

  offset = p - start + allocation_size;
  if (arena->used + offset > arena->size) {
    fprintf(stderr,"Ran out of memory (while trying to allocate %zu bytes)\n", allocation_size);
    free(arena->base);
    exit(1);
  }
  arena->used += offset;
  // printf("size: %d, n: %d\n", size, n);
  memset((void *)p, 0, allocation_size);
  return (void *)p;
}
