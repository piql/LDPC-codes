/* CHANNEL.C - Procedures and variables regarding channels. */

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
#include <string.h>

#include "channel.h"


/* PARSE A COMMAND-LINE SPECIFICATION OF A CHANNEL.  Takes a pointer to an
   argument list and an argument count; returns the number of arguments that 
   make up a channel specification at this point in the command line.  Returns
   zero if the argument list does not start with a channel specification.
   Returns -1 if there seems to be a channel specification here, but it's 
   invalid.
 */

int channel_parse
( char **argv,		/* Pointer to argument list */
  int argc,		/* Number of arguments in list */
  channel_type *channel, /* (output parameter) Type of channel */
  double *channel_data /* (output parameter) BSC: Error probability, AWGN: Noise standard deviation, AWLN: Width of noise distribution */
)
{ 
  char junk;

  if (argc==0) return 0;

  if (strcmp(argv[0],"bsc")==0  || strcmp(argv[0],"BSC")==0)
  { 
    *channel = BSC;
    if (argc<2 || sscanf(argv[1],"%lf%c",channel_data,&junk)!=1
     || *channel_data<=0 || *channel_data>=1)
    { return -1;
    }
    else
    { return 2;
    }
  }
  else if (strcmp(argv[0],"awgn")==0 || strcmp(argv[0],"AWGN")==0)
  { 
    *channel = AWGN;
    if (argc<2 || sscanf(argv[1],"%lf%c",channel_data,&junk)!=1
     || *channel_data<=0)
    { return -1;
    }
    else
    { return 2;
    }
  }
  else if (strcmp(argv[0],"awln")==0 || strcmp(argv[0],"AWLN")==0)
  {
    *channel = AWLN;
    if (argc<2 || sscanf(argv[1],"%lf%c",channel_data,&junk)!=1
     || *channel_data<=0)
    { return -1;
    }
    else
    { return 2;
    }
  }
  else if (strcmp(argv[0],"misc")==0 || strcmp(argv[0],"MISC")==0)
  {
    *channel = MISC;
    return 2;
  }
  else
  { 
    return 0;
  }
}


/* PRINT USAGE MESSAGE REGARDING CHANNEL SPECIFICATIONS. */

void channel_usage(void)
{
  fprintf(stderr,
    "Channel: bsc error-probability | awgn standard-deviation | awln width | misc 0.0\n");
}
