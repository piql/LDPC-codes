/* CHANNEL.H - Declarations regarding channels. */

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


/* TYPES OF CHANNEL */

typedef enum { BSC, AWGN, AWLN , MISC} channel_type;


/* PROCEDURES TO DO WITH CHANNELS. */

int  channel_parse (char **, int, channel_type *channel, double *channel_data);
void channel_usage (void);
