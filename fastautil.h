// Author:  Aaron D. Vose
// License: This is software is released into the public domain by the author in March 2013.

#ifndef FASTAUTIL_H
#define FASTAUTIL_H


// Information struct about the input fasta file this
// turns out to be a handy way of interacting with it.
typedef struct st_fastainfo {
  char       *queries;      // Pointer to start of input queries
  char       *cquery;       // Current input query
  char       *equery;       // One byte past last valid byte of query
} fastainfo_t;


////////////////////////////////////////////////////////////
// Interface
////////////////////////////////////////////////////////////
#ifndef FASTAUTIL_C
extern fastainfo_t* fasta_map(char *fn);
extern void         fasta_unmap(fastainfo_t *fi);
extern char*        fasta_getseq(fastainfo_t *fi, int *len);
extern int          fasta_seqlen(fastainfo_t *fi, char *seq);
extern char*        fasta_species(fastainfo_t *fi, char *seq);
#endif


#endif
