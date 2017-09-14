// Author:  Aaron D. Vose
// License: This is software is released into the public domain by the author in March 2013.

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>


#define FASTAUTIL_C
#include "fastautil.h"


// Finds the length of the query sequence seq (bytes)
int fasta_seqlen(fastainfo_t *fi, char *seq)
{
  char *p,*nl,*ns;

  // Look forward until end of sequence.
  for(p=seq,nl=NULL,ns=NULL; p < fi->equery; p++) {
    if( !nl ) {
      // If we haven't found a newline character yet, we are still
      // in the sequence id line of the fasta format for this sequence.
      if( *p == '\n' ) {
	// Record position of first newline (end of id line)
	nl = p;
      }
    } else {
      // Now, the terminator is not a newline, as some fasta files
      // have sequence lines of AAs shortened to 80 chars or so,
      // separated by newlines.  Now, the terminator should be '>'.
      if( *p == '>' ) {
	ns = p;
	break;
      }
    }
  }

  // Return the diff of start and end
  if( ns ) {
    return ((int)(ns-seq));
  } else {
    return ((int)(p-seq));
  }
}


// Returns the species name of seq
char* fasta_species(fastainfo_t *fi, char *seq)
{
  static char  buf[2048];
  int          i,j,sns,len;

  // Extract species name
  len = fasta_seqlen(fi, seq);
  for(i=1,j=0,sns=0; (i<len) && (seq[i] != '\n'); i++) {
    if( seq[i] == '[' ) {
      // Start of species name section
      sns = 1;
      j = 0;
      continue;
    }
    if( seq[i] == ']' ) {
      // End of species name section
      sns = 0;
      continue;
    }
    if( sns ) {
      buf[j++] = seq[i];
    }
  }
  if( j ) {
    // Found a name; terminate the string and return
    buf[j] = '\0';
    return buf;
  } else {
    // No species name found; insert a stub instead
    buf[0] = '\0';
    return buf;
  }
}


// Returns a pointer to the start of the next sequence
char* fasta_getseq(fastainfo_t *fi, int *len)
{
  char *c;
  int   sl;
  
  if( fi->cquery >= fi->equery ) {
    // No more queries, return NULL
    return NULL;
  } else {
    c = fi->cquery;
    // Advance to the next sequence
    sl = fasta_seqlen(fi,fi->cquery);
    fi->cquery += sl;
    // Return the current sequence
    if( len ) {
      (*len) = sl;
    }
    return c;
  }
}


// Unmaps and frees the passed fasta info object
void fasta_unmap(fastainfo_t *fi)
{
  if( fi ) {
    if( fi->queries && fi->equery ) {
      // Unmap pages
      munmap(fi,fi->equery-fi->queries);
    }
    // Free the allocated fasta info object itself
    free(fi);
  }
}


// Sets up the global FastaInfo struct for later use
fastainfo_t* fasta_map(char *fn)
{
  fastainfo_t *fi;
  struct stat  statbf;
  int          f;
  
  // Allocate a new fasta into struct
  if( !(fi=malloc(sizeof(fastainfo_t))) ) {
    fprintf(stderr,"fasta_map(): Failed to allocate a new fasta info struct.\n");
    exit(1);
  }

  // Open and memory map the input sequence file
  if( (f = open(fn,O_RDONLY)) < 0 ) {
    fprintf(stderr,"fasta_map(): Could not open() fasta file \"%s\".\n",fn);
    exit(1);
  }
  if( fstat(f, &statbf) < 0 ) {
    close(f);
    fprintf(stderr,"fasta_map(): Could not fstat() opened fasta file.\n");
    exit(1);
  }
  fi->queries = mmap(NULL,statbf.st_size,PROT_READ,MAP_PRIVATE,f,0);
  close(f);
  if( fi->queries == MAP_FAILED ) {
    fprintf(stderr,"fasta_map(): Could not mmap() opened fasta file.\n");
    exit(1);
  }
  
  // Fast file is mapped; setup pointers to iterate through the sequences
  fi->cquery = fi->queries;
  fi->equery = fi->queries+statbf.st_size;

  // Return the new mapping
  return fi;
}
