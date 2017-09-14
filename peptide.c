// Author:  Aaron D. Vose
// License: This is software is released into the public domain by the author in March 2013.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fastautil.h"
#include "peptide.h"


double *IsotopicMass=NULL;


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


// Returns the number of characters consumed from p
int parse_atom_and_count(char *p, char *atom, int *c)
{
  int n,rv;

  // Sanity check and copy the first atom
  if( Atom2Index[(int)*p] == -1 ) {
    fprintf(stderr,"parse_atom_and_count(): Unknown atom '%c' from \"%s\".\n",*p,p);
    exit(1);
  }
  *atom = *p;
  
  // Find the atom count
  *c = 0;
  rv = sscanf(p+1,"%d%n",c,&n);
  if( (rv != 1) && (rv != 2) ) {
    // No count could be read, so default to 1
    *c = 1;
    n  = 0;
  } else {
    // Sanity check the atom count
    if( !*c ) {
      fprintf(stderr,"parse_atom_and_count(): Atom count is 0 from start of \"%s\".\n",p+1);
      exit(1);
    }
  }

  // Return the number of characters "consumed" from p
  return 1+n;
}


// Returns the weight of the amino acid aa
double AminoAcidWeight(char aa)
{
  double  weight;
  char   *p,*atoms,atom;
  int     n,l;

  // Find out the atomic composition of the amino acid
  atoms = AminoAcid2Atoms[(int)aa];

  // Sum the weight for each atom in the set
  weight = 0.0;
  for(p=atoms; *p; p+=l) {
    l = parse_atom_and_count(p,&atom,&n);
    weight += IsotopicMass[Atom2Index[(int)atom]] * n;
  }

  // Return the molecular weight of the amino acid
  return weight;
}


// Returns the weight of the protein/peptide
double PeptideWeight(char *sequence, int len)
{
  double weight=0.0;
  int    i;

  // Just sum the weight of each amino acid in the sequence
  for(i=0; i<len; i++) {
    if( sequence[i] != '\n' ) {
      weight += AminoAcidWeight(sequence[i]);
    }
  }

  // Take off some weight due to the fact that the bonded amino acid
  // chain is not the same as a the sum of the individual amino acids.
  weight -= (IsotopicMass[Atom2Index[(int)'O']]+IsotopicMass[Atom2Index[(int)'H']]*2) * (len-1);
  weight += IsotopicMass[Atom2Index[(int)'O']]+IsotopicMass[Atom2Index[(int)'H']]*2;

  // Retrun the total weight of the sequence
  return weight;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


// Prints useage information and then exits
void error_usage()
{
  fprintf(stderr,
	  "usage:\n\tpeptide [-m|-a] <fastafile>\n\n\n"
	  "\tWhere:\n\n"
	  "\t-m           specifies monoisotopic weight.\n"
	  "\t-a           specifies averageisotopic weight.\n"
	  "\t<fastafile>  specifies the input data file.\n\n"
	  "\tResults are written to STDOUT.\n\n");
  exit(1);
}


int main(int argc, char **argv)
{
  fastainfo_t *fi;
  char        *sequence,*s;
  int          slen;
  double       weight;

  // Process command line args
  if( argc != 3 ) {
    error_usage();
  }
  if( !strcmp(argv[1],"-m") ) {
    // Use Mono Isotopic Mass
    IsotopicMass = MonoIsotopicMass;
  } else if( !strcmp(argv[1],"-a") ) {
    // Average Isotopic Mass
    IsotopicMass = AvgIsotopicMass;
  } else {
    error_usage();
  }

  // Open the input fasta file
  fi = fasta_map(argv[2]);

  // Process each amino acid sequence in the file
  while( (sequence=fasta_getseq(fi, &slen)) != NULL ) {
    // Print out the first "word"/ID of the fasta sequence
    for(s=sequence; ((s-sequence)<slen) && (*s != ' ') && (*s != '\t') && (*s != '\n'); s++) {
      printf("%c",*s);
    }
    printf("\t");
    // The sequence has a fasta ID line header; skip past it.
    for(s=sequence; ((s-sequence)<slen) && (*s != '\n'); s++);
    if( ((s-sequence)<(slen+1)) && (*s == '\n') ) {
      slen     -= s-sequence+1;
      sequence  = s+1;
    }
    // Compute the molecular weight of the "protein"
    weight = PeptideWeight(sequence,slen);
    printf("%lf\n",weight);
  }

  // Cleanup and return success
  fasta_unmap(fi);
  return 0;
}
