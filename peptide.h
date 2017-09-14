// Author:  Aaron D. Vose
// License: This is software is released into the public domain by the author in March 2013.

#ifndef PEPTIDE_H
#define PEPTIDE_H


#include "types.h"


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


// List of the amino acids and their atomic makeup
#define AA_A "C3H7NO2"
#define AA_R "C6H14N4O2"
#define AA_N "C4H8N2O3"
#define AA_D "C4H7NO4"
#define AA_C "C3H7NO2S"
#define AA_E "C5H9NO4"
#define AA_Q "C5H10N2O3"
#define AA_G "C2H5NO2"
#define AA_H "C6H9N3O2"
#define AA_I "C6H13NO2"
#define AA_L "C6H13NO2"
#define AA_K "C6H14N2O2"
#define AA_M "C5H11NO2S"
#define AA_F "C9H11NO2"
#define AA_P "C5H9NO2"
#define AA_S "C3H7NO3"
#define AA_T "C4H9NO3"
#define AA_W "C11H12N2O2"
#define AA_Y "C9H11NO3"
#define AA_V "C5H11NO2"
#define AA_U "C3H7NO2e"    // I use 'e' as short for 'Se'
#define AA_O "C12H21N3O3"

// Translates an amino-acid alphabet character into a set of atoms in string format
char *AminoAcid2Atoms[256] = 
  {
    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,   NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,   NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,   NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,   NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
    NULL, AA_A, NULL, AA_C, AA_D, AA_E, AA_F, AA_G,   AA_H, AA_I, NULL, AA_K, AA_L, AA_M, AA_N, AA_O,
    AA_P, AA_Q, AA_R, AA_S, AA_T, AA_U, AA_V, AA_W,   NULL, AA_Y, NULL, NULL, NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,   NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,   NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
    
    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,   NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,   NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,   NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,   NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,   NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,   NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,   NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,   NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
  };


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


// From: "ATOMIC WEIGHTS OF THE ELEMENTS: REVIEW 2000" (IUPAC Technical Report)
//       J. R. de LAETER et al.
//
// "monoisotopic mass": 
//  C: 12.00000, H: 1.007825, O: 15.994915, N: 14.003074, P: 30.973761, S: 31.972070, Se: 79.9165221
//
// "average mass":
//  C: 12.01078, H: 1.00794,  O: 15.99943,  N: 14.00672,  P: 30.973761, S: 32.06550,  Se: 78.96
//
// For indicies into the mass tables below:
// C='C'=0,  H='H'=1,  O='O'=2,  N='N'=3,  P='P'=4,  S='S'=5,  Se='e'=6.


// Translates an ASCII atom to an index into the isotope weight table(s).
s8b_t Atom2Index[256] = 
  {
    -1, -1, -1, -1, -1, -1, -1, -1,  -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,  -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,  -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,  -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1,  0, -1, -1, -1, -1,   1, -1, -1, -1, -1, -1,  3,  2,
     4, -1, -1,  5, -1, -1, -1, -1,  -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1,  6, -1, -1,  -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,  -1, -1, -1, -1, -1, -1, -1, -1,

    -1, -1, -1, -1, -1, -1, -1, -1,  -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,  -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,  -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,  -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,  -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,  -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,  -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,  -1, -1, -1, -1, -1, -1, -1, -1
  };


// Can be indexed into with an atom index to get the "monoisotopic mass"
double MonoIsotopicMass[] = {12.00000, 1.007825, 15.994915, 14.003074, 30.973761, 31.972070, 79.9165221};


// Can be indexed into with an atom index to get the "average isotopic mass"
double AvgIsotopicMass[] = {12.01078, 1.00794, 15.99943, 14.00672, 30.973761, 32.06550, 78.96};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


#endif
