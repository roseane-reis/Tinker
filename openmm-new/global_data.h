#ifndef GLOBAL_DATA_H_
#define GLOBAL_DATA_H_

void setNullTerminator (char* string, int maxLength, char* buffer);

extern struct openmm_st  {
   void* ommHandle;
   char* cudaPrecision;
   char* ommPlatform;
   char* cudaDevice;
} openmm__;


extern "C" void set_openmm_data_ (void** ommHandle, char* cudaPrecision,
                      char* ommPlatform, char* cudaDevice);

extern struct sizes_st {
  int maxatm,maxtyp,maxclass,maxval,maxref,maxgrp,maxres,maxfix;
} sizes__;

extern struct units_st {
   double avogadro;
   double lightspd;
   double boltzmann;
   double gasconst;
   double elemchg;
   double vacperm;
   double emass;
   double planck;
   double joule;
   double ekcal;
   double bohr;
   double hartree;
   double evolt;
   double efreq;
   double coulomb;
   double debye;
   double prescon;
} units__;


extern "C" void set_sizes_data_ (int* maxatm, int* maxtyp, int* maxclass, int* maxval,
                      int* maxref, int* maxgrp, int* maxres, int* maxfix) ;

extern "C" void  set_units_data_ (double* avogadro, double* lightspd, double* boltzmann,
                      double* gasconst, double* elemchg, double* vacperm,
                      double* emass, double* planck, double* joule,
                      double* ekcal, double* bohr, double* hartree,
                      double* evolt, double* efreq, double* coulomb,
                      double* debye, double* prescon);

#include "global_data_decl.hh"

#endif
