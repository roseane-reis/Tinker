#include "global_data.h"


openmm_st openmm__;

void set_openmm_data_ (void** ommHandle, char* cudaPrecision,
                      char* ommPlatform, char* cudaDevice) {

   openmm__.ommHandle = *ommHandle;
   setNullTerminator (cudaPrecision, 6, openmm__.cudaPrecision);
   setNullTerminator (ommPlatform, 9, openmm__.ommPlatform);
   setNullTerminator (cudaDevice, 16, openmm__.cudaDevice);
}

sizes_st sizes__;

void set_sizes_data_ (int* maxatm, int* maxtyp, int* maxclass, int* maxval,
                      int* maxref, int* maxgrp, int* maxres, int* maxfix) {
   sizes__.maxatm = *maxatm;
   sizes__.maxtyp = *maxtyp;
   sizes__.maxclass = *maxclass;
   sizes__.maxval = *maxval;
   sizes__.maxref = *maxref;
   sizes__.maxgrp = *maxgrp;
   sizes__.maxres = *maxres;
   sizes__.maxfix = *maxfix;
}

units_st units__;

void set_units_data_ (double* avogadro, double* lightspd, double* boltzmann,
                      double* gasconst, double* elemchg, double* vacperm,
                      double* emass, double* planck, double* joule,
                      double* ekcal, double* bohr, double* hartree,
                      double* evolt, double* efreq, double* coulomb,
                      double* debye, double* prescon) {

   units__.avogadro = *avogadro;
   units__.lightspd = *lightspd;
   units__.boltzmann = *boltzmann;
   units__.gasconst = *gasconst;
   units__.elemchg = *elemchg;
   units__.vacperm = *vacperm;
   units__.emass = *emass;
   units__.planck = *planck;
   units__.joule = *joule;
   units__.ekcal = *ekcal;
   units__.bohr = *bohr;
   units__.hartree = *hartree;
   units__.evolt = *evolt;
   units__.efreq = *efreq;
   units__.coulomb = *coulomb;
   units__.debye = *debye;
   units__.prescon = *prescon;
}

#include "global_data_def.hh"
