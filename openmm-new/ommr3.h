void setupHippoNonbondedForce (OpenMM_System* system, FILE* log) {
   OpenMM_HippoNonbondedForce* hippoForce;
   hippoForce = OpenMM_HippoNonbondedForce_create();
   OpenMM_System_addForce(system, (OpenMM_Force*) hippoForce);
   // OpenMM_Force_setForceGroup((OpenMM_Force*) hippoForce, 1);

   if (*chgtrn__.nct != *atoms__.n || *repel__.nrep != *atoms__.n || *chgpen__.ncp != *atoms__.n) {
      // FIXME: print out error
   }

   double dipoleConversion = OpenMM_NmPerAngstrom;
   double quadrupoleConversion = OpenMM_NmPerAngstrom*OpenMM_NmPerAngstrom;

   // multipole
   OpenMM_DoubleArray* dipoles = OpenMM_DoubleArray_create(3);
   OpenMM_DoubleArray* quadrupoles = OpenMM_DoubleArray_create(9);

   for (int ii = 0; ii < *atoms__.n; ++ii) {

      OpenMM_HippoNonbondedForce_ParticleAxisTypes axisType = OpenMM_HippoNonbondedForce_NoAxisType;
      int atomz = 0;
      int atomx = 0;
      int atomy = 0;
      double charge = 0.0;
      for (int jj = 0; jj < 3; ++jj)
         OpenMM_DoubleArray_set(dipoles, jj, 0.0);
      for (int jj = 0; jj < 9; ++jj)
         OpenMM_DoubleArray_set(quadrupoles, jj, 0.0);

      if (potent__.use_repuls || true) {
         // set axis type for this atom
         char* axisPtr = mpole__.polaxe + ii*8;
         if (strncasecmp (axisPtr, "Z-then-X", 8) == 0) {
            axisType = OpenMM_HippoNonbondedForce_ZThenX;
         } else if (strncasecmp(axisPtr, "Bisector", 8) == 0) {
            axisType = OpenMM_HippoNonbondedForce_Bisector;
         } else if (strncasecmp(axisPtr, "Z-Bisect", 8) == 0) {
            axisType = OpenMM_HippoNonbondedForce_ZBisect;
         } else if (strncasecmp(axisPtr, "3-Fold", 6) == 0) {
            axisType = OpenMM_HippoNonbondedForce_ThreeFold;
         } else if (strncasecmp(axisPtr, "Z-Only", 6) == 0) {
            axisType = OpenMM_HippoNonbondedForce_ZOnly;
         } else if (strncasecmp(axisPtr, "None", 4) == 0 || strncasecmp( axisPtr, "    ",   4 ) == 0 ) {
            axisType = OpenMM_HippoNonbondedForce_NoAxisType;
         } else {
            // FIXME print error message
         }
         atomz = *(mpole__.zaxis+ii)-1;
         atomx = *(mpole__.xaxis+ii)-1;
         atomy = *(mpole__.yaxis+ii)-1;

         // set permanent charge for this atom
         charge = *(mpole__.pole + ii*mpole__.maxpole);
         // set permanent dipole and quadrupole for this atom
         double* polePtr = mpole__.pole + ii*mpole__.maxpole + 1;
         for (int jj = 0; jj < 3; jj++) {
            OpenMM_DoubleArray_set (dipoles, jj, (*(polePtr))*dipoleConversion);
            polePtr++;
         }
         for (int jj = 0; jj < 9; jj++) {
            OpenMM_DoubleArray_set (quadrupoles, jj, (*(polePtr))*quadrupoleConversion);
            polePtr++;
         }  
      }
      // Charge Transfer parameters
      double ct_alpha = 1.0;
      double ct_chgval = 0.0;

      if (potent__.use_chgtrn) {
         ct_alpha = chgtrn__.dmpct[ii] / OpenMM_NmPerAngstrom; // Tinker unit: 1/angstrom
         ct_chgval = chgtrn__.chgct[ii] * 4.184; // Tinker unit: kcal/mol
      }

      // Dispersion 
      double csixval = 0.0;
      int dispclass = 1;

      if (potent__.use_disp) {
         csixval = disp__.csix[ii] * sqrt(4.184)*OpenMM_NmPerAngstrom*OpenMM_NmPerAngstrom*OpenMM_NmPerAngstrom; // Tinker unit: sqrt(ang^6*kcal/mol)
         dispclass = atomid__.class_[ii];
      }

      // Repulsion
      double sizpr = 0.0;
      double dmppr = 1.0;
      double elepr = 1;
      if (potent__.use_repuls) {
         // repulsion parameters for this atom
         sizpr = repel__.sizpr[ii] * sqrt(4.184 * 0.1); // Tinker unit: kcal/mol
         dmppr = repel__.dmppr[ii] / OpenMM_NmPerAngstrom; // Tinker unit: 1/angstrom
         elepr = repel__.elepr[ii]; // Tinker unit: electron
      }

      double pcore = 1.0;
      double pval = -1.0;
      double palpha = 1000.0;
      double polarity = 0.0;

      if (mplpot__.use_chgpen) {
         pcore = chgpen__.pcore[ii]; // Tinker unit: electron
         pval = chgpen__.pval[ii]; // Tinker unit: electron
         palpha = chgpen__.palpha[ii] / OpenMM_NmPerAngstrom; // Tinker unit: 1/angstrom
         polarity = polar__.polarity[ii]*OpenMM_NmPerAngstrom*OpenMM_NmPerAngstrom*OpenMM_NmPerAngstrom;

      }
    OpenMM_HippoNonbondedForce_addParticle(hippoForce, charge, dipoles, quadrupoles,
                                           pcore, palpha,
                                           ct_chgval, ct_alpha, csixval,
                                           sizpr, elepr, dmppr,
                                           polarity, axisType,
                                           atomz, atomx, atomy);

   }

   OpenMM_DoubleArray_destroy (dipoles);
   OpenMM_DoubleArray_destroy (quadrupoles);

   double pme_cutoffdistance = 0.001;
   double dpme_cutoffdistance = 0.001;
   if (limits__.use_ewald) {
      OpenMM_HippoNonbondedForce_setNonbondedMethod(hippoForce, OpenMM_HippoNonbondedForce_PME);

      pme_cutoffdistance = *limits__.ewaldcut*OpenMM_NmPerAngstrom;

      OpenMM_HippoNonbondedForce_setPMEParameters(hippoForce, *ewald__.aewald/OpenMM_NmPerAngstrom, 
                                                  *pme__.nfft1, *pme__.nfft2, *pme__.nfft3);

      // OpenMM_HippoNonbondedForce_setPMEOrder (hippoForce, pme__.bsorder);
      OpenMM_HippoNonbondedForce_setEwaldErrorTolerance (hippoForce, 1.0e-04);
   }

   // Dispersion PME
   if (limits__.use_dewald) {
       OpenMM_HippoNonbondedForce_setDPMEParameters(hippoForce, *ewald__.adewald/OpenMM_NmPerAngstrom, 
                                                  *pme__.ndfft1, *pme__.ndfft2, *pme__.ndfft3);

       dpme_cutoffdistance = *limits__.dewaldcut*OpenMM_NmPerAngstrom;
   }
   
   double cutoffdistance = *limits__.ewaldcut*OpenMM_NmPerAngstrom;
   double taperdistance = *limits__.mpoletaper*OpenMM_NmPerAngstrom;

   OpenMM_HippoNonbondedForce_setCutoffDistance (hippoForce, cutoffdistance);
                                   
   OpenMM_HippoNonbondedForce_setSwitchingDistance(hippoForce, taperdistance);

   OpenMM_Boolean useCorrection;
   useCorrection = OpenMM_False;
   if (dsppot__.use_dcorr)  {
    useCorrection = OpenMM_True;

 // OpenMM_NonbondedForce_setUseDispersionCorrection(hippoForce, useCorrection);
    OpenMM_NonbondedForce* dispCorrectionForce;
    dispCorrectionForce = OpenMM_NonbondedForce_create();
    OpenMM_System_addForce(system, (OpenMM_Force*) dispCorrectionForce);
    OpenMM_NonbondedForce_setUseDispersionCorrection(dispCorrectionForce, useCorrection);
   }

   int maxn13 = 3*sizes__.maxval;
   int maxn14 = 9*sizes__.maxval;
   int maxn15 = 27*sizes__.maxval;
   int maxp11 = polgrp__.maxp11;
   int nn;

   double multipoleMultipoleScale = 1.0; // the factor by which to scale the Coulomb interaction between fixed multipoles
   double dipoleMultipoleScale = 1.0;    // the factor by which to scale the Coulomb interaction between an induced dipole and a fixed multipole
   double dipoleDipoleScale = 1.0;       // the factor by which to scale the Coulomb interaction between induced dipoles
   double dispersionScale = 1.0;         // the factor by which to scale the dispersion interaction
   double repulsionScale = 1.0;          // the factor by which to scale the Pauli repulsion
   double chargeTransferScale = 1.0;     // the factor by which to scale the charge transfer interaction
   
   for (int ii = 0; ii < *atoms__.n; ++ii) {

      int ptrpol = ii * maxp11;
      bool samegrp = false;
      int testAtom;
      int ptratm;

      // covalent12 scale factors
      for (int j = 0; j < couple__.n12[ii]; ++j) {
        int k = couple__.i12[j] - 1;

        multipoleMultipoleScale = *mplpot__.m2scale;
        dipoleDipoleScale = *polpot__.w2scale;
        dispersionScale = *dsppot__.dsp2scale;
        repulsionScale = *reppot__.r2scale;
        chargeTransferScale = *mplpot__.m2scale;
        
        // check if atoms belong to the same group
        samegrp = false;
        for (int l = 0; l < polgrp__.np11[ii]; ++l) {
          testAtom = polgrp__.ip11[l + ptrpol] - 1;

          if (k == testAtom) {
            samegrp = true;
            break;  
          }
        }
        if (samegrp) {
          // intra-group
          dipoleMultipoleScale = *polpot__.p2iscale;           
        } else {
          // inter-group
          dipoleMultipoleScale = *polpot__.p2scale;
        }

        OpenMM_HippoNonbondedForce_addException(hippoForce, ii, k, 
          multipoleMultipoleScale, dipoleMultipoleScale, 
          dipoleDipoleScale, dispersionScale, 
          repulsionScale, chargeTransferScale, OpenMM_False);
      }

      // covalent13 scale factors
      
      nn = couple__.n13[ii];
      ptratm = ii * maxn13;

      for (int j = 0; j < nn; ++j) {
        int k = couple__.i13[j + ptratm] - 1;

        multipoleMultipoleScale = *mplpot__.m3scale;
        dipoleDipoleScale = *polpot__.w3scale;
        dispersionScale = *dsppot__.dsp3scale;
        repulsionScale = *reppot__.r3scale;
        chargeTransferScale = *mplpot__.m3scale;

        // check if atoms belong to the same group
        samegrp = false;
        for (int l = 0; l < polgrp__.np11[ii]; ++l) {
          testAtom = polgrp__.ip11[l + ptrpol] - 1;

          if (k == testAtom) {
            samegrp = true;
            break;  
          }
        }
        if (samegrp) {
          // intra-group
          dipoleMultipoleScale = *polpot__.p3iscale;           
        } else {
          // inter-group
          dipoleMultipoleScale = *polpot__.p3scale;
        }

        OpenMM_HippoNonbondedForce_addException(hippoForce, ii, k, 
          multipoleMultipoleScale, dipoleMultipoleScale, 
          dipoleDipoleScale, dispersionScale, 
          repulsionScale, chargeTransferScale, OpenMM_False);
      }

      // covalent14 scale factors
      samegrp = false;
      nn = couple__.n14[ii];
      ptratm = ii * maxn14;

      for (int j = 0; j < nn; ++j) {
        int k = couple__.i14[j + ptratm] - 1;

        multipoleMultipoleScale = *mplpot__.m4scale;
        dipoleDipoleScale = *polpot__.w4scale;
        dispersionScale = *dsppot__.dsp4scale;
        repulsionScale = *reppot__.r4scale;
        chargeTransferScale = *mplpot__.m4scale;

        // check if atoms belong to the same group
        samegrp = false;
        for (int l = 0; l < polgrp__.np11[ii]; ++l) {
          testAtom = polgrp__.ip11[l + ptrpol] - 1;

          if (k == testAtom) {
            samegrp = true;
            break;  
          }
        }

        if (samegrp) {
          // intra-group
          dipoleMultipoleScale = *polpot__.p4iscale;           
        } else {
          // inter-group
          dipoleMultipoleScale = *polpot__.p4scale;
        }

        OpenMM_HippoNonbondedForce_addException(hippoForce, ii, k, 
          multipoleMultipoleScale, dipoleMultipoleScale, 
          dipoleDipoleScale, dispersionScale, 
          repulsionScale, chargeTransferScale, OpenMM_False);
      }

      // covalent15 scale factors
      samegrp = false;
      nn = couple__.n15[ii];
      ptratm = ii * maxn15;

      for (int j = 0; j < nn; ++j) {
        int k = couple__.i15[j + ptratm] - 1;

        multipoleMultipoleScale = *mplpot__.m5scale;
        dispersionScale = *dsppot__.dsp5scale;
        repulsionScale = *reppot__.r5scale;
        chargeTransferScale = *mplpot__.m5scale;

        // check if atoms belong to the same group
        samegrp = false;
        for (int l = 0; l < polgrp__.np11[ii]; ++l) {
          testAtom = polgrp__.ip11[l + ptrpol] - 1;

          if (k == testAtom) {
            samegrp = true;
            break;
          }
        }

        if (samegrp) {
          // intra-group
          dipoleMultipoleScale = *polpot__.p5iscale;           
        } else {
          // inter-group
          dipoleMultipoleScale = *polpot__.p5scale;
        }

        OpenMM_HippoNonbondedForce_addException(hippoForce, ii, k, 
          multipoleMultipoleScale, dipoleMultipoleScale, 
          dipoleDipoleScale, dispersionScale, 
          repulsionScale, chargeTransferScale, OpenMM_False);
      }
   }
      // set up PBC
   if (limits__.use_ewald) {
      setDefaultPeriodicBoxVectors (system, log);
   }
}