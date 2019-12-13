void setDefaultPeriodicBoxVectors (OpenMM_System* system, FILE* log) {

   OpenMM_Vec3 a;
   OpenMM_Vec3 b;
   OpenMM_Vec3 c;

   if (bound__.use_bounds) {

      a.x = boxes__.lvec[0] * OpenMM_NmPerAngstrom;
      b.x = boxes__.lvec[1] * OpenMM_NmPerAngstrom;
      c.x = boxes__.lvec[2] * OpenMM_NmPerAngstrom;

      a.y = boxes__.lvec[3] * OpenMM_NmPerAngstrom;
      b.y = boxes__.lvec[4] * OpenMM_NmPerAngstrom;
      c.y = boxes__.lvec[5] * OpenMM_NmPerAngstrom;

      a.z = boxes__.lvec[6] * OpenMM_NmPerAngstrom;
      b.z = boxes__.lvec[7] * OpenMM_NmPerAngstrom;
      c.z = boxes__.lvec[8] * OpenMM_NmPerAngstrom;

      OpenMM_System_setDefaultPeriodicBoxVectors (system, &a, &b, &c);
   }
}

void printDefaultPeriodicBoxVectors (OpenMM_System* system, FILE* log) {

   OpenMM_Vec3 a;
   OpenMM_Vec3 b;
   OpenMM_Vec3 c;

   OpenMM_System_getDefaultPeriodicBoxVectors (system, &a, &b, &c);
   a.x = a.x / OpenMM_NmPerAngstrom;
   a.y = a.y / OpenMM_NmPerAngstrom;
   a.z = a.z / OpenMM_NmPerAngstrom;
   b.x = b.x / OpenMM_NmPerAngstrom;
   b.y = b.y / OpenMM_NmPerAngstrom;
   b.z = b.z / OpenMM_NmPerAngstrom;
   c.x = c.x / OpenMM_NmPerAngstrom;
   c.y = c.y / OpenMM_NmPerAngstrom;
   c.z = c.z / OpenMM_NmPerAngstrom;

   if (bound__.use_bounds) {
      (void) fprintf (log, "\n Box Vectors:  %12.4f %12.4f %12.4f",
                         a.x, a.y, a.z);
      (void) fprintf (log, "\n    (Ang)      %12.4f %12.4f %12.4f",
                         b.x, b.y, b.z);
      (void) fprintf (log, "\n               %12.4f %12.4f %12.4f\n",
                         c.x, c.y, c.z);
   }
}

void setupAmoebaChargeForce (OpenMM_System* system, FILE* log) {

   if (*chgpot__.c2scale != 0.0 || *chgpot__.c3scale != 0.0 || *chgpot__.c5scale != 1.0) return;

   setDefaultPeriodicBoxVectors (system, log);

   char buffer[128];
   OpenMM_NonbondedForce* coulombForce;
   coulombForce = OpenMM_NonbondedForce_create ();

   for (int ii = 0; ii < *atoms__.n; ++ii) {
      // set vdw params sig = 1 and eps = 0
      OpenMM_NonbondedForce_addParticle (coulombForce, charge__.pchg[ii], 1.0, 0.0);
   }

   OpenMM_BondArray* bondArray;
   bondArray = OpenMM_BondArray_create (0);
   for (int ii = 0; ii < *atoms__.n; ++ii) {
      for (int jj = 0; jj < *(couple__.n12 + ii); ++jj) {
         int atomj = *(couple__.i12 + sizes__.maxval*ii + jj) - 1;
         if (ii < atomj) {
            OpenMM_BondArray_append (bondArray, ii, atomj);
         }
      }
   }
   OpenMM_NonbondedForce_createExceptionsFromBonds (coulombForce, bondArray, *chgpot__.c4scale, 1.0);
   OpenMM_BondArray_destroy (bondArray);

   double cutoffdistance = 0.0;
   if (bound__.use_bounds) {
      if (limits__.use_ewald) {
         OpenMM_NonbondedForce_setNonbondedMethod (coulombForce, OpenMM_NonbondedForce_PME);
         cutoffdistance = *limits__.ewaldcut;
         OpenMM_NonbondedForce_setPMEParameters (coulombForce, 10*(*ewald__.aewald),
                                                 *pme__.nfft1, *pme__.nfft2, *pme__.nfft3);
      } else {
         // OpenMM_NonbondedForce_setReactionFieldDielectric (coulombForce, chgpot__.dielec);
         // OpenMM_NonbondedForce_setNonbondedMethod (coulombForce, OpenMM_NonbondedForce_CutoffPeriodic);
         // cutoffdistance = limits__.chgcut;
         // if (limits__.chgtaper < limits__.chgcut) {
         //    OpenMM_NonbondedForce_setUseSwitchingFunction (coulombForce, OpenMM_True);
         //    OpenMM_NonbondedForce_setSwitchingDistance (coulombForce, limits__.chgtaper * OpenMM_NmPerAngstrom);
         // }
         fprintf(stderr, " EXIT -- Non-PME is not supported.\n");
         exit (-1);
      }
   } else {
      // OpenMM_NonbondedForce_setNonbondedMethod (coulombForce, OpenMM_NonbondedForce_CutoffNonPeriodic);
      // cutoffdistance = limits__.chgcut;
      fprintf(stderr, " EXIT -- Nonperiodic box is not supported.\n");
      exit (-1);
   }
   cutoffdistance *= OpenMM_NmPerAngstrom;
   OpenMM_NonbondedForce_setCutoffDistance (coulombForce, cutoffdistance);
   // Turn off vdw switching
   OpenMM_NonbondedForce_setUseSwitchingFunction (coulombForce, OpenMM_False);

   OpenMM_Force_setForceGroup ((OpenMM_Force*) coulombForce, 1);
   OpenMM_System_addForce (system, (OpenMM_Force*) coulombForce);
}

void loadCovalentArray (int numberToLoad, int* valuesToLoad,
                               OpenMM_IntArray* covaletMap) {

   int ii;
   OpenMM_IntArray_resize (covaletMap, numberToLoad);
   for (ii = 0; ii < numberToLoad; ii++) {
      OpenMM_IntArray_set (covaletMap, ii, *(valuesToLoad +ii) - 1);
   }
}

void setupAmoebaMultipoleForce (OpenMM_System* system, FILE* log) {

   char buffer[128];
   char* axisPtr;
   int ii, jj;
   int errorReport;
   double* polePtr;

   OpenMM_AmoebaMultipoleForce_MultipoleAxisTypes axisType;
   OpenMM_DoubleArray* dipoles;
   OpenMM_DoubleArray* quadrupoles;

   OpenMM_IntArray* covalentMap;
   OpenMM_IntArray* gridDimensions;

   OpenMM_DoubleArray* exptCoefficients;

   double dipoleConversion = OpenMM_NmPerAngstrom;
   double quadrupoleConversion = OpenMM_NmPerAngstrom*OpenMM_NmPerAngstrom;
   double polarityConversion = OpenMM_NmPerAngstrom*OpenMM_NmPerAngstrom
                                        *OpenMM_NmPerAngstrom;
   double dampingFactorConversion = sqrt(OpenMM_NmPerAngstrom);

   OpenMM_AmoebaMultipoleForce* amoebaMultipoleForce;
   amoebaMultipoleForce = OpenMM_AmoebaMultipoleForce_create();
   OpenMM_System_addForce(system, (OpenMM_Force*) amoebaMultipoleForce);
   OpenMM_Force_setForceGroup((OpenMM_Force*) amoebaMultipoleForce, 1);

   dipoles = OpenMM_DoubleArray_create(3);
   quadrupoles = OpenMM_DoubleArray_create(9);
   errorReport = 0;
   for (ii = 0; ii < *atoms__.n; ii++) {

      axisPtr = mpole__.polaxe + ii*8;
      if (strncasecmp (axisPtr, "Z-then-X", 8) == 0) {
         axisType = OpenMM_AmoebaMultipoleForce_ZThenX;
      } else if (strncasecmp(axisPtr, "Bisector", 8) == 0) {
         axisType = OpenMM_AmoebaMultipoleForce_Bisector;
      } else if (strncasecmp(axisPtr, "Z-Bisect", 8) == 0) {
         axisType = OpenMM_AmoebaMultipoleForce_ZBisect;
      } else if (strncasecmp(axisPtr, "3-Fold", 6) == 0) {
         axisType = OpenMM_AmoebaMultipoleForce_ThreeFold;
      } else if (strncasecmp(axisPtr, "Z-Only", 6) == 0) {
         axisType = OpenMM_AmoebaMultipoleForce_ZOnly;
      } else if (strncasecmp(axisPtr, "None", 4) == 0
                   || strncasecmp( axisPtr, "    ",   4 ) == 0 ) {
         axisType = OpenMM_AmoebaMultipoleForce_NoAxisType;
      } else {
         errorReport++;
         setNullTerminator (axisPtr, 8, buffer);
         (void) fprintf (stderr,
    "setupAmoebaMultipoleForce: Axis Type=%s for Atom %7d Not Supported\n",
                 buffer, ii );
         if (errorReport > 20) {
            (void) fflush (stderr);
            exit (-1);
         }
      }

      polePtr = mpole__.pole + ii*mpole__.maxpole + 1;
      for (jj = 0; jj < 3; jj++) {
         OpenMM_DoubleArray_set (dipoles, jj,
                                   (*(polePtr))*dipoleConversion);
         polePtr++;
      }
      for (jj = 0; jj < 9; jj++) {
         OpenMM_DoubleArray_set (quadrupoles, jj,
                                   (*(polePtr))*quadrupoleConversion);
         polePtr++;
      }

      polePtr = mpole__.pole + ii*mpole__.maxpole;
      OpenMM_AmoebaMultipoleForce_addMultipole (amoebaMultipoleForce,
                                *polePtr, dipoles, quadrupoles, axisType,
                                *(mpole__.zaxis+ii)-1,
                                *(mpole__.xaxis+ii)-1,
                                *(mpole__.yaxis+ii)-1,
                                polar__.thole[ii],
                                polar__.pdamp[ii]*dampingFactorConversion,
                                polar__.polarity[ii]*polarityConversion);
   }

   if (errorReport) {
      exit (-1);
   }

   if (limits__.use_ewald) {

      double ewaldTolerance = 1.0e-04;
      OpenMM_AmoebaMultipoleForce_setNonbondedMethod (amoebaMultipoleForce,
                                   OpenMM_AmoebaMultipoleForce_PME);
      OpenMM_AmoebaMultipoleForce_setCutoffDistance (amoebaMultipoleForce,
                                   *limits__.ewaldcut*OpenMM_NmPerAngstrom);
      OpenMM_AmoebaMultipoleForce_setAEwald (amoebaMultipoleForce,
                                   *ewald__.aewald/OpenMM_NmPerAngstrom);
      // OpenMM_AmoebaMultipoleForce_setPmeBSplineOrder (amoebaMultipoleForce,
      //                           pme__.bsorder);

      gridDimensions = OpenMM_IntArray_create (3);

      OpenMM_IntArray_set (gridDimensions, 0, *pme__.nfft1);
      OpenMM_IntArray_set (gridDimensions, 1, *pme__.nfft2);
      OpenMM_IntArray_set (gridDimensions, 2, *pme__.nfft3);

      OpenMM_AmoebaMultipoleForce_setPmeGridDimensions (amoebaMultipoleForce,
                                                        gridDimensions);
      OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance
                                   (amoebaMultipoleForce, ewaldTolerance);
      OpenMM_IntArray_destroy (gridDimensions);
      setDefaultPeriodicBoxVectors (system, log);

   } else {
      OpenMM_AmoebaMultipoleForce_setNonbondedMethod (amoebaMultipoleForce,
                                   OpenMM_AmoebaMultipoleForce_NoCutoff);
   }

   if (strncasecmp (polpot__.poltyp, "DIRECT", 6) == 0) {
      OpenMM_AmoebaMultipoleForce_setPolarizationType (amoebaMultipoleForce,
                                   OpenMM_AmoebaMultipoleForce_Direct);
   } else if (strncasecmp (polpot__.poltyp, "OPT", 3) == 0) {
      if (polopt__.copt[4] != 0.0) {
         exptCoefficients = OpenMM_DoubleArray_create (5);
         OpenMM_DoubleArray_set (exptCoefficients, 0, polopt__.copt[0]);
         OpenMM_DoubleArray_set (exptCoefficients, 1, polopt__.copt[1]);
         OpenMM_DoubleArray_set (exptCoefficients, 2, polopt__.copt[2]);
         OpenMM_DoubleArray_set (exptCoefficients, 3, polopt__.copt[3]);
         OpenMM_DoubleArray_set (exptCoefficients, 4, polopt__.copt[4]);
      } else if (polopt__.copt[3] != 0.0) {
         exptCoefficients = OpenMM_DoubleArray_create (4);
         OpenMM_DoubleArray_set (exptCoefficients, 0, polopt__.copt[0]);
         OpenMM_DoubleArray_set (exptCoefficients, 1, polopt__.copt[1]);
         OpenMM_DoubleArray_set (exptCoefficients, 2, polopt__.copt[2]);
         OpenMM_DoubleArray_set (exptCoefficients, 3, polopt__.copt[3]);
      } else if (polopt__.copt[2] != 0.0) {
         exptCoefficients = OpenMM_DoubleArray_create (3);
         OpenMM_DoubleArray_set (exptCoefficients, 0, polopt__.copt[0]);
         OpenMM_DoubleArray_set (exptCoefficients, 1, polopt__.copt[1]);
         OpenMM_DoubleArray_set (exptCoefficients, 2, polopt__.copt[2]);
      } else if (polopt__.copt[1] != 0.0) {
         exptCoefficients = OpenMM_DoubleArray_create (2);
         OpenMM_DoubleArray_set (exptCoefficients, 0, polopt__.copt[0]);
         OpenMM_DoubleArray_set (exptCoefficients, 1, polopt__.copt[1]);
      }
      OpenMM_AmoebaMultipoleForce_setExtrapolationCoefficients
                                  (amoebaMultipoleForce,exptCoefficients);
      OpenMM_AmoebaMultipoleForce_setPolarizationType (amoebaMultipoleForce,
                                   OpenMM_AmoebaMultipoleForce_Extrapolated);
   } else {
      OpenMM_AmoebaMultipoleForce_setPolarizationType (amoebaMultipoleForce,
                                   OpenMM_AmoebaMultipoleForce_Mutual);
   }

   int PolType_out = OpenMM_AmoebaMultipoleForce_getPolarizationType
                                    (amoebaMultipoleForce);

   OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations
                                    (amoebaMultipoleForce, 500);
   OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon
                                    (amoebaMultipoleForce, *polpot__.poleps);

   covalentMap = OpenMM_IntArray_create (0);
   for (ii = 0; ii < *atoms__.n; ii++) {

      loadCovalentArray (*(couple__.n12 + ii),
                         (couple__.i12 + sizes__.maxval*ii), covalentMap);
      OpenMM_AmoebaMultipoleForce_setCovalentMap (amoebaMultipoleForce, ii,
                         OpenMM_AmoebaMultipoleForce_Covalent12,
                         covalentMap);

      loadCovalentArray (*(couple__.n13 + ii),
                         (couple__.i13 + 3*sizes__.maxval*ii), covalentMap);
      OpenMM_AmoebaMultipoleForce_setCovalentMap (amoebaMultipoleForce, ii,
                         OpenMM_AmoebaMultipoleForce_Covalent13,
                         covalentMap);

      loadCovalentArray (*(couple__.n14 + ii),
                         (couple__.i14 + 9*sizes__.maxval*ii), covalentMap);
      OpenMM_AmoebaMultipoleForce_setCovalentMap (amoebaMultipoleForce, ii,
                         OpenMM_AmoebaMultipoleForce_Covalent14,
                         covalentMap);

      loadCovalentArray (*(couple__.n15 + ii),
                         (couple__.i15 + 27*sizes__.maxval*ii), covalentMap);
      OpenMM_AmoebaMultipoleForce_setCovalentMap (amoebaMultipoleForce, ii,
                         OpenMM_AmoebaMultipoleForce_Covalent15,
                         covalentMap);

      loadCovalentArray (*(polgrp__.np11 + ii),
                         (polgrp__.ip11 + polgrp__.maxp11*ii), covalentMap);
      OpenMM_AmoebaMultipoleForce_setCovalentMap (amoebaMultipoleForce, ii,
                         OpenMM_AmoebaMultipoleForce_PolarizationCovalent11,
                         covalentMap);

      loadCovalentArray (*(polgrp__.np12 + ii),
                         (polgrp__.ip12 + polgrp__.maxp12*ii), covalentMap);
      OpenMM_AmoebaMultipoleForce_setCovalentMap (amoebaMultipoleForce, ii,
                         OpenMM_AmoebaMultipoleForce_PolarizationCovalent12,
                         covalentMap);

      loadCovalentArray (*(polgrp__.np13 + ii),
                         (polgrp__.ip13 + polgrp__.maxp13*ii), covalentMap);
      OpenMM_AmoebaMultipoleForce_setCovalentMap (amoebaMultipoleForce, ii,
                         OpenMM_AmoebaMultipoleForce_PolarizationCovalent13,
                         covalentMap);

      loadCovalentArray (*(polgrp__.np14 + ii),
                         (polgrp__.ip14 + polgrp__.maxp14*ii), covalentMap);
      OpenMM_AmoebaMultipoleForce_setCovalentMap( amoebaMultipoleForce, ii,
                         OpenMM_AmoebaMultipoleForce_PolarizationCovalent14,
                         covalentMap);
   }

   OpenMM_DoubleArray_destroy (dipoles);
   OpenMM_DoubleArray_destroy (quadrupoles);

   OpenMM_IntArray_destroy (covalentMap);
}

void setupAmoebaWcaDispersionForce (OpenMM_System* system, FILE* log) {

   int ii;
   double cdispTotal = 0.0;
   double epso = 0.1100;
   double epsh = 0.0135;
   double rmino = 1.7025;
   double rminh = 1.3275;
   double awater = 0.033428;
   double slevy = 1.0;
   double dispoff = 0.26;
   double shctd = 0.81;

   OpenMM_AmoebaWcaDispersionForce*  amoebaWcaDispersionForce;
   amoebaWcaDispersionForce = OpenMM_AmoebaWcaDispersionForce_create ();
   OpenMM_System_addForce (system, (OpenMM_Force*) amoebaWcaDispersionForce);

   for (ii = 0; ii < *atoms__.n; ii++) {
      cdispTotal += nonpol__.cdisp[ii];
      OpenMM_AmoebaWcaDispersionForce_addParticle (amoebaWcaDispersionForce,
              OpenMM_NmPerAngstrom*kvdws__.rad[atomid__.class_[ii]-1],
              OpenMM_KJPerKcal*kvdws__.eps[atomid__.class_[ii]-1]);
   }

   OpenMM_AmoebaWcaDispersionForce_setEpso (amoebaWcaDispersionForce,
                                            epso*OpenMM_KJPerKcal);
   OpenMM_AmoebaWcaDispersionForce_setEpsh (amoebaWcaDispersionForce,
                                            epsh*OpenMM_KJPerKcal);
   OpenMM_AmoebaWcaDispersionForce_setRmino (amoebaWcaDispersionForce,
                                             rmino*OpenMM_NmPerAngstrom);
   OpenMM_AmoebaWcaDispersionForce_setRminh (amoebaWcaDispersionForce,
                                             rminh*OpenMM_NmPerAngstrom);
   OpenMM_AmoebaWcaDispersionForce_setDispoff (amoebaWcaDispersionForce,
                                               dispoff*OpenMM_NmPerAngstrom);
   OpenMM_AmoebaWcaDispersionForce_setAwater (amoebaWcaDispersionForce,
      awater/(OpenMM_NmPerAngstrom*OpenMM_NmPerAngstrom*OpenMM_NmPerAngstrom));
   OpenMM_AmoebaWcaDispersionForce_setSlevy (amoebaWcaDispersionForce, slevy);
   OpenMM_AmoebaWcaDispersionForce_setShctd (amoebaWcaDispersionForce, shctd);
}

