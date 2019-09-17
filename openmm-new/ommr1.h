struct ConstraintMap {

   int* constraintOffset;  // offset into constraint list
   int* constraintCount;   // number of constraints for atom i
   int* constraintList;    // list of constraints sorted by atom index

   // For constraint involving atom i and j with i < j,
   //     constraintList[offset+kk] = j
   // where offset=constraintOffset[i], and 0 <= kk < constraintCount[i]
   // Note: one constraintOffset or constraintCount could be eliminated
   // since constraintCount[i] = constraintOffset[i+1] - constraintOffset[i]
};

static void freeConstraintMap (struct ConstraintMap* map) {

   free (map->constraintOffset);
   free (map->constraintCount);
   free (map->constraintList);
}

static void mapConstraints (struct ConstraintMap* map, FILE* log) {

   int ii, jj;
   int p1, p2;
   int offset, count;

   int numberOfParticles = atoms__.n;

   int* constraintCount = (int*) malloc (sizeof(int)*numberOfParticles);
   int* constraintOffset = (int*) malloc (sizeof(int)*numberOfParticles);

   memset (constraintCount, 0, sizeof(int)*numberOfParticles);
   memset (constraintOffset, 0, sizeof(int)*numberOfParticles);

   // count number of constraints for each particle where that particle
   // has the smaller index of the two constrained particles

   for (ii = 0; ii < freeze__.nrat; ii++) {
      p1 = *(freeze__.irat+2*ii) - 1;
      p2 = *(freeze__.irat + 2*ii +1) - 1;
      if (p1 > p2) {
         p1 = p2;
      }
      constraintCount[p1]++;
   }

   // set the offset value

   constraintOffset[0] = 0;
   for (ii = 1; ii < numberOfParticles; ii++){
      constraintOffset[ii] = constraintOffset[ii-1] + constraintCount[ii-1];
   }

   // allocate constraint list and load

   int* constraintList = (int*)  malloc (sizeof(int)*freeze__.nrat);
   memset (constraintCount, 0, sizeof(int)*numberOfParticles);
   for (ii = 0; ii < freeze__.nrat; ii++) {
      p1 = *(freeze__.irat+2*ii) - 1;
      p2 = *(freeze__.irat + 2*ii +1) - 1;
      if (p1 > p2) {
         int p3 = p2;
         p2 = p1;
         p1 = p3;
      }
      offset = constraintOffset[p1];
      count = constraintCount[p1];
      constraintCount[p1]++;
      constraintList[offset+count] = p2;
   }

   if (log && 0) {
      for (ii = 0; ii < numberOfParticles; ii++) {
         offset = constraintOffset[ii];
         count = constraintCount[ii];
         (void) fprintf (stderr, "%5d Offset=%5d count=%5d: ",
                         ii, offset, count);
         for (jj = 0; jj < count; jj++) {
              (void) fprintf (stderr, "%5d ", constraintList[offset+jj] );
         }
         (void) fprintf (stderr, "\n ");
      }
   }

   map->constraintCount = constraintCount;
   map->constraintOffset = constraintOffset;
   map->constraintList = constraintList;
}

static int checkForConstraint (struct ConstraintMap* map,
                               int p1, int p2, FILE* log) {

   int ii;
   int offset;
   int match = 0;

   if (p1 > p2) {
      int p3 = p2;
      p2 = p1;
      p1 = p3;
   }

   offset = map->constraintOffset[p1];
   for (ii = 0; ii < map->constraintCount[p1] && match == 0; ii++) {
      if (map->constraintList[offset+ii] == p2) {
         match = 1;
      }
   }
   return match;
}

static void setupAmoebaBondForce (OpenMM_System* system,
                                  int removeConstrainedBonds, FILE* log) {

   int ii;
   int match;
   int* bondPtr;
   double kParameterConversion;

   struct ConstraintMap map;
   if (removeConstrainedBonds) {
      mapConstraints (&map, log);
   }

   OpenMM_AmoebaBondForce* amoebaBondForce;
   amoebaBondForce = OpenMM_AmoebaBondForce_create ();

   kParameterConversion = OpenMM_KJPerKcal
                             / (OpenMM_NmPerAngstrom*OpenMM_NmPerAngstrom);
   bondPtr = bndstr__.ibnd;
   for (ii = 0; ii < bndstr__.nbond; ii++) {

      match = removeConstrainedBonds ? checkForConstraint (&map, (*bondPtr)-1,
                                               *(bondPtr+1)-1, log ) : 0;

      if (match == 0) {
         OpenMM_AmoebaBondForce_addBond (amoebaBondForce, (*bondPtr)-1,
                   *(bondPtr+1)-1, bndstr__.bl[ii]*OpenMM_NmPerAngstrom,
                   kParameterConversion*bndpot__.bndunit*bndstr__.bk[ii]);
      }
      bondPtr += 2;
   }

   OpenMM_AmoebaBondForce_setAmoebaGlobalBondCubic (amoebaBondForce,
                     bndpot__.cbnd/OpenMM_NmPerAngstrom);
   OpenMM_AmoebaBondForce_setAmoebaGlobalBondQuartic (amoebaBondForce,
                 bndpot__.qbnd/(OpenMM_NmPerAngstrom*OpenMM_NmPerAngstrom));

   // if (OpenMM_AmoebaBondForce_getNumBonds (amoebaBondForce) > 0 ) {
        OpenMM_System_addForce (system, (OpenMM_Force*) amoebaBondForce);
   // }

   if (removeConstrainedBonds) {
       freeConstraintMap (&map);
   }

}

static void setupAmoebaAngleForce (OpenMM_System* system,
                                   int removeConstrainedAngles, FILE* log) {

   int ii;
   int* angleIndexPtr;
   int match;
   char* angleTypPtr;

   // Tinker harmonic and in-plane angles are separate terms in OpenMM

   struct ConstraintMap map;
   if (removeConstrainedAngles) {
      mapConstraints (&map, log);
   }

   OpenMM_AmoebaAngleForce* amoebaAngleForce;
   amoebaAngleForce = OpenMM_AmoebaAngleForce_create ();
   OpenMM_System_addForce (system, (OpenMM_Force*) amoebaAngleForce);

   angleIndexPtr = angbnd__.iang;
   angleTypPtr = angpot__.angtyp;
   for (ii = 0; ii < angbnd__.nangle; ii++) {
      if (strncasecmp ("HARMONIC", angleTypPtr, 8) == 0 ) {
         match = removeConstrainedAngles ? checkForConstraint (&map,
                    *(angleIndexPtr)-1, *(angleIndexPtr+2)-1, log) : 0;

            match = 0;
         if (match == 0) {
            OpenMM_AmoebaAngleForce_addAngle (amoebaAngleForce,
                   *(angleIndexPtr)-1, (*(angleIndexPtr+1))-1,
                   (*(angleIndexPtr+2))-1, angbnd__.anat[ii],
                   OpenMM_KJPerKcal*angpot__.angunit*angbnd__.ak[ii]);
         }
      }
      angleIndexPtr += 4;
      angleTypPtr += 8;
   }

   OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleCubic (amoebaAngleForce,
                                                      angpot__.cang);
   OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleQuartic (amoebaAngleForce,
                                                        angpot__.qang);
   OpenMM_AmoebaAngleForce_setAmoebaGlobalAnglePentic (amoebaAngleForce,
                                                       angpot__.pang);
   OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleSextic (amoebaAngleForce,
                                                       angpot__.sang);

   if (removeConstrainedAngles) {
       freeConstraintMap (&map);
   }

}

static void setupAmoebaInPlaneAngleForce (OpenMM_System* system,
                                  int removeConstrainedBonds, FILE* log) {

   int ii;
   int* angleIndexPtr;
   char* angleTypPtr;
   int match;

   // Tinker harmonic and in-plane angles are separate terms in OpenMM

   struct ConstraintMap map;
   if (removeConstrainedBonds) {
      mapConstraints (&map, log);
   }

   OpenMM_AmoebaInPlaneAngleForce* amoebaInPlaneAngleForce;
   amoebaInPlaneAngleForce = OpenMM_AmoebaInPlaneAngleForce_create ();
   OpenMM_System_addForce (system, (OpenMM_Force*) amoebaInPlaneAngleForce);

   angleIndexPtr = angbnd__.iang;
   angleTypPtr = angpot__.angtyp;
   for (ii = 0; ii < angbnd__.nangle; ii++) {
      if (strncasecmp ("IN-PLANE", angleTypPtr, 8 ) == 0) {
         match = removeConstrainedBonds ? checkForConstraint (&map,
                     *(angleIndexPtr)-1, (*(angleIndexPtr+2))-1, log) : 0;
         if (match == 0) {
            OpenMM_AmoebaInPlaneAngleForce_addAngle (amoebaInPlaneAngleForce,
                      *(angleIndexPtr)-1, (*(angleIndexPtr+1))-1,
                      (*(angleIndexPtr+2))-1, (*(angleIndexPtr+3))-1,
                      angbnd__.anat[ii],
                      OpenMM_KJPerKcal*angpot__.angunit*angbnd__.ak[ii]);
         }
      }
      angleIndexPtr += 4;
      angleTypPtr += 8;
   }

   OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAngleCubic
                       (amoebaInPlaneAngleForce, angpot__.cang);
   OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAngleQuartic
                       (amoebaInPlaneAngleForce, angpot__.qang);
   OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAnglePentic
                       (amoebaInPlaneAngleForce, angpot__.pang);
   OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAngleSextic
                       (amoebaInPlaneAngleForce, angpot__.sang);

   if (removeConstrainedBonds) {
      freeConstraintMap (&map);
   }
}

static void setupAmoebaStretchBendForce (OpenMM_System* system,
                                 int removeConstrainedBonds, FILE* log) {

   int ii, index, abIndex, cbIndex;
   int* angleIndexPtr;
   double bondLengthAB;
   double bondLengthCB;
   double* bondLengthPtr;
   int match;

   struct ConstraintMap map;
   if (removeConstrainedBonds) {
       mapConstraints (&map, log);
   }

   OpenMM_AmoebaStretchBendForce* amoebaStretchBendForce;
   amoebaStretchBendForce = OpenMM_AmoebaStretchBendForce_create ();
   OpenMM_System_addForce (system, (OpenMM_Force*) amoebaStretchBendForce);

   for (ii = 0; ii < strbnd__.nstrbnd; ii++) {

      index = *(strbnd__.isb + 3*ii) - 1;

      abIndex = *(strbnd__.isb + 3*ii + 1);
      cbIndex = *(strbnd__.isb + 3*ii + 2);

      if (abIndex != 0) {
         bondLengthPtr = bndstr__.bl + abIndex - 1;
         bondLengthAB = (*bondLengthPtr)*OpenMM_NmPerAngstrom;
      } else {
         bondLengthAB = -1.0;
      }

      if (cbIndex != 0) {
         bondLengthPtr = bndstr__.bl + cbIndex - 1;
         bondLengthCB = (*bondLengthPtr)*OpenMM_NmPerAngstrom;
      } else {
         bondLengthCB = -1.0;
      }

      angleIndexPtr = angbnd__.iang + 4*index;

      match = removeConstrainedBonds ? checkForConstraint (&map,
                      *(angleIndexPtr)-1, (*(angleIndexPtr+2))-1, log) : 0;
      if (match == 0) {
         OpenMM_AmoebaStretchBendForce_addStretchBend
                    (amoebaStretchBendForce, *(angleIndexPtr)-1,
                    (*(angleIndexPtr+1))-1, (*(angleIndexPtr+2))-1,
                    bondLengthAB, bondLengthCB,
                    OpenMM_RadiansPerDegree*(*(angbnd__.anat +index)),
                    (OpenMM_KJPerKcal/
         OpenMM_NmPerAngstrom)*angpot__.stbnunit*(*(strbnd__.sbk+2*ii)),
                    (OpenMM_KJPerKcal/
         OpenMM_NmPerAngstrom)*angpot__.stbnunit*(*(strbnd__.sbk+2*ii+1)) );
      }
   }

   if (removeConstrainedBonds) {
      freeConstraintMap (&map);
   }
}

static void setupAmoebaUreyBradleyForce (OpenMM_System* system,
                                 int removeConstrainedBonds, FILE* log) {

   int ii;
   int* angleIndexPtr;
   double kParameterConversion;
   int match;

   struct ConstraintMap map;
   if (removeConstrainedBonds) {
      mapConstraints (&map, log);
   }

   OpenMM_HarmonicBondForce* harmonicBondForce;
   harmonicBondForce = OpenMM_HarmonicBondForce_create ();
   OpenMM_System_addForce (system, (OpenMM_Force*) harmonicBondForce);

   angleIndexPtr = urey__.iury;
   kParameterConversion = urypot__.ureyunit*OpenMM_KJPerKcal/
                              (OpenMM_NmPerAngstrom*OpenMM_NmPerAngstrom)*2;
   for (ii = 0; ii < urey__.nurey; ii++) {
      match = removeConstrainedBonds ? checkForConstraint (&map,
                      *(angleIndexPtr)-1, (*(angleIndexPtr+2))-1, log ) : 0;
      if (match == 0) {
          OpenMM_HarmonicBondForce_addBond (harmonicBondForce,
                              *(angleIndexPtr)-1, (*(angleIndexPtr+2))-1,
                              OpenMM_NmPerAngstrom*urey__.ul[ii],
                              kParameterConversion*urey__.uk[ii]);
      }
      angleIndexPtr += 3;
   }

   // OpenMM_AmoebaBondForce_setAmoebaGlobalBondCubic (amoebaBondForce,
   //                                                  urypot__.cury);
   // OpenMM_AmoebaBondForce_setAmoebaGlobalBondQuartic (amoebaBondForce,
   //                                                    urypot__.qury);

   if (removeConstrainedBonds) {
      freeConstraintMap (&map);
   }
}

static void setupAmoebaOutOfPlaneBendForce (OpenMM_System* system, FILE* log) {

   int ii, index;
   int* angleIndexPtr;

   OpenMM_AmoebaOutOfPlaneBendForce* amoebaOutOfPlaneBendForce;
   amoebaOutOfPlaneBendForce = OpenMM_AmoebaOutOfPlaneBendForce_create ();
   OpenMM_System_addForce (system, (OpenMM_Force*) amoebaOutOfPlaneBendForce);

   for (ii = 0; ii < opbend__.nopbend; ii++) {
      index = *(opbend__.iopb + ii) - 1;
      angleIndexPtr = angbnd__.iang + 4*index;
      OpenMM_AmoebaOutOfPlaneBendForce_addOutOfPlaneBend
                  (amoebaOutOfPlaneBendForce, *(angleIndexPtr)-1,
                  (*(angleIndexPtr+1))-1, (*(angleIndexPtr+2))-1,
                  (*(angleIndexPtr+3))-1,
                  OpenMM_KJPerKcal*angpot__.opbunit*(*(opbend__.opbk +ii)));
   }

   OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendCubic
                   (amoebaOutOfPlaneBendForce, angpot__.cang);
   OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendQuartic
                   (amoebaOutOfPlaneBendForce, angpot__.qang);
   OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendPentic
                   (amoebaOutOfPlaneBendForce, angpot__.pang);
   OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendSextic
                   (amoebaOutOfPlaneBendForce, angpot__.sang);
}

static void setupAmoebaImproperTorsionForce (OpenMM_System* system, FILE* log) {

   int ii;
   double torsunit;
   int* torsIndexPtr;
   double* torsPtr;

   OpenMM_DoubleArray* torsion1;
   OpenMM_DoubleArray* torsion2;
   OpenMM_DoubleArray* torsion3;

   OpenMM_PeriodicTorsionForce* amoebaTorsionForce;

   torsion1 = OpenMM_DoubleArray_create(2);
   torsion2 = OpenMM_DoubleArray_create(2);
   torsion3 = OpenMM_DoubleArray_create(2);

   amoebaTorsionForce = OpenMM_PeriodicTorsionForce_create ();
   OpenMM_System_addForce (system, (OpenMM_Force*) amoebaTorsionForce);

   torsunit = OpenMM_KJPerKcal*torpot__.itorunit;
   torsIndexPtr = imptor__.iitors;
   for (ii = 0; ii < imptor__.nitors; ii++) {

      torsPtr = imptor__.itors1 + ii*4;
      OpenMM_DoubleArray_set (torsion1, 0, torsunit*(*torsPtr));
      OpenMM_DoubleArray_set (torsion1, 1, acos((*(torsPtr+2))));

      torsPtr = imptor__.itors2 + ii*4;
      OpenMM_DoubleArray_set (torsion2, 0, torsunit*(*torsPtr));
      OpenMM_DoubleArray_set (torsion2, 1, acos((*(torsPtr+2))));

      torsPtr = imptor__.itors3 + ii*4;
      OpenMM_DoubleArray_set (torsion3, 0, torsunit*(*torsPtr));
      OpenMM_DoubleArray_set (torsion3, 1, acos((*(torsPtr+2))));

      OpenMM_PeriodicTorsionForce_addTorsion (amoebaTorsionForce,
                (*torsIndexPtr) - 1, (*(torsIndexPtr+1)) - 1,
                (*(torsIndexPtr+2)) - 1, (*(torsIndexPtr+3)) - 1, 1,
                OpenMM_DoubleArray_get (torsion1,1),
                OpenMM_DoubleArray_get (torsion1,0));

      OpenMM_PeriodicTorsionForce_addTorsion (amoebaTorsionForce,
                (*torsIndexPtr) - 1, (*(torsIndexPtr+1)) - 1,
                (*(torsIndexPtr+2)) - 1, (*(torsIndexPtr+3)) - 1, 2,
                OpenMM_DoubleArray_get(torsion2,1),
                OpenMM_DoubleArray_get(torsion2,0));

      OpenMM_PeriodicTorsionForce_addTorsion (amoebaTorsionForce,
                (*torsIndexPtr) - 1, (*(torsIndexPtr+1)) - 1,
                (*(torsIndexPtr+2)) - 1, (*(torsIndexPtr+3)) - 1, 3,
                OpenMM_DoubleArray_get(torsion3,1),
                OpenMM_DoubleArray_get(torsion3,0));

      torsIndexPtr += 4;
   }

   OpenMM_DoubleArray_destroy (torsion1);
   OpenMM_DoubleArray_destroy (torsion2);
   OpenMM_DoubleArray_destroy (torsion3);
}

static void setupAmoebaTorsionForce (OpenMM_System* system, FILE* log) {

   int ii;
   double torsunit;
   int* torsIndexPtr;
   double* torsPtr;

   OpenMM_DoubleArray* torsion1;
   OpenMM_DoubleArray* torsion2;
   OpenMM_DoubleArray* torsion3;
   OpenMM_DoubleArray* torsion4;
   OpenMM_DoubleArray* torsion5;
   OpenMM_DoubleArray* torsion6;

   OpenMM_PeriodicTorsionForce* amoebaTorsionForce;

   torsion1 = OpenMM_DoubleArray_create(2);
   torsion2 = OpenMM_DoubleArray_create(2);
   torsion3 = OpenMM_DoubleArray_create(2);
   torsion4 = OpenMM_DoubleArray_create(2);
   torsion5 = OpenMM_DoubleArray_create(2);
   torsion6 = OpenMM_DoubleArray_create(2);

   amoebaTorsionForce = OpenMM_PeriodicTorsionForce_create ();
   OpenMM_System_addForce (system, (OpenMM_Force*) amoebaTorsionForce);

   torsunit = OpenMM_KJPerKcal*torpot__.torsunit;
   torsIndexPtr = tors__.itors;
   for (ii = 0; ii < tors__.ntors; ii++) {

      torsPtr = tors__.tors1 + ii*4;
      OpenMM_DoubleArray_set (torsion1, 0, torsunit*(*torsPtr));
      OpenMM_DoubleArray_set (torsion1, 1, acos((*(torsPtr+2))));

      torsPtr = tors__.tors2 + ii*4;
      OpenMM_DoubleArray_set (torsion2, 0, torsunit*(*torsPtr));
      OpenMM_DoubleArray_set (torsion2, 1, acos((*(torsPtr+2))));

      torsPtr = tors__.tors3 + ii*4;
      OpenMM_DoubleArray_set (torsion3, 0, torsunit*(*torsPtr));
      OpenMM_DoubleArray_set (torsion3, 1, acos((*(torsPtr+2))));

      torsPtr = tors__.tors4 + ii*4;
      OpenMM_DoubleArray_set (torsion4, 0, torsunit*(*torsPtr));
      OpenMM_DoubleArray_set (torsion4, 1, acos((*(torsPtr+2))));

      torsPtr = tors__.tors5 + ii*4;
      OpenMM_DoubleArray_set (torsion5, 0, torsunit*(*torsPtr));
      OpenMM_DoubleArray_set (torsion5, 1, acos((*(torsPtr+2))));

      torsPtr = tors__.tors6 + ii*4;
      OpenMM_DoubleArray_set (torsion6, 0, torsunit*(*torsPtr));
      OpenMM_DoubleArray_set (torsion6, 1, acos((*(torsPtr+2))));

      OpenMM_PeriodicTorsionForce_addTorsion (amoebaTorsionForce,
                (*torsIndexPtr) - 1, (*(torsIndexPtr+1)) - 1,
                (*(torsIndexPtr+2)) - 1, (*(torsIndexPtr+3)) - 1, 1,
                OpenMM_DoubleArray_get (torsion1,1),
                OpenMM_DoubleArray_get (torsion1,0));

      OpenMM_PeriodicTorsionForce_addTorsion (amoebaTorsionForce,
                (*torsIndexPtr) - 1, (*(torsIndexPtr+1)) - 1,
                (*(torsIndexPtr+2)) - 1, (*(torsIndexPtr+3)) - 1, 2,
                OpenMM_DoubleArray_get(torsion2,1),
                OpenMM_DoubleArray_get(torsion2,0));

      OpenMM_PeriodicTorsionForce_addTorsion (amoebaTorsionForce,
                (*torsIndexPtr) - 1, (*(torsIndexPtr+1)) - 1,
                (*(torsIndexPtr+2)) - 1, (*(torsIndexPtr+3)) - 1, 3,
                OpenMM_DoubleArray_get(torsion3,1),
                OpenMM_DoubleArray_get(torsion3,0));

      OpenMM_PeriodicTorsionForce_addTorsion (amoebaTorsionForce,
                (*torsIndexPtr) - 1, (*(torsIndexPtr+1)) - 1,
                (*(torsIndexPtr+2)) - 1, (*(torsIndexPtr+3)) - 1, 4,
                OpenMM_DoubleArray_get(torsion4,1),
                OpenMM_DoubleArray_get(torsion4,0));

      OpenMM_PeriodicTorsionForce_addTorsion (amoebaTorsionForce,
                (*torsIndexPtr) - 1, (*(torsIndexPtr+1)) - 1,
                (*(torsIndexPtr+2)) - 1, (*(torsIndexPtr+3)) - 1, 5,
                OpenMM_DoubleArray_get(torsion5,1),
                OpenMM_DoubleArray_get(torsion5,0));

      OpenMM_PeriodicTorsionForce_addTorsion (amoebaTorsionForce,
                (*torsIndexPtr) - 1, (*(torsIndexPtr+1)) - 1,
                (*(torsIndexPtr+2)) - 1, (*(torsIndexPtr+3)) - 1, 6,
                OpenMM_DoubleArray_get(torsion6,1),
                OpenMM_DoubleArray_get(torsion6,0));

      torsIndexPtr += 4;
   }

   OpenMM_DoubleArray_destroy (torsion1);
   OpenMM_DoubleArray_destroy (torsion2);
   OpenMM_DoubleArray_destroy (torsion3);
   OpenMM_DoubleArray_destroy (torsion4);
   OpenMM_DoubleArray_destroy (torsion5);
   OpenMM_DoubleArray_destroy (torsion6);
}

static void setupAmoebaPiTorsionForce (OpenMM_System* system, FILE* log) {

   int ii;
   int* piTorsIndexPtr;

   OpenMM_AmoebaPiTorsionForce* amoebaPiTorsionForce;
   amoebaPiTorsionForce = OpenMM_AmoebaPiTorsionForce_create ();
   OpenMM_System_addForce (system, (OpenMM_Force*) amoebaPiTorsionForce);

   piTorsIndexPtr = pitors__.ipit;
   for (ii = 0; ii < pitors__.npitors; ii++) {
      OpenMM_AmoebaPiTorsionForce_addPiTorsion (amoebaPiTorsionForce,
                    (*piTorsIndexPtr) -1, (*(piTorsIndexPtr+1))-1,
                    (*(piTorsIndexPtr+2))-1, (*(piTorsIndexPtr+3))-1,
                    (*(piTorsIndexPtr+4))-1, (*(piTorsIndexPtr+5))-1,
                    OpenMM_KJPerKcal*torpot__.ptorunit*pitors__.kpit[ii]);
      piTorsIndexPtr += 6;
   }
}

static int getChiralIndex (int atomB, int atomC, int atomD) {

   int ii, j, m, k;
   int chiralAtom;

   // test for chirality at the central torsion-torsion site

   chiralAtom = -1;
   if (*(couple__.n12 + atomC) == 4) {
      j = 0;
      for (ii = 0; ii < 4; ii++) {
         m = *(couple__.i12 + sizes__.maxval*atomC + ii) - 1;
         if (m != atomB && m != atomD) {
            if (j == 0) {
               j = m;
            } else {
               k = m;
            }
         }
      }
      if (atoms__.type[j] > atoms__.type[k])  chiralAtom = j;
      if (atoms__.type[k] > atoms__.type[j])  chiralAtom = k;
      if (atomid__.atomic[j] > atomid__.atomic[k])  chiralAtom = j;
      if (atomid__.atomic[k] > atomid__.atomic[j])  chiralAtom = k;
   }
   return chiralAtom;
}

static void setupAmoebaTorsionTorsionForce (OpenMM_System* system, FILE* log) {

   int ii, jj, index, count;
   int ia, ib, ic, id, ie, ichiral;
   int gridIndex;
   int xIndex, yIndex;
   int* ibitorPtr;
   int* ittPtr;
   int maxntt, maxtgrd, maxtgrd2;
   int addIndex;
   int numberOfTorsionTorsions;
   OpenMM_DoubleArray* values;
   OpenMM_3D_DoubleArray* grid;

   OpenMM_AmoebaTorsionTorsionForce* amoebaTorsionTorsionForce;
   amoebaTorsionTorsionForce = OpenMM_AmoebaTorsionTorsionForce_create ();
   OpenMM_System_addForce (system, (OpenMM_Force*) amoebaTorsionTorsionForce);

   ittPtr = tortor__.itt;
   for (ii = 0; ii < tortor__.ntortor; ii++) {
      index = *(ittPtr) - 1;
      gridIndex = *(ittPtr+1) - 1;
      ibitorPtr = bitor__.ibitor + index*5;
      if (*(ittPtr+2) == 1) {
         count = 0;
         ia = *(ibitorPtr + count++) - 1;
         ib = *(ibitorPtr + count++) - 1;
         ic = *(ibitorPtr + count++) - 1;
         id = *(ibitorPtr + count++) - 1;
         ie = *(ibitorPtr + count++) - 1;
      } else {
         count = 4;
         ia = *(ibitorPtr + count--) - 1;
         ib = *(ibitorPtr + count--) - 1;
         ic = *(ibitorPtr + count--) - 1;
         id = *(ibitorPtr + count--) - 1;
         ie = *(ibitorPtr + count--) - 1;
      }
      ichiral = getChiralIndex (ib, ic, id);
      OpenMM_AmoebaTorsionTorsionForce_addTorsionTorsion
                      (amoebaTorsionTorsionForce, ia, ib, ic, id, ie,
                       ichiral, gridIndex );
      ittPtr += 3;
   }

   numberOfTorsionTorsions  =
                 OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsions
                    (amoebaTorsionTorsionForce);

   maxntt = ktrtor__.maxntt;
   maxtgrd = ktrtor__.maxtgrd;
   maxtgrd2 = maxtgrd*maxtgrd;
   values = OpenMM_DoubleArray_create(6);

   if (numberOfTorsionTorsions) {

      int kk = 0;
      for (int i = 0; i < maxntt; i++) {
         char char21[21];
         for (int j = 0; j < 20; j++) {
            char21[j] = ktrtor__.ktt[i].s20[j];
         }
         char21[20] = '\n';
         if (char21[0] != ' ') {
            kk++;
         }
      }

      for (ii = 0; ii < kk; ii++) {
         grid = OpenMM_3D_DoubleArray_create (*(ktrtor__.tnx+ii),
                                              *(ktrtor__.tny+ii), 6);
         xIndex = 0;
         yIndex = 0;

         for (jj = 0; jj < *(ktrtor__.tnx+ii)*(*(ktrtor__.tny+ii)); jj++) {
            addIndex  = 0;
            OpenMM_DoubleArray_set (values, addIndex++,
                *(ktrtor__.ttx + maxtgrd*ii+xIndex));
            OpenMM_DoubleArray_set (values, addIndex++,
                *(ktrtor__.tty + maxtgrd*ii+yIndex));

            OpenMM_DoubleArray_set (values, addIndex++,
                OpenMM_KJPerKcal*(*(ktrtor__.tbf  + maxtgrd2*ii + jj)));
            OpenMM_DoubleArray_set (values, addIndex++,
                OpenMM_KJPerKcal*(*(ktrtor__.tbx  + maxtgrd2*ii + jj)));
            OpenMM_DoubleArray_set (values, addIndex++,
                OpenMM_KJPerKcal*(*(ktrtor__.tby  + maxtgrd2*ii + jj)));
            OpenMM_DoubleArray_set (values, addIndex++,
                OpenMM_KJPerKcal*(*(ktrtor__.tbxy + maxtgrd2*ii + jj)));

            OpenMM_3D_DoubleArray_set (grid, yIndex, xIndex, values);

            xIndex++;
            if (xIndex == *(ktrtor__.tnx+ii)) {
               xIndex = 0;
               yIndex++;
            }
         }

         OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid
                 (amoebaTorsionTorsionForce, ii, grid);
         OpenMM_3D_DoubleArray_destroy (grid);
      }
   }

   OpenMM_DoubleArray_destroy (values);
}


