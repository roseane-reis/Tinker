#ifndef SETUP_FORCES_H_
#define SETUP_FORCES_H_

#include "OpenMMCWrapper.h"
#include "AmoebaOpenMMCWrapper.h"
#include "OpenMM.h"
using namespace OpenMM;

/*
 *  #############################################################
 *  ConstraintMap  --  list of distance constraints for each atom
 *  #############################################################
 */

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

void setupSystemParticles (OpenMM_System* , FILE* );
void setupCMMotionRemover (OpenMM_System* , FILE* );
void freeConstraintMap (ConstraintMap*);
void mapConstraints (ConstraintMap* map, FILE* );
int checkForConstraint (ConstraintMap*, int, int, FILE*);

void setupAmoebaBondForce (OpenMM_System*, int, FILE*);
void setupAmoebaAngleForce(OpenMM_System*, int, FILE*);
void setupAmoebaInPlaneAngleForce(OpenMM_System*, int, FILE*);
void setupAmoebaStrechBendForce(OpenMM_System*, int, FILE*);
void setupAmoebaUreyBradleyForce(OpenMM_System*, int, FILE*);
void setupAmoebaOutOfPlaneBendForce (OpenMM_System*, FILE*);
void setupAmoebaImproperTorsionForce(OpenMM_System*, FILE*);
void setupAmoebaTorsionForce(OpenMM_System*, FILE*);
void setupAmoebaPiTorsionForce(OpenMM_System*, FILE*);
void setupAmoebaTorsionTorsionForce(OpenMM_System*, FILE*);

#endif
