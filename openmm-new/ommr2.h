double getObcShct (int atmnum) {

   double shct = 0.80;
   if (atmnum == 1)  shct = 0.85;
   if (atmnum == 6)  shct = 0.72;
   if (atmnum == 7)  shct = 0.79;
   if (atmnum == 8)  shct = 0.85;
   if (atmnum == 9)  shct = 0.88;
   if (atmnum == 15)  shct = 0.86;
   if (atmnum == 16)  shct = 0.96;
   if (atmnum == 26)  shct = 0.88;
   return shct;
}

void setupAmoebaGeneralizedKirkwoodForce (OpenMM_System* system,
                                         int includeCavityTerm, FILE* log) {

   int ii;
   int useGrycuk;
   int useObc;
   double* polePtr;
   double shct;
   char buffer[80];

   // check Born radius type; force use of "Grycuk" for now

   setNullTerminator (solute__.borntyp, 8, buffer);
   useGrycuk = 1;
   if (strncasecmp (buffer, "GRYCUK", 6 ) != 0) {
      if (log) {
         (void) fprintf (log, "setupAmoebaGeneralizedKirkwoodForce: Born type=%s, forcing Born type to 'Grycuk'.\n", buffer);
      } else {
         (void) fprintf (stderr, "setupAmoebaGeneralizedKirkwoodForce: Born type=%s, forcing Born type to 'Grycuk'.\n", buffer);
      }
   }

   OpenMM_AmoebaGeneralizedKirkwoodForce* amoebaGeneralizedKirkwoodForce;
   amoebaGeneralizedKirkwoodForce =
       OpenMM_AmoebaGeneralizedKirkwoodForce_create ();
   OpenMM_System_addForce (system, (OpenMM_Force*)
                           amoebaGeneralizedKirkwoodForce);

   OpenMM_AmoebaGeneralizedKirkwoodForce_setSolventDielectric
                          (amoebaGeneralizedKirkwoodForce, 78.3);
   OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric
                          (amoebaGeneralizedKirkwoodForce, *chgpot__.dielec);
   OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm
                          (amoebaGeneralizedKirkwoodForce, includeCavityTerm);
   // OpenMM_AmoebaGeneralizedGeneralizedKirkwoodForce_setDielectricOffset
   //                     (amoebaGeneralizedGeneralizedKirkwoodForce,
   //                      solute__.doffset*OpenMM_NmPerAngstrom);
   // OpenMM_AmoebaGeneralizedGeneralizedKirkwoodForce_setProbeRadius
   //                     (OpenMM_AmoebaGeneralizedGeneralizedKirkwoodForce,
   //                      0.14);
   // OpenMM_AmoebaGeneralizedGeneralizedKirkwoodForce_setSurfaceAreaFactor
   //                     (amoebaGeneralizedGeneralizedKirkwoodForce,
   //                      surfaceAreaFactor);

   for (ii = 0; ii < *atoms__.n; ii++) {
      polePtr = mpole__.pole + ii*mpole__.maxpole;
      if (useGrycuk) {
         shct = solute__.shct[ii];
      } else {
         shct = getObcShct (atomid__.atomic[ii]);
      }
      OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle
                         (amoebaGeneralizedKirkwoodForce, *(polePtr),
                          OpenMM_NmPerAngstrom*solute__.rsolv[ii], shct);
   }
}

/*
 *    ############################################################
 *          Setup Geometric Restraint Potential Energy Terms
 *    ############################################################
 */

#define RADIANS_PER_DEGREE 0.0174532925

void setupPositionalRestraints (OpenMM_System* system, FILE* log) {

   double convert = OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom*OpenMM_NmPerAngstrom);
   double nmPerAng = OpenMM_NmPerAngstrom;

   OpenMM_CustomExternalForce* force;

   if (boxes__.orthogonal) {
      force = OpenMM_CustomExternalForce_create("k*(max(0.0,periodicdistance(x, y, z, x0, y0, z0)-range))^2");
   } else if (boxes__.monoclinic) {
      force = OpenMM_CustomExternalForce_create("k*(max(0.0,periodicdistance(x, y, z, x0, y0, z0)-range))^2");
      OpenMM_CustomExternalForce_addGlobalParameter (force, "beta_sin", *boxes__.beta_sin);
      OpenMM_CustomExternalForce_addGlobalParameter (force, "beta_cos", *boxes__.beta_cos);
   } else if (boxes__.triclinic) {
      force = OpenMM_CustomExternalForce_create("k*(max(0.0,periodicdistance(x, y, z, x0, y0, z0)-range))^2");
      OpenMM_CustomExternalForce_addGlobalParameter (force, "beta_sin", *boxes__.beta_sin);
      OpenMM_CustomExternalForce_addGlobalParameter (force, "beta_cos", *boxes__.beta_cos);
      OpenMM_CustomExternalForce_addGlobalParameter (force, "gamma_term", *boxes__.gamma_term);
      OpenMM_CustomExternalForce_addGlobalParameter (force, "beta_term", *boxes__.beta_term);
      OpenMM_CustomExternalForce_addGlobalParameter (force, "gamma_sin", *boxes__.gamma_sin);
      OpenMM_CustomExternalForce_addGlobalParameter (force, "gamma_cos", *boxes__.gamma_cos);
   } else if (boxes__.octahedron) {
      force = OpenMM_CustomExternalForce_create("k*(max(0.0,periodicdistance(x, y, z, x0, y0, z0)-range))^2");
      OpenMM_CustomExternalForce_addGlobalParameter (force, "box34", *boxes__.box34*nmPerAng);
      OpenMM_CustomExternalForce_addGlobalParameter (force, "xbox", *boxes__.xbox*nmPerAng);
      OpenMM_CustomExternalForce_addGlobalParameter (force, "xbox2", *boxes__.xbox2*nmPerAng);
      OpenMM_CustomExternalForce_addGlobalParameter (force, "ybox", *boxes__.ybox*nmPerAng);
      OpenMM_CustomExternalForce_addGlobalParameter (force, "ybox2", *boxes__.ybox2*nmPerAng);
      OpenMM_CustomExternalForce_addGlobalParameter (force, "zbox", *boxes__.zbox*nmPerAng);
      OpenMM_CustomExternalForce_addGlobalParameter (force, "zbox2", *boxes__.zbox2*nmPerAng);
   }

   OpenMM_CustomExternalForce_addPerParticleParameter (force, "k");
   OpenMM_CustomExternalForce_addPerParticleParameter (force, "range");
   OpenMM_CustomExternalForce_addPerParticleParameter (force, "x0");
   OpenMM_CustomExternalForce_addPerParticleParameter (force, "y0");
   OpenMM_CustomExternalForce_addPerParticleParameter (force, "z0");
   OpenMM_CustomExternalForce_addPerParticleParameter (force, "use_bounds");
   OpenMM_CustomExternalForce_addGlobalParameter (force, "xcell", *cell__.xcell*nmPerAng);
   OpenMM_CustomExternalForce_addGlobalParameter (force, "xcell2", *cell__.xcell2*nmPerAng);
   OpenMM_CustomExternalForce_addGlobalParameter (force, "ycell", *cell__.ycell*nmPerAng);
   OpenMM_CustomExternalForce_addGlobalParameter (force, "ycell2", *cell__.ycell2*nmPerAng);
   OpenMM_CustomExternalForce_addGlobalParameter (force, "zcell", *cell__.zcell*nmPerAng);
   OpenMM_CustomExternalForce_addGlobalParameter (force, "zcell2", *cell__.zcell2*nmPerAng);

   for (int i = 0; i < *restrn__.npfix; ++i) {
      double use_bounds = 0.0;
      if (bound__.use_bounds)
         use_bounds = 1.0;
      OpenMM_DoubleArray* positionalParameters = OpenMM_DoubleArray_create(0);
      OpenMM_DoubleArray_append(positionalParameters, restrn__.pfix[i*2]*convert);
      OpenMM_DoubleArray_append(positionalParameters, restrn__.pfix[i*2 + 1] * nmPerAng);
      OpenMM_DoubleArray_append(positionalParameters, restrn__.xpfix[i] * nmPerAng);
      OpenMM_DoubleArray_append(positionalParameters, restrn__.ypfix[i] * nmPerAng);
      OpenMM_DoubleArray_append(positionalParameters, restrn__.zpfix[i] * nmPerAng);
      OpenMM_DoubleArray_append(positionalParameters, use_bounds);
      OpenMM_CustomExternalForce_addParticle(force, restrn__.ipfix[i]-1, positionalParameters);
      OpenMM_DoubleArray_destroy(positionalParameters);
   }
   OpenMM_System_addForce(system, (OpenMM_Force*) force);
}

void setupDistanceRestraints (OpenMM_System* system, FILE* log) {

   const double convert = OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom*OpenMM_NmPerAngstrom);
   const double nmPerAng = OpenMM_NmPerAngstrom;

   OpenMM_CustomCompoundBondForce* force;

   if (boxes__.orthogonal) {
      force = OpenMM_CustomCompoundBondForce_create(2,"k*(max(max(0.0,distMin-r),max(0.0,r-distMax)))^2;\
         r = sqrt(xr^2+yr^2+zr^2);\
         xr = xr_pre-use_bounds*(xcell*ceil(max(0.0,(abs(xr_pre)-xcell2))/xcell)*(step(xr_pre)*2-1));\
         yr = yr_pre-use_bounds*(ycell*ceil(max(0.0,(abs(yr_pre)-ycell2))/ycell)*(step(yr_pre)*2-1));\
         zr = zr_pre-use_bounds*(zcell*ceil(max(0.0,(abs(zr_pre)-zcell2))/zcell)*(step(zr_pre)*2-1));\
         xr_pre=x1-x2;\
         yr_pre=y1-y2;\
         zr_pre=z1-z2");
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "xcell", *cell__.xcell*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "xcell2", *cell__.xcell2*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "ycell", *cell__.ycell*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "ycell2", *cell__.ycell2*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "zcell", *cell__.zcell*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "zcell2", *cell__.zcell2*nmPerAng);
   } else if (boxes__.monoclinic) {
      force = OpenMM_CustomCompoundBondForce_create(2, "k*(max(max(0.0,distMin-r),max(0.0,r-distMax)))^2;\
         r = sqrt(xr^2+yr^2+zr^2);\
         xr = xr_converted + zr_converted*beta_cos;\
         zr = zr_converted*beta_sin;\
         yr = yr_pre-use_bounds*(ycell*ceil(max(0,(abs(yr_pre)-ycell2))/ycell)*(step(yr_pre)*2-1));\
         xr_converted = xr_pre-use_bounds*(xcell*ceil(max(0,(abs(xr_pre)-xcell2))/xcell)*(step(xr_pre)*2-1));\
         zr_converted = zr_pre-use_bounds*(zcell*ceil(max(0,(abs(zr_pre)-zcell2))/zcell)*(step(zr_pre)*2-1));\
         xr_pre = x1-x2 - zr_pre*beta_cos;\
         yr_pre = y1-y2;\
         zr_pre = (z1-z2)/beta_sin");
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "beta_sin", *boxes__.beta_sin);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "beta_cos", *boxes__.beta_cos);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "xcell", *cell__.xcell*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "xcell2", *cell__.xcell2*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "ycell", *cell__.ycell*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "ycell2", *cell__.ycell2*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "zcell", *cell__.zcell*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "zcell2", *cell__.zcell2*nmPerAng);
   } else if (boxes__.triclinic) {
      force = OpenMM_CustomCompoundBondForce_create(2, "k*(max(max(0.0,distMin-r),max(0.0,r-distMax)))^2;\
         r = sqrt(xr^2+yr^2+zr^2);\
         xr = xr_converted + yr_converted*gamma_cos + zr_converted*beta_cos;\
         yr = yr_converted*gamma_sin + zr_converted*beta_term;\
         zr = zr_converted*gamma_term;\
         xr_converted = xr_pre-use_bounds*(xcell*ceil(max(0,(abs(xr_pre)-xcell2))/xcell)*(step(xr_pre)*2-1));\
         yr_converted = yr_pre-use_bounds*(ycell*ceil(max(0,(abs(yr_pre)-ycell2))/ycell)*(step(yr_pre)*2-1));\
         zr_converted = zr_pre-use_bounds*(zcell*ceil(max(0,(abs(zr_pre)-zcell2))/zcell)*(step(zr_pre)*2-1));\
         xr_pre = x1-x2 - yr_pre*gamma_cos - zr_pre*beta_cos;\
         yr_pre = (y1-y2 - zr_pre*beta_term)/gamma_sin;\
         zr_pre = (z1-z2)/gamma_term");
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "beta_sin", *boxes__.beta_sin);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "beta_cos", *boxes__.beta_cos);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "gamma_term", *boxes__.gamma_term);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "beta_term", *boxes__.beta_term);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "gamma_sin", *boxes__.gamma_sin);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "gamma_cos", *boxes__.gamma_cos);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "xcell", *cell__.xcell*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "xcell2", *cell__.xcell2*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "ycell", *cell__.ycell*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "ycell2", *cell__.ycell2*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "zcell", *cell__.zcell*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "zcell2", *cell__.zcell2*nmPerAng);
   } else if (boxes__.octahedron) {
      force = OpenMM_CustomCompoundBondForce_create(2, "k*(max(max(0.0,distMin-r),max(0.0,r-distMax)))^2;\
         r = sqrt(xr^2+yr^2+zr^2);\
         xr = xr_converted - absDist*xbox2 * (step(xr_converted)*2-1);\
         yr = yr_converted - absDist*ybox2 * (step(yr_converted)*2-1);\
         zr = zr_converted - absDist*zbox2 * (step(zr_converted)*2-1);\
         absDist = step(abs(xr_converted)+abs(yr_converted)+abs(zr_converted)-box34);\
         xr_converted = xr_pre-use_bounds*(xbox*ceil(max(0,(abs(xr_pre)-xbox2))/xbox)*(step(xr_pre)*2-1));\
         yr_converted = yr_pre-use_bounds*(ybox*ceil(max(0,(abs(yr_pre)-ybox2))/ybox)*(step(yr_pre)*2-1));\
         zr_converted = zr_pre-use_bounds*(zbox*ceil(max(0,(abs(zr_pre)-zbox2))/zbox)*(step(zr_pre)*2-1));\
         xr_pre = x1-x2;\
         yr_pre = y1-y2;\
         zr_pre = z1-z2");
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "box34", *boxes__.box34*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "xbox", *boxes__.xbox*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "xbox2", *boxes__.xbox2*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "ybox", *boxes__.ybox*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "ybox2", *boxes__.ybox2*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "zbox", *boxes__.zbox*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "zbox2", *boxes__.zbox2*nmPerAng);
   }

   OpenMM_CustomCompoundBondForce_addPerBondParameter (force, "k");
   OpenMM_CustomCompoundBondForce_addPerBondParameter (force, "distMin");
   OpenMM_CustomCompoundBondForce_addPerBondParameter (force, "distMax");
   OpenMM_CustomCompoundBondForce_addPerBondParameter (force, "use_bounds");

   for (int i = 0; i < *restrn__.ndfix; ++i) {
      double use_bounds = 0.0;
      if (bound__.use_bounds && molcul__.molcule[restrn__.idfix[i*2]-1]
                               != molcul__.molcule[restrn__.idfix[i*2+1]-1]) {
         use_bounds = 1.0;
      }
      OpenMM_IntArray* particles = OpenMM_IntArray_create(0);
      OpenMM_IntArray_append (particles, restrn__.idfix[i*2]-1);
      OpenMM_IntArray_append (particles, restrn__.idfix[i*2+1]-1);
      OpenMM_DoubleArray* distanceParameters = OpenMM_DoubleArray_create(0);
      OpenMM_DoubleArray_append (distanceParameters, restrn__.dfix[i*3]*convert);
      OpenMM_DoubleArray_append (distanceParameters, restrn__.dfix[i*3 + 1]*nmPerAng);
      OpenMM_DoubleArray_append (distanceParameters, restrn__.dfix[i*3 + 2]*nmPerAng);
      OpenMM_DoubleArray_append (distanceParameters, use_bounds);
      OpenMM_CustomCompoundBondForce_addBond (force, particles, distanceParameters);
      OpenMM_IntArray_destroy (particles);
      OpenMM_DoubleArray_destroy (distanceParameters);
   }
   OpenMM_CustomCompoundBondForce_setUsesPeriodicBoundaryConditions (force, OpenMM_True);
   OpenMM_System_addForce (system, (OpenMM_Force*) force);
}

void setupAngleRestraints (OpenMM_System* system, FILE* log) {

   double convert = OpenMM_KJPerKcal / RADIANS_PER_DEGREE / RADIANS_PER_DEGREE;

   OpenMM_CustomAngleForce* force =
      OpenMM_CustomAngleForce_create("k*(max(max(0.0,thetaMin-theta),max(0.0,theta-thetaMax)))^2");

   OpenMM_CustomAngleForce_addPerAngleParameter (force, "k");
   OpenMM_CustomAngleForce_addPerAngleParameter (force, "thetaMin");
   OpenMM_CustomAngleForce_addPerAngleParameter (force, "thetaMax");

   for (int i = 0; i < *restrn__.nafix; ++i) {
      OpenMM_DoubleArray* AngleParameters = OpenMM_DoubleArray_create(0);
      OpenMM_DoubleArray_append (AngleParameters, restrn__.afix[i*2]*convert);
      OpenMM_DoubleArray_append (AngleParameters, restrn__.afix[i*2 + 1]*RADIANS_PER_DEGREE);
      OpenMM_DoubleArray_append (AngleParameters, restrn__.afix[i*2 + 2]*RADIANS_PER_DEGREE);
      OpenMM_CustomAngleForce_addAngle (force, restrn__.iafix[i*3]-1, restrn__.iafix[i*3+1]-1,
                                       restrn__.iafix[i*3+2]-1, AngleParameters);
      OpenMM_DoubleArray_destroy (AngleParameters);
   }
   OpenMM_System_addForce (system, (OpenMM_Force*) force);
}

void setupTorsionRestraints (OpenMM_System* system, FILE* log) {

   double convert;
   convert = OpenMM_KJPerKcal / RADIANS_PER_DEGREE / RADIANS_PER_DEGREE;

   OpenMM_CustomTorsionForce* force =
         OpenMM_CustomTorsionForce_create ("k*max(\
            (step(thetaMin-theta)*(min(min(abs(theta-thetaMin),abs(theta-thetaMin-6.28318530718)), abs(theta-thetaMin+6.28318530718)))),\
            (step(theta-thetaMax)*(min(min(abs(theta-thetaMax),abs(theta-thetaMax-6.28318530718)), abs(theta-thetaMax+6.28318530718))))\
               )^2");

   OpenMM_CustomTorsionForce_addPerTorsionParameter (force, "k");
   OpenMM_CustomTorsionForce_addPerTorsionParameter (force, "thetaMin");
   OpenMM_CustomTorsionForce_addPerTorsionParameter (force, "thetaMax");

   for (int i = 0; i < *restrn__.ntfix; i++) {
      float thetaMin = restrn__.tfix[i*3+1];
      float thetaMax = restrn__.tfix[i*3+2];
      if (thetaMin > 180.0f)
         thetaMin -= 360.0f;
      else if (thetaMin < -180.0f)
         thetaMin += 360.0f;
      if (thetaMax > 180.0f)
         thetaMax -= 360.0f;
      else if (thetaMax < -180.0f)
         thetaMax += 360.0f;

      OpenMM_DoubleArray* lowerTorsionParameters = OpenMM_DoubleArray_create (0);
      OpenMM_DoubleArray_append (lowerTorsionParameters, restrn__.tfix[i*3]*convert);
      OpenMM_DoubleArray_append (lowerTorsionParameters, thetaMin*RADIANS_PER_DEGREE);
      OpenMM_DoubleArray_append (lowerTorsionParameters, thetaMax*RADIANS_PER_DEGREE);
      OpenMM_CustomTorsionForce_addTorsion (force, restrn__.itfix[i*4]-1,
                                            restrn__.itfix[i*4+1]-1,
                                            restrn__.itfix[i*4+2]-1,
                                            restrn__.itfix[i*4+3]-1,
                                            lowerTorsionParameters);
      OpenMM_DoubleArray_destroy (lowerTorsionParameters);
   }
   OpenMM_System_addForce (system, (OpenMM_Force*) force);
}

OpenMM_IntArray* getGroup (int* kgrp, int* igrp, int idx) {

   int start = igrp[2 * idx] - 1;
   int end = igrp[2 * idx + 1] - 1;
   OpenMM_IntArray* group = OpenMM_IntArray_create (0);
   for (int i = start; i < end + 1; i++) {
      OpenMM_IntArray_append (group, kgrp[i] - 1);
   }
   return group;
}

OpenMM_DoubleArray* getWeights (int* kgrp, int* igrp, int idx,
                                       OpenMM_System* system) {

   int start = igrp[2 * idx] - 1;
   int end = igrp[2 * idx + 1] - 1;
   OpenMM_DoubleArray* weights = OpenMM_DoubleArray_create (0);
   for (int i = start; i < end + 1; i++){
      OpenMM_DoubleArray_append (weights, OpenMM_System_getParticleMass(system,
                                 kgrp[i] - 1));
   }
   return weights;
}

void setupCentroidRestraints (OpenMM_System* system, FILE* log) {

   double convert;
   convert = OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom*OpenMM_NmPerAngstrom);

   // In the expression below, u and l are the upper and lower threshold

   OpenMM_CustomCentroidBondForce* force =
      OpenMM_CustomCentroidBondForce_create (2,
                 "step(distance(g1,g2)-u)*k*(distance(g1,g2)-u)^2+step(l-distance(g1,g2))*k*(distance(g1,g2)-l)^2");
   OpenMM_CustomCentroidBondForce_addPerBondParameter (force, "k");
   OpenMM_CustomCentroidBondForce_addPerBondParameter (force, "l");
   OpenMM_CustomCentroidBondForce_addPerBondParameter (force, "u");

   for (int j = 1; j < *group__.ngrp + 1; j++) {
      OpenMM_IntArray* igroup = getGroup(group__.kgrp, group__.igrp, j);
      OpenMM_DoubleArray* iweight = getWeights(group__.kgrp, group__.igrp, j, system);
      OpenMM_CustomCentroidBondForce_addGroup (force, igroup, iweight);
      OpenMM_IntArray_destroy (igroup);
      OpenMM_DoubleArray_destroy (iweight);
   }

   for (int i = 0; i < *restrn__.ngfix; i++) {
      OpenMM_IntArray* bondGroups = OpenMM_IntArray_create (0);
      OpenMM_IntArray_append (bondGroups, restrn__.igfix[2*i] - 1);
      OpenMM_IntArray_append (bondGroups, restrn__.igfix[2*i + 1] - 1);

      OpenMM_DoubleArray* bondParameters = OpenMM_DoubleArray_create (0);
      OpenMM_DoubleArray_append (bondParameters, restrn__.gfix[3*i]
                                    *convert);
      OpenMM_DoubleArray_append (bondParameters, restrn__.gfix[3*i + 1]
                                    *OpenMM_NmPerAngstrom);
      OpenMM_DoubleArray_append (bondParameters, restrn__.gfix[3*i + 2]
                                    *OpenMM_NmPerAngstrom);

      OpenMM_CustomCentroidBondForce_addBond (force, bondGroups,
                                              bondParameters);
      OpenMM_IntArray_destroy (bondGroups);
      OpenMM_DoubleArray_destroy (bondParameters);
   }
   OpenMM_CustomCentroidBondForce_setUsesPeriodicBoundaryConditions (force, OpenMM_True);
   OpenMM_System_addForce (system, (OpenMM_Force*) force);
}
/*
 *    ############################################################
 *          Create Standalone Thermotat and Barostat Methods
 *    ############################################################
 */

void setupAndersenThermostat (OpenMM_System* system, FILE* log) {

   OpenMM_AndersenThermostat* andersenThermostat;

   andersenThermostat = OpenMM_AndersenThermostat_create (*bath__.kelvin,
                                                          *bath__.collide);
   OpenMM_System_addForce (system, (OpenMM_Force*) andersenThermostat);

   if (log) {
      (void) fprintf (log, "\n Andersen Thermostat:\n" );
      (void) fprintf (log, "\n Target Temperature   %15.2f K",
                      OpenMM_AndersenThermostat_getDefaultTemperature
                      (andersenThermostat));
      (void) fprintf (log, "\n Collision Frequency  %15.7e ps^(-1)",
                      OpenMM_AndersenThermostat_getDefaultCollisionFrequency
                      (andersenThermostat));
      (void) fprintf (log, "\n Random Number Seed                 %d\n",
                      OpenMM_AndersenThermostat_getRandomNumberSeed
                      (andersenThermostat));
   }
}

void setupMonteCarloBarostat (OpenMM_System* system, FILE* log) {

   OpenMM_MonteCarloBarostat* monteCarloBarostat;

   int frequency = *bath__.voltrial;
   monteCarloBarostat = OpenMM_MonteCarloBarostat_create
                            (*bath__.atmsph*1.01325, *bath__.kelvin, frequency);
   OpenMM_System_addForce (system, (OpenMM_Force*) monteCarloBarostat);

   if (log) {
      (void) fprintf (log, "\n MonteCarlo Barostat :\n");
      (void) fprintf (log, "\n Target Temperature   %15.2f K",
                      OpenMM_MonteCarloBarostat_getDefaultTemperature
                      (monteCarloBarostat));

      //
      // Change needed for latest Stanford OpenMM update (M. Harger)
      //
      // (void) fprintf (log, "\n Target Temperature   %15.2f K",
      //                 OpenMM_MonteCarloBarostat_getDefaultTemperature
      //                 (monteCarloBarostat));

      (void) fprintf (log, "\n Target Pressure      %15.4f atm",
                      OpenMM_MonteCarloBarostat_getDefaultPressure
                      (monteCarloBarostat) / 1.01325);
      (void) fprintf (log, "\n Trial Frequency                   %d",
                      OpenMM_MonteCarloBarostat_getFrequency
                      (monteCarloBarostat));
      (void) fprintf (log, "\n Random Number Seed                 %d\n",
                      OpenMM_MonteCarloBarostat_getRandomNumberSeed
                      (monteCarloBarostat));
   }
}

/*
 *    ############################################################
 *         Setup Positions, Velocities and Rattle Constraints
 *    ############################################################
 */

void setupPositions (OpenMM_Vec3Array* initialPosInNm, FILE* log) {

   int ii;
   for (ii = 0; ii < *atoms__.n; ii++) {
      OpenMM_Vec3 posInNm;
      posInNm.x = atoms__.x[ii]*OpenMM_NmPerAngstrom;
      posInNm.y = atoms__.y[ii]*OpenMM_NmPerAngstrom;
      posInNm.z = atoms__.z[ii]*OpenMM_NmPerAngstrom;
      OpenMM_Vec3Array_append (initialPosInNm, posInNm);
   }
}

void setupVelocities (OpenMM_Vec3Array* initialVelInNm, FILE* log) {

   int ii;
   for (ii = 0; ii < *atoms__.n; ii++) {
      OpenMM_Vec3 velInNm;
      int offset;
      offset = 3*ii;
      velInNm.x = moldyn__.v[offset]*OpenMM_NmPerAngstrom;
      velInNm.y = moldyn__.v[offset+1]*OpenMM_NmPerAngstrom;
      velInNm.z = moldyn__.v[offset+2]*OpenMM_NmPerAngstrom;
      OpenMM_Vec3Array_append (initialVelInNm, velInNm);
   }
}

void setupConstraints (OpenMM_System* system, FILE* log) {

   int ii;
   for (ii = 0; ii < *freeze__.nrat; ii++) {
      OpenMM_System_addConstraint (system, *(freeze__.irat+2*ii) -1,
                                   *(freeze__.irat + 2*ii +1)-1,
                          (*(freeze__.krat +ii))*OpenMM_NmPerAngstrom);
   }
}

/*
 *    ############################################################
 *           Set Calculation Platform to Reference or CUDA
 *    ############################################################
 */

#include "gpu_cards.h"

static OpenMM_Platform* getReferencePlatform (FILE* log) {

   OpenMM_Platform* platform = OpenMM_Platform_getPlatformByName ("Reference");
   if (platform == NULL) {
      if (log) {
         (void) fprintf (log, "Reference Platform Unavailable\n");
      }
      return platform;
   }
   return platform;
}

static OpenMM_Platform* getCUDAPlatform (FILE* log) {

   char buffer[8];
   char* deviceId;
   int device_number;
   bool device_key = false;
   OpenMM_Platform* platform = OpenMM_Platform_getPlatformByName ("CUDA");
   if (platform == NULL) {
      if (log) {
         (void) fprintf (log, "\n CUDA Platform Unavailable\n");
      }
      return platform;
   }

   if (openmm__.cudadevice[0] != 0) {
      deviceId = &openmm__.cudadevice[0];
      device_key = true;
   } else {
      device_number = findBestCUDACard();
      if (device_number < 0) {
         deviceId = NULL;
      } else {
         deviceId = buffer;
         sprintf(deviceId, "%d", device_number);
      }
   }

   if (device_key) {
      OpenMM_Platform_setPropertyDefaultValue (platform, "cudadeviceIndex",
                                               deviceId);
      if (log) {
         (void) fprintf (log, "\n Platform CUDA :  Setting Device ID to %s from CUDA-DEVICE keyword\n", deviceId);
      }
   } else if (log && inform__.verbose) {
      (void) fprintf (log, "\n Platform CUDA :  Setting Device ID to %s \n", deviceId);
   }

   if (strncmp(openmm__.cudaprecision,"DOUBLE",6) == 0) {
      OpenMM_Platform_setPropertyDefaultValue (platform, "Precision",
                                               "double" );
   } else if (strncmp(openmm__.cudaprecision,"SINGLE",6) == 0) {
      OpenMM_Platform_setPropertyDefaultValue (platform, "Precision",
                                               "single" );
   } else {
      OpenMM_Platform_setPropertyDefaultValue (platform, "Precision",
                                               "mixed" );
   }

   if (log) {
      (void) fprintf (log, "\n Platform CUDA :  Setting Precision to %s via CUDA-PRECISION\n", openmm__.cudaprecision);
   }

   return platform;
}
