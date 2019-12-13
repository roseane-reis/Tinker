extern "C" {

extern void fatal_ ();

void openmm_init_ (void** ommHandle, double* dt) { 

   int ii;
   int mdMode = 0;
   int removeConstrainedCovalentIxns = 0;
   char buffer[128];
   char buffer2[128];
   FILE* log = stderr;

   // Check if we can run the simulation with OpenMM
   /*
   if (boxes__.octahedron) {
      printf("\n Truncated octahedral boxes are not supported in OpenMM.\n");
      fatal_(); //This will never return
   }
   */

   // allocate space for opaque handle to hold OpenMM objects
   // such as system, integrator, context, etc.
   OpenMMData* omm = (OpenMMData*) malloc(sizeof(struct OpenMMData_s));

   // These are temporary OpenMM objects used and discarded here

   OpenMM_Vec3Array*       initialPosInNm;
   OpenMM_Vec3Array*       initialVelInNm;
   OpenMM_StringArray*     pluginList;
   OpenMM_Platform*        platform;

   // Load all OpenMM plugin libraries from their default location;
   // Call the plugin loading routine twice to fix an issue with MacOS
   // where the first library in the alphabetical list gets skipped

   pluginList = OpenMM_Platform_loadPluginsFromDirectory
                    (OpenMM_Platform_getDefaultPluginsDirectory());
   pluginList = OpenMM_Platform_loadPluginsFromDirectory
                    (OpenMM_Platform_getDefaultPluginsDirectory());

   if (*inform__.verbose) {
      (void) fprintf (stderr, "\n Default OpenMM Plugin Directory :  %s\n\n",
                          OpenMM_Platform_getDefaultPluginsDirectory());
      for (ii = 0; ii < OpenMM_StringArray_getSize(pluginList); ii++) {
         (void) fprintf (stderr, " Plugin Library :  %s\n",
                             OpenMM_StringArray_get(pluginList, ii));
      }
   }

   OpenMM_StringArray_destroy (pluginList);
   (void) fflush (NULL);

   // create System and Force objects within the System; retain a reference
   // to each force object so we can fill in the forces; mote the OpenMM
   // System takes ownership of the force objects, don't delete them yourself

   omm->system = OpenMM_System_create ();
   setupSystemParticles (omm->system, log);
   setupCMMotionRemover (omm->system, log);
   
   if (*potent__.use_bond) {
      setupAmoebaBondForce (omm->system,
                            removeConstrainedCovalentIxns, log);
   }

   if (*potent__.use_angle) {
      setupAmoebaAngleForce (omm->system,
                             removeConstrainedCovalentIxns, log);
      setupAmoebaInPlaneAngleForce (omm->system,
                                    removeConstrainedCovalentIxns, log);
   }
   
   if (*potent__.use_strbnd) {
      setupAmoebaStretchBendForce (omm->system,
                                   removeConstrainedCovalentIxns, log);
   }

   if (*potent__.use_urey) {
      setupAmoebaUreyBradleyForce (omm->system,
                                   removeConstrainedCovalentIxns, log);
   }

   if (*potent__.use_opbend) {
      setupAmoebaOutOfPlaneBendForce (omm->system, log);
   }

   if (*potent__.use_improp) {
      setupAmoebaImproperTorsionForce (omm->system, log);
   }

   if (*potent__.use_tors) {
      setupAmoebaTorsionForce (omm->system, log);
   }

   if (*potent__.use_pitors) {
      setupAmoebaPiTorsionForce (omm->system, log);
   }

   if (*potent__.use_tortor) {
      setupAmoebaTorsionTorsionForce (omm->system, log);
   }
   
   if (*potent__.use_charge) {
      setupAmoebaChargeForce (omm->system, log);
   }
   
   if (*potent__.use_chgtrn || *potent__.use_repuls) {
      setupHippoNonbondedForce (omm->system, log);
   }
   
   if (*potent__.use_mpole && ! *mplpot__.use_chgpen) {
      setupAmoebaMultipoleForce (omm->system, log);
      
      if (*potent__.use_solv) {
         setupAmoebaGeneralizedKirkwoodForce (omm->system, 1, log);
      }
      
   }
  
   if (*potent__.use_solv) {
      setupAmoebaWcaDispersionForce (omm->system, log);
   }
    
   if (*bath__.isobaric && *bath__.atmsph > 0.0) {
      mdMode = 4;
      setupMonteCarloBarostat (omm->system, log);
   }

   
   // setup of constraints, restraints, positions and velocities
   setupConstraints (omm->system, log);

   setupTorsionRestraints (omm->system, log);
   setupDistanceRestraints (omm->system, log);
   setupPositionalRestraints (omm->system, log);
   (void) fflush (NULL);
   setupAngleRestraints (omm->system, log);
   setupCentroidRestraints (omm->system, log);

   

   initialPosInNm = OpenMM_Vec3Array_create (0);
   setupPositions (initialPosInNm, log);

   initialVelInNm = OpenMM_Vec3Array_create (0);
   setupVelocities (initialVelInNm, log);

   

   // choose an Integrator, and a Context connecting the System with the
   // Integrator; let the Context choose the best available Platform;
   // initialize the configuration from default positions collected above

   setNullTerminator (mdstuf__.integrate, 11, buffer);
   setNullTerminator (bath__.thermostat, 11, buffer2);

   if (strncasecmp (buffer, "RESPA", 5) == 0) {
      omm->integrator = (OpenMM_Integrator*) OpenMM_CustomIntegrator_create (*dt);
      OpenMM_CustomIntegrator* integrator = (OpenMM_CustomIntegrator*) omm->integrator;
      OpenMM_CustomIntegrator_addUpdateContextState (integrator);

      int n = int (round ((*dt) / *mdstuf__.arespa));
      char n_char[16] = {0};
      sprintf (n_char, "%d", n);

      char e1[64] = {0};
      strcat (e1, "v+0.5*(dt/");
      strcat (e1, n_char);
      strcat (e1, ")*f0/m");

      char e11[64] = {0};
      strcat (e11, n_char);
      strcat (e11, "*(x-x1)/dt+");
      strcat (e11, e1);

      char e2[64] = {0};
      strcat (e2, "x+(dt/");
      strcat (e2, n_char);
      strcat (e2, ")*v");

      OpenMM_CustomIntegrator_addPerDofVariable (integrator, "x1", 0.0);

      OpenMM_CustomIntegrator_addComputePerDof (integrator, "v", "v+0.5*dt*f1/m");
      for (int i = 0; i < n; i++) {
         OpenMM_CustomIntegrator_addComputePerDof (integrator, "v", e1);
         OpenMM_CustomIntegrator_addComputePerDof (integrator, "x", e2);
         OpenMM_CustomIntegrator_addComputePerDof (integrator, "x1", "x");
         OpenMM_CustomIntegrator_addConstrainPositions (integrator);
         OpenMM_CustomIntegrator_addComputePerDof (integrator, "v", e11);
         OpenMM_CustomIntegrator_addConstrainVelocities (integrator);
      }
      OpenMM_CustomIntegrator_addComputePerDof (integrator, "v", "v+0.5*dt*f1/m");
      OpenMM_CustomIntegrator_addConstrainVelocities (integrator);

      // setup and computation for the Bussi-Parrinello thermostat

      int n_constraints = OpenMM_System_getNumConstraints (omm->system) + 3;
      int adjlen = (n_constraints + 2) * 16;
      char* adjustment = (char*) malloc (adjlen * sizeof(char));
      char* e3 = (char*) malloc (adjlen * sizeof(char));
      for (int i = 0; i < adjlen; ++i) {
         adjustment[i] = 0;
         e3[i] = 0;
      }

      for (int i = 0; i < n_constraints; i++) {
         strcat (adjustment, "gaussian^2+");
      }
      strcat (adjustment, "gaussian^2");

      strcat (e3, "s-(");
      strcat (e3, adjustment);
      strcat (e3, ")");

      if (strncasecmp (buffer2, "BUSSI", 5) == 0) {
         OpenMM_CustomIntegrator_addGlobalVariable (integrator, "df", *mdstuf__.nfree);
         OpenMM_CustomIntegrator_addGlobalVariable (integrator, "gasconst",
             units__.gasconst * units__.joule);
         OpenMM_CustomIntegrator_addGlobalVariable (integrator, "tautemp", *bath__.tautemp);
         OpenMM_CustomIntegrator_addGlobalVariable (integrator, "kelvin", *bath__.kelvin);
         OpenMM_CustomIntegrator_addGlobalVariable (integrator, "ke", 0.0);
         OpenMM_CustomIntegrator_addGlobalVariable (integrator, "temp", 0.0);
         OpenMM_CustomIntegrator_addGlobalVariable (integrator, "scale", 0.0);
         OpenMM_CustomIntegrator_addGlobalVariable (integrator, "c", exp (-(*dt)/ *bath__.tautemp));
         OpenMM_CustomIntegrator_addGlobalVariable (integrator, "d", 0.0);
         OpenMM_CustomIntegrator_addGlobalVariable (integrator, "r", 0.0);
         OpenMM_CustomIntegrator_addGlobalVariable (integrator, "s", 0.0);
         OpenMM_CustomIntegrator_addPerDofVariable (integrator, "si", 0.0);

         OpenMM_CustomIntegrator_addComputeSum (integrator, "ke", "m*v*v/2.0");
         OpenMM_CustomIntegrator_addComputeGlobal (integrator, "temp",
             "2.0*ke/(df*gasconst)");
         OpenMM_CustomIntegrator_addComputeGlobal (integrator, "temp",
             "temp + step(0.0000000001-temp) * 0.1");
         OpenMM_CustomIntegrator_addComputeGlobal (integrator, "c", "exp(-dt/tautemp)");
         OpenMM_CustomIntegrator_addComputeGlobal (integrator, "d", "(1.0-c)*(kelvin/temp)/df");
         OpenMM_CustomIntegrator_addComputeGlobal (integrator, "r", "gaussian");
         OpenMM_CustomIntegrator_addComputePerDof (integrator, "si", "gaussian");
         OpenMM_CustomIntegrator_addComputeSum (integrator, "s", "si*si");
         OpenMM_CustomIntegrator_addComputeGlobal (integrator, "s", e3);
         OpenMM_CustomIntegrator_addComputeGlobal (integrator, "scale",
              "sqrt(c+(s+r*r)*d+2.0*r*sqrt(c*d))");
         OpenMM_CustomIntegrator_addComputeGlobal (integrator, "scale",
             "scale-2.0*scale*step(-r-sqrt(c/d))");
         OpenMM_CustomIntegrator_addComputePerDof (integrator, "v", "scale*v");

      } else if (strncasecmp (buffer2, "ANDERSEN", 8) == 0) {
         setupAndersenThermostat (omm->system, log);
      }

      free (adjustment);
      free (e3);

   } else if (strncasecmp (buffer, "VERLET", 6) == 0) {
      omm->integrator = (OpenMM_Integrator*)OpenMM_VerletIntegrator_create (*dt);
      if (strncasecmp (buffer2, "ANDERSEN", 8) == 0) {
         setupAndersenThermostat (omm->system, log);
      }

   } else if (strncasecmp (buffer, "STOCHASTIC", 10) == 0) {

      if (log) {
         (void) fprintf (log, "\n Stochastic Integrator:\n");
         (void) fprintf (log, "\n Temperature          %15.4f K",
                             bath__.kelvin );
         (void) fprintf (log, "\n Friction             %15.4f ps^(-1)",
                             stodyn__.friction );
         (void) fprintf (log, "\n TimeStep             %15.4f ps\n",
                             *dt );
         (void) fflush (log);
      }
      omm->integrator = (OpenMM_Integrator*)OpenMM_LangevinIntegrator_create
                             (*bath__.kelvin, *stodyn__.friction, *dt);

   } else {
      (void) fprintf (stderr, "\n Integrator %s is Not Supported\n", buffer);
      (void) fflush (stderr);
      exit (-1);
   }

   // choose either the reference or the CUDA platform

   // platform = getReferencePlatform (log);
   platform = getCUDAPlatform (log);
   if (platform == NULL) {
      exit (-1);
   }

   // modification of context creation to avoid bug on large systems

   omm->context = OpenMM_Context_create_2 (omm->system, omm->integrator,
                                           platform);
   /*
   OpenMM_PropertyArray* properties = OpenMM_PropertyArray_create ();
   OpenMM_PropertyArray_add (properties, "DisablePmeStream", "true");
   omm->context = OpenMM_Context_create_3 (omm->system, omm->integrator,
                                           platform, properties);
   OpenMM_PropertyArray_destroy (properties);
   */
   OpenMM_Platform_setPropertyValue(platform, omm->context, "DisablePmeStream", "true");


   if (*inform__.debug) {
      (void) fprintf (log, "\n OpenMMDataHandle:  %x\n", (void*)(omm));
      (void) fprintf (log, "\n Integrator:  %x\n", (void*)(omm->integrator));
   }

   OpenMM_Context_setPositions (omm->context, initialPosInNm);
   OpenMM_Context_setVelocities (omm->context, initialVelInNm);

   if (*inform__.verbose) {
      int arraySz;
      int maxPrint;
      double x1, x2, x3;
      double v1, v2, v3;
      arraySz = OpenMM_Vec3Array_getSize (initialPosInNm);
      maxPrint = 5;
      (void) fprintf (log, "\n Initial Positions and Velocities :\n\n");
      for (ii = 0; ii < arraySz; ii++) {
         x1 = (*OpenMM_Vec3Array_get(initialPosInNm, ii)).x
                   / OpenMM_NmPerAngstrom;
         x2 = (*OpenMM_Vec3Array_get(initialPosInNm, ii)).y
                   / OpenMM_NmPerAngstrom;
         x3 = (*OpenMM_Vec3Array_get(initialPosInNm, ii)).z
                   / OpenMM_NmPerAngstrom;
         v1 = (*OpenMM_Vec3Array_get(initialVelInNm, ii)).x
                   / OpenMM_NmPerAngstrom;
         v2 = (*OpenMM_Vec3Array_get(initialVelInNm, ii)).y
                   / OpenMM_NmPerAngstrom;
         v3 = (*OpenMM_Vec3Array_get(initialVelInNm, ii)).z
                  / OpenMM_NmPerAngstrom;
         (void) fprintf (log, "%7d   POS    %17.8e %17.8e %17.8e\n",
                         ii+1, x1, x2, x3);
         (void) fprintf (log, "%7d   VEL    %17.8e %17.8e %17.8e\n",
                         ii+1, v1, v2, v3);
         if (ii == maxPrint-1 && ii < (arraySz-maxPrint-1))
             ii = arraySz - maxPrint - 1;
      }
   }

   *ommHandle = (void*) omm;
}

/*
 *    ############################################################
 *       Update Tinker Data Structures; Call mdstat and mdsave
 *    ############################################################
 *
 *    @param omm          handle containing OpenMM data structures
 *    @param dt           simulation time step in ps
 *    @param istep        current step in MD loop
 *    @param callMdStat   if nonzero, call Tinker mdstat routine
 *    @param callMdSave   if nonzero, call Tinker mdsave routine
 */

void openmm_update_ (void** omm, double* dt, int* istep,
                     int* callMdStat, int* callMdSave) {

   OpenMM_State* state;
   const OpenMM_Vec3Array* positionArray;
   const OpenMM_Vec3Array* velocityArray;
   const OpenMM_Vec3Array* forceArray;
   OpenMM_Vec3 aBox;
   OpenMM_Vec3 bBox;
   OpenMM_Vec3 cBox;
   int infoMask;
   int ii;
   double amass;
   double positionConvert;
   double velocityConvert;
   double forceConvert;

   static const double gasconst = units__.gasconst;

   double totalEnergy, potentialEnergy, kineticEnergy;

   OpenMMData* openMMDataHandle;

   openMMDataHandle = (OpenMMData*) (*omm);

   infoMask = OpenMM_State_Positions;
   infoMask += OpenMM_State_Velocities;
   infoMask += OpenMM_State_Forces;
   infoMask += OpenMM_State_Energy;

   // state object is created here and must be destroyed below

   state = OpenMM_Context_getState (openMMDataHandle->context, infoMask, 0);
   OpenMM_State_getPeriodicBoxVectors (state, &aBox, &bBox, &cBox);

   *(boxes__.xbox) = aBox.x / OpenMM_NmPerAngstrom;
   *(boxes__.ybox) = bBox.y / OpenMM_NmPerAngstrom;
   *(boxes__.zbox) = cBox.z / OpenMM_NmPerAngstrom;
   *(boxes__.xbox2) = 0.5 * (*(boxes__.xbox));
   *(boxes__.ybox2) = 0.5 * (*(boxes__.ybox));
   *(boxes__.zbox2) = 0.5 * (*(boxes__.zbox));

   // fprintf (stderr, "openmm_update_ %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e \n",
   // *(boxes__.xbox), *(boxes__.ybox), *(boxes__.zbox),
   // aBox.x, bBox.y, cBox.z );

   // load positions/velocities and energies

   positionConvert = 1.0 / OpenMM_NmPerAngstrom;
   velocityConvert = 1.0 / OpenMM_NmPerAngstrom;
   forceConvert = 10.0;

   positionArray = OpenMM_State_getPositions (state);
   velocityArray = OpenMM_State_getVelocities (state);
   forceArray = OpenMM_State_getForces (state);

   for (ii = 0; ii < *atoms__.n; ii++) {
      atoms__.x[ii] = (*OpenMM_Vec3Array_get(positionArray, ii)).x
                          * positionConvert;
      atoms__.y[ii] = (*OpenMM_Vec3Array_get(positionArray, ii)).y
                          * positionConvert;
      atoms__.z[ii] = (*OpenMM_Vec3Array_get(positionArray, ii)).z
                          * positionConvert;
   }

   for (ii = 0; ii < *atoms__.n; ii++) {
      int offset = ii*3;
      moldyn__.v[offset] = (*OpenMM_Vec3Array_get(velocityArray, ii)).x
                                * velocityConvert;
      moldyn__.v[offset+1] = (*OpenMM_Vec3Array_get(velocityArray, ii)).y
                                  * velocityConvert;
      moldyn__.v[offset+2] = (*OpenMM_Vec3Array_get(velocityArray, ii)).z
                                  * velocityConvert;
   }

   for (ii = 0; ii < *atoms__.n; ii++) {
      int offset = ii*3;
      amass = 1.0 / atomid__.mass[ii];
      moldyn__.a[offset] = (*OpenMM_Vec3Array_get(forceArray, ii)).x
                                * amass * forceConvert;
      moldyn__.a[offset+1] = (*OpenMM_Vec3Array_get(forceArray, ii)).y
                                  * amass * forceConvert;
      moldyn__.a[offset+2] = (*OpenMM_Vec3Array_get(forceArray, ii)).z
                                  * amass * forceConvert;
   }

   for (ii = 0; ii < *atoms__.n; ii++) {
      int offset = ii*3;
      moldyn__.aalt[offset] = 0.0;
      moldyn__.aalt[offset+1] = 0.0;
      moldyn__.aalt[offset+2] = 0.0;
   }

   potentialEnergy = (OpenMM_State_getPotentialEnergy(state))
                         * OpenMM_KcalPerKJ;
   kineticEnergy = (OpenMM_State_getKineticEnergy(state))
                       * OpenMM_KcalPerKJ;
   totalEnergy = potentialEnergy + kineticEnergy;

   if (inform__.debug) {

      double eksum, temp, pres;
      double ekin[3][3];

      kinetic_ (&eksum, &ekin, &temp);

      (void) fprintf (stderr, "\n State: E=%15.7e [%15.7e %15.7e]\n",
                      totalEnergy, potentialEnergy, kineticEnergy);
      (void) fprintf (stderr, "        t=%15.7e ps T=%12.5e InitT=%12.3f\n\n",
                      (*dt)*(*istep), temp, bath__.kelvin );

      for (ii = 0; ii < *atoms__.n; ii++) {
         (void) fprintf (stderr, "%7d   POS   %17.8e %17.8e %17.8e\n", ii+1,
             (*OpenMM_Vec3Array_get(positionArray,ii)).x*positionConvert,
             (*OpenMM_Vec3Array_get(positionArray,ii)).y*positionConvert,
             (*OpenMM_Vec3Array_get(positionArray,ii)).z*positionConvert);
         (void) fprintf (stderr, "%7d   VEL   %17.8e %17.8e %17.8e\n", ii+1,
             (*OpenMM_Vec3Array_get(velocityArray,ii)).x*velocityConvert,
             (*OpenMM_Vec3Array_get(velocityArray,ii)).y*velocityConvert,
             (*OpenMM_Vec3Array_get(velocityArray,ii)).z*velocityConvert);
         (void) fprintf (stderr, "%7d   FRC   %17.8e %17.8e %17.8e\n", ii+1,
             (*OpenMM_Vec3Array_get(forceArray,ii)).x*forceConvert,
             (*OpenMM_Vec3Array_get(forceArray,ii)).y*forceConvert,
             (*OpenMM_Vec3Array_get(forceArray,ii)).z*forceConvert);
      }
      (void) fflush (stderr);
   }

   // make calls to mdstat and/or mdsave if flags are set
   //
   // note: call to "bounds" below enforces periodic boundaries on output
   // coordinates as in canonical Tinker, but can break center-of-mass
   // restraints and perhaps other things computed via OpenMM, so it is
   // advised to allow the molecules to diffuse out of the periodic box,
   // and then "wrap" the MD trajectory after the simulation if desired

   if (*callMdStat || *callMdSave) {
      double eksum, temp, pres;
      double ekin[3][3];
      kinetic_ (&eksum, &ekin, &temp);

      if (*callMdStat) {
         pres = 0.0;
         if (bath__.isobaric) {
            lattice_ ();
         }
         mdstat_ (istep,dt,&totalEnergy,&potentialEnergy,&eksum,&temp,&pres);
      }
      if (*callMdSave) {
         bounds_ ();
         mdsave_ (istep,dt,&potentialEnergy,&eksum);
      }
   }
   OpenMM_State_destroy (state);
}

void openmm_take_steps_ (void** omm, int* numSteps) {

   OpenMMData* openMMDataHandle = (OpenMMData*) (*omm);
   OpenMM_Integrator_step (openMMDataHandle->integrator, *numSteps);
}

void openmm_cleanup_ (void** omm) {

   // clean up top-level heap allocated objects we are done with

   OpenMMData* openMMDataHandle = (OpenMMData*) (*omm);
   OpenMM_Context_destroy (openMMDataHandle->context);
   OpenMM_Integrator_destroy (openMMDataHandle->integrator);
   OpenMM_System_destroy (openMMDataHandle->system);
   free (openMMDataHandle);
}

}

static void zeroTinkerForce (double* tinkerForce) {

   int ii;
   for (ii = 0; ii < *atoms__.n; ii++){
      *(tinkerForce + 3*ii + 0) = 0.0;
      *(tinkerForce + 3*ii + 1) = 0.0;
      *(tinkerForce + 3*ii + 2) = 0.0;
   }
}

static void zeroVec3Force (OpenMM_Vec3Array* tinkerForce) {

   int ii;
   for (ii = 0; ii < *atoms__.n; ii++) {
      OpenMM_Vec3 force;
      force.x = force.y = force.z = 0.0;
      OpenMM_Vec3Array_set (tinkerForce, ii, force);
   }
}

static void setTinker1DArray (int size, double* tinkerArray, double value) {

   int ii;
   for (ii = 0; ii < size; ii++) {
      tinkerArray[ii] = value;
   }
}

static void loadTinkerForce (double* tinkerForce, int add,
                             OpenMM_Vec3Array* arrayToLoad) {

   int ii;
   if (add == 0) {
      for (ii = 0; ii < *atoms__.n; ii++) {
         OpenMM_Vec3 force;
         force.x = *(tinkerForce + 3*ii);
         force.y = *(tinkerForce + 3*ii + 1);
         force.z = *(tinkerForce + 3*ii + 2);
         OpenMM_Vec3Array_set (arrayToLoad, ii, force);
      }
   } else {
      double factor = (double) add;
      for (ii = 0; ii < *atoms__.n; ii++) {
         OpenMM_Vec3 force;
         const OpenMM_Vec3* currentForce = OpenMM_Vec3Array_get
                                               (arrayToLoad, ii);
         force.x = currentForce->x + factor * (*(tinkerForce + 3*ii));
         force.y = currentForce->y + factor * (*(tinkerForce + 3*ii + 1));
         force.z = currentForce->z + factor * (*(tinkerForce + 3*ii + 2));
         OpenMM_Vec3Array_set (arrayToLoad, ii, force);
      }
   }
}

static int usingImplicitSolvent (void) {

   int implicitSolventActive = -1;
   char solvatationType[16];
   char bornType[16];

   setNullTerminator (solute__.solvtyp, 8, solvatationType);
   setNullTerminator (solute__.borntyp, 8, bornType);

   // return <0 if parameter/option combination is unsupported
   //         0 if explicit solvent (Ewald is in use)
   //         1 if GK implicit solvent via OBC
   //         2 if GK implicit solvent via Grycuk

   if (strncasecmp (solvatationType, "GK", 2) == 0 &&
            ((strncasecmp (bornType, "OBC", 3) == 0) ||
             (strncasecmp (bornType, "GRYCUK", 6) == 0))) {
      if (limits__.use_ewald) {
         implicitSolventActive = -2;
      } else {
         if ((strncasecmp( bornType, "OBC", 3) == 0)) {
             implicitSolventActive = 1;
         } else {
             implicitSolventActive = 2;
         }
      }
   } else if (limits__.use_ewald) {
      implicitSolventActive = 0;
   }
   return implicitSolventActive;
}

/*
 *    ############################################################
 *           Compare Tinker and OpenMM Energies and Forces
 *    ############################################################
 */

extern "C" {

void openmm_test_ (void) {

   OpenMM_Vec3Array* initialPosInNm;
   OpenMM_Platform* platform;
   OpenMM_System* system;
   OpenMM_Context* context;
   OpenMM_Integrator* integrator;
   OpenMM_State* state;

   const OpenMM_Vec3Array* posArrayInNm;
   const OpenMM_Vec3Array* openMMForces;

   int infoMask;
   int ii, jj;
   int countActiveForces;
   char const* testName;
   double conversion, delta, dot;
   double tinkerNorm, openMMNorm;
   double openMMPotentialEnergy;
   OpenMM_Vec3Array* tinkerForce;
   double tinkerEnergy;
   double maxEDelta, maxFDelta;
   double minDot, avgDot;
   int maxFDeltaIndex, minDotIndex;
   int implicitSolventActive;

   int removeConstrainedCovalentIxns = 0;
   FILE* log = stderr;

   if (log) {
      (void) fprintf (log, "\n Testing Tinker vs OpenMM Energy & Force :\n");
   }

   implicitSolventActive = usingImplicitSolvent ();

   tinkerForce = OpenMM_Vec3Array_create (*atoms__.n);

   // Create a System and Force objects within the System. Retain a reference
   // to each force object so we can fill in the forces. Note: the OpenMM
   // System takes ownership of the force objects; don't delete them yourself.

   testName = NULL;
   system = OpenMM_System_create ();
   setupSystemParticles (system, log);

   countActiveForces = 0;
   if (*potent__.use_bond)  countActiveForces++;
   if (*potent__.use_angle)  countActiveForces++;
   if (*potent__.use_strbnd)  countActiveForces++;
   if (*potent__.use_urey)  countActiveForces++;
   if (*potent__.use_angang)  countActiveForces++;
   if (*potent__.use_opbend)  countActiveForces++;
   if (*potent__.use_opdist)  countActiveForces++;
   if (*potent__.use_improp)  countActiveForces++;
   if (*potent__.use_imptor)  countActiveForces++;
   if (*potent__.use_tors)  countActiveForces++;
   if (*potent__.use_pitors)  countActiveForces++;
   if (*potent__.use_strtor)  countActiveForces++;
   if (*potent__.use_angtor)  countActiveForces++;
   if (*potent__.use_tortor)  countActiveForces++;
   if (*potent__.use_vdw)  countActiveForces++;
   if (*potent__.use_charge)  countActiveForces++;
   if (*potent__.use_chgdpl)  countActiveForces++;
   if (*potent__.use_dipole)  countActiveForces++;
   if (*potent__.use_chgtrn)  countActiveForces++;
   if (*potent__.use_repuls)  countActiveForces++;
   if (*potent__.use_disp)  countActiveForces++;
   if (*potent__.use_mpole)  countActiveForces++;
   if (*potent__.use_polar)  countActiveForces++;
   if (*potent__.use_rxnfld)  countActiveForces++;
   if (*potent__.use_solv)  countActiveForces++;
   if (*potent__.use_metal)  countActiveForces++;
   if (*potent__.use_geom)  countActiveForces++;

   if (log) {
      (void) fprintf (log, "\n Potential Terms Used in Tinker :\n" );
      (void) fprintf (log, "\n    Bond=    %d", abs(*potent__.use_bond));
      (void) fprintf (log, "    Angle=   %d", abs(*potent__.use_angle));
      (void) fprintf (log, "    StrBnd=  %d", abs(*potent__.use_strbnd));
      (void) fprintf (log, "    Urey=    %d", abs(*potent__.use_urey));
      (void) fprintf (log, "    AngAng=  %d", abs(*potent__.use_angang));
      (void) fprintf (log, "\n    OPBend=  %d", abs(*potent__.use_opbend));
      (void) fprintf (log, "    OPDist=  %d", abs(*potent__.use_opdist));
      (void) fprintf (log, "    ImProp=  %d", abs(*potent__.use_improp));
      (void) fprintf (log, "    ImpTor=  %d", abs(*potent__.use_imptor));
      (void) fprintf (log, "    Tors=    %d", abs(*potent__.use_tors));
      (void) fprintf (log, "\n    PiTors=  %d", abs(*potent__.use_pitors));
      (void) fprintf (log, "    StrTor=  %d", abs(*potent__.use_strtor));
      (void) fprintf (log, "    AngTor=  %d", abs(*potent__.use_angtor));
      (void) fprintf (log, "    TorTor=  %d", abs(*potent__.use_tortor));
      (void) fprintf (log, "    Vdw=     %d", abs(*potent__.use_vdw));
      (void) fprintf (log, "\n    Charge=  %d", abs(*potent__.use_charge));
      (void) fprintf (log, "    ChgDpl=  %d", abs(*potent__.use_chgdpl));
      (void) fprintf (log, "    Dipole=  %d", abs(*potent__.use_dipole));
      (void) fprintf (log, "    MPole=   %d", abs(*potent__.use_mpole));
      (void) fprintf (log, "    Polar=   %d", abs(*potent__.use_polar));
      (void) fprintf (log, "\n    RxnFld=  %d", abs(*potent__.use_rxnfld));
      (void) fprintf (log, "    Solv=    %d", abs(*potent__.use_solv));
      (void) fprintf (log, "    LigFld=  %d", abs(*potent__.use_metal));
      (void) fprintf (log, "    Restrn=  %d", abs(*potent__.use_geom));
      (void) fprintf (log, "    Extra=   %d\n", abs(*potent__.use_extra));
      (void) fprintf (log, "\n    Repel=   %d\n", abs(*potent__.use_repuls));
      (void) fprintf (log, "    Disp=    %d\n", abs(*potent__.use_disp));  
      (void) fprintf (log, "    Chgtrn=   %d\n", abs(*potent__.use_chgtrn));
      (void) fprintf (log, "    Chgpen=    %d\n", abs(*mplpot__.use_chgpen));    
   }

   if (countActiveForces > 1) {

      if (*potent__.use_bond) {
         setupAmoebaBondForce (system, removeConstrainedCovalentIxns, log);
      }
      if (*potent__.use_angle) {
         setupAmoebaAngleForce (system, removeConstrainedCovalentIxns, log);
         setupAmoebaInPlaneAngleForce (system, removeConstrainedCovalentIxns,
                                       log);
      }
      if (*potent__.use_strbnd) {
         setupAmoebaStretchBendForce (system, removeConstrainedCovalentIxns,
                                      log);
      }
      if (*potent__.use_urey) {
         setupAmoebaUreyBradleyForce (system, removeConstrainedCovalentIxns,
                                      log);
      }
      if (*potent__.use_opbend) {
         setupAmoebaOutOfPlaneBendForce (system, log);
      }
      if (*potent__.use_imptor) {
         setupAmoebaImproperTorsionForce (system, log);
      }
      if (*potent__.use_tors) {
         setupAmoebaTorsionForce (system, log);
      }
      if (*potent__.use_pitors) {
         setupAmoebaPiTorsionForce (system, log);
      }
      if (*potent__.use_tortor) {
         setupAmoebaTorsionTorsionForce (system, log);
      }
      
      if (*potent__.use_charge) {
         setupAmoebaChargeForce (system, log);
      }
      
      if (*potent__.use_chgtrn || *potent__.use_repuls || *potent__.use_disp) {
         setupHippoNonbondedForce (system, log);
      }
      
      if (*potent__.use_mpole && ! *mplpot__.use_chgpen) {
         setupAmoebaMultipoleForce (system, log);
      }
      
      if (*potent__.use_solv) {
         setupAmoebaWcaDispersionForce (system, log);
         setupAmoebaGeneralizedKirkwoodForce (system, 1, log);
      }
      
      if (*potent__.use_geom) {
         setupTorsionRestraints (system, log);
         setupDistanceRestraints (system, log);
         setupPositionalRestraints (system, log);
         setupAngleRestraints (system, log);
         setupCentroidRestraints (system, log);
      }

      loadTinkerForce (deriv__.desum, 0, tinkerForce);
      tinkerEnergy = *energi__.esum;
      testName = "PotentialsTest";

      if (log) {
         (void) fprintf (log,
                    "\n Potential Energy Components from Tinker :\n\n"
                    "    EB=  %15.7e   EA=  %15.7e   EBA= %15.7e\n"
                    "    EUB= %15.7e   EAA= %15.7e   EOPB=%15.7e\n"
                    "    EOPD=%15.7e   EID= %15.7e   EIT= %15.7e\n"
                    "    ET=  %15.7e   EPT= %15.7e   EBT= %15.7e\n"
                    "    EAT= %15.7e   ETT= %15.7e   EV=  %15.7e\n"
                    "    EC=  %15.7e   ECD= %15.7e   ED=  %15.7e\n"
                    "    EM=  %15.7e   EP=  %15.7e   ER=  %15.7e\n"
                    "    ES=  %15.7e   ELF= %15.7e   EG=  %15.7e\n"
                    "    EDSP=%15.7e   ECT= %15.7e   EX=  %15.7e\n",
                    *energi__.eb,   *energi__.ea,  *energi__.eba,
                    *energi__.eub,  *energi__.eaa, *energi__.eopb,
                    *energi__.eopd, *energi__.eid, *energi__.eit,
                    *energi__.et,   *energi__.ept, *energi__.ebt,
                    *energi__.eat,  *energi__.ett, *energi__.ev,
                    *energi__.ec,   *energi__.ecd, *energi__.ed,
                    *energi__.em,   *energi__.ep,  *energi__.er,
                    *energi__.es,   *energi__.elf, *energi__.eg,
                    *energi__.edsp, *energi__.ect, *energi__.ex );
                    
         (void) fflush (log);
      }

   } else if (*potent__.use_bond) {

      setupAmoebaBondForce (system, removeConstrainedCovalentIxns, log);
      loadTinkerForce (deriv__.deb, 0, tinkerForce);
      tinkerEnergy = *energi__.eb;
      testName = "AmoebaHarmonicBondTest";

   } else if (*potent__.use_angle) {

      // note Tinker angle = OpenMM (Angle + InPlaneAngle)

      setupAmoebaAngleForce (system, removeConstrainedCovalentIxns, log);
      setupAmoebaInPlaneAngleForce (system, removeConstrainedCovalentIxns,
                                    log);
      loadTinkerForce (deriv__.dea, 0, tinkerForce);
      tinkerEnergy = *energi__.ea;
      testName = "AmoebaHarmonicAngleTest";

   } else if (*potent__.use_strbnd) {

      setupAmoebaStretchBendForce (system, removeConstrainedCovalentIxns,
                                   log);
      loadTinkerForce (deriv__.deba, 0, tinkerForce);
      tinkerEnergy = *energi__.eba;
      testName = "AmoebaStretchBendTest";

   } else if (*potent__.use_urey) {

      setupAmoebaUreyBradleyForce (system, removeConstrainedCovalentIxns,
                                   log);
      loadTinkerForce (deriv__.deub, 0, tinkerForce);
      tinkerEnergy = *energi__.eub;
      testName = "AmoebaUreyBradleyForceTest";

   } else if (*potent__.use_opbend) {

      setupAmoebaOutOfPlaneBendForce (system, log);
      loadTinkerForce (deriv__.deopb, 0, tinkerForce);
      tinkerEnergy = *energi__.eopb;
      testName = "AmoebaOutOfPlaneBendTest";

   } else if (*potent__.use_imptor) {

      setupAmoebaImproperTorsionForce (system, log);
      loadTinkerForce (deriv__.deit, 0, tinkerForce);
      tinkerEnergy = *energi__.eit;
      testName = "AmoebaImproperTorsionForce";

   } else if (*potent__.use_tors) {

      setupAmoebaTorsionForce (system, log);
      loadTinkerForce (deriv__.det, 0, tinkerForce);
      tinkerEnergy = *energi__.et;
      testName = "AmoebaTorsionTest";

   } else if (*potent__.use_pitors) {

      setupAmoebaPiTorsionForce (system, log);
      loadTinkerForce (deriv__.dept, 0, tinkerForce);
      tinkerEnergy = *energi__.ept;
      testName = "AmoebaPiTorsionTest";

   } else if (*potent__.use_geom) {

      setupTorsionRestraints (system, log);
      setupDistanceRestraints (system, log);
      setupPositionalRestraints (system, log);
      setupAngleRestraints (system, log);
      setupCentroidRestraints (system, log);
      loadTinkerForce (deriv__.deg, 0, tinkerForce);
      tinkerEnergy = *energi__.eg;
      testName = "AmoebaRestraintTest";

   } else if (*potent__.use_charge) {

      setupAmoebaChargeForce (system, log);
      loadTinkerForce (deriv__.dec, 0, tinkerForce);
      tinkerEnergy = *energi__.ec;
      testName = "AmoebaChargeTest";

   } else if (*potent__.use_chgtrn || *potent__.use_repuls || *potent__.use_disp) {

      setupHippoNonbondedForce (system, log);
      loadTinkerForce (deriv__.dect, 0, tinkerForce);
      tinkerEnergy = *energi__.ect;
      testName = "HippoChargeTransferTest";

   } else if (*potent__.use_mpole && ! *potent__.use_solv && ! *mplpot__.use_chgpen) {

      if (log) {
         (void) fprintf (log, "ImplicitSolventActive=%d\n",
                         implicitSolventActive);
      }
      if (implicitSolventActive == 0) {
         setupAmoebaBondForce (system, removeConstrainedCovalentIxns, log);
         loadTinkerForce (deriv__.deb, 0, tinkerForce);
         tinkerEnergy = *energi__.eb;
      } else {
         zeroVec3Force (tinkerForce);
         tinkerEnergy = 0.0;
      }
      setupAmoebaMultipoleForce (system, log);
      loadTinkerForce (deriv__.dem, 1, tinkerForce);
      loadTinkerForce (deriv__.dep, 1, tinkerForce);
      tinkerEnergy += *energi__.em + *energi__.ep;

      if (strncasecmp (polpot__.poltyp, "DIRECT", 6) == 0) {
         testName = "AmoebaMultipoleDirectTest";
      } else {
         testName = "AmoebaMultipoleMutualTest";
      }

   } else if (*potent__.use_solv && implicitSolventActive > 0 &&
              ! *potent__.use_mpole) {

      // to get Tinker WCA, zero deriv__.des, then call ewca1

      zeroTinkerForce (deriv__.des);
      ewca1_ (&tinkerEnergy);
      loadTinkerForce (deriv__.des, 0, tinkerForce);
      setupAmoebaWcaDispersionForce (system, log);
      testName = "AmoebaWcaDispersionTest";

   } else if (implicitSolventActive > 0 && *potent__.use_mpole) {

      // generalized Kirkwood; OpenMM should have cavity term turned off

      double ecav, edisp;
      if (log) {
         char buffer[128];
         setNullTerminator (solute__.borntyp, 8, buffer);
         (void) fprintf (log, "Born radius type=%s ", buffer);
         setNullTerminator (solute__.solvtyp, 8, buffer);
         (void) fprintf (log, "Solvation type=%s\n", buffer);
      }
      setupAmoebaMultipoleForce (system, log);
      setupAmoebaGeneralizedKirkwoodForce (system, 1, log);
      zeroTinkerForce (deriv__.des);
      setTinker1DArray (*atoms__.n, solute__.drb, 0.0);
      setTinker1DArray (*atoms__.n, solute__.drbp, 0.0);

      *energi__.es = 0.0;
      born_ ();
      empole1_ ();
      egk1_ ();

      loadTinkerForce (deriv__.des, 0, tinkerForce);
      loadTinkerForce (deriv__.dem, 1, tinkerForce);
      loadTinkerForce (deriv__.dep, 1, tinkerForce);
      tinkerEnergy = *energi__.es + *energi__.em + *energi__.ep;

      enp1_ (&ecav,&edisp);

      if (log) {
         (void) fprintf (log, "Energies total=%15.7e Gk=%15.7e (cavity=%15.7e dispersion=%15.7e) Multipole=%15.7e\n",
                         tinkerEnergy, *energi__.es, ecav, edisp,
                         *energi__.em + *energi__.ep );

         //for (ii = 0; ii < atoms__.n; ii++) {
         //    fprintf (stderr, "Fs %5d [%15.7e %15.7e %15.7e] [%15.7e %15.7e%15.7e] [%15.7e %15.7e %15.7e]\n",
         //             ii, *(deriv__.des + 3*ii), *(deriv__.des + 3*ii+1),
         //             *(deriv__.des + 3*ii+2), *(deriv__.dem + 3*ii),
         //             *(deriv__.dem + 3*ii+1), *(deriv__.dem + 3*ii+2),
         //             *(deriv__.dep + 3*ii), *(deriv__.dep + 3*ii+1),
         //             *(deriv__.dep + 3*ii+2));
         //}

      }
      if (strncasecmp (polpot__.poltyp, "DIRECT", 6) == 0) {
         testName = "AmoebaKirkwoodDirectTest";
      } else {
         testName = "AmoebaKirkwoodMutualTest";
      }
   }

   if (log) {
      if (testName) {
         (void) fprintf (log, "\n Test Option :  %s\n", testName);
         (void) fflush (NULL);
      } else {
         (void) fprintf (log, "\n Test Option not Recognized; Exiting\n");
         (void) fflush (log);
         exit (-1);
      }
   }

   initialPosInNm = OpenMM_Vec3Array_create (0);
   setupPositions (initialPosInNm, log);

   integrator = (OpenMM_Integrator*) OpenMM_VerletIntegrator_create (0.001);

   // choose either the reference or the CUDA platform

   // platform = getReferencePlatform (log);
   platform = getCUDAPlatform (log);
   if (platform == NULL) {
      exit (-1);
   }

   // modification of context creation to avoid bug on large systems

   context = OpenMM_Context_create_2 (system, integrator, platform);
   // OpenMM_PropertyArray* properties = OpenMM_PropertyArray_create ();
   // OpenMM_PropertyArray_add (properties, "DisablePmeStream", "true");
   // omm->context = OpenMM_Context_create_3 (omm->system, omm->integrator,
   //                                         platform, properties);
   // OpenMM_PropertyArray_destroy (properties);

   // Zhi's new fix to avoid using "PropertyArray"
   OpenMM_Platform_setPropertyValue(platform, context, "DisablePmeStream", "true");

   OpenMM_Context_setPositions (context, initialPosInNm);

   infoMask = OpenMM_State_Positions;
   infoMask += OpenMM_State_Forces;
   infoMask += OpenMM_State_Energy;
   
   state = OpenMM_Context_getState (context, infoMask, 0);
   
   openMMForces = OpenMM_State_getForces (state);
   openMMPotentialEnergy = OpenMM_State_getPotentialEnergy(state)
                              / OpenMM_KJPerKcal;

   conversion = -OpenMM_NmPerAngstrom / OpenMM_KJPerKcal;
   
   setDefaultPeriodicBoxVectors (system, log);
   if (log) {
      printDefaultPeriodicBoxVectors (system, log);
      (void) fflush (log);
   }
   
   // find differences in Tinker and OpenMM energies and forces

   maxFDelta = 0.0;
   maxFDeltaIndex = -1;
   minDot = 1.0;
   avgDot = 0.0;
   minDotIndex = -1;

   int maxPrint;
   maxPrint = 5;

   if (log) {

      (void) fprintf (log, "\n Tinker vs OpenMM Energy Values :\n");
      (void) fprintf (log, "\n Tinker Potential Energy   %18.4f",
                      tinkerEnergy);
      (void) fprintf (log, "\n OpenMM Potential Energy   %18.4f\n",
                      openMMPotentialEnergy);
      maxEDelta = fabs(tinkerEnergy - openMMPotentialEnergy);

      maxFDelta = 0.0;
      for (ii = 0; ii < *atoms__.n; ii++) {
         double relxNrm;
         double dot;
         OpenMM_Vec3 force;
         const OpenMM_Vec3* tinkerF;
         force.x = conversion*(*OpenMM_Vec3Array_get(openMMForces, ii)).x;
         force.y = conversion*(*OpenMM_Vec3Array_get(openMMForces, ii)).y;
         force.z = conversion*(*OpenMM_Vec3Array_get(openMMForces, ii)).z;
         tinkerF = OpenMM_Vec3Array_get(tinkerForce, ii);

         tinkerNorm = sqrt(tinkerF->x*tinkerF->x + tinkerF->y*tinkerF->y
                               + tinkerF->z*tinkerF->z);
         openMMNorm = sqrt(force.x*force.x + force.y*force.y
                               + force.z*force.z );

         delta = sqrt((tinkerF->x-force.x)*(tinkerF->x-force.x)
                    + (tinkerF->y-force.y)*(tinkerF->y-force.y)
                    + (tinkerF->z-force.z)*(tinkerF->z-force.z));
         dot = ((tinkerNorm > 0.0) && (openMMNorm > 0.0)) ?
                (tinkerF->x*force.x + tinkerF->y*force.y
                     + tinkerF->z*force.z) / (tinkerNorm*openMMNorm) : 0.0;

         if (delta > maxFDelta) {
            maxFDelta = delta;
            maxFDeltaIndex = ii + 1;
         }
         if (dot < minDot) {
            minDot = dot;
            minDotIndex = ii + 1;
         }
         avgDot += dot;

         if (ii == 0) {
            (void) fprintf (log, "\n Tinker vs OpenMM Force Values :\n\n");
         }
         if (ii < maxPrint || *atoms__.n - ii - 1 < maxPrint) {
            (void) fprintf (stderr, "%6d   Tinker %17.8e %17.8e %17.8e",
                            ii+1, tinkerF->x, tinkerF->y, tinkerF->z);
            (void) fprintf (stderr, "\n         OpenMM %17.8e %17.8e %17.8e\n",
                            force.x, force.y, force.z);
         }
      }

      if (*atoms__.n) {
         avgDot /= (double)(*atoms__.n);
      }

      (void) fprintf (log, "\n Summary of Tinker vs OpenMM Comparison :\n");
      (void) fprintf (log, "\n EnergyDelta                    %19.8e\n",
                      maxEDelta);
      (void) fprintf (log, " MaxForceDelta at   %11d %19.8e\n",
                      maxFDeltaIndex, maxFDelta);
      (void) fprintf (log, " MinDotProduct at   %11d %19.8e\n",
                      minDotIndex, minDot);
      (void) fprintf (log, " AvgDotProduct                  %19.8e\n",
                      avgDot);
   }

   OpenMM_Vec3Array_destroy (tinkerForce);
   OpenMM_State_destroy (state);
   OpenMM_Context_destroy (context);
   OpenMM_Integrator_destroy (integrator);
   OpenMM_System_destroy (system);
}

void openmm_bar_energy_ (void** ommHandle, double* energyInKcal) {
   OpenMMData_s* omm = (OpenMMData_s*) (*ommHandle);

   // copy periodic box from Tinker to OpenMM
   OpenMM_Vec3 a;
   OpenMM_Vec3 b;
   OpenMM_Vec3 c;
   a.x = boxes__.lvec[0] * OpenMM_NmPerAngstrom;
   a.y = boxes__.lvec[3] * OpenMM_NmPerAngstrom;
   a.z = boxes__.lvec[6] * OpenMM_NmPerAngstrom;
   b.x = boxes__.lvec[1] * OpenMM_NmPerAngstrom;
   b.y = boxes__.lvec[4] * OpenMM_NmPerAngstrom;
   b.z = boxes__.lvec[7] * OpenMM_NmPerAngstrom;
   c.x = boxes__.lvec[2] * OpenMM_NmPerAngstrom;
   c.y = boxes__.lvec[5] * OpenMM_NmPerAngstrom;
   c.z = boxes__.lvec[8] * OpenMM_NmPerAngstrom;
   OpenMM_Context_setPeriodicBoxVectors (omm->context, &a, &b, &c);

   // copy coordinates from Tinker to OpenMM
   // never call OpenMM_Vec3Array_destroy (posInMNm);
   static OpenMM_Vec3Array* posInNm = OpenMM_Vec3Array_create (*atoms__.n);
   for (int ii = 0; ii < *atoms__.n; ++ii) {
      OpenMM_Vec3 r;
      r.x = atoms__.x[ii] * OpenMM_NmPerAngstrom;
      r.y = atoms__.y[ii] * OpenMM_NmPerAngstrom;
      r.z = atoms__.z[ii] * OpenMM_NmPerAngstrom;
      OpenMM_Vec3Array_set (posInNm, ii, r);
   }
   OpenMM_Context_setPositions (omm->context, posInNm);

   // get OpenMM state
   int infoMask = 0;
   infoMask += OpenMM_State_Energy;
   OpenMM_State* state = OpenMM_Context_getState (omm->context, infoMask, 0);
   *energyInKcal = OpenMM_State_getPotentialEnergy (state) * OpenMM_KcalPerKJ;
}

}