/*
 *    ################################################################
 *    ##                     COPYRIGHT (C) 2015                     ##
 *    ##  by Mark Friedrichs, Lee-Ping Wang, Kailong Mao, Chao Lu,  ##
 *    ##       Zhi Wang, Matthew Harger and Jay William Ponder      ##
 *    ##                     All Rights Reserved                    ##
 *    ################################################################
 *
 *    ############################################################
 *    ##                                                        ##
 *    ##  ommstuff.cpp  --  Tinker interface to the OpenMM API  ##
 *    ##                                                        ##
 *    ############################################################
 */

/*
 *    ############################################################
 *                     System-Level Include Files
 *    ############################################################
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

const int MAX_STRING = 240;

/*
 *    ############################################################
 *           Tinker Routines Called Directly from Interface
 *        (to convert C to C++, must enclose in extern "C" {})
 *    ############################################################
 */

extern "C" {

   void born_ ();
   void bounds_ ();
   void egk1_ ();
   void empole1_ ();
   void enp1_ (double*, double*);
   void ewca1_ (double*);
   void kinetic_ (double*, double(*)[3][3], double*);
   void lattice_ ();
   void mdsave_ (int*, double*, double*, double*);
   void mdstat_ (int*, double*, double*, double*, double*, double*, double*);

}

/*
 *    ############################################################
 *                    OpenMM Wrapper Include Files
 *    ############################################################
 */

#include "OpenMMCWrapper.h"
#include "AmoebaOpenMMCWrapper.h"
#include "OpenMM.h"
using namespace OpenMM;

typedef struct OpenMMData_s OpenMMData;

struct OpenMMData_s {
   OpenMM_System* system;
   OpenMM_Context* context;
   OpenMM_Integrator* integrator;
};

typedef struct {
   char s20[20];
} char20;

/*
 *    ############################################################
 *            C++ Data Structures Corresponding to Tinker
 *    ############################################################
 */

struct {
   int nangle;
   int* iang;
   double* ak;
   double* anat;
   double* afld;
} angbnd__;

struct {
   double angunit;
   double stbnunit;
   double aaunit;
   double opbunit;
   double opdunit;
   double cang;
   double qang;
   double pang;
   double sang;
   double copb;
   double qopb;
   double popb;
   double sopb;
   double copd;
   double qopd;
   double popd;
   double sopd;
   char opbtyp[MAX_STRING];
   char* angtyp;
} angpot__;

static struct {
   int nangtor;
   int* iat;
   double* kant;
} angtor__;

struct {
   int* tag;
   int* classs;   // variable "class" not allowed in C++; add an extra "s"
   int* atomic;
   int* valence;
   double* mass;
   char* name;
   char* story;
} atomid__;

struct {
   int n;
   int* type;
   double* x;
   double* y;
   double* z;
} atoms__;

struct {
   int maxnose;
   int voltrial;
   double kelvin;
   double atmsph;
   double tautemp;
   double taupres;
   double compress;
   double collide;
   double eta;
   double volmove;
   double vbar;
   double qbar;
   double gbar;
   double* vnh;
   double* qnh;
   double* gnh;
   int isothermal;
   int isobaric;
   int anisotrop;
   char volscale[MAX_STRING];
   char barostat[MAX_STRING];
   char thermostat[MAX_STRING];
} bath__;

struct {
   int nbitor;
   int* ibitor;
} bitor__;

struct {
   double cbnd;
   double qbnd;
   double bndunit;
   char bndtyp[MAX_STRING];
} bndpot__;

struct {
   int nbond;
   int* ibnd;
   double* bk;
   double* bl;
} bndstr__;

struct {
   double polycut;
   double polycut2;
   int use_bounds;
   int use_replica;
   int use_polymer;
} bound__;

struct {
   double* xbox;
   double* ybox;
   double* zbox;
   double* alpha;
   double* beta;
   double* gamma;
   double* xbox2;
   double* ybox2;
   double* zbox2;
   double box34;
   double volbox;
   double beta_sin;
   double beta_cos;
   double gamma_sin;
   double gamma_cos;
   double beta_term;
   double gamma_term;
   double* lvec;
   double* recip;
   int orthogonal;
   int monoclinic;
   int triclinic;
   int octahedron;
   char spacegrp[MAX_STRING];
} boxes__;

struct {
   int ncell;
   int* icell;
   double xcell;
   double ycell;
   double zcell;
   double xcell2;
   double ycell2;
   double zcell2;
} cell__;

struct {
   int nion;
   int* iion;
   int* jion;
   int* kion;
   double* pchg;
} charge__;

struct {
   double electric;
   double dielec;
   double ebuffer;
   double c2scale;
   double c3scale;
   double c4scale;
   double c5scale;
   int neutnbr;
   int neutcut;
} chgpot__;

struct {
   int* n12;
   int* n13;
   int* n14;
   int* n15;
   int* i12;
   int* i13;
   int* i14;
   int* i15;
} couple__;

struct {
   double* desum;
   double* deb;
   double* dea;
   double* deba;
   double* deub;
   double* deaa;
   double* deopb;
   double* deopd;
   double* deid;
   double* deit;
   double* det;
   double* dept;
   double* debt;
   double* deat;
   double* dett;
   double* dev;
   double* der;
   double* dedsp;
   double* dec;
   double* decd;
   double* ded;
   double* dem;
   double* dep;
   double* dect;
   double* derxf;
   double* des;
   double* delf;
   double* deg;
   double* dex;
} deriv__;

struct {
   double* esum;
   double* eb;
   double* ea;
   double* eba;
   double* eub;
   double* eaa;
   double* eopb;
   double* eopd;
   double* eid;
   double* eit;
   double* et;
   double* ept;
   double* ebt;
   double* eat;
   double* ett;
   double* ev;
   double* er;
   double* edsp;
   double* ec;
   double* ecd;
   double* ed;
   double* em;
   double* ep;
   double* ect;
   double* erxf;
   double* es;
   double* elf;
   double* eg;
   double* ex;
} energi__;

struct {
   double aewald;
   double aeewald;
   double apewald;
   double adewald;
   char boundary[MAX_STRING];
} ewald__;

struct {
   int nrat;
   int nratx;
   int* iratx;
   int* kratx;
   int* irat;
   double rateps;
   double* krat;
   int use_rattle;
   int* ratimage;
} freeze__;

struct {
   int ngrp;
   int* kgrp;
   int* grplist;
   int* igrp;
   double* grpmass;
   double* wgrp;
   int use_group;
   int use_intra;
   int use_inter;
} group__;

struct {
   int nitors;
   int* iitors;
   double* itors1;
   double* itors2;
   double* itors3;
} imptor__;

struct {
   int maxask;
   int digits;
   int iprint;
   int iwrite;
   int isend;
   int silent;
   int verbose;
   int debug;
   int holdup;
   int abort;
} inform__;

struct {
   int maxntt;
   int maxtgrd;
   int maxtgrd2;
   int* tnx;
   int* tny;
   double* ttx;
   double* tty;
   double* tbf;
   double* tbx;
   double* tby;
   double* tbxy;
   char20* ktt;
} ktrtor__;

struct {
   int maxnvp;
   double* radpr;
   double* epspr;
   char* kvpr;
} kvdwpr__;

struct {
   double* rad;
   double* eps;
   double* rad4;
   double* eps4;
   double* reduct;
} kvdws__;

struct {
   double vdwcut;
   double repcut;
   double dispcut;
   double chgcut;
   double dplcut;
   double mpolecut;
   double ctrncut;
   double vdwtaper;
   double reptaper;
   double disptaper;
   double chgtaper;
   double dpltaper;
   double mpoletaper;
   double ctrntaper;
   double ewaldcut;
   double dewaldcut;
   double usolvcut;
   int use_ewald;
   int use_dewald;
   int use_lights;
   int use_list;
   int use_vlist;
   int use_dlist;
   int use_clist;
   int use_mlist;
   int use_ulist;
} limits__;

struct {
   int nfree;
   int irest;
   int bmnmix;
   double arespa;
   int dorest;
   char integrate[MAX_STRING];
} mdstuf__;

struct {
   int nmol;
   int* imol;
   int* kmol;
   int* molcule;
   double totmass;
   double* molmass;
} molcul__;

struct {
   double* v;
   double* a;
   double* aalt;
} moldyn__;

struct {
   double m2scale;
   double m3scale;
   double m4scale;
   double m5scale;
   int use_chgpen;
} mplpot__;

struct {
   int maxpole;
   int npole;
   int* ipole;
   int* polsiz;
   int* pollist;
   int* zaxis;
   int* xaxis;
   int* yaxis;
   double* pole;
   double* rpole;
   double* spole;
   double* srpole;
   char* polaxe;
} mpole__;

struct {
   int nmut;
   int vcouple;
   int* imut;
   int* type0;
   int* class0;
   int* type1;
   int* class1;
   double lambda;
   double tlambda;
   double vlambda;
   double elambda;
   double scexp;
   double scalpha;
   int* mut;
} mutant__;

struct {
   double epso;
   double epsh;
   double rmino;
   double rminh;
   double awater;
   double slevy;
   double solvprs;
   double surften;
   double spcut;
   double spoff;
   double stcut;
   double stoff;
   double* rcav;
   double* rdisp;
   double* cdisp;
} nonpol__;

struct {
   int nopbend;
   int* iopb;
   double* opbk;
} opbend__;

struct {
   void* ommHandle;
   char cudaPrecision[MAX_STRING];
   char ommPlatform[MAX_STRING];
   char cudaDevice[MAX_STRING];
} openmm__;

struct {
   int npitors;
   int* ipit;
   double* kpit;
} pitors__;

struct {
   int nfft1;
   int nfft2;
   int nfft3;
   int nefft1;
   int nefft2;
   int nefft3;
   int ndfft1;
   int ndfft2;
   int ndfft3;
   int bsorder;
   int bseorder;
   int bsporder;
   int bsdorder;
   int* igrid;
   double* bsmod1;
   double* bsmod2;
   double* bsmod3;
   double* bsbuild;
   double*** thetai1;
   double*** thetai2;
   double*** thetai3;
   double**** qgrid;
   double*** qfac;
} pme__;

struct {
   int npolar;
   int* ipolar;
   double* polarity;
   double* thole;
   double* dirdamp;
   double* pdamp;
   double* udir;
   double* udirp;
   double* udirs;
   double* udirps;
   double* uind;
   double* uinp;
   double* uinds;
   double* uinps;
   double* uexact;
   int* douind;
} polar__;

struct {
   int maxp11;
   int maxp12;
   int maxp13;
   int maxp14;
   int* np11;
   int* np12;
   int* np13;
   int* np14;
   int* ip11;
   int* ip12;
   int* ip13;
   int* ip14;
} polgrp__;

struct {
   int maxopt;
   int optorder;
   int optlevel;
   double* copt;
   double* copm;
   double* uopt;
   double* uoptp;
   double* uopts;
   double* uoptps;
   double* fopt;
   double* foptp;
} polopt__;

struct {
   int politer;
   double poleps;
   double p2scale;
   double p3scale;
   double p4scale;
   double p5scale;
   double p2iscale;
   double p3iscale;
   double p4iscale;
   double p5iscale;
   double d1scale;
   double d2scale;
   double d3scale;
   double d4scale;
   double u1scale;
   double u2scale;
   double u3scale;
   double u4scale;
   double w2scale;
   double w3scale;
   double w4scale;
   double w5scale;
   double udiag;
   int use_thole;
   char poltyp[MAX_STRING];
} polpot__;

struct {
   int use_bond;
   int use_angle;
   int use_strbnd;
   int use_urey;
   int use_angang;
   int use_opbend;
   int use_opdist;
   int use_improp;
   int use_imptor;
   int use_tors;
   int use_pitors;
   int use_strtor;
   int use_angtor;
   int use_tortor;
   int use_vdw;
   int use_repuls;
   int use_disp;
   int use_charge;
   int use_chgdpl;
   int use_dipole;
   int use_mpole;
   int use_polar;
   int use_chgtrn;
   int use_rxnfld;
   int use_solv;
   int use_metal;
   int use_geom;
   int use_extra;
   int use_born;
   int use_orbit;
} potent__;

struct {
   int npfix;
   int ndfix;
   int nafix;
   int ntfix;
   int ngfix;
   int nchir;
   int* ipfix;
   int* kpfix;
   int* idfix;
   int* iafix;
   int* itfix;
   int* igfix;
   int* ichir;
   double depth;
   double width;
   double rwall;
   double* xpfix;
   double* ypfix;
   double* zpfix;
   double* pfix;
   double* dfix;
   double* afix;
   double* tfix;
   double* gfix;
   double* chir;
   int use_basin;
   int use_wall;
} restrn__;

struct {
   int maxatm;
   int maxtyp;
   int maxclass;
   int maxval;
   int maxref;
   int maxgrp;
   int maxres;
   int maxfix;
} sizes__;

struct {
   double doffset;
   double p1;
   double p2;
   double p3;
   double p4;
   double p5;
   double* rsolv;
   double* asolv;
   double* rborn;
   double* drb;
   double* drbp;
   double* drobc;
   double* gpol;
   double* shct;
   double* aobc;
   double* bobc;
   double* gobc;
   double* vsolv;
   double* wace;
   double* s2ace;
   double* uace;
   char solvtyp[MAX_STRING];
   char borntyp[MAX_STRING];
} solute__;

struct {
   double friction;
   double* fgamma;
   int use_sdarea;
} stodyn__;

struct {
   int nstrbnd;
   int* isb;
   double* sbk;
} strbnd__;

static struct {
   int nstrtor;
   int* ist;
   double* kst;
} strtor__;

struct {
   double idihunit;
   double itorunit;
   double torsunit;
   double ptorunit;
   double storunit;
   double atorunit;
   double ttorunit;
} torpot__;

struct {
   int ntors;
   int* itors;
   double* tors1;
   double* tors2;
   double* tors3;
   double* tors4;
   double* tors5;
   double* tors6;
} tors__;

struct {
   int ntortor;
   int* itt;
} tortor__;

struct {
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

struct {
   int nurey;
   int* iury;
   double* uk;
   double* ul;
} urey__;

struct {
   double cury;
   double qury;
   double ureyunit;
} urypot__;

struct {
   int* nuse;
   int* iuse;
   int* use;
} usage__;

struct {
   int nvdw;
   int* ivdw;
   int* jvdw;
   int* ired;
   double* kred;
   double* xred;
   double* yred;
   double* zred;
   double* radmin;
   double* epsilon;
   double* radmin4;
   double* epsilon4;
   double* radhbnd;
   double* epshbnd;
} vdw__;

struct {
   int maxgauss;
   int ngauss;
   double* igauss;
   double abuck;
   double bbuck;
   double cbuck;
   double ghal;
   double dhal;
   double v2scale;
   double v3scale;
   double v4scale;
   double v5scale;
   int use_vcorr;
   char vdwindex[MAX_STRING];
   char radtyp[MAX_STRING];
   char radsiz[MAX_STRING];
   char gausstyp[MAX_STRING];
   char radrule[MAX_STRING];
   char epsrule[MAX_STRING];
   char vdwtyp[MAX_STRING];
} vdwpot__;

static void setNullTerminator (char* string, int maxLength, char* buffer) {

   int count;
   int ptr_i;
   char* ptr_c;
   ptr_c = string;
   ptr_i = (int) (*ptr_c);
   count = 0;

   // add contents of string to buffer until a non-character is found
   // or end of the string is reached, then add NULL to end of buffer

   while (ptr_i > 33 && ptr_i < 126 && count < maxLength) {
      buffer[count++] = (*ptr_c);
      ptr_c++;
      ptr_i = (int) (*ptr_c);
   }

   buffer[count] = '\0';

   return;
}

/*
 *    ############################################################
 *            Copy Tinker Fortran Data into C++ Structures
 *    ############################################################
 */

extern "C" {

void set_angbnd_data_ (int* nangle, int* iang, double* ak,
                       double* anat, double* afld) {

   angbnd__.nangle = *nangle;
   angbnd__.iang = iang;
   angbnd__.ak = ak;
   angbnd__.anat = anat;
   angbnd__.afld = afld;
}

void set_angpot_data_ (double* angunit, double* stbnunit, double* aaunit,
                       double* opbunit, double* opdunit, double* cang,
                       double* qang, double* pang, double* sang, double* copb,
                       double* qopb, double* popb, double* sopb, double* copd,
                       double* qopd, double* popd, double* sopd, char* opbtyp,
                       char* angtyp) {

   angpot__.angunit = *angunit;
   angpot__.stbnunit = *stbnunit;
   angpot__.aaunit = *aaunit;
   angpot__.opbunit = *opbunit;
   angpot__.opdunit = *opdunit;
   angpot__.cang = *cang;
   angpot__.qang = *qang;
   angpot__.pang = *pang;
   angpot__.sang = *sang;
   angpot__.copb = *copb;
   angpot__.qopb = *qopb;
   angpot__.popb = *popb;
   angpot__.sopb = *sopb;
   angpot__.copd = *copd;
   angpot__.qopd = *qopd;
   angpot__.popd = *popd;
   angpot__.sopd = *sopd;
   setNullTerminator (opbtyp, 8, angpot__.opbtyp);
   angpot__.angtyp = angtyp;
}

void set_angtor_data_ (int* nangtor, int* iat, double* kant) {

   angtor__.nangtor = *nangtor;
   angtor__.iat = iat;
   angtor__.kant = kant;
}

void set_atomid_data_ (int* tag, int* classs, int* atomic, int* valence,
                       double* mass, char* name, char* story) {

   atomid__.tag = tag;
   atomid__.classs = classs;
   atomid__.atomic = atomic;
   atomid__.valence = valence;
   atomid__.mass = mass;
   atomid__.name = name;
   atomid__.story = story;
}

void set_atoms_data_ (int* n, int* type, double* x, double* y, double* z) {

   atoms__.n = *n;
   atoms__.type = type;
   atoms__.x = x;
   atoms__.y = y;
   atoms__.z = z;
}

void set_bath_data_ (int* maxnose, int* voltrial, double* kelvin,
                     double* atmsph, double* tautemp, double* taupres,
                     double* compress, double* collide, double* eta,
                     double* volmove, double* vbar, double* qbar,
                     double* gbar, double* vnh, double* qnh, double* gnh,
                     int* isothermal, int* isobaric, int* anisotrop,
                     char* volscale, char* barostat, char* thermostat) {

   bath__.maxnose = *maxnose;
   bath__.voltrial = *voltrial;
   bath__.kelvin = *kelvin;
   bath__.atmsph = *atmsph;
   bath__.tautemp = *tautemp;
   bath__.taupres = *taupres;
   bath__.compress = *compress;
   bath__.collide = *collide;
   bath__.eta = *eta;
   bath__.volmove = *volmove;
   bath__.vbar = *vbar;
   bath__.qbar = *qbar;
   bath__.gbar = *gbar;
   bath__.vnh = vnh;
   bath__.qnh = qnh;
   bath__.gnh = gnh;
   bath__.isothermal = *isothermal;
   bath__.isobaric = *isobaric;
   bath__.anisotrop = *anisotrop;
   setNullTerminator (volscale, 9, bath__.volscale);
   setNullTerminator (barostat, 11, bath__.barostat);
   setNullTerminator (thermostat, 11, bath__.thermostat);
}

void set_bitor_data_ (int* nbitor, int* ibitor) {

   bitor__.nbitor = *nbitor;
   bitor__.ibitor = ibitor;
}

void set_bndpot_data_ (double* cbnd, double* qbnd, double* bndunit,
                       char* bndtyp) {

   bndpot__.cbnd = *cbnd;
   bndpot__.qbnd = *qbnd;
   bndpot__.bndunit = *bndunit;
   setNullTerminator (bndtyp, 8, bndpot__.bndtyp);
}

void set_bndstr_data_ (int* nbond, int* ibnd, double* bk, double* bl) {

   bndstr__.nbond = *nbond;
   bndstr__.ibnd = ibnd;
   bndstr__.bk = bk;
   bndstr__.bl = bl;
}

void set_bound_data_ (double* polycut, double* polycut2, int* use_bounds,
                      int* use_replica, int* use_polymer) {

   bound__.polycut = *polycut;
   bound__.polycut2 = *polycut2;
   bound__.use_bounds = *use_bounds;
   bound__.use_replica = *use_replica;
   bound__.use_polymer = *use_polymer;
}

void set_boxes_data_ (double* xbox, double* ybox, double* zbox,
                      double* alpha, double* beta, double* gamma,
                      double* xbox2, double* ybox2, double* zbox2,
                      double* box34, double* volbox, double* beta_sin,
                      double* beta_cos, double* gamma_sin, double* gamma_cos,
                      double* beta_term, double* gamma_term, double* lvec,
                      double* recip, int* orthogonal, int* monoclinic,
                      int* triclinic, int* octahedron, char* spacegrp) {

   boxes__.xbox = xbox;
   boxes__.ybox = ybox;
   boxes__.zbox = zbox;
   boxes__.alpha = alpha;
   boxes__.beta = beta;
   boxes__.gamma = gamma;
   boxes__.xbox2 = xbox2;
   boxes__.ybox2 = ybox2;
   boxes__.zbox2 = zbox2;
   boxes__.box34 = *box34;
   boxes__.volbox = *volbox;
   boxes__.beta_sin = *beta_sin;
   boxes__.beta_cos = *beta_cos;
   boxes__.gamma_sin = *gamma_sin;
   boxes__.gamma_cos = *gamma_cos;
   boxes__.beta_term = *beta_term;
   boxes__.gamma_term = *gamma_term;
   boxes__.lvec = lvec;
   boxes__.recip = recip;
   boxes__.orthogonal = *orthogonal;
   boxes__.monoclinic = *monoclinic;
   boxes__.triclinic = *triclinic;
   boxes__.octahedron = *octahedron;
   setNullTerminator (spacegrp, 10, boxes__.spacegrp);
}

void set_cell_data_ (int* ncell, int* icell, double* xcell, double* ycell,
                     double* zcell, double* xcell2, double* ycell2,
                     double* zcell2) {

   cell__.ncell = *ncell;
   cell__.icell = icell;
   cell__.xcell = *xcell;
   cell__.ycell = *ycell;
   cell__.zcell = *zcell;
   cell__.xcell2 = *xcell2;
   cell__.ycell2 = *ycell2;
   cell__.zcell2 = *zcell2;
}

void set_charge_data_ (int* nion, int* iion, int* jion,
                       int* kion, double* pchg) {

   charge__.nion = *nion;
   charge__.iion = iion;
   charge__.jion = jion;
   charge__.kion = kion;
   charge__.pchg = pchg;
}

void set_chgpot_data_ (double* electric, double* dielec, double* ebuffer,
                       double* c2scale, double* c3scale, double* c4scale,
                       double* c5scale, int* neutnbr, int* neutcut) {

   chgpot__.electric = *electric;
   chgpot__.dielec = *dielec;
   chgpot__.ebuffer = *ebuffer;
   chgpot__.c2scale = *c2scale;
   chgpot__.c3scale = *c3scale;
   chgpot__.c4scale = *c4scale;
   chgpot__.c5scale = *c5scale;
   chgpot__.neutnbr = *neutnbr;
   chgpot__.neutcut = *neutcut;
}

void set_couple_data_ (int* n12, int* n13, int* n14, int* n15,
                       int* i12, int* i13, int* i14, int* i15) {

   couple__.n12 = n12;
   couple__.n13 = n13;
   couple__.n14 = n14;
   couple__.n15 = n15;
   couple__.i12 = i12;
   couple__.i13 = i13;
   couple__.i14 = i14;
   couple__.i15 = i15;
}

void set_deriv_data_ (double* desum, double* deb, double* dea,
                      double* deba, double* deub, double* deaa,
                      double* deopb, double* deopd, double* deid,
                      double* deit, double* det, double* dept,
                      double* debt, double* deat, double* dett,
                      double* dev, double* der, double* dedsp,
                      double* dec, double* decd, double* ded,
                      double* dem, double* dep, double* dect,
                      double* derxf, double* des, double* delf,
                      double* deg, double* dex) {

   deriv__.desum = desum;
   deriv__.deb = deb;
   deriv__.dea = dea;
   deriv__.deba = deba;
   deriv__.deub = deub;
   deriv__.deaa = deaa;
   deriv__.deopb = deopb;
   deriv__.deopd = deopd;
   deriv__.deid = deid;
   deriv__.deit = deit;
   deriv__.det = det;
   deriv__.dept = dept;
   deriv__.debt = debt;
   deriv__.deat = deat;
   deriv__.dett = dett;
   deriv__.dev = dev;
   deriv__.der = der;
   deriv__.dedsp = dedsp;
   deriv__.dec = dec;
   deriv__.decd = decd;
   deriv__.ded = ded;
   deriv__.dem = dem;
   deriv__.dep = dep;
   deriv__.dect = dect;
   deriv__.derxf = derxf;
   deriv__.des = des;
   deriv__.delf = delf;
   deriv__.deg = deg;
   deriv__.dex = dex;
}

void set_energi_data_ (double* esum, double* eb, double* ea,
                       double* eba, double* eub, double* eaa,
                       double* eopb, double* eopd, double* eid,
                       double* eit, double* et, double* ept,
                       double* ebt, double* eat, double* ett,
                       double* ev, double* er, double* edsp,
                       double* ec, double* ecd, double* ed,
                       double* em, double* ep, double* ect,
                       double* erxf, double* es, double* elf,
                       double* eg, double* ex) {

   energi__.esum = esum;
   energi__.eb = eb;
   energi__.ea = ea;
   energi__.eba = eba;
   energi__.eub = eub;
   energi__.eaa = eaa;
   energi__.eopb = eopb;
   energi__.eopd = eopd;
   energi__.eid = eid;
   energi__.eit = eit;
   energi__.et = et;
   energi__.ept = ept;
   energi__.ebt = ebt;
   energi__.eat = eat;
   energi__.ett = ett;
   energi__.ev = ev;
   energi__.er = er;
   energi__.edsp = edsp;
   energi__.ec = ec;
   energi__.ecd = ecd;
   energi__.ed = ed;
   energi__.em = em;
   energi__.ep = ep;
   energi__.ect = ect;
   energi__.erxf = erxf;
   energi__.es = es;
   energi__.elf = elf;
   energi__.eg = eg;
   energi__.ex = ex;
}

void set_ewald_data_ (double* aewald, double* aeewald, double* apewald,
                      double* adewald, char* boundary) {

   ewald__.aewald = *aewald;
   ewald__.aeewald = *aeewald;
   ewald__.apewald = *apewald;
   ewald__.adewald = *adewald;
   setNullTerminator (boundary, 7, ewald__.boundary);
}

void set_freeze_data_ (int* nrat, int* nratx, int* iratx, int* kratx,
                       int* irat, double* rateps, double* krat,
                       int* use_rattle, int* ratimage) {

   freeze__.nrat = *nrat;
   freeze__.nratx = *nratx;
   freeze__.iratx = iratx;
   freeze__.kratx = kratx;
   freeze__.irat = irat;
   freeze__.rateps = *rateps;
   freeze__.krat = krat;
   freeze__.use_rattle = *use_rattle;
   freeze__.ratimage = ratimage;
}

void set_group_data_ (int* ngrp, int* kgrp, int* grplist, int* igrp,
                      double* grpmass, double* wgrp, int* use_group,
                      int* use_intra, int* use_inter) {

   group__.ngrp = *ngrp;
   group__.kgrp = kgrp;
   group__.grplist = grplist;
   group__.igrp = igrp;
   group__.grpmass = grpmass;
   group__.wgrp = wgrp;
   group__.use_group = *use_group;
   group__.use_intra = *use_intra;
   group__.use_inter = *use_inter;
}

void set_imptor_data_ (int* nitors, int* iitors, double* itors1,
                       double* itors2, double* itors3) {

   imptor__.nitors = *nitors;
   imptor__.iitors = iitors;
   imptor__.itors1 = itors1;
   imptor__.itors2 = itors2;
   imptor__.itors3 = itors3;
}

void set_inform_data_ (int* maxask, int* digits, int* iprint, int* iwrite,
                       int* isend, int* silent, int* verbose, int* debug,
                       int* holdup, int* abort) {

   inform__.maxask = *maxask;
   inform__.digits = *digits;
   inform__.iprint = *iprint;
   inform__.iwrite = *iwrite;
   inform__.isend = *isend;
   inform__.silent = *silent;
   inform__.verbose = *verbose;
   inform__.debug = *debug;
   inform__.holdup = *holdup;
   inform__.abort = *abort;
}

void set_ktrtor_data_ (int* maxntt, int* maxtgrd, int* maxtgrd2,
                       int* tnx, int* tny, double* ttx, double* tty,
                       double* tbf, double* tbx, double* tby,
                       double* tbxy, char20* ktt) {

   ktrtor__.maxntt = *maxntt;
   ktrtor__.maxtgrd = *maxtgrd;
   ktrtor__.maxtgrd2 = *maxtgrd2;
   ktrtor__.tnx = tnx;
   ktrtor__.tny = tny;
   ktrtor__.ttx = ttx;
   ktrtor__.tty = tty;
   ktrtor__.tbf = tbf;
   ktrtor__.tbx = tbx;
   ktrtor__.tby = tby;
   ktrtor__.tbxy = tbxy;
   ktrtor__.ktt = ktt;
}

void set_kvdwpr_data_ (int* maxnvp, double* radpr, double* epspr,
                       char* kvpr) {

   kvdwpr__.maxnvp = *maxnvp;
   kvdwpr__.radpr = radpr;
   kvdwpr__.epspr = epspr;
   kvdwpr__.kvpr = kvpr;
}

void set_kvdws_data_ (double* rad, double* eps, double* rad4,
                      double* eps4, double* reduct) {

   kvdws__.rad = rad;
   kvdws__.eps = eps;
   kvdws__.rad4 = rad4;
   kvdws__.eps4 = eps4;
   kvdws__.reduct = reduct;
}

void set_limits_data_ (double* vdwcut, double* repcut, double* dispcut,
                       double* chgcut, double* dplcut, double* mpolecut,
                       double* ctrncut, double* vdwtaper, double* reptaper,
                       double* disptaper, double* chgtaper, double* dpltaper,
                       double* mpoletaper, double* ctrntaper, double* ewaldcut,
                       double* dewaldcut, double* usolvcut, int* use_ewald,
                       int* use_dewald, int* use_lights, int* use_list,
                       int* use_vlist, int* use_dlist, int* use_clist,
                       int* use_mlist, int* use_ulist) {

   limits__.vdwcut = *vdwcut;
   limits__.repcut = *repcut;
   limits__.dispcut = *dispcut;
   limits__.chgcut = *chgcut;
   limits__.dplcut = *dplcut;
   limits__.mpolecut = *mpolecut;
   limits__.ctrncut = *ctrncut;
   limits__.vdwtaper = *vdwtaper;
   limits__.reptaper = *reptaper;
   limits__.disptaper = *disptaper;
   limits__.chgtaper = *chgtaper;
   limits__.dpltaper = *dpltaper;
   limits__.mpoletaper = *mpoletaper;
   limits__.ctrntaper = *ctrntaper;
   limits__.ewaldcut = *ewaldcut;
   limits__.dewaldcut = *dewaldcut;
   limits__.usolvcut = *usolvcut;
   limits__.use_ewald = *use_ewald;
   limits__.use_dewald = *use_dewald;
   limits__.use_lights = *use_lights;
   limits__.use_list = *use_list;
   limits__.use_vlist = *use_vlist;
   limits__.use_dlist = *use_dlist;
   limits__.use_clist = *use_clist;
   limits__.use_mlist = *use_mlist;
   limits__.use_ulist = *use_ulist;
}

void set_mdstuf_data_ (int* nfree, int* irest, int* bmnmix,
                       double* arespa, int* dorest, char* integrate) {

   mdstuf__.nfree = *nfree;
   mdstuf__.irest = *irest;
   mdstuf__.bmnmix = *bmnmix;
   mdstuf__.arespa = *arespa;
   mdstuf__.dorest = *dorest;
   setNullTerminator (integrate, 11, mdstuf__.integrate);
}

void set_molcul_data_ (int* nmol, int* imol, int* kmol, int* molcule,
                       double* totmass, double* molmass) {

   molcul__.nmol = *nmol;
   molcul__.imol = imol;
   molcul__.kmol = kmol;
   molcul__.molcule = molcule;
   molcul__.totmass = *totmass;
   molcul__.molmass = molmass;
}

void set_moldyn_data_ (double* v, double* a, double* aalt) {

   moldyn__.v = v;
   moldyn__.a = a;
   moldyn__.aalt = aalt;
}

void set_mplpot_data_ (double* m2scale, double* m3scale, double* m4scale,
                       double* m5scale, int* use_chgpen) {

   mplpot__.m2scale = *m2scale;
   mplpot__.m3scale = *m3scale;
   mplpot__.m4scale = *m4scale;
   mplpot__.m5scale = *m5scale;
   mplpot__.use_chgpen = *use_chgpen;
}

void set_mpole_data_ (int* maxpole, int* npole, int* ipole, int* polsiz,
                      int* pollist, int* zaxis, int* xaxis, int* yaxis,
                      double* pole, double* rpole, double* spole,
                      double* srpole, char* polaxe ) {

   mpole__.maxpole = *maxpole;
   mpole__.npole = *npole;
   mpole__.ipole = ipole;
   mpole__.polsiz = polsiz;
   mpole__.pollist = pollist;
   mpole__.zaxis = zaxis;
   mpole__.xaxis = xaxis;
   mpole__.yaxis = yaxis;
   mpole__.pole = pole;
   mpole__.rpole = rpole;
   mpole__.spole = spole;
   mpole__.srpole = srpole;
   mpole__.polaxe = polaxe;
}

void set_mutant_data_ (int* nmut, int* vcouple, int* imut, int* type0,
                       int* class0, int* type1, int* class1, double* lambda,
                       double* tlambda, double* vlambda, double* elambda,
                       double* scexp, double* scalpha, int* mut) {

   mutant__.nmut = *nmut;
   mutant__.vcouple = *vcouple;
   mutant__.imut = imut;
   mutant__.type0 = type0;
   mutant__.class0 = class0;
   mutant__.type1 = type1;
   mutant__.class1 = class1;
   mutant__.lambda = *lambda;
   mutant__.tlambda = *tlambda;
   mutant__.vlambda = *vlambda;
   mutant__.elambda = *elambda;
   mutant__.scexp = *scexp;
   mutant__.scalpha = *scalpha;
   mutant__.mut = mut;
}

void set_nonpol_data_ (double* epso, double* epsh, double* rmino,
                       double* rminh, double* awater, double* slevy,
                       double* solvprs, double* surften, double* spcut,
                       double* spoff, double* stcut, double* stoff,
                       double* rcav, double* rdisp, double* cdisp) {

   nonpol__.epso = *epso;
   nonpol__.epsh = *epsh;
   nonpol__.rmino = *rmino;
   nonpol__.rminh = *rminh;
   nonpol__.awater = *awater;
   nonpol__.slevy = *slevy;
   nonpol__.solvprs = *solvprs;
   nonpol__.surften = *surften;
   nonpol__.spcut = *spcut;
   nonpol__.spoff = *spoff;
   nonpol__.stcut = *stcut;
   nonpol__.stoff = *stoff;
   nonpol__.rcav = rcav;
   nonpol__.rdisp = rdisp;
   nonpol__.cdisp = cdisp;
}

void set_opbend_data_ (int* nopbend, int* iopb, double* opbk) {

   opbend__.nopbend = *nopbend;
   opbend__.iopb = iopb;
   opbend__.opbk = opbk;
}

void set_openmm_data_ (void** ommHandle, char* cudaPrecision,
                      char* ommPlatform, char* cudaDevice) {

   openmm__.ommHandle = *ommHandle;
   setNullTerminator (cudaPrecision, 6, openmm__.cudaPrecision);
   setNullTerminator (ommPlatform, 9, openmm__.ommPlatform);
   setNullTerminator (cudaDevice, 16, openmm__.cudaDevice);
}

void set_pitors_data_ (int* npitors, int* ipit, double* kpit) {

   pitors__.npitors = *npitors;
   pitors__.ipit = ipit;
   pitors__.kpit = kpit;
}

void set_pme_data_ (int* nfft1, int* nfft2, int* nfft3, int* nefft1,
                    int* nefft2, int* nefft3, int* ndfft1, int* ndfft2,
                    int* ndfft3, int* bsorder, int* bseorder, int* bsporder,
                    int* bsdorder, int* igrid, double* bsmod1, double* bsmod2,
                    double* bsmod3, double* bsbuild, double*** thetai1,
                    double*** thetai2, double*** thetai3, double**** qgrid,
                    double*** qfac) {

   pme__.nfft1 = *nfft1;
   pme__.nfft2 = *nfft2;
   pme__.nfft3 = *nfft3;
   pme__.nefft1 = *nefft1;
   pme__.nefft2 = *nefft2;
   pme__.nefft3 = *nefft3;
   pme__.ndfft1 = *ndfft1;
   pme__.ndfft2 = *ndfft2;
   pme__.ndfft3 = *ndfft3;
   pme__.bsorder = *bsorder;
   pme__.bseorder = *bseorder;
   pme__.bsporder = *bsporder;
   pme__.bsdorder = *bsdorder;
   pme__.igrid = igrid;
   pme__.bsmod1 = bsmod1;
   pme__.bsmod2 = bsmod2;
   pme__.bsmod3 = bsmod3;
   pme__.bsbuild = bsbuild;
   pme__.thetai1 = thetai1;
   pme__.thetai2 = thetai2;
   pme__.thetai3 = thetai3;
   pme__.qgrid = qgrid;
   pme__.qfac = qfac;
}

void set_polar_data_ (int* npolar, int* ipolar, double* polarity,
                      double* thole, double* dirdamp, double* pdamp,
                      double* udir,double* udirp, double* udirs,
                      double* udirps,double* uind, double* uinp,
                      double* uinds,double* uinps, double* uexact,
                      int* douind) {

   polar__.npolar = *npolar;
   polar__.ipolar = ipolar;
   polar__.polarity = polarity;
   polar__.thole = thole;
   polar__.dirdamp = dirdamp;
   polar__.pdamp = pdamp;
   polar__.udir = udir;
   polar__.udirp = udirp;
   polar__.udirs = udirs;
   polar__.udirps = udirps;
   polar__.uind = uind;
   polar__.uinp = uinp;
   polar__.uinds = uinds;
   polar__.uinps = uinps;
   polar__.uexact = uexact;
   polar__.douind = douind;
}

void set_polgrp_data_ (int* maxp11, int* maxp12, int* maxp13, int* maxp14,
                       int* np11, int* np12, int* np13, int* np14,
                       int* ip11, int* ip12, int* ip13, int* ip14) {

   polgrp__.maxp11 = *maxp11;
   polgrp__.maxp12 = *maxp12;
   polgrp__.maxp13 = *maxp13;
   polgrp__.maxp14 = *maxp14;
   polgrp__.np11 = np11;
   polgrp__.np12 = np12;
   polgrp__.np13 = np13;
   polgrp__.np14 = np14;
   polgrp__.ip11 = ip11;
   polgrp__.ip12 = ip12;
   polgrp__.ip13 = ip13;
   polgrp__.ip14 = ip14;
}

void set_polopt_data_ (int* maxopt, int* optorder, int* optlevel,
                       double* copt, double* copm, double* uopt,
                       double* uoptp, double* uopts, double* uoptps,
                       double* fopt, double* foptp) {

   polopt__.maxopt = *maxopt;
   polopt__.optorder = *optorder;
   polopt__.optlevel = *optlevel;
   polopt__.copt = copt;
   polopt__.copm = copm;
   polopt__.uopt = uopt;
   polopt__.uoptp = uoptp;
   polopt__.uopts = uopts;
   polopt__.uoptps = uoptps;
   polopt__.fopt = fopt;
   polopt__.foptp = foptp;
}

void set_polpot_data_ (int* politer, double* poleps, double* p2scale,
                       double* p3scale, double* p4scale, double* p5scale,
                       double* p2iscale, double* p3iscale, double* p4iscale,
                       double* p5iscale, double* d1scale, double* d2scale,
                       double* d3scale, double* d4scale, double* u1scale,
                       double* u2scale, double* u3scale, double* u4scale,
                       double* w2scale, double* w3scale, double* w4scale,
                       double* w5scale, double* udiag, int* use_thole,
                       char* poltyp) {

   polpot__.politer = *politer;
   polpot__.poleps = *poleps;
   polpot__.p2scale = *p2scale;
   polpot__.p3scale = *p3scale;
   polpot__.p4scale = *p4scale;
   polpot__.p2iscale = *p2iscale;
   polpot__.p3iscale = *p3iscale;
   polpot__.p4iscale = *p4iscale;
   polpot__.p5iscale = *p5iscale;
   polpot__.p5scale = *p5scale;
   polpot__.d1scale = *d1scale;
   polpot__.d2scale = *d2scale;
   polpot__.d3scale = *d3scale;
   polpot__.d4scale = *d4scale;
   polpot__.u1scale = *u1scale;
   polpot__.u2scale = *u2scale;
   polpot__.u3scale = *u3scale;
   polpot__.u4scale = *u4scale;
   polpot__.w2scale = *w2scale;
   polpot__.w3scale = *w3scale;
   polpot__.w4scale = *w4scale;
   polpot__.w5scale = *w5scale;
   polpot__.use_thole = *use_thole;
   polpot__.udiag = *udiag;
   setNullTerminator (poltyp, 6, polpot__.poltyp);
}

void set_potent_data_ (int* use_bond, int* use_angle, int* use_strbnd,
                       int* use_urey, int* use_angang, int* use_opbend,
                       int* use_opdist, int* use_improp, int* use_imptor,
                       int* use_tors, int* use_pitors, int* use_strtor,
                       int* use_angtor, int* use_tortor, int* use_vdw,
                       int* use_repuls, int* use_disp, int* use_charge,
                       int* use_chgdpl, int* use_dipole, int* use_mpole,
                       int* use_polar, int* use_chgtrn, int* use_rxnfld,
                       int* use_solv, int* use_metal, int* use_geom,
                       int* use_extra, int* use_born, int* use_orbit) {

   potent__.use_bond = *use_bond;
   potent__.use_angle = *use_angle;
   potent__.use_urey = *use_urey;
   potent__.use_strbnd = *use_strbnd;
   potent__.use_angang = *use_angang;
   potent__.use_opbend = *use_opbend;
   potent__.use_opdist = *use_opdist;
   potent__.use_improp = *use_improp;
   potent__.use_imptor = *use_imptor;
   potent__.use_tors = *use_tors;
   potent__.use_pitors = *use_pitors;
   potent__.use_strtor = *use_strtor;
   potent__.use_angtor = *use_angtor;
   potent__.use_tortor = *use_tortor;
   potent__.use_vdw = *use_vdw;
   potent__.use_repuls = *use_repuls;
   potent__.use_disp = *use_disp;
   potent__.use_charge = *use_charge;
   potent__.use_chgdpl = *use_chgdpl;
   potent__.use_dipole = *use_dipole;
   potent__.use_mpole = *use_mpole;
   potent__.use_polar = *use_polar;
   potent__.use_chgtrn = *use_chgtrn;
   potent__.use_rxnfld = *use_rxnfld;
   potent__.use_solv = *use_solv;
   potent__.use_metal = *use_metal;
   potent__.use_geom = *use_geom;
   potent__.use_extra = *use_extra;
   potent__.use_born = *use_born;
   potent__.use_orbit = *use_orbit;
}

void set_restrn_data_ (int* npfix, int* ndfix, int* nafix, int* ntfix,
                       int* ngfix, int* nchir, int* ipfix, int* kpfix,
                       int* idfix, int* iafix, int* itfix, int* igfix,
                       int* ichir, double* depth, double* width,
                       double* rwall, double* xpfix, double* ypfix,
                       double* zpfix, double* pfix, double* dfix,
                       double* afix, double* tfix, double* gfix,
                       double* chir, int* use_basin, int* use_wall) {

   restrn__.npfix = *npfix;
   restrn__.ndfix = *ndfix;
   restrn__.nafix = *nafix;
   restrn__.ntfix = *ntfix;
   restrn__.ngfix = *ngfix;
   restrn__.nchir = *nchir;
   restrn__.ipfix = ipfix;
   restrn__.kpfix = kpfix;
   restrn__.idfix = idfix;
   restrn__.iafix = iafix;
   restrn__.itfix = itfix;
   restrn__.igfix = igfix;
   restrn__.ichir = ichir;
   restrn__.depth = *depth;
   restrn__.width = *width;
   restrn__.rwall = *rwall;
   restrn__.xpfix = xpfix;
   restrn__.ypfix = ypfix;
   restrn__.zpfix = zpfix;
   restrn__.pfix = pfix;
   restrn__.dfix = dfix;
   restrn__.afix = afix;
   restrn__.tfix = tfix;
   restrn__.gfix = gfix;
   restrn__.chir = chir;
   restrn__.use_basin = *use_basin;
   restrn__.use_wall = *use_wall;
}

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

void set_solute_data_ (double* doffset, double* p1, double* p2, double* p3,
                       double* p4, double* p5, double* rsolv, double* asolv,
                       double* rborn, double* drb, double* drbp, double* drobc,
                       double* gpol, double* shct, double* aobc, double* bobc,
                       double* gobc, double* vsolv, double* wace,
                       double* s2ace, double* uace, char* solvtyp,
                       char* borntyp) {

   solute__.doffset = *doffset;
   solute__.p1 = *p1;
   solute__.p2 = *p2;
   solute__.p3 = *p3;
   solute__.p4 = *p4;
   solute__.p5 = *p5;
   solute__.rsolv = rsolv;
   solute__.asolv = asolv;
   solute__.rborn = rborn;
   solute__.drb = drb;
   solute__.drbp = drbp;
   solute__.drobc = drobc;
   solute__.gpol = gpol;
   solute__.shct = shct;
   solute__.aobc = aobc;
   solute__.bobc = bobc;
   solute__.gobc = gobc;
   solute__.vsolv = vsolv;
   solute__.wace = wace;
   solute__.s2ace = s2ace;
   solute__.uace = uace;
   setNullTerminator (solvtyp, 8, solute__.solvtyp);
   setNullTerminator (borntyp, 8, solute__.borntyp);
}

void set_stodyn_data_ (double* friction, double* fgamma, int* use_sdarea) {

   stodyn__.friction = *friction;
   stodyn__.fgamma = fgamma;
   stodyn__.use_sdarea = *use_sdarea;
}

void set_strbnd_data_ (int* nstrbnd, int* isb, double* sbk) {

   strbnd__.nstrbnd = *nstrbnd;
   strbnd__.isb = isb;
   strbnd__.sbk = sbk;
}

void set_strtor_data_ (int* nstrtor, int* ist, double* kst) {

   strtor__.nstrtor = *nstrtor;
   strtor__.ist = ist;
   strtor__.kst = kst;
}

void set_torpot_data_ (double* idihunit, double* itorunit, double* torsunit,
                       double* ptorunit, double* storunit, double* atorunit,
                       double* ttorunit) {

   torpot__.idihunit = *idihunit;
   torpot__.itorunit = *itorunit;
   torpot__.torsunit = *torsunit;
   torpot__.ptorunit = *ptorunit;
   torpot__.storunit = *storunit;
   torpot__.atorunit = *atorunit;
   torpot__.ttorunit = *ttorunit;
}

void set_tors_data_ (int* ntors, int* itors, double* tors1,
                     double* tors2, double* tors3, double* tors4,
                     double* tors5, double* tors6) {

   tors__.ntors = *ntors;
   tors__.itors = itors;
   tors__.tors1 = tors1;
   tors__.tors2 = tors2;
   tors__.tors3 = tors3;
   tors__.tors4 = tors4;
   tors__.tors5 = tors5;
   tors__.tors6 = tors6;
}

void set_tortor_data_ (int* ntortor, int* itt) {

   tortor__.ntortor = *ntortor;
   tortor__.itt = itt;
}

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

void set_urey_data_ (int* nurey, int* iury, double* uk, double* ul) {

   urey__.nurey = *nurey;
   urey__.iury = iury;
   urey__.uk = uk;
   urey__.ul = ul;
}

void set_urypot_data_ (double* cury, double* qury, double* ureyunit) {

   urypot__.cury = *cury;
   urypot__.qury = *qury;
   urypot__.ureyunit = *ureyunit;
}

void set_usage_data_ (int* nuse, int* iuse, int* use) {

   usage__.nuse = nuse;
   usage__.iuse = iuse;
   usage__.use = use;
}

void set_vdw_data_ (int* nvdw, int* ivdw, int* jvdw, int* ired,
                    double* kred, double* xred, double* yred, double* zred,
                    double* radmin, double* epsilon, double* radmin4,
                    double* epsilon4, double* radhbnd, double* epshbnd) {

   vdw__.nvdw = *nvdw;
   vdw__.ivdw = ivdw;
   vdw__.jvdw = jvdw;
   vdw__.ired = ired;
   vdw__.kred = kred;
   vdw__.xred = xred;
   vdw__.yred = yred;
   vdw__.zred = zred;
   vdw__.radmin = radmin;
   vdw__.epsilon = epsilon;
   vdw__.radmin4 = radmin4;
   vdw__.epsilon4 = epsilon4;
   vdw__.radhbnd = radhbnd;
   vdw__.epshbnd = epshbnd;
}

void set_vdwpot_data_ (int* maxgauss, int* ngauss, double* igauss,
                       double* abuck, double* bbuck, double* cbuck,
                       double* ghal, double* dhal, double* v2scale,
                       double* v3scale, double* v4scale, double* v5scale,
                       int* use_vcorr, char* vdwindex, char* radtyp,
                       char* radsiz, char* gausstyp, char* radrule,
                       char* epsrule, char* vdwtyp) {

   vdwpot__.maxgauss = *maxgauss;
   vdwpot__.ngauss = *ngauss;
   vdwpot__.igauss = igauss;
   vdwpot__.abuck = *abuck;
   vdwpot__.bbuck = *bbuck;
   vdwpot__.cbuck = *cbuck;
   vdwpot__.ghal = *ghal;
   vdwpot__.dhal = *dhal;
   vdwpot__.v2scale = *v2scale;
   vdwpot__.v3scale = *v3scale;
   vdwpot__.v4scale = *v4scale;
   vdwpot__.v5scale = *v5scale;
   vdwpot__.use_vcorr = *use_vcorr;
   setNullTerminator (vdwindex, 5, vdwpot__.vdwindex);
   setNullTerminator (radtyp, 5, vdwpot__.radtyp);
   setNullTerminator (radsiz, 8, vdwpot__.radsiz);
   setNullTerminator (gausstyp, 8, vdwpot__.gausstyp);
   setNullTerminator (radrule, 10, vdwpot__.radrule);
   setNullTerminator (epsrule, 10, vdwpot__.epsrule);
   setNullTerminator (vdwtyp, 13, vdwpot__.vdwtyp);
}

}

/*
 *  ###############################################################
 *  setupSystemParticles  --  add individual atoms to OpenMM System
 *  ###############################################################
 */

static void setupSystemParticles (OpenMM_System* system, FILE* log) {

   int ii;
   for (ii = 0; ii < atoms__.n; ii++) {
      OpenMM_System_addParticle (system, atomid__.mass[ii]);
   }
}

/*
 *  ####################################################################
 *  setupCMMotionRemover  --  frequency of center-of-mass motion removal
 *  ####################################################################
 */

static void setupCMMotionRemover (OpenMM_System* system, FILE* log) {

   int frequency = mdstuf__.irest > 0 ? mdstuf__.irest : 100;
   OpenMM_CMMotionRemover* cMMotionRemover;
   cMMotionRemover = OpenMM_CMMotionRemover_create (frequency);
   OpenMM_System_addForce (system, (OpenMM_Force*) cMMotionRemover);

}

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

/*
 *  ############################################################
 *  freeConstraintMap  --  removal of an existing constraint map
 *  ############################################################
 */

static void freeConstraintMap (struct ConstraintMap* map) {

   free (map->constraintOffset);
   free (map->constraintCount);
   free (map->constraintList);
}

/*
 *  ##############################################################
 *  mapConstraints  --  generate lists of Shake/Rattle constraints
 *  ##############################################################
 */

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

/*
 *  ###################################################################
 *  checkForConstraint  --  check of a constraint between pair of atoms
 *  ###################################################################
 */

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

/*
 *  ########################################################
 *  setupAmoebaBondForce  --  setup AMOEBA bond stretch term
 *  ########################################################
 */

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
   OpenMM_System_addForce (system, (OpenMM_Force*) amoebaBondForce);

   if (removeConstrainedBonds) {
       freeConstraintMap (&map);
   }
}

/*
 *  #############################################################
 *  openmm_init  --  initialization of the OpenMM data structures
 *  #############################################################
 *
 */

extern "C" void openmm_init_ (void** ommHandle, double* dt) {

   int ii;
   int mdMode = 0;
   int removeConstrainedCovalentIxns = 0;
   char buffer[128];
   char buffer2[128];
   FILE* log = stderr;

   // allocate space for opaque handle to hold OpenMM objects
   // such as system, integrator, context, etc.

   OpenMMData* omm = (OpenMMData*) malloc(sizeof(struct OpenMMData_s));

   // temporary OpenMM objects used and discarded here

   OpenMM_Vec3Array*       initialPosInNm;
   OpenMM_Vec3Array*       initialVelInNm;
   OpenMM_StringArray*     pluginList;
   OpenMM_Platform*        platform;

   // load the OpenMM plugin libraries from their default location;
   // call the plugin loading routine twice to fix an issue with MacOS
   // where the first library in the alphabetical list gets skipped

   pluginList = OpenMM_Platform_loadPluginsFromDirectory
                    (OpenMM_Platform_getDefaultPluginsDirectory());
   pluginList = OpenMM_Platform_loadPluginsFromDirectory
                    (OpenMM_Platform_getDefaultPluginsDirectory());

   if (inform__.verbose) {
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
   // to each force object so we can fill in the forces; note the OpenMM
   // System takes ownership of force objects, so don't delete them yourself

   omm->system = OpenMM_System_create ();
   setupSystemParticles (omm->system, log);
   setupCMMotionRemover (omm->system, log);

   if (potent__.use_bond) {
      setupAmoebaBondForce (omm->system,
                            removeConstrainedCovalentIxns, log);
   }

/*
   if (potent__.use_angle) {
      setupAmoebaAngleForce (omm->system,
                             removeConstrainedCovalentIxns, log);
      setupAmoebaInPlaneAngleForce (omm->system,
                                    removeConstrainedCovalentIxns, log);
   }

   if (potent__.use_strbnd) {
      setupAmoebaStretchBendForce (omm->system,
                                   removeConstrainedCovalentIxns, log);
   }

   if (potent__.use_urey) {
      setupAmoebaUreyBradleyForce (omm->system,
                                   removeConstrainedCovalentIxns, log);
   }

   if (potent__.use_opbend) {
      setupAmoebaOutOfPlaneBendForce (omm->system, log);
   }

   if (potent__.use_improp) {
      setupAmoebaImproperTorsionForce (omm->system, log);
   }

   if (potent__.use_tors) {
      setupAmoebaTorsionForce (omm->system, log);
   }

   if (potent__.use_strtor) {
      setupAmoebaStretchTorsionForce (omm->system, log);
   }

   if (potent__.use_angtor) {
      setupAmoebaAngleTorsionForce (omm->system, log);
   }

   if (potent__.use_pitors) {
      setupAmoebaPiTorsionForce (omm->system, log);
   }

   if (potent__.use_tortor) {
      setupAmoebaTorsionTorsionForce (omm->system, log);
   }

   if (potent__.use_vdw) {
      setupAmoebaVdwForce (omm->system, log);
   }

   if (potent__.use_charge) {
      setupAmoebaChargeForce (omm->system, log);
   }

   if (potent__.use_mpole) {
      setupAmoebaMultipoleForce (omm->system, log);
   }
*/

}

extern "C" void openmm_take_steps_ (void** omm, int* numSteps) {}

extern "C" void openmm_update_(void** omm, double* dt, int* istep,int* callMdStat, int* callMdSave) {}

extern "C" void openmm_cleanup_ (void** omm) {}

extern "C" void openmm_test_ () {}

extern "C" void openmm_bar_energy_ (void** ommHandle, double* energyInKcal) {}

#include "ommz1.h"
#include "ommz2.h"
#include "ommz3.h"
#include "ommz4.h"
#include "ommz5.h"

#include "ommj1.h"
#include "ommj2.h"
#include "ommj3.h"
#include "ommj4.h"
#include "ommj5.h"

#include "ommr1.h"
#include "ommr2.h"
#include "ommr3.h"
#include "ommr4.h"
#include "ommr5.h"
