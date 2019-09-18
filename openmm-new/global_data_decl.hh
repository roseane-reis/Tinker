extern struct action_st {
int* neb;
int* nea;
int* neba;
int* neub;
int* neaa;
int* neopb;
int* neopd;
int* neid;
int* neit;
int* net;
int* nept;
int* nebt;
int* neat;
int* nett;
int* nev;
int* ner;
int* nedsp;
int* nec;
int* necd;
int* ned;
int* nem;
int* nep;
int* nect;
int* new_;
int* nerxf;
int* nes;
int* nelf;
int* neg;
int* nex;
} action__;

extern "C" void set_action_data_ (int* neb, int* nea, int* neba, int* neub, int* neaa, int* neopb, int* neopd, int* neid, int* neit, int* net, int* nept, int* nebt, int* neat, int* nett, int* nev, int* ner, int* nedsp, int* nec, int* necd, int* ned, int* nem, int* nep, int* nect, int* new_, int* nerxf, int* nes, int* nelf, int* neg, int* nex);

extern struct align_st {
int* nfit;
int* ifit;
double* wfit;
} align__;

extern "C" void set_align_data_ (int* nfit, int* ifit, double* wfit);

extern struct analyz_st {
double* aesum;
double* aeb;
double* aea;
double* aeba;
double* aeub;
double* aeaa;
double* aeopb;
double* aeopd;
double* aeid;
double* aeit;
double* aet;
double* aept;
double* aebt;
double* aeat;
double* aett;
double* aev;
double* aer;
double* aedsp;
double* aec;
double* aecd;
double* aed;
double* aem;
double* aep;
double* aect;
double* aerxf;
double* aes;
double* aelf;
double* aeg;
double* aex;
} analyz__;

extern "C" void set_analyz_data_ (double* aesum, double* aeb, double* aea, double* aeba, double* aeub, double* aeaa, double* aeopb, double* aeopd, double* aeid, double* aeit, double* aet, double* aept, double* aebt, double* aeat, double* aett, double* aev, double* aer, double* aedsp, double* aec, double* aecd, double* aed, double* aem, double* aep, double* aect, double* aerxf, double* aes, double* aelf, double* aeg, double* aex);

extern struct angang_st {
int* nangang;
int* iaa;
double* kaa;
} angang__;

extern "C" void set_angang_data_ (int* nangang, int* iaa, double* kaa);

extern struct angbnd_st {
int* nangle;
int* iang;
double* ak;
double* anat;
double* afld;
} angbnd__;

extern "C" void set_angbnd_data_ (int* nangle, int* iang, double* ak, double* anat, double* afld);

extern struct angpot_st {
double* angunit;
double* stbnunit;
double* aaunit;
double* opbunit;
double* opdunit;
double* cang;
double* qang;
double* pang;
double* sang;
double* copb;
double* qopb;
double* popb;
double* sopb;
double* copd;
double* qopd;
double* popd;
double* sopd;
char* opbtyp;
char* angtyp;
} angpot__;

extern "C" void set_angpot_data_ (double* angunit, double* stbnunit, double* aaunit, double* opbunit, double* opdunit, double* cang, double* qang, double* pang, double* sang, double* copb, double* qopb, double* popb, double* sopb, double* copd, double* qopd, double* popd, double* sopd, char* opbtyp, char* angtyp);

extern struct angtor_st {
int* nangtor;
int* iat;
double* kant;
} angtor__;

extern "C" void set_angtor_data_ (int* nangtor, int* iat, double* kant);

extern struct argue_st {
int maxarg;
int* narg;
int* listarg;
char* arg;
} argue__;

extern "C" void set_argue_data_ (int* maxarg, int* narg, int* listarg, char* arg);

extern struct ascii_st {
int null;
int tab;
int linefeed;
int formfeed;
int carriage;
int escape;
int space;
int exclamation;
int quote;
int pound;
int dollar;
int percent;
int ampersand;
int apostrophe;
int asterisk;
int plus;
int comma;
int minus;
int period;
int frontslash;
int colon;
int semicolon;
int equal;
int question;
int atsign;
int backslash;
int caret;
int underbar;
int vertical;
int tilde;
} ascii__;

extern "C" void set_ascii_data_ (int* null, int* tab, int* linefeed, int* formfeed, int* carriage, int* escape, int* space, int* exclamation, int* quote, int* pound, int* dollar, int* percent, int* ampersand, int* apostrophe, int* asterisk, int* plus, int* comma, int* minus, int* period, int* frontslash, int* colon, int* semicolon, int* equal, int* question, int* atsign, int* backslash, int* caret, int* underbar, int* vertical, int* tilde);

extern struct atmlst_st {
int* bndlist;
int* anglist;
} atmlst__;

extern "C" void set_atmlst_data_ (int* bndlist, int* anglist);

extern struct atomid_st {
int* tag;
int* class_;
int* atomic;
int* valence;
double* mass;
char* name;
char* story;
} atomid__;

extern "C" void set_atomid_data_ (int* tag, int* class_, int* atomic, int* valence, double* mass, char* name, char* story);

extern struct atoms_st {
int* n;
int* type;
double* x;
double* y;
double* z;
} atoms__;

extern "C" void set_atoms_data_ (int* n, int* type, double* x, double* y, double* z);

extern struct bath_st {
int maxnose;
int* voltrial;
double* kelvin;
double* atmsph;
double* tautemp;
double* taupres;
double* compress;
double* collide;
double* eta;
double* volmove;
double* vbar;
double* qbar;
double* gbar;
double* vnh;
double* qnh;
double* gnh;
int* isothermal;
int* isobaric;
int* anisotrop;
char* volscale;
char* barostat;
char* thermostat;
} bath__;

extern "C" void set_bath_data_ (int* maxnose, int* voltrial, double* kelvin, double* atmsph, double* tautemp, double* taupres, double* compress, double* collide, double* eta, double* volmove, double* vbar, double* qbar, double* gbar, double* vnh, double* qnh, double* gnh, int* isothermal, int* isobaric, int* anisotrop, char* volscale, char* barostat, char* thermostat);

extern struct bitor_st {
int* nbitor;
int* ibitor;
} bitor__;

extern "C" void set_bitor_data_ (int* nbitor, int* ibitor);

extern struct bndpot_st {
double* cbnd;
double* qbnd;
double* bndunit;
char* bndtyp;
} bndpot__;

extern "C" void set_bndpot_data_ (double* cbnd, double* qbnd, double* bndunit, char* bndtyp);

extern struct bndstr_st {
int* nbond;
int* ibnd;
double* bk;
double* bl;
} bndstr__;

extern "C" void set_bndstr_data_ (int* nbond, int* ibnd, double* bk, double* bl);

extern struct bound_st {
double* polycut;
double* polycut2;
int* use_bounds;
int* use_replica;
int* use_polymer;
} bound__;

extern "C" void set_bound_data_ (double* polycut, double* polycut2, int* use_bounds, int* use_replica, int* use_polymer);

extern struct boxes_st {
double* xbox;
double* ybox;
double* zbox;
double* alpha;
double* beta;
double* gamma;
double* xbox2;
double* ybox2;
double* zbox2;
double* box34;
double* volbox;
double* beta_sin;
double* beta_cos;
double* gamma_sin;
double* gamma_cos;
double* beta_term;
double* gamma_term;
double* lvec;
double* recip;
int* orthogonal;
int* monoclinic;
int* triclinic;
int* octahedron;
char* spacegrp;
} boxes__;

extern "C" void set_boxes_data_ (double* xbox, double* ybox, double* zbox, double* alpha, double* beta, double* gamma, double* xbox2, double* ybox2, double* zbox2, double* box34, double* volbox, double* beta_sin, double* beta_cos, double* gamma_sin, double* gamma_cos, double* beta_term, double* gamma_term, double* lvec, double* recip, int* orthogonal, int* monoclinic, int* triclinic, int* octahedron, char* spacegrp);

extern struct cell_st {
int* ncell;
int* icell;
double* xcell;
double* ycell;
double* zcell;
double* xcell2;
double* ycell2;
double* zcell2;
} cell__;

extern "C" void set_cell_data_ (int* ncell, int* icell, double* xcell, double* ycell, double* zcell, double* xcell2, double* ycell2, double* zcell2);

extern struct charge_st {
int* nion;
int* iion;
int* jion;
int* kion;
double* pchg;
} charge__;

extern "C" void set_charge_data_ (int* nion, int* iion, int* jion, int* kion, double* pchg);

extern struct chgpen_st {
int* ncp;
double* pcore;
double* pval;
double* palpha;
} chgpen__;

extern "C" void set_chgpen_data_ (int* ncp, double* pcore, double* pval, double* palpha);

extern struct chgpot_st {
double* electric;
double* dielec;
double* ebuffer;
double* c2scale;
double* c3scale;
double* c4scale;
double* c5scale;
int* neutnbr;
int* neutcut;
} chgpot__;

extern "C" void set_chgpot_data_ (double* electric, double* dielec, double* ebuffer, double* c2scale, double* c3scale, double* c4scale, double* c5scale, int* neutnbr, int* neutcut);

extern struct chgtrn_st {
int* nct;
double* chgct;
double* dmpct;
} chgtrn__;

extern "C" void set_chgtrn_data_ (int* nct, double* chgct, double* dmpct);

extern struct chrono_st {
double* twall;
double* tcpu;
} chrono__;

extern "C" void set_chrono_data_ (double* twall, double* tcpu);

extern struct chunks_st {
int* nchunk;
int* nchk1;
int* nchk2;
int* nchk3;
int* ngrd1;
int* ngrd2;
int* ngrd3;
int* nlpts;
int* nrpts;
int* grdoff;
int* pmetable;
} chunks__;

extern "C" void set_chunks_data_ (int* nchunk, int* nchk1, int* nchk2, int* nchk3, int* ngrd1, int* ngrd2, int* ngrd3, int* nlpts, int* nrpts, int* grdoff, int* pmetable);

extern struct couple_st {
int* n12;
int* n13;
int* n14;
int* n15;
int* i12;
int* i13;
int* i14;
int* i15;
} couple__;

extern "C" void set_couple_data_ (int* n12, int* n13, int* n14, int* n15, int* i12, int* i13, int* i14, int* i15);

extern struct ctrpot_st {
char* ctrntyp;
} ctrpot__;

extern "C" void set_ctrpot_data_ (char* ctrntyp);

extern struct deriv_st {
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

extern "C" void set_deriv_data_ (double* desum, double* deb, double* dea, double* deba, double* deub, double* deaa, double* deopb, double* deopd, double* deid, double* deit, double* det, double* dept, double* debt, double* deat, double* dett, double* dev, double* der, double* dedsp, double* dec, double* decd, double* ded, double* dem, double* dep, double* dect, double* derxf, double* des, double* delf, double* deg, double* dex);

extern struct dipole_st {
int* ndipole;
int* idpl;
double* bdpl;
double* sdpl;
} dipole__;

extern "C" void set_dipole_data_ (int* ndipole, int* idpl, double* bdpl, double* sdpl);

extern struct disgeo_st {
double* vdwmax;
double* compact;
double* pathmax;
double* dbnd;
double* georad;
int* use_invert;
int* use_anneal;
} disgeo__;

extern "C" void set_disgeo_data_ (double* vdwmax, double* compact, double* pathmax, double* dbnd, double* georad, int* use_invert, int* use_anneal);

extern struct disp_st {
int* ndisp;
int* idisp;
double* csixpr;
double* csix;
double* adisp;
} disp__;

extern "C" void set_disp_data_ (int* ndisp, int* idisp, double* csixpr, double* csix, double* adisp);

extern struct dma_st {
double* mp;
double* dpx;
double* dpy;
double* dpz;
double* q20;
double* q21c;
double* q21s;
double* q22c;
double* q22s;
} dma__;

extern "C" void set_dma_data_ (double* mp, double* dpx, double* dpy, double* dpz, double* q20, double* q21c, double* q21s, double* q22c, double* q22s);

extern struct domega_st {
double* tesum;
double* teb;
double* tea;
double* teba;
double* teub;
double* teaa;
double* teopb;
double* teopd;
double* teid;
double* teit;
double* tet;
double* tept;
double* tebt;
double* teat;
double* tett;
double* tev;
double* ter;
double* tedsp;
double* tec;
double* tecd;
double* ted;
double* tem;
double* tep;
double* tect;
double* terxf;
double* tes;
double* telf;
double* teg;
double* tex;
} domega__;

extern "C" void set_domega_data_ (double* tesum, double* teb, double* tea, double* teba, double* teub, double* teaa, double* teopb, double* teopd, double* teid, double* teit, double* tet, double* tept, double* tebt, double* teat, double* tett, double* tev, double* ter, double* tedsp, double* tec, double* tecd, double* ted, double* tem, double* tep, double* tect, double* terxf, double* tes, double* telf, double* teg, double* tex);

extern struct dsppot_st {
double* dsp2scale;
double* dsp3scale;
double* dsp4scale;
double* dsp5scale;
int* use_dcorr;
} dsppot__;

extern "C" void set_dsppot_data_ (double* dsp2scale, double* dsp3scale, double* dsp4scale, double* dsp5scale, int* use_dcorr);

extern struct energi_st {
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

extern "C" void set_energi_data_ (double* esum, double* eb, double* ea, double* eba, double* eub, double* eaa, double* eopb, double* eopd, double* eid, double* eit, double* et, double* ept, double* ebt, double* eat, double* ett, double* ev, double* er, double* edsp, double* ec, double* ecd, double* ed, double* em, double* ep, double* ect, double* erxf, double* es, double* elf, double* eg, double* ex);

extern struct ewald_st {
double* aewald;
double* aeewald;
double* apewald;
double* adewald;
char* boundary;
} ewald__;

extern "C" void set_ewald_data_ (double* aewald, double* aeewald, double* apewald, double* adewald, char* boundary);

extern struct faces_st {
int* maxcls;
int* maxtt;
int* maxt;
int* maxp;
int* maxv;
int* maxen;
int* maxfn;
int* maxc;
int* maxep;
int* maxfs;
int* maxcy;
int* mxcyep;
int* maxfp;
int* mxfpcy;
int* na;
double* pr;
double* ar;
double* axyz;
int* skip;
int* nosurf;
int* afree;
int* abur;
int* cls;
int* clst;
int* acls;
int* ntt;
int* ttfe;
int* ttle;
int* enext;
int* tta;
int* ttbur;
int* ttfree;
int* nt;
int* tfe;
int* ta;
double* tr;
double* t;
double* tax;
int* tfree;
int* np;
int* pa;
double* p;
int* nv;
int* va;
int* vp;
double* vxyz;
int* nen;
int* nfn;
int* env;
int* fnen;
int* nc;
int* ca;
int* ct;
double* cr;
double* c;
int* nep;
int* epc;
int* epv;
int* afe;
int* ale;
int* epnext;
int* nfs;
int* fsen;
int* fsep;
int* ncy;
int* cynep;
int* cyep;
int* nfp;
int* fpa;
int* fpncy;
int* fpcy;
} faces__;

extern "C" void set_faces_data_ (int* maxcls, int* maxtt, int* maxt, int* maxp, int* maxv, int* maxen, int* maxfn, int* maxc, int* maxep, int* maxfs, int* maxcy, int* mxcyep, int* maxfp, int* mxfpcy, int* na, double* pr, double* ar, double* axyz, int* skip, int* nosurf, int* afree, int* abur, int* cls, int* clst, int* acls, int* ntt, int* ttfe, int* ttle, int* enext, int* tta, int* ttbur, int* ttfree, int* nt, int* tfe, int* ta, double* tr, double* t, double* tax, int* tfree, int* np, int* pa, double* p, int* nv, int* va, int* vp, double* vxyz, int* nen, int* nfn, int* env, int* fnen, int* nc, int* ca, int* ct, double* cr, double* c, int* nep, int* epc, int* epv, int* afe, int* ale, int* epnext, int* nfs, int* fsen, int* fsep, int* ncy, int* cynep, int* cyep, int* nfp, int* fpa, int* fpncy, int* fpcy);

extern struct fft_st {
int maxprime;
int* iprime;
unsigned long long* planf;
unsigned long long* planb;
double* ffttable;
char* ffttyp;
} fft__;

extern "C" void set_fft_data_ (int* maxprime, int* iprime, unsigned long long* planf, unsigned long long* planb, double* ffttable, char* ffttyp);

extern struct fields_st {
int maxbio;
int* biotyp;
char* forcefield;
} fields__;

extern "C" void set_fields_data_ (int* maxbio, int* biotyp, char* forcefield);

extern struct files_st {
int* nprior;
int* ldir;
int* leng;
char* filename;
char* outfile;
} files__;

extern "C" void set_files_data_ (int* nprior, int* ldir, int* leng, char* filename, char* outfile);

extern struct fracs_st {
double* xfrac;
double* yfrac;
double* zfrac;
} fracs__;

extern "C" void set_fracs_data_ (double* xfrac, double* yfrac, double* zfrac);

extern struct freeze_st {
int* nrat;
int* nratx;
int* iratx;
int* kratx;
int* irat;
double* rateps;
double* krat;
int* use_rattle;
int* ratimage;
} freeze__;

extern "C" void set_freeze_data_ (int* nrat, int* nratx, int* iratx, int* kratx, int* irat, double* rateps, double* krat, int* use_rattle, int* ratimage);

extern struct gkstuf_st {
double* gkc;
double* gkr;
} gkstuf__;

extern "C" void set_gkstuf_data_ (double* gkc, double* gkr);

extern struct group_st {
int* ngrp;
int* kgrp;
int* grplist;
int* igrp;
double* grpmass;
double* wgrp;
int* use_group;
int* use_intra;
int* use_inter;
} group__;

extern "C" void set_group_data_ (int* ngrp, int* kgrp, int* grplist, int* igrp, double* grpmass, double* wgrp, int* use_group, int* use_intra, int* use_inter);

extern struct hescut_st {
double* hesscut;
} hescut__;

extern "C" void set_hescut_data_ (double* hesscut);

extern struct hessn_st {
double* hessx;
double* hessy;
double* hessz;
} hessn__;

extern "C" void set_hessn_data_ (double* hessx, double* hessy, double* hessz);

extern struct hpmf_st {
double rcarbon;
double rwater;
double acsurf;
double safact;
double tslope;
double toffset;
double hpmfcut;
double h1;
double h2;
double h3;
double c1;
double c2;
double c3;
double w1;
double w2;
double w3;
int* npmf;
int* ipmf;
double* rpmf;
double* acsa;
} hpmf__;

extern "C" void set_hpmf_data_ (double* rcarbon, double* rwater, double* acsurf, double* safact, double* tslope, double* toffset, double* hpmfcut, double* h1, double* h2, double* h3, double* c1, double* c2, double* c3, double* w1, double* w2, double* w3, int* npmf, int* ipmf, double* rpmf, double* acsa);

extern struct ielscf_st {
int* nfree_aux;
double* tautemp_aux;
double* kelvin_aux;
double* uaux;
double* upaux;
double* vaux;
double* vpaux;
double* aaux;
double* apaux;
int* use_ielscf;
} ielscf__;

extern "C" void set_ielscf_data_ (int* nfree_aux, double* tautemp_aux, double* kelvin_aux, double* uaux, double* upaux, double* vaux, double* vpaux, double* aaux, double* apaux, int* use_ielscf);

extern struct improp_st {
int* niprop;
int* iiprop;
double* kprop;
double* vprop;
} improp__;

extern "C" void set_improp_data_ (int* niprop, int* iiprop, double* kprop, double* vprop);

extern struct imptor_st {
int* nitors;
int* iitors;
double* itors1;
double* itors2;
double* itors3;
} imptor__;

extern "C" void set_imptor_data_ (int* nitors, int* iitors, double* itors1, double* itors2, double* itors3);

extern struct inform_st {
int maxask;
int* digits;
int* iprint;
int* iwrite;
int* isend;
int* silent;
int* verbose;
int* debug;
int* holdup;
int* abort;
} inform__;

extern "C" void set_inform_data_ (int* maxask, int* digits, int* iprint, int* iwrite, int* isend, int* silent, int* verbose, int* debug, int* holdup, int* abort);

extern struct inter_st {
double* einter;
} inter__;

extern "C" void set_inter_data_ (double* einter);

extern struct iounit_st {
int* input;
int* iout;
} iounit__;

extern "C" void set_iounit_data_ (int* input, int* iout);

extern struct kanang_st {
double* anan;
} kanang__;

extern "C" void set_kanang_data_ (double* anan);

extern struct kangs_st {
int maxna;
int maxna5;
int maxna4;
int maxna3;
int maxnap;
int maxnaf;
double* acon;
double* acon5;
double* acon4;
double* acon3;
double* aconp;
double* aconf;
double* ang;
double* ang5;
double* ang4;
double* ang3;
double* angp;
double* angf;
char* ka;
char* ka5;
char* ka4;
char* ka3;
char* kap;
char* kaf;
} kangs__;

extern "C" void set_kangs_data_ (int* maxna, int* maxna5, int* maxna4, int* maxna3, int* maxnap, int* maxnaf, double* acon, double* acon5, double* acon4, double* acon3, double* aconp, double* aconf, double* ang, double* ang5, double* ang4, double* ang3, double* angp, double* angf, char* ka, char* ka5, char* ka4, char* ka3, char* kap, char* kaf);

extern struct kantor_st {
int maxnat;
double* atcon;
char* kat;
} kantor__;

extern "C" void set_kantor_data_ (int* maxnat, double* atcon, char* kat);

extern struct katoms_st {
int* atmcls;
int* atmnum;
int* ligand;
double* weight;
char* symbol;
char* describe;
} katoms__;

extern "C" void set_katoms_data_ (int* atmcls, int* atmnum, int* ligand, double* weight, char* symbol, char* describe);

extern struct kbonds_st {
int maxnb;
int maxnb5;
int maxnb4;
int maxnb3;
int maxnel;
double* bcon;
double* bcon5;
double* bcon4;
double* bcon3;
double* blen;
double* blen5;
double* blen4;
double* blen3;
double* dlen;
char* kb;
char* kb5;
char* kb4;
char* kb3;
char* kel;
} kbonds__;

extern "C" void set_kbonds_data_ (int* maxnb, int* maxnb5, int* maxnb4, int* maxnb3, int* maxnel, double* bcon, double* bcon5, double* bcon4, double* bcon3, double* blen, double* blen5, double* blen4, double* blen3, double* dlen, char* kb, char* kb5, char* kb4, char* kb3, char* kel);

extern struct kchrge_st {
double* chg;
} kchrge__;

extern "C" void set_kchrge_data_ (double* chg);

extern struct kcpen_st {
double* cpele;
double* cpalp;
} kcpen__;

extern "C" void set_kcpen_data_ (double* cpele, double* cpalp);

extern struct kctrn_st {
double* ctchg;
double* ctdmp;
} kctrn__;

extern "C" void set_kctrn_data_ (double* ctchg, double* ctdmp);

extern struct kdipol_st {
int maxnd;
int maxnd5;
int maxnd4;
int maxnd3;
double* dpl;
double* dpl5;
double* dpl4;
double* dpl3;
double* pos;
double* pos5;
double* pos4;
double* pos3;
char* kd;
char* kd5;
char* kd4;
char* kd3;
} kdipol__;

extern "C" void set_kdipol_data_ (int* maxnd, int* maxnd5, int* maxnd4, int* maxnd3, double* dpl, double* dpl5, double* dpl4, double* dpl3, double* pos, double* pos5, double* pos4, double* pos3, char* kd, char* kd5, char* kd4, char* kd3);

extern struct kdsp_st {
double* dspsix;
double* dspdmp;
} kdsp__;

extern "C" void set_kdsp_data_ (double* dspsix, double* dspdmp);

extern struct keys_st {
int maxkey;
int* nkey;
char* keyline;
} keys__;

extern "C" void set_keys_data_ (int* maxkey, int* nkey, char* keyline);

extern struct khbond_st {
int maxnhb;
double* radhb;
double* epshb;
char* khb;
} khbond__;

extern "C" void set_khbond_data_ (int* maxnhb, double* radhb, double* epshb, char* khb);

extern struct kiprop_st {
int maxndi;
double* dcon;
double* tdi;
char* kdi;
} kiprop__;

extern "C" void set_kiprop_data_ (int* maxndi, double* dcon, double* tdi, char* kdi);

extern struct kitors_st {
int maxnti;
double* ti1;
double* ti2;
double* ti3;
char* kti;
} kitors__;

extern "C" void set_kitors_data_ (int* maxnti, double* ti1, double* ti2, double* ti3, char* kti);

extern struct kmulti_st {
int maxnmp;
double* multip;
char* mpaxis;
char* kmp;
} kmulti__;

extern "C" void set_kmulti_data_ (int* maxnmp, double* multip, char* mpaxis, char* kmp);

extern struct kopbnd_st {
int maxnopb;
double* opbn;
char* kopb;
} kopbnd__;

extern "C" void set_kopbnd_data_ (int* maxnopb, double* opbn, char* kopb);

extern struct kopdst_st {
int maxnopd;
double* opds;
char* kopd;
} kopdst__;

extern "C" void set_kopdst_data_ (int* maxnopd, double* opds, char* kopd);

extern struct korbs_st {
int maxnpi;
int maxnpi5;
int maxnpi4;
double* sslope;
double* sslope5;
double* sslope4;
double* tslope;
double* tslope5;
double* tslope4;
double* electron;
double* ionize;
double* repulse;
char* kpi;
char* kpi5;
char* kpi4;
} korbs__;

extern "C" void set_korbs_data_ (int* maxnpi, int* maxnpi5, int* maxnpi4, double* sslope, double* sslope5, double* sslope4, double* tslope, double* tslope5, double* tslope4, double* electron, double* ionize, double* repulse, char* kpi, char* kpi5, char* kpi4);

extern struct kpitor_st {
int maxnpt;
double* ptcon;
char* kpt;
} kpitor__;

extern "C" void set_kpitor_data_ (int* maxnpt, double* ptcon, char* kpt);

extern struct kpolr_st {
int* pgrp;
double* polr;
double* athl;
double* ddir;
} kpolr__;

extern "C" void set_kpolr_data_ (int* pgrp, double* polr, double* athl, double* ddir);

extern struct krepl_st {
double* prsiz;
double* prdmp;
double* prele;
} krepl__;

extern "C" void set_krepl_data_ (double* prsiz, double* prdmp, double* prele);

extern struct kstbnd_st {
int maxnsb;
double* stbn;
char* ksb;
} kstbnd__;

extern "C" void set_kstbnd_data_ (int* maxnsb, double* stbn, char* ksb);

extern struct ksttor_st {
int maxnbt;
double* btcon;
char* kbt;
} ksttor__;

extern "C" void set_ksttor_data_ (int* maxnbt, double* btcon, char* kbt);

extern struct ktorsn_st {
int maxnt;
int maxnt5;
int maxnt4;
double* t1;
double* t2;
double* t3;
double* t4;
double* t5;
double* t6;
double* t15;
double* t25;
double* t35;
double* t45;
double* t55;
double* t65;
double* t14;
double* t24;
double* t34;
double* t44;
double* t54;
double* t64;
char* kt;
char* kt5;
char* kt4;
} ktorsn__;

extern "C" void set_ktorsn_data_ (int* maxnt, int* maxnt5, int* maxnt4, double* t1, double* t2, double* t3, double* t4, double* t5, double* t6, double* t15, double* t25, double* t35, double* t45, double* t55, double* t65, double* t14, double* t24, double* t34, double* t44, double* t54, double* t64, char* kt, char* kt5, char* kt4);

extern struct ktrtor_st {
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
char* ktt;
} ktrtor__;

extern "C" void set_ktrtor_data_ (int* maxntt, int* maxtgrd, int* maxtgrd2, int* tnx, int* tny, double* ttx, double* tty, double* tbf, double* tbx, double* tby, double* tbxy, char* ktt);

extern struct kurybr_st {
int maxnu;
double* ucon;
double* dst13;
char* ku;
} kurybr__;

extern "C" void set_kurybr_data_ (int* maxnu, double* ucon, double* dst13, char* ku);

extern struct kvdwpr_st {
int maxnvp;
double* radpr;
double* epspr;
char* kvpr;
} kvdwpr__;

extern "C" void set_kvdwpr_data_ (int* maxnvp, double* radpr, double* epspr, char* kvpr);

extern struct kvdws_st {
double* rad;
double* eps;
double* rad4;
double* eps4;
double* reduct;
} kvdws__;

extern "C" void set_kvdws_data_ (double* rad, double* eps, double* rad4, double* eps4, double* reduct);

extern struct light_st {
int* nlight;
int* kbx;
int* kby;
int* kbz;
int* kex;
int* key;
int* kez;
int* locx;
int* locy;
int* locz;
int* rgx;
int* rgy;
int* rgz;
} light__;

extern "C" void set_light_data_ (int* nlight, int* kbx, int* kby, int* kbz, int* kex, int* key, int* kez, int* locx, int* locy, int* locz, int* rgx, int* rgy, int* rgz);

extern struct limits_st {
double* vdwcut;
double* repcut;
double* dispcut;
double* chgcut;
double* dplcut;
double* mpolecut;
double* ctrncut;
double* vdwtaper;
double* reptaper;
double* disptaper;
double* chgtaper;
double* dpltaper;
double* mpoletaper;
double* ctrntaper;
double* ewaldcut;
double* dewaldcut;
double* usolvcut;
int* use_ewald;
int* use_dewald;
int* use_lights;
int* use_list;
int* use_vlist;
int* use_dlist;
int* use_clist;
int* use_mlist;
int* use_ulist;
} limits__;

extern "C" void set_limits_data_ (double* vdwcut, double* repcut, double* dispcut, double* chgcut, double* dplcut, double* mpolecut, double* ctrncut, double* vdwtaper, double* reptaper, double* disptaper, double* chgtaper, double* dpltaper, double* mpoletaper, double* ctrntaper, double* ewaldcut, double* dewaldcut, double* usolvcut, int* use_ewald, int* use_dewald, int* use_lights, int* use_list, int* use_vlist, int* use_dlist, int* use_clist, int* use_mlist, int* use_ulist);

extern struct linmin_st {
int* intmax;
double* stpmin;
double* stpmax;
double* cappa;
double* slpmax;
double* angmax;
} linmin__;

extern "C" void set_linmin_data_ (int* intmax, double* stpmin, double* stpmax, double* cappa, double* slpmax, double* angmax);

extern struct math_st {
double pi;
double elog;
double radian;
double logten;
double twosix;
double sqrtpi;
double sqrttwo;
double sqrtthree;
} math__;

extern "C" void set_math_data_ (double* pi, double* elog, double* radian, double* logten, double* twosix, double* sqrtpi, double* sqrttwo, double* sqrtthree);

extern struct mdstuf_st {
int* nfree;
int* irest;
int* bmnmix;
double* arespa;
int* dorest;
char* integrate;
} mdstuf__;

extern "C" void set_mdstuf_data_ (int* nfree, int* irest, int* bmnmix, double* arespa, int* dorest, char* integrate);

extern struct merck_st {
int* nlignes;
int* bt_1;
int* eqclass;
int* crd;
int* val;
int* pilp;
int* mltb;
int* arom;
int* lin;
int* sbmb;
int* mmffarom;
int* mmffaromc;
int* mmffaroma;
double* rad0;
double* paulel;
double* r0ref;
double* kbref;
double* mmff_kb;
double* mmff_kb1;
double* mmff_b0;
double* mmff_b1;
double* mmff_ka;
double* mmff_ka1;
double* mmff_ka2;
double* mmff_ka3;
double* mmff_ka4;
double* mmff_ka5;
double* mmff_ka6;
double* mmff_ka7;
double* mmff_ka8;
double* mmff_ang0;
double* mmff_ang1;
double* mmff_ang2;
double* mmff_ang3;
double* mmff_ang4;
double* mmff_ang5;
double* mmff_ang6;
double* mmff_ang7;
double* mmff_ang8;
double* stbn_abc;
double* stbn_cba;
double* stbn_abc1;
double* stbn_cba1;
double* stbn_abc2;
double* stbn_cba2;
double* stbn_abc3;
double* stbn_cba3;
double* stbn_abc4;
double* stbn_cba4;
double* stbn_abc5;
double* stbn_cba5;
double* stbn_abc6;
double* stbn_cba6;
double* stbn_abc7;
double* stbn_cba7;
double* stbn_abc8;
double* stbn_cba8;
double* stbn_abc9;
double* stbn_cba9;
double* stbn_abc10;
double* stbn_cba10;
double* stbn_abc11;
double* stbn_cba11;
double* defstbn_abc;
double* defstbn_cba;
double* t1_1;
double* t2_1;
double* t3_1;
double* t1_2;
double* t2_2;
double* t3_2;
char* kt_1;
char* kt_2;
double* g;
double* alph;
double* nn;
char* da;
double* bci;
double* bci_1;
double* pbci;
double* fcadj;
} merck__;

extern "C" void set_merck_data_ (int* nlignes, int* bt_1, int* eqclass, int* crd, int* val, int* pilp, int* mltb, int* arom, int* lin, int* sbmb, int* mmffarom, int* mmffaromc, int* mmffaroma, double* rad0, double* paulel, double* r0ref, double* kbref, double* mmff_kb, double* mmff_kb1, double* mmff_b0, double* mmff_b1, double* mmff_ka, double* mmff_ka1, double* mmff_ka2, double* mmff_ka3, double* mmff_ka4, double* mmff_ka5, double* mmff_ka6, double* mmff_ka7, double* mmff_ka8, double* mmff_ang0, double* mmff_ang1, double* mmff_ang2, double* mmff_ang3, double* mmff_ang4, double* mmff_ang5, double* mmff_ang6, double* mmff_ang7, double* mmff_ang8, double* stbn_abc, double* stbn_cba, double* stbn_abc1, double* stbn_cba1, double* stbn_abc2, double* stbn_cba2, double* stbn_abc3, double* stbn_cba3, double* stbn_abc4, double* stbn_cba4, double* stbn_abc5, double* stbn_cba5, double* stbn_abc6, double* stbn_cba6, double* stbn_abc7, double* stbn_cba7, double* stbn_abc8, double* stbn_cba8, double* stbn_abc9, double* stbn_cba9, double* stbn_abc10, double* stbn_cba10, double* stbn_abc11, double* stbn_cba11, double* defstbn_abc, double* defstbn_cba, double* t1_1, double* t2_1, double* t3_1, double* t1_2, double* t2_2, double* t3_2, char* kt_1, char* kt_2, double* g, double* alph, double* nn, char* da, double* bci, double* bci_1, double* pbci, double* fcadj);

extern struct minima_st {
int* maxiter;
int* nextiter;
double* fctmin;
double* hguess;
} minima__;

extern "C" void set_minima_data_ (int* maxiter, int* nextiter, double* fctmin, double* hguess);

extern struct molcul_st {
int* nmol;
int* imol;
int* kmol;
int* molcule;
double* totmass;
double* molmass;
} molcul__;

extern "C" void set_molcul_data_ (int* nmol, int* imol, int* kmol, int* molcule, double* totmass, double* molmass);

extern struct moldyn_st {
double* v;
double* a;
double* aalt;
} moldyn__;

extern "C" void set_moldyn_data_ (double* v, double* a, double* aalt);

extern struct moment_st {
double* netchg;
double* netdpl;
double* netqdp;
double* xdpl;
double* ydpl;
double* zdpl;
double* xxqdp;
double* xyqdp;
double* xzqdp;
double* yxqdp;
double* yyqdp;
double* yzqdp;
double* zxqdp;
double* zyqdp;
double* zzqdp;
} moment__;

extern "C" void set_moment_data_ (double* netchg, double* netdpl, double* netqdp, double* xdpl, double* ydpl, double* zdpl, double* xxqdp, double* xyqdp, double* xzqdp, double* yxqdp, double* yyqdp, double* yzqdp, double* zxqdp, double* zyqdp, double* zzqdp);

extern struct mplpot_st {
double* m2scale;
double* m3scale;
double* m4scale;
double* m5scale;
int* use_chgpen;
char* pentyp;
} mplpot__;

extern "C" void set_mplpot_data_ (double* m2scale, double* m3scale, double* m4scale, double* m5scale, int* use_chgpen, char* pentyp);

extern struct mpole_st {
int maxpole;
int* npole;
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

extern "C" void set_mpole_data_ (int* maxpole, int* npole, int* ipole, int* polsiz, int* pollist, int* zaxis, int* xaxis, int* yaxis, double* pole, double* rpole, double* spole, double* srpole, char* polaxe);

extern struct mrecip_st {
double* vmxx;
double* vmyy;
double* vmzz;
double* vmxy;
double* vmxz;
double* vmyz;
double* cmp;
double* fmp;
double* cphi;
double* fphi;
} mrecip__;

extern "C" void set_mrecip_data_ (double* vmxx, double* vmyy, double* vmzz, double* vmxy, double* vmxz, double* vmyz, double* cmp, double* fmp, double* cphi, double* fphi);

extern struct mutant_st {
int* nmut;
int* vcouple;
int* imut;
int* type0;
int* class0;
int* type1;
int* class1;
double* lambda;
double* tlambda;
double* vlambda;
double* elambda;
double* scexp;
double* scalpha;
int* mut;
} mutant__;

extern "C" void set_mutant_data_ (int* nmut, int* vcouple, int* imut, int* type0, int* class0, int* type1, int* class1, double* lambda, double* tlambda, double* vlambda, double* elambda, double* scexp, double* scalpha, int* mut);

extern struct neigh_st {
int* maxvlst;
int* maxelst;
int* maxulst;
int* nvlst;
int* vlst;
int* nelst;
int* elst;
int* nulst;
int* ulst;
double* lbuffer;
double* pbuffer;
double* lbuf2;
double* pbuf2;
double* vbuf2;
double* vbufx;
double* dbuf2;
double* dbufx;
double* cbuf2;
double* cbufx;
double* mbuf2;
double* mbufx;
double* ubuf2;
double* ubufx;
double* xvold;
double* yvold;
double* zvold;
double* xeold;
double* yeold;
double* zeold;
double* xuold;
double* yuold;
double* zuold;
int* dovlst;
int* dodlst;
int* doclst;
int* domlst;
int* doulst;
} neigh__;

extern "C" void set_neigh_data_ (int* maxvlst, int* maxelst, int* maxulst, int* nvlst, int* vlst, int* nelst, int* elst, int* nulst, int* ulst, double* lbuffer, double* pbuffer, double* lbuf2, double* pbuf2, double* vbuf2, double* vbufx, double* dbuf2, double* dbufx, double* cbuf2, double* cbufx, double* mbuf2, double* mbufx, double* ubuf2, double* ubufx, double* xvold, double* yvold, double* zvold, double* xeold, double* yeold, double* zeold, double* xuold, double* yuold, double* zuold, int* dovlst, int* dodlst, int* doclst, int* domlst, int* doulst);

extern struct nonpol_st {
double epso;
double epsh;
double rmino;
double rminh;
double awater;
double slevy;
double* solvprs;
double* surften;
double* spcut;
double* spoff;
double* stcut;
double* stoff;
double* rcav;
double* rdisp;
double* cdisp;
} nonpol__;

extern "C" void set_nonpol_data_ (double* epso, double* epsh, double* rmino, double* rminh, double* awater, double* slevy, double* solvprs, double* surften, double* spcut, double* spoff, double* stcut, double* stoff, double* rcav, double* rdisp, double* cdisp);

extern struct nucleo_st {
int* pucker;
double* glyco;
double* bkbone;
int* dblhlx;
int* deoxy;
char* hlxform;
} nucleo__;

extern "C" void set_nucleo_data_ (int* pucker, double* glyco, double* bkbone, int* dblhlx, int* deoxy, char* hlxform);

extern struct omega_st {
int* nomega;
int* iomega;
int* zline;
double* dihed;
} omega__;

extern "C" void set_omega_data_ (int* nomega, int* iomega, int* zline, double* dihed);

extern struct opbend_st {
int* nopbend;
int* iopb;
double* opbk;
} opbend__;

extern "C" void set_opbend_data_ (int* nopbend, int* iopb, double* opbk);

extern struct opdist_st {
int* nopdist;
int* iopd;
double* opdk;
} opdist__;

extern "C" void set_opdist_data_ (int* nopdist, int* iopd, double* opdk);

extern struct openmm_st {
unsigned long long* ommhandle;
char* cudaprecision;
char* ommplatform;
char* cudadevice;
} openmm__;

extern "C" void set_openmm_data_ (unsigned long long* ommhandle, char* cudaprecision, char* ommplatform, char* cudadevice);

extern struct openmp_st {
int* nproc;
int* nthread;
} openmp__;

extern "C" void set_openmp_data_ (int* nproc, int* nthread);

extern struct orbits_st {
double* qorb;
double* worb;
double* emorb;
} orbits__;

extern "C" void set_orbits_data_ (double* qorb, double* worb, double* emorb);

extern struct output_st {
int* archive;
int* noversion;
int* overwrite;
int* cyclesave;
int* velsave;
int* frcsave;
int* uindsave;
char* coordtype;
} output__;

extern "C" void set_output_data_ (int* archive, int* noversion, int* overwrite, int* cyclesave, int* velsave, int* frcsave, int* uindsave, char* coordtype);

extern struct params_st {
int maxprm;
int* nprm;
char* prmline;
} params__;

extern "C" void set_params_data_ (int* maxprm, int* nprm, char* prmline);

extern struct paths_st {
double* pnorm;
double* acoeff;
double* pc0;
double* pc1;
double* pvect;
double* pstep;
double* pzet;
double* gc;
} paths__;

extern "C" void set_paths_data_ (double* pnorm, double* acoeff, double* pc0, double* pc1, double* pvect, double* pstep, double* pzet, double* gc);

extern struct pbstuf_st {
int maxion;
int* ionn;
int* dime;
int* ionq;
double* pbe;
double* pdie;
double* sdie;
double* srad;
double* swin;
double* sdens;
double* smin;
double* grid;
double* gcent;
double* cgrid;
double* cgcent;
double* fgrid;
double* fgcent;
double* ionr;
double* ionc;
double* apbe;
double* pbr;
double* pbep;
double* pbfp;
double* pbtp;
double* pbeuind;
double* pbeuinp;
char* pbtyp;
char* pbsoln;
char* bcfl;
char* chgm;
char* srfm;
} pbstuf__;

extern "C" void set_pbstuf_data_ (int* maxion, int* ionn, int* dime, int* ionq, double* pbe, double* pdie, double* sdie, double* srad, double* swin, double* sdens, double* smin, double* grid, double* gcent, double* cgrid, double* cgcent, double* fgrid, double* fgcent, double* ionr, double* ionc, double* apbe, double* pbr, double* pbep, double* pbfp, double* pbtp, double* pbeuind, double* pbeuinp, char* pbtyp, char* pbsoln, char* bcfl, char* chgm, char* srfm);

extern struct pdb_st {
int* npdb;
int* nres;
int* resnum;
int* resatm;
int* npdb12;
int* ipdb12;
int* pdblist;
double* xpdb;
double* ypdb;
double* zpdb;
char* altsym;
char* pdbres;
char* pdbatm;
char* pdbtyp;
char* chnsym;
char* instyp;
} pdb__;

extern "C" void set_pdb_data_ (int* npdb, int* nres, int* resnum, int* resatm, int* npdb12, int* ipdb12, int* pdblist, double* xpdb, double* ypdb, double* zpdb, char* altsym, char* pdbres, char* pdbatm, char* pdbtyp, char* chnsym, char* instyp);

extern struct phipsi_st {
int* chiral;
int* disulf;
double* phi;
double* psi;
double* omega;
double* chi;
} phipsi__;

extern "C" void set_phipsi_data_ (int* chiral, int* disulf, double* phi, double* psi, double* omega, double* chi);

extern struct piorbs_st {
int* norbit;
int* nconj;
int* reorbit;
int* nbpi;
int* ntpi;
int* iorbit;
int* iconj;
int* kconj;
int* piperp;
int* ibpi;
int* itpi;
double* pbpl;
double* pnpl;
int* listpi;
} piorbs__;

extern "C" void set_piorbs_data_ (int* norbit, int* nconj, int* reorbit, int* nbpi, int* ntpi, int* iorbit, int* iconj, int* kconj, int* piperp, int* ibpi, int* itpi, double* pbpl, double* pnpl, int* listpi);

extern struct pistuf_st {
double* bkpi;
double* blpi;
double* kslope;
double* lslope;
double* torsp2;
} pistuf__;

extern "C" void set_pistuf_data_ (double* bkpi, double* blpi, double* kslope, double* lslope, double* torsp2);

extern struct pitors_st {
int* npitors;
int* ipit;
double* kpit;
} pitors__;

extern "C" void set_pitors_data_ (int* npitors, int* ipit, double* kpit);

extern struct pme_st {
int* nfft1;
int* nfft2;
int* nfft3;
int* nefft1;
int* nefft2;
int* nefft3;
int* ndfft1;
int* ndfft2;
int* ndfft3;
int* bsorder;
int* bseorder;
int* bsporder;
int* bsdorder;
int* igrid;
double* bsmod1;
double* bsmod2;
double* bsmod3;
double* bsbuild;
double* thetai1;
double* thetai2;
double* thetai3;
double* qgrid;
double* qfac;
} pme__;

extern "C" void set_pme_data_ (int* nfft1, int* nfft2, int* nfft3, int* nefft1, int* nefft2, int* nefft3, int* ndfft1, int* ndfft2, int* ndfft3, int* bsorder, int* bseorder, int* bsporder, int* bsdorder, int* igrid, double* bsmod1, double* bsmod2, double* bsmod3, double* bsbuild, double* thetai1, double* thetai2, double* thetai3, double* qgrid, double* qfac);

extern struct polar_st {
int* npolar;
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

extern "C" void set_polar_data_ (int* npolar, int* ipolar, double* polarity, double* thole, double* dirdamp, double* pdamp, double* udir, double* udirp, double* udirs, double* udirps, double* uind, double* uinp, double* uinds, double* uinps, double* uexact, int* douind);

extern struct polgrp_st {
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

extern "C" void set_polgrp_data_ (int* maxp11, int* maxp12, int* maxp13, int* maxp14, int* np11, int* np12, int* np13, int* np14, int* ip11, int* ip12, int* ip13, int* ip14);

extern struct polopt_st {
int maxopt;
int* optorder;
int* optlevel;
double* copt;
double* copm;
double* uopt;
double* uoptp;
double* uopts;
double* uoptps;
double* fopt;
double* foptp;
} polopt__;

extern "C" void set_polopt_data_ (int* maxopt, int* optorder, int* optlevel, double* copt, double* copm, double* uopt, double* uoptp, double* uopts, double* uoptps, double* fopt, double* foptp);

extern struct polpcg_st {
int* mindex;
double* pcgpeek;
double* minv;
int* pcgprec;
int* pcgguess;
} polpcg__;

extern "C" void set_polpcg_data_ (int* mindex, double* pcgpeek, double* minv, int* pcgprec, int* pcgguess);

extern struct polpot_st {
int* politer;
double* poleps;
double* p2scale;
double* p3scale;
double* p4scale;
double* p5scale;
double* p2iscale;
double* p3iscale;
double* p4iscale;
double* p5iscale;
double* d1scale;
double* d2scale;
double* d3scale;
double* d4scale;
double* u1scale;
double* u2scale;
double* u3scale;
double* u4scale;
double* w2scale;
double* w3scale;
double* w4scale;
double* w5scale;
double* udiag;
int* dpequal;
int* use_thole;
int* use_dirdamp;
char* poltyp;
} polpot__;

extern "C" void set_polpot_data_ (int* politer, double* poleps, double* p2scale, double* p3scale, double* p4scale, double* p5scale, double* p2iscale, double* p3iscale, double* p4iscale, double* p5iscale, double* d1scale, double* d2scale, double* d3scale, double* d4scale, double* u1scale, double* u2scale, double* u3scale, double* u4scale, double* w2scale, double* w3scale, double* w4scale, double* w5scale, double* udiag, int* dpequal, int* use_thole, int* use_dirdamp, char* poltyp);

extern struct poltcg_st {
int* tcgorder;
int* tcgnab;
double* tcgpeek;
double* uad;
double* uap;
double* ubd;
double* ubp;
int* tcgguess;
} poltcg__;

extern "C" void set_poltcg_data_ (int* tcgorder, int* tcgnab, double* tcgpeek, double* uad, double* uap, double* ubd, double* ubp, int* tcgguess);

extern struct potent_st {
int* use_bond;
int* use_angle;
int* use_strbnd;
int* use_urey;
int* use_angang;
int* use_opbend;
int* use_opdist;
int* use_improp;
int* use_imptor;
int* use_tors;
int* use_pitors;
int* use_strtor;
int* use_angtor;
int* use_tortor;
int* use_vdw;
int* use_repuls;
int* use_disp;
int* use_charge;
int* use_chgdpl;
int* use_dipole;
int* use_mpole;
int* use_polar;
int* use_chgtrn;
int* use_rxnfld;
int* use_solv;
int* use_metal;
int* use_geom;
int* use_extra;
int* use_born;
int* use_orbit;
} potent__;

extern "C" void set_potent_data_ (int* use_bond, int* use_angle, int* use_strbnd, int* use_urey, int* use_angang, int* use_opbend, int* use_opdist, int* use_improp, int* use_imptor, int* use_tors, int* use_pitors, int* use_strtor, int* use_angtor, int* use_tortor, int* use_vdw, int* use_repuls, int* use_disp, int* use_charge, int* use_chgdpl, int* use_dipole, int* use_mpole, int* use_polar, int* use_chgtrn, int* use_rxnfld, int* use_solv, int* use_metal, int* use_geom, int* use_extra, int* use_born, int* use_orbit);

extern struct potfit_st {
int* nconf;
int* namax;
int* ngatm;
int* nfatm;
int* npgrid;
int* ipgrid;
double* resp;
double* xdpl0;
double* ydpl0;
double* zdpl0;
double* xxqdp0;
double* xyqdp0;
double* xzqdp0;
double* yyqdp0;
double* yzqdp0;
double* zzqdp0;
double* fit0;
double* fchg;
double* fpol;
double* pgrid;
double* epot;
int* use_dpl;
int* use_qdp;
int* fit_mpl;
int* fit_dpl;
int* fit_qdp;
int* fitchg;
int* fitpol;
int* gatm;
int* fatm;
} potfit__;

extern "C" void set_potfit_data_ (int* nconf, int* namax, int* ngatm, int* nfatm, int* npgrid, int* ipgrid, double* resp, double* xdpl0, double* ydpl0, double* zdpl0, double* xxqdp0, double* xyqdp0, double* xzqdp0, double* yyqdp0, double* yzqdp0, double* zzqdp0, double* fit0, double* fchg, double* fpol, double* pgrid, double* epot, int* use_dpl, int* use_qdp, int* fit_mpl, int* fit_dpl, int* fit_qdp, int* fitchg, int* fitpol, int* gatm, int* fatm);

extern struct ptable_st {
int maxele;
double* atmass;
double* vdwrad;
double* covrad;
char* elemnt;
} ptable__;

extern "C" void set_ptable_data_ (int* maxele, double* atmass, double* vdwrad, double* covrad, char* elemnt);

extern struct qmstuf_st {
int* ngatom;
double* egau;
double* gx;
double* gy;
double* gz;
double* gfreq;
double* gforce;
double* gh;
} qmstuf__;

extern "C" void set_qmstuf_data_ (int* ngatom, double* egau, double* gx, double* gy, double* gz, double* gfreq, double* gforce, double* gh);

extern struct refer_st {
int* nref;
int* refltitle;
int* refleng;
int* reftyp;
int* n12ref;
int* i12ref;
double* xboxref;
double* yboxref;
double* zboxref;
double* alpharef;
double* betaref;
double* gammaref;
double* xref;
double* yref;
double* zref;
char* refnam;
char* reffile;
char* reftitle;
} refer__;

extern "C" void set_refer_data_ (int* nref, int* refltitle, int* refleng, int* reftyp, int* n12ref, int* i12ref, double* xboxref, double* yboxref, double* zboxref, double* alpharef, double* betaref, double* gammaref, double* xref, double* yref, double* zref, char* refnam, char* reffile, char* reftitle);

extern struct repel_st {
int* nrep;
double* sizpr;
double* dmppr;
double* elepr;
} repel__;

extern "C" void set_repel_data_ (int* nrep, double* sizpr, double* dmppr, double* elepr);

extern struct reppot_st {
double* r2scale;
double* r3scale;
double* r4scale;
double* r5scale;
} reppot__;

extern "C" void set_reppot_data_ (double* r2scale, double* r3scale, double* r4scale, double* r5scale);

extern struct resdue_st {
int maxamino;
int maxnuc;
int* ntyp;
int* catyp;
int* ctyp;
int* hntyp;
int* otyp;
int* hatyp;
int* cbtyp;
int* nntyp;
int* cantyp;
int* cntyp;
int* hnntyp;
int* ontyp;
int* hantyp;
int* nctyp;
int* cactyp;
int* cctyp;
int* hnctyp;
int* octyp;
int* hactyp;
int* o5typ;
int* c5typ;
int* h51typ;
int* h52typ;
int* c4typ;
int* h4typ;
int* o4typ;
int* c1typ;
int* h1typ;
int* c3typ;
int* h3typ;
int* c2typ;
int* h21typ;
int* o2typ;
int* h22typ;
int* o3typ;
int* ptyp;
int* optyp;
int* h5ttyp;
int* h3ttyp;
char* amino1;
char* nuclz1;
char* amino;
char* nuclz;
} resdue__;

extern "C" void set_resdue_data_ (int* maxamino, int* maxnuc, int* ntyp, int* catyp, int* ctyp, int* hntyp, int* otyp, int* hatyp, int* cbtyp, int* nntyp, int* cantyp, int* cntyp, int* hnntyp, int* ontyp, int* hantyp, int* nctyp, int* cactyp, int* cctyp, int* hnctyp, int* octyp, int* hactyp, int* o5typ, int* c5typ, int* h51typ, int* h52typ, int* c4typ, int* h4typ, int* o4typ, int* c1typ, int* h1typ, int* c3typ, int* h3typ, int* c2typ, int* h21typ, int* o2typ, int* h22typ, int* o3typ, int* ptyp, int* optyp, int* h5ttyp, int* h3ttyp, char* amino1, char* nuclz1, char* amino, char* nuclz);

extern struct restrn_st {
int* npfix;
int* ndfix;
int* nafix;
int* ntfix;
int* ngfix;
int* nchir;
int* ipfix;
int* kpfix;
int* idfix;
int* iafix;
int* itfix;
int* igfix;
int* ichir;
double* depth;
double* width;
double* rwall;
double* xpfix;
double* ypfix;
double* zpfix;
double* pfix;
double* dfix;
double* afix;
double* tfix;
double* gfix;
double* chir;
int* use_basin;
int* use_wall;
} restrn__;

extern "C" void set_restrn_data_ (int* npfix, int* ndfix, int* nafix, int* ntfix, int* ngfix, int* nchir, int* ipfix, int* kpfix, int* idfix, int* iafix, int* itfix, int* igfix, int* ichir, double* depth, double* width, double* rwall, double* xpfix, double* ypfix, double* zpfix, double* pfix, double* dfix, double* afix, double* tfix, double* gfix, double* chir, int* use_basin, int* use_wall);

extern struct rgddyn_st {
double* xcmo;
double* ycmo;
double* zcmo;
double* vcm;
double* wcm;
double* lm;
double* vc;
double* wc;
int* linear;
} rgddyn__;

extern "C" void set_rgddyn_data_ (double* xcmo, double* ycmo, double* zcmo, double* vcm, double* wcm, double* lm, double* vc, double* wc, int* linear);

extern struct rigid_st {
double* xrb;
double* yrb;
double* zrb;
double* rbc;
int* use_rigid;
} rigid__;

extern "C" void set_rigid_data_ (double* xrb, double* yrb, double* zrb, double* rbc, int* use_rigid);

extern struct ring_st {
int* nring3;
int* nring4;
int* nring5;
int* nring6;
int* nring7;
int* iring3;
int* iring4;
int* iring5;
int* iring6;
int* iring7;
} ring__;

extern "C" void set_ring_data_ (int* nring3, int* nring4, int* nring5, int* nring6, int* nring7, int* iring3, int* iring4, int* iring5, int* iring6, int* iring7);

extern struct rotbnd_st {
int* nrot;
int* rot;
int* use_short;
} rotbnd__;

extern "C" void set_rotbnd_data_ (int* nrot, int* rot, int* use_short);

extern struct rxnfld_st {
int* ijk;
double* b1;
double* b2;
} rxnfld__;

extern "C" void set_rxnfld_data_ (int* ijk, double* b1, double* b2);

extern struct rxnpot_st {
int* rfterms;
double* rfsize;
double* rfbulkd;
} rxnpot__;

extern "C" void set_rxnpot_data_ (int* rfterms, double* rfsize, double* rfbulkd);

extern struct scales_st {
double* scale;
int* set_scale;
} scales__;

extern "C" void set_scales_data_ (double* scale, int* set_scale);

extern struct sequen_st {
int* nseq;
int* nchain;
int* ichain;
int* seqtyp;
char* chnnam;
char* seq;
char* chntyp;
} sequen__;

extern "C" void set_sequen_data_ (int* nseq, int* nchain, int* ichain, int* seqtyp, char* chnnam, char* seq, char* chntyp);

extern struct shunt_st {
double* off;
double* off2;
double* cut;
double* cut2;
double* c0;
double* c1;
double* c2;
double* c3;
double* c4;
double* c5;
double* f0;
double* f1;
double* f2;
double* f3;
double* f4;
double* f5;
double* f6;
double* f7;
} shunt__;

extern "C" void set_shunt_data_ (double* off, double* off2, double* cut, double* cut2, double* c0, double* c1, double* c2, double* c3, double* c4, double* c5, double* f0, double* f1, double* f2, double* f3, double* f4, double* f5, double* f6, double* f7);

extern struct sizes_st {
int maxatm;
int maxtyp;
int maxclass;
int maxval;
int maxref;
int maxgrp;
int maxres;
int maxfix;
} sizes__;

extern "C" void set_sizes_data_ (int* maxatm, int* maxtyp, int* maxclass, int* maxval, int* maxref, int* maxgrp, int* maxres, int* maxfix);

extern struct socket_st {
int* skttyp;
int* cstep;
double* cdt;
double* cenergy;
int* sktstart;
int* sktstop;
int* use_socket;
} socket__;

extern "C" void set_socket_data_ (int* skttyp, int* cstep, double* cdt, double* cenergy, int* sktstart, int* sktstop, int* use_socket);

extern struct solute_st {
double* doffset;
double* p1;
double* p2;
double* p3;
double* p4;
double* p5;
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
char* solvtyp;
char* borntyp;
} solute__;

extern "C" void set_solute_data_ (double* doffset, double* p1, double* p2, double* p3, double* p4, double* p5, double* rsolv, double* asolv, double* rborn, double* drb, double* drbp, double* drobc, double* gpol, double* shct, double* aobc, double* bobc, double* gobc, double* vsolv, double* wace, double* s2ace, double* uace, char* solvtyp, char* borntyp);

extern struct stodyn_st {
double* friction;
double* fgamma;
int* use_sdarea;
} stodyn__;

extern "C" void set_stodyn_data_ (double* friction, double* fgamma, int* use_sdarea);

extern struct strbnd_st {
int* nstrbnd;
int* isb;
double* sbk;
} strbnd__;

extern "C" void set_strbnd_data_ (int* nstrbnd, int* isb, double* sbk);

extern struct strtor_st {
int* nstrtor;
int* ist;
double* kst;
} strtor__;

extern "C" void set_strtor_data_ (int* nstrtor, int* ist, double* kst);

extern struct syntrn_st {
double* tpath;
double* ppath;
double* xmin1;
double* xmin2;
double* xm;
} syntrn__;

extern "C" void set_syntrn_data_ (double* tpath, double* ppath, double* xmin1, double* xmin2, double* xm);

extern struct tarray_st {
int* ntpair;
int* tindex;
double* tdipdip;
} tarray__;

extern "C" void set_tarray_data_ (int* ntpair, int* tindex, double* tdipdip);

extern struct titles_st {
int* ltitle;
char* title;
} titles__;

extern "C" void set_titles_data_ (int* ltitle, char* title);

extern struct torpot_st {
double* idihunit;
double* itorunit;
double* torsunit;
double* ptorunit;
double* storunit;
double* atorunit;
double* ttorunit;
} torpot__;

extern "C" void set_torpot_data_ (double* idihunit, double* itorunit, double* torsunit, double* ptorunit, double* storunit, double* atorunit, double* ttorunit);

extern struct tors_st {
int* ntors;
int* itors;
double* tors1;
double* tors2;
double* tors3;
double* tors4;
double* tors5;
double* tors6;
} tors__;

extern "C" void set_tors_data_ (int* ntors, int* itors, double* tors1, double* tors2, double* tors3, double* tors4, double* tors5, double* tors6);

extern struct tortor_st {
int* ntortor;
int* itt;
} tortor__;

extern "C" void set_tortor_data_ (int* ntortor, int* itt);

extern struct tree_st {
int maxpss;
int* nlevel;
double* etree;
double* ilevel;
} tree__;

extern "C" void set_tree_data_ (int* maxpss, int* nlevel, double* etree, double* ilevel);

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

extern "C" void set_units_data_ (double* avogadro, double* lightspd, double* boltzmann, double* gasconst, double* elemchg, double* vacperm, double* emass, double* planck, double* joule, double* ekcal, double* bohr, double* hartree, double* evolt, double* efreq, double* coulomb, double* debye, double* prescon);

extern struct uprior_st {
int maxpred;
int* nualt;
int* maxualt;
double* gear;
double* aspc;
double* bpred;
double* bpredp;
double* bpreds;
double* bpredps;
double* udalt;
double* upalt;
double* usalt;
double* upsalt;
int* use_pred;
char* polpred;
} uprior__;

extern "C" void set_uprior_data_ (int* maxpred, int* nualt, int* maxualt, double* gear, double* aspc, double* bpred, double* bpredp, double* bpreds, double* bpredps, double* udalt, double* upalt, double* usalt, double* upsalt, int* use_pred, char* polpred);

extern struct urey_st {
int* nurey;
int* iury;
double* uk;
double* ul;
} urey__;

extern "C" void set_urey_data_ (int* nurey, int* iury, double* uk, double* ul);

extern struct urypot_st {
double* cury;
double* qury;
double* ureyunit;
} urypot__;

extern "C" void set_urypot_data_ (double* cury, double* qury, double* ureyunit);

extern struct usage_st {
int* nuse;
int* iuse;
int* use;
} usage__;

extern "C" void set_usage_data_ (int* nuse, int* iuse, int* use);

extern struct valfit_st {
int* fit_bond;
int* fit_angle;
int* fit_strbnd;
int* fit_urey;
int* fit_opbend;
int* fit_tors;
int* fit_struct;
int* fit_force;
} valfit__;

extern "C" void set_valfit_data_ (int* fit_bond, int* fit_angle, int* fit_strbnd, int* fit_urey, int* fit_opbend, int* fit_tors, int* fit_struct, int* fit_force);

extern struct vdw_st {
int* nvdw;
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

extern "C" void set_vdw_data_ (int* nvdw, int* ivdw, int* jvdw, int* ired, double* kred, double* xred, double* yred, double* zred, double* radmin, double* epsilon, double* radmin4, double* epsilon4, double* radhbnd, double* epshbnd);

extern struct vdwpot_st {
int maxgauss;
int* ngauss;
double* igauss;
double* abuck;
double* bbuck;
double* cbuck;
double* ghal;
double* dhal;
double* v2scale;
double* v3scale;
double* v4scale;
double* v5scale;
int* use_vcorr;
char* vdwindex;
char* radtyp;
char* radsiz;
char* gausstyp;
char* radrule;
char* epsrule;
char* vdwtyp;
} vdwpot__;

extern "C" void set_vdwpot_data_ (int* maxgauss, int* ngauss, double* igauss, double* abuck, double* bbuck, double* cbuck, double* ghal, double* dhal, double* v2scale, double* v3scale, double* v4scale, double* v5scale, int* use_vcorr, char* vdwindex, char* radtyp, char* radsiz, char* gausstyp, char* radrule, char* epsrule, char* vdwtyp);

extern struct vibs_st {
double* phi;
double* phik;
double* pwork;
} vibs__;

extern "C" void set_vibs_data_ (double* phi, double* phik, double* pwork);

extern struct virial_st {
double* vir;
int* use_virial;
} virial__;

extern "C" void set_virial_data_ (double* vir, int* use_virial);

extern struct warp_st {
double* deform;
double* difft;
double* diffv;
double* diffc;
double* m2;
int* use_smooth;
int* use_dem;
int* use_gda;
int* use_tophat;
int* use_stophat;
} warp__;

extern "C" void set_warp_data_ (double* deform, double* difft, double* diffv, double* diffc, double* m2, int* use_smooth, int* use_dem, int* use_gda, int* use_tophat, int* use_stophat);

extern struct xtals_st {
int maxlsq;
int maxrsd;
int* nxtal;
int* nvary;
int* ivary;
int* iresid;
int* vary;
double* e0_lattice;
char* vartyp;
char* rsdtyp;
} xtals__;

extern "C" void set_xtals_data_ (int* maxlsq, int* maxrsd, int* nxtal, int* nvary, int* ivary, int* iresid, int* vary, double* e0_lattice, char* vartyp, char* rsdtyp);

extern struct zclose_st {
int* nadd;
int* ndel;
int* iadd;
int* idel;
} zclose__;

extern "C" void set_zclose_data_ (int* nadd, int* ndel, int* iadd, int* idel);

extern struct zcoord_st {
int* iz;
double* zbond;
double* zang;
double* ztors;
} zcoord__;

extern "C" void set_zcoord_data_ (int* iz, double* zbond, double* zang, double* ztors);

