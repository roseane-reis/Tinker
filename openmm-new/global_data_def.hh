action_st action__;

void set_action_data_ (int* neb, int* nea, int* neba, int* neub, int* neaa, int* neopb, int* neopd, int* neid, int* neit, int* net, int* nept, int* nebt, int* neat, int* nett, int* nev, int* ner, int* nedsp, int* nec, int* necd, int* ned, int* nem, int* nep, int* nect, int* new_, int* nerxf, int* nes, int* nelf, int* neg, int* nex) {
action__.neb = neb;
action__.nea = nea;
action__.neba = neba;
action__.neub = neub;
action__.neaa = neaa;
action__.neopb = neopb;
action__.neopd = neopd;
action__.neid = neid;
action__.neit = neit;
action__.net = net;
action__.nept = nept;
action__.nebt = nebt;
action__.neat = neat;
action__.nett = nett;
action__.nev = nev;
action__.ner = ner;
action__.nedsp = nedsp;
action__.nec = nec;
action__.necd = necd;
action__.ned = ned;
action__.nem = nem;
action__.nep = nep;
action__.nect = nect;
action__.new_ = new_;
action__.nerxf = nerxf;
action__.nes = nes;
action__.nelf = nelf;
action__.neg = neg;
action__.nex = nex;
}

align_st align__;

void set_align_data_ (int* nfit, int* ifit, double* wfit) {
align__.nfit = nfit;
align__.ifit = ifit;
align__.wfit = wfit;
}

analyz_st analyz__;

void set_analyz_data_ (double* aesum, double* aeb, double* aea, double* aeba, double* aeub, double* aeaa, double* aeopb, double* aeopd, double* aeid, double* aeit, double* aet, double* aept, double* aebt, double* aeat, double* aett, double* aev, double* aer, double* aedsp, double* aec, double* aecd, double* aed, double* aem, double* aep, double* aect, double* aerxf, double* aes, double* aelf, double* aeg, double* aex) {
analyz__.aesum = aesum;
analyz__.aeb = aeb;
analyz__.aea = aea;
analyz__.aeba = aeba;
analyz__.aeub = aeub;
analyz__.aeaa = aeaa;
analyz__.aeopb = aeopb;
analyz__.aeopd = aeopd;
analyz__.aeid = aeid;
analyz__.aeit = aeit;
analyz__.aet = aet;
analyz__.aept = aept;
analyz__.aebt = aebt;
analyz__.aeat = aeat;
analyz__.aett = aett;
analyz__.aev = aev;
analyz__.aer = aer;
analyz__.aedsp = aedsp;
analyz__.aec = aec;
analyz__.aecd = aecd;
analyz__.aed = aed;
analyz__.aem = aem;
analyz__.aep = aep;
analyz__.aect = aect;
analyz__.aerxf = aerxf;
analyz__.aes = aes;
analyz__.aelf = aelf;
analyz__.aeg = aeg;
analyz__.aex = aex;
}

angang_st angang__;

void set_angang_data_ (int* nangang, int* iaa, double* kaa) {
angang__.nangang = nangang;
angang__.iaa = iaa;
angang__.kaa = kaa;
}

angbnd_st angbnd__;

void set_angbnd_data_ (int* nangle, int* iang, double* ak, double* anat, double* afld) {
angbnd__.nangle = nangle;
angbnd__.iang = iang;
angbnd__.ak = ak;
angbnd__.anat = anat;
angbnd__.afld = afld;
}

angpot_st angpot__;

void set_angpot_data_ (double* angunit, double* stbnunit, double* aaunit, double* opbunit, double* opdunit, double* cang, double* qang, double* pang, double* sang, double* copb, double* qopb, double* popb, double* sopb, double* copd, double* qopd, double* popd, double* sopd, char* opbtyp, char* angtyp) {
angpot__.angunit = angunit;
angpot__.stbnunit = stbnunit;
angpot__.aaunit = aaunit;
angpot__.opbunit = opbunit;
angpot__.opdunit = opdunit;
angpot__.cang = cang;
angpot__.qang = qang;
angpot__.pang = pang;
angpot__.sang = sang;
angpot__.copb = copb;
angpot__.qopb = qopb;
angpot__.popb = popb;
angpot__.sopb = sopb;
angpot__.copd = copd;
angpot__.qopd = qopd;
angpot__.popd = popd;
angpot__.sopd = sopd;
angpot__.opbtyp = opbtyp;
angpot__.angtyp = angtyp;
}

angtor_st angtor__;

void set_angtor_data_ (int* nangtor, int* iat, double* kant) {
angtor__.nangtor = nangtor;
angtor__.iat = iat;
angtor__.kant = kant;
}

argue_st argue__;

void set_argue_data_ (int* maxarg, int* narg, int* listarg, char* arg) {
argue__.maxarg = *maxarg;
argue__.narg = narg;
argue__.listarg = listarg;
argue__.arg = arg;
}

ascii_st ascii__;

void set_ascii_data_ (int* null, int* tab, int* linefeed, int* formfeed, int* carriage, int* escape, int* space, int* exclamation, int* quote, int* pound, int* dollar, int* percent, int* ampersand, int* apostrophe, int* asterisk, int* plus, int* comma, int* minus, int* period, int* frontslash, int* colon, int* semicolon, int* equal, int* question, int* atsign, int* backslash, int* caret, int* underbar, int* vertical, int* tilde) {
ascii__.null = *null;
ascii__.tab = *tab;
ascii__.linefeed = *linefeed;
ascii__.formfeed = *formfeed;
ascii__.carriage = *carriage;
ascii__.escape = *escape;
ascii__.space = *space;
ascii__.exclamation = *exclamation;
ascii__.quote = *quote;
ascii__.pound = *pound;
ascii__.dollar = *dollar;
ascii__.percent = *percent;
ascii__.ampersand = *ampersand;
ascii__.apostrophe = *apostrophe;
ascii__.asterisk = *asterisk;
ascii__.plus = *plus;
ascii__.comma = *comma;
ascii__.minus = *minus;
ascii__.period = *period;
ascii__.frontslash = *frontslash;
ascii__.colon = *colon;
ascii__.semicolon = *semicolon;
ascii__.equal = *equal;
ascii__.question = *question;
ascii__.atsign = *atsign;
ascii__.backslash = *backslash;
ascii__.caret = *caret;
ascii__.underbar = *underbar;
ascii__.vertical = *vertical;
ascii__.tilde = *tilde;
}

atmlst_st atmlst__;

void set_atmlst_data_ (int* bndlist, int* anglist) {
atmlst__.bndlist = bndlist;
atmlst__.anglist = anglist;
}

atomid_st atomid__;

void set_atomid_data_ (int* tag, int* class_, int* atomic, int* valence, double* mass, char* name, char* story) {
atomid__.tag = tag;
atomid__.class_ = class_;
atomid__.atomic = atomic;
atomid__.valence = valence;
atomid__.mass = mass;
atomid__.name = name;
atomid__.story = story;
}

atoms_st atoms__;

void set_atoms_data_ (int* n, int* type, double* x, double* y, double* z) {
atoms__.n = n;
atoms__.type = type;
atoms__.x = x;
atoms__.y = y;
atoms__.z = z;
}

bath_st bath__;

void set_bath_data_ (int* maxnose, int* voltrial, double* kelvin, double* atmsph, double* tautemp, double* taupres, double* compress, double* collide, double* eta, double* volmove, double* vbar, double* qbar, double* gbar, double* vnh, double* qnh, double* gnh, int* isothermal, int* isobaric, int* anisotrop, char* volscale, char* barostat, char* thermostat) {
bath__.maxnose = *maxnose;
bath__.voltrial = voltrial;
bath__.kelvin = kelvin;
bath__.atmsph = atmsph;
bath__.tautemp = tautemp;
bath__.taupres = taupres;
bath__.compress = compress;
bath__.collide = collide;
bath__.eta = eta;
bath__.volmove = volmove;
bath__.vbar = vbar;
bath__.qbar = qbar;
bath__.gbar = gbar;
bath__.vnh = vnh;
bath__.qnh = qnh;
bath__.gnh = gnh;
bath__.isothermal = isothermal;
bath__.isobaric = isobaric;
bath__.anisotrop = anisotrop;
bath__.volscale = volscale;
bath__.barostat = barostat;
bath__.thermostat = thermostat;
}

bitor_st bitor__;

void set_bitor_data_ (int* nbitor, int* ibitor) {
bitor__.nbitor = nbitor;
bitor__.ibitor = ibitor;
}

bndpot_st bndpot__;

void set_bndpot_data_ (double* cbnd, double* qbnd, double* bndunit, char* bndtyp) {
bndpot__.cbnd = cbnd;
bndpot__.qbnd = qbnd;
bndpot__.bndunit = bndunit;
bndpot__.bndtyp = bndtyp;
}

bndstr_st bndstr__;

void set_bndstr_data_ (int* nbond, int* ibnd, double* bk, double* bl) {
bndstr__.nbond = nbond;
bndstr__.ibnd = ibnd;
bndstr__.bk = bk;
bndstr__.bl = bl;
}

bound_st bound__;

void set_bound_data_ (double* polycut, double* polycut2, int* use_bounds, int* use_replica, int* use_polymer) {
bound__.polycut = polycut;
bound__.polycut2 = polycut2;
bound__.use_bounds = use_bounds;
bound__.use_replica = use_replica;
bound__.use_polymer = use_polymer;
}

boxes_st boxes__;

void set_boxes_data_ (double* xbox, double* ybox, double* zbox, double* alpha, double* beta, double* gamma, double* xbox2, double* ybox2, double* zbox2, double* box34, double* volbox, double* beta_sin, double* beta_cos, double* gamma_sin, double* gamma_cos, double* beta_term, double* gamma_term, double* lvec, double* recip, int* orthogonal, int* monoclinic, int* triclinic, int* octahedron, char* spacegrp) {
boxes__.xbox = xbox;
boxes__.ybox = ybox;
boxes__.zbox = zbox;
boxes__.alpha = alpha;
boxes__.beta = beta;
boxes__.gamma = gamma;
boxes__.xbox2 = xbox2;
boxes__.ybox2 = ybox2;
boxes__.zbox2 = zbox2;
boxes__.box34 = box34;
boxes__.volbox = volbox;
boxes__.beta_sin = beta_sin;
boxes__.beta_cos = beta_cos;
boxes__.gamma_sin = gamma_sin;
boxes__.gamma_cos = gamma_cos;
boxes__.beta_term = beta_term;
boxes__.gamma_term = gamma_term;
boxes__.lvec = lvec;
boxes__.recip = recip;
boxes__.orthogonal = orthogonal;
boxes__.monoclinic = monoclinic;
boxes__.triclinic = triclinic;
boxes__.octahedron = octahedron;
boxes__.spacegrp = spacegrp;
}

cell_st cell__;

void set_cell_data_ (int* ncell, int* icell, double* xcell, double* ycell, double* zcell, double* xcell2, double* ycell2, double* zcell2) {
cell__.ncell = ncell;
cell__.icell = icell;
cell__.xcell = xcell;
cell__.ycell = ycell;
cell__.zcell = zcell;
cell__.xcell2 = xcell2;
cell__.ycell2 = ycell2;
cell__.zcell2 = zcell2;
}

charge_st charge__;

void set_charge_data_ (int* nion, int* iion, int* jion, int* kion, double* pchg) {
charge__.nion = nion;
charge__.iion = iion;
charge__.jion = jion;
charge__.kion = kion;
charge__.pchg = pchg;
}

chgpen_st chgpen__;

void set_chgpen_data_ (int* ncp, double* pcore, double* pval, double* palpha) {
chgpen__.ncp = ncp;
chgpen__.pcore = pcore;
chgpen__.pval = pval;
chgpen__.palpha = palpha;
}

chgpot_st chgpot__;

void set_chgpot_data_ (double* electric, double* dielec, double* ebuffer, double* c2scale, double* c3scale, double* c4scale, double* c5scale, int* neutnbr, int* neutcut) {
chgpot__.electric = electric;
chgpot__.dielec = dielec;
chgpot__.ebuffer = ebuffer;
chgpot__.c2scale = c2scale;
chgpot__.c3scale = c3scale;
chgpot__.c4scale = c4scale;
chgpot__.c5scale = c5scale;
chgpot__.neutnbr = neutnbr;
chgpot__.neutcut = neutcut;
}

chgtrn_st chgtrn__;

void set_chgtrn_data_ (int* nct, double* chgct, double* dmpct) {
chgtrn__.nct = nct;
chgtrn__.chgct = chgct;
chgtrn__.dmpct = dmpct;
}

chrono_st chrono__;

void set_chrono_data_ (double* twall, double* tcpu) {
chrono__.twall = twall;
chrono__.tcpu = tcpu;
}

chunks_st chunks__;

void set_chunks_data_ (int* nchunk, int* nchk1, int* nchk2, int* nchk3, int* ngrd1, int* ngrd2, int* ngrd3, int* nlpts, int* nrpts, int* grdoff, int* pmetable) {
chunks__.nchunk = nchunk;
chunks__.nchk1 = nchk1;
chunks__.nchk2 = nchk2;
chunks__.nchk3 = nchk3;
chunks__.ngrd1 = ngrd1;
chunks__.ngrd2 = ngrd2;
chunks__.ngrd3 = ngrd3;
chunks__.nlpts = nlpts;
chunks__.nrpts = nrpts;
chunks__.grdoff = grdoff;
chunks__.pmetable = pmetable;
}

couple_st couple__;

void set_couple_data_ (int* n12, int* n13, int* n14, int* n15, int* i12, int* i13, int* i14, int* i15) {
couple__.n12 = n12;
couple__.n13 = n13;
couple__.n14 = n14;
couple__.n15 = n15;
couple__.i12 = i12;
couple__.i13 = i13;
couple__.i14 = i14;
couple__.i15 = i15;
}

ctrpot_st ctrpot__;

void set_ctrpot_data_ (char* ctrntyp) {
ctrpot__.ctrntyp = ctrntyp;
}

deriv_st deriv__;

void set_deriv_data_ (double* desum, double* deb, double* dea, double* deba, double* deub, double* deaa, double* deopb, double* deopd, double* deid, double* deit, double* det, double* dept, double* debt, double* deat, double* dett, double* dev, double* der, double* dedsp, double* dec, double* decd, double* ded, double* dem, double* dep, double* dect, double* derxf, double* des, double* delf, double* deg, double* dex) {
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

dipole_st dipole__;

void set_dipole_data_ (int* ndipole, int* idpl, double* bdpl, double* sdpl) {
dipole__.ndipole = ndipole;
dipole__.idpl = idpl;
dipole__.bdpl = bdpl;
dipole__.sdpl = sdpl;
}

disgeo_st disgeo__;

void set_disgeo_data_ (double* vdwmax, double* compact, double* pathmax, double* dbnd, double* georad, int* use_invert, int* use_anneal) {
disgeo__.vdwmax = vdwmax;
disgeo__.compact = compact;
disgeo__.pathmax = pathmax;
disgeo__.dbnd = dbnd;
disgeo__.georad = georad;
disgeo__.use_invert = use_invert;
disgeo__.use_anneal = use_anneal;
}

disp_st disp__;

void set_disp_data_ (int* ndisp, int* idisp, double* csixpr, double* csix, double* adisp) {
disp__.ndisp = ndisp;
disp__.idisp = idisp;
disp__.csixpr = csixpr;
disp__.csix = csix;
disp__.adisp = adisp;
}

dma_st dma__;

void set_dma_data_ (double* mp, double* dpx, double* dpy, double* dpz, double* q20, double* q21c, double* q21s, double* q22c, double* q22s) {
dma__.mp = mp;
dma__.dpx = dpx;
dma__.dpy = dpy;
dma__.dpz = dpz;
dma__.q20 = q20;
dma__.q21c = q21c;
dma__.q21s = q21s;
dma__.q22c = q22c;
dma__.q22s = q22s;
}

domega_st domega__;

void set_domega_data_ (double* tesum, double* teb, double* tea, double* teba, double* teub, double* teaa, double* teopb, double* teopd, double* teid, double* teit, double* tet, double* tept, double* tebt, double* teat, double* tett, double* tev, double* ter, double* tedsp, double* tec, double* tecd, double* ted, double* tem, double* tep, double* tect, double* terxf, double* tes, double* telf, double* teg, double* tex) {
domega__.tesum = tesum;
domega__.teb = teb;
domega__.tea = tea;
domega__.teba = teba;
domega__.teub = teub;
domega__.teaa = teaa;
domega__.teopb = teopb;
domega__.teopd = teopd;
domega__.teid = teid;
domega__.teit = teit;
domega__.tet = tet;
domega__.tept = tept;
domega__.tebt = tebt;
domega__.teat = teat;
domega__.tett = tett;
domega__.tev = tev;
domega__.ter = ter;
domega__.tedsp = tedsp;
domega__.tec = tec;
domega__.tecd = tecd;
domega__.ted = ted;
domega__.tem = tem;
domega__.tep = tep;
domega__.tect = tect;
domega__.terxf = terxf;
domega__.tes = tes;
domega__.telf = telf;
domega__.teg = teg;
domega__.tex = tex;
}

dsppot_st dsppot__;

void set_dsppot_data_ (double* dsp2scale, double* dsp3scale, double* dsp4scale, double* dsp5scale, int* use_dcorr) {
dsppot__.dsp2scale = dsp2scale;
dsppot__.dsp3scale = dsp3scale;
dsppot__.dsp4scale = dsp4scale;
dsppot__.dsp5scale = dsp5scale;
dsppot__.use_dcorr = use_dcorr;
}

energi_st energi__;

void set_energi_data_ (double* esum, double* eb, double* ea, double* eba, double* eub, double* eaa, double* eopb, double* eopd, double* eid, double* eit, double* et, double* ept, double* ebt, double* eat, double* ett, double* ev, double* er, double* edsp, double* ec, double* ecd, double* ed, double* em, double* ep, double* ect, double* erxf, double* es, double* elf, double* eg, double* ex) {
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

ewald_st ewald__;

void set_ewald_data_ (double* aewald, double* aeewald, double* apewald, double* adewald, char* boundary) {
ewald__.aewald = aewald;
ewald__.aeewald = aeewald;
ewald__.apewald = apewald;
ewald__.adewald = adewald;
ewald__.boundary = boundary;
}

faces_st faces__;

void set_faces_data_ (int* maxcls, int* maxtt, int* maxt, int* maxp, int* maxv, int* maxen, int* maxfn, int* maxc, int* maxeq, int* maxfs, int* maxfq, int* maxcy, int* mxcyeq, int* mxfqcy, int* na, double* pr, double* ar, double* axyz, int* skip, int* nosurf, int* afree, int* abur, int* cls, int* clst, int* acls, int* ntt, int* ttfe, int* ttle, int* enext, int* tta, int* ttbur, int* ttfree, int* nt, int* tfe, int* ta, double* tr, double* t, double* tax, int* tfree, int* np, int* pa, double* p, int* nv, int* va, int* vp, double* vxyz, int* nen, int* nfn, int* env, int* fnen, int* nc, int* ca, int* ct, double* cr, double* c, int* neq, int* eqc, int* eqv, int* afe, int* ale, int* eqnext, int* nfs, int* fsen, int* fseq, int* ncy, int* cyneq, int* cyeq, int* nfq, int* fqa, int* fqncy, int* fqcy) {
faces__.maxcls = maxcls;
faces__.maxtt = maxtt;
faces__.maxt = maxt;
faces__.maxp = maxp;
faces__.maxv = maxv;
faces__.maxen = maxen;
faces__.maxfn = maxfn;
faces__.maxc = maxc;
faces__.maxeq = maxeq;
faces__.maxfs = maxfs;
faces__.maxfq = maxfq;
faces__.maxcy = maxcy;
faces__.mxcyeq = mxcyeq;
faces__.mxfqcy = mxfqcy;
faces__.na = na;
faces__.pr = pr;
faces__.ar = ar;
faces__.axyz = axyz;
faces__.skip = skip;
faces__.nosurf = nosurf;
faces__.afree = afree;
faces__.abur = abur;
faces__.cls = cls;
faces__.clst = clst;
faces__.acls = acls;
faces__.ntt = ntt;
faces__.ttfe = ttfe;
faces__.ttle = ttle;
faces__.enext = enext;
faces__.tta = tta;
faces__.ttbur = ttbur;
faces__.ttfree = ttfree;
faces__.nt = nt;
faces__.tfe = tfe;
faces__.ta = ta;
faces__.tr = tr;
faces__.t = t;
faces__.tax = tax;
faces__.tfree = tfree;
faces__.np = np;
faces__.pa = pa;
faces__.p = p;
faces__.nv = nv;
faces__.va = va;
faces__.vp = vp;
faces__.vxyz = vxyz;
faces__.nen = nen;
faces__.nfn = nfn;
faces__.env = env;
faces__.fnen = fnen;
faces__.nc = nc;
faces__.ca = ca;
faces__.ct = ct;
faces__.cr = cr;
faces__.c = c;
faces__.neq = neq;
faces__.eqc = eqc;
faces__.eqv = eqv;
faces__.afe = afe;
faces__.ale = ale;
faces__.eqnext = eqnext;
faces__.nfs = nfs;
faces__.fsen = fsen;
faces__.fseq = fseq;
faces__.ncy = ncy;
faces__.cyneq = cyneq;
faces__.cyeq = cyeq;
faces__.nfq = nfq;
faces__.fqa = fqa;
faces__.fqncy = fqncy;
faces__.fqcy = fqcy;
}

fft_st fft__;

void set_fft_data_ (int* maxprime, int* iprime, unsigned long long* planf, unsigned long long* planb, double* ffttable, char* ffttyp) {
fft__.maxprime = *maxprime;
fft__.iprime = iprime;
fft__.planf = planf;
fft__.planb = planb;
fft__.ffttable = ffttable;
fft__.ffttyp = ffttyp;
}

fields_st fields__;

void set_fields_data_ (int* maxbio, int* biotyp, char* forcefield) {
fields__.maxbio = *maxbio;
fields__.biotyp = biotyp;
fields__.forcefield = forcefield;
}

files_st files__;

void set_files_data_ (int* nprior, int* ldir, int* leng, char* filename, char* outfile) {
files__.nprior = nprior;
files__.ldir = ldir;
files__.leng = leng;
files__.filename = filename;
files__.outfile = outfile;
}

fracs_st fracs__;

void set_fracs_data_ (double* xfrac, double* yfrac, double* zfrac) {
fracs__.xfrac = xfrac;
fracs__.yfrac = yfrac;
fracs__.zfrac = zfrac;
}

freeze_st freeze__;

void set_freeze_data_ (int* nrat, int* nratx, int* iratx, int* kratx, int* irat, double* rateps, double* krat, int* use_rattle, int* ratimage) {
freeze__.nrat = nrat;
freeze__.nratx = nratx;
freeze__.iratx = iratx;
freeze__.kratx = kratx;
freeze__.irat = irat;
freeze__.rateps = rateps;
freeze__.krat = krat;
freeze__.use_rattle = use_rattle;
freeze__.ratimage = ratimage;
}

gkstuf_st gkstuf__;

void set_gkstuf_data_ (double* gkc, double* gkr) {
gkstuf__.gkc = gkc;
gkstuf__.gkr = gkr;
}

group_st group__;

void set_group_data_ (int* ngrp, int* kgrp, int* grplist, int* igrp, double* grpmass, double* wgrp, int* use_group, int* use_intra, int* use_inter) {
group__.ngrp = ngrp;
group__.kgrp = kgrp;
group__.grplist = grplist;
group__.igrp = igrp;
group__.grpmass = grpmass;
group__.wgrp = wgrp;
group__.use_group = use_group;
group__.use_intra = use_intra;
group__.use_inter = use_inter;
}

hescut_st hescut__;

void set_hescut_data_ (double* hesscut) {
hescut__.hesscut = hesscut;
}

hessn_st hessn__;

void set_hessn_data_ (double* hessx, double* hessy, double* hessz) {
hessn__.hessx = hessx;
hessn__.hessy = hessy;
hessn__.hessz = hessz;
}

hpmf_st hpmf__;

void set_hpmf_data_ (double* rcarbon, double* rwater, double* acsurf, double* safact, double* tslope, double* toffset, double* hpmfcut, double* hd1, double* hd2, double* hd3, double* hc1, double* hc2, double* hc3, double* hw1, double* hw2, double* hw3, int* npmf, int* ipmf, double* rpmf, double* acsa) {
hpmf__.rcarbon = *rcarbon;
hpmf__.rwater = *rwater;
hpmf__.acsurf = *acsurf;
hpmf__.safact = *safact;
hpmf__.tslope = *tslope;
hpmf__.toffset = *toffset;
hpmf__.hpmfcut = *hpmfcut;
hpmf__.hd1 = *hd1;
hpmf__.hd2 = *hd2;
hpmf__.hd3 = *hd3;
hpmf__.hc1 = *hc1;
hpmf__.hc2 = *hc2;
hpmf__.hc3 = *hc3;
hpmf__.hw1 = *hw1;
hpmf__.hw2 = *hw2;
hpmf__.hw3 = *hw3;
hpmf__.npmf = npmf;
hpmf__.ipmf = ipmf;
hpmf__.rpmf = rpmf;
hpmf__.acsa = acsa;
}

ielscf_st ielscf__;

void set_ielscf_data_ (int* nfree_aux, double* tautemp_aux, double* kelvin_aux, double* uaux, double* upaux, double* vaux, double* vpaux, double* aaux, double* apaux, int* use_ielscf) {
ielscf__.nfree_aux = nfree_aux;
ielscf__.tautemp_aux = tautemp_aux;
ielscf__.kelvin_aux = kelvin_aux;
ielscf__.uaux = uaux;
ielscf__.upaux = upaux;
ielscf__.vaux = vaux;
ielscf__.vpaux = vpaux;
ielscf__.aaux = aaux;
ielscf__.apaux = apaux;
ielscf__.use_ielscf = use_ielscf;
}

improp_st improp__;

void set_improp_data_ (int* niprop, int* iiprop, double* kprop, double* vprop) {
improp__.niprop = niprop;
improp__.iiprop = iiprop;
improp__.kprop = kprop;
improp__.vprop = vprop;
}

imptor_st imptor__;

void set_imptor_data_ (int* nitors, int* iitors, double* itors1, double* itors2, double* itors3) {
imptor__.nitors = nitors;
imptor__.iitors = iitors;
imptor__.itors1 = itors1;
imptor__.itors2 = itors2;
imptor__.itors3 = itors3;
}

inform_st inform__;

void set_inform_data_ (int* maxask, int* digits, int* iprint, int* iwrite, int* isend, int* silent, int* verbose, int* debug, int* holdup, int* abort) {
inform__.maxask = *maxask;
inform__.digits = digits;
inform__.iprint = iprint;
inform__.iwrite = iwrite;
inform__.isend = isend;
inform__.silent = silent;
inform__.verbose = verbose;
inform__.debug = debug;
inform__.holdup = holdup;
inform__.abort = abort;
}

inter_st inter__;

void set_inter_data_ (double* einter) {
inter__.einter = einter;
}

iounit_st iounit__;

void set_iounit_data_ (int* input, int* iout) {
iounit__.input = input;
iounit__.iout = iout;
}

kanang_st kanang__;

void set_kanang_data_ (double* anan) {
kanang__.anan = anan;
}

kangs_st kangs__;

void set_kangs_data_ (int* maxna, int* maxna5, int* maxna4, int* maxna3, int* maxnap, int* maxnaf, double* acon, double* acon5, double* acon4, double* acon3, double* aconp, double* aconf, double* ang, double* ang5, double* ang4, double* ang3, double* angp, double* angf, char* ka, char* ka5, char* ka4, char* ka3, char* kap, char* kaf) {
kangs__.maxna = *maxna;
kangs__.maxna5 = *maxna5;
kangs__.maxna4 = *maxna4;
kangs__.maxna3 = *maxna3;
kangs__.maxnap = *maxnap;
kangs__.maxnaf = *maxnaf;
kangs__.acon = acon;
kangs__.acon5 = acon5;
kangs__.acon4 = acon4;
kangs__.acon3 = acon3;
kangs__.aconp = aconp;
kangs__.aconf = aconf;
kangs__.ang = ang;
kangs__.ang5 = ang5;
kangs__.ang4 = ang4;
kangs__.ang3 = ang3;
kangs__.angp = angp;
kangs__.angf = angf;
kangs__.ka = ka;
kangs__.ka5 = ka5;
kangs__.ka4 = ka4;
kangs__.ka3 = ka3;
kangs__.kap = kap;
kangs__.kaf = kaf;
}

kantor_st kantor__;

void set_kantor_data_ (int* maxnat, double* atcon, char* kat) {
kantor__.maxnat = *maxnat;
kantor__.atcon = atcon;
kantor__.kat = kat;
}

katoms_st katoms__;

void set_katoms_data_ (int* atmcls, int* atmnum, int* ligand, double* weight, char* symbol, char* describe) {
katoms__.atmcls = atmcls;
katoms__.atmnum = atmnum;
katoms__.ligand = ligand;
katoms__.weight = weight;
katoms__.symbol = symbol;
katoms__.describe = describe;
}

kbonds_st kbonds__;

void set_kbonds_data_ (int* maxnb, int* maxnb5, int* maxnb4, int* maxnb3, int* maxnel, double* bcon, double* bcon5, double* bcon4, double* bcon3, double* blen, double* blen5, double* blen4, double* blen3, double* dlen, char* kb, char* kb5, char* kb4, char* kb3, char* kel) {
kbonds__.maxnb = *maxnb;
kbonds__.maxnb5 = *maxnb5;
kbonds__.maxnb4 = *maxnb4;
kbonds__.maxnb3 = *maxnb3;
kbonds__.maxnel = *maxnel;
kbonds__.bcon = bcon;
kbonds__.bcon5 = bcon5;
kbonds__.bcon4 = bcon4;
kbonds__.bcon3 = bcon3;
kbonds__.blen = blen;
kbonds__.blen5 = blen5;
kbonds__.blen4 = blen4;
kbonds__.blen3 = blen3;
kbonds__.dlen = dlen;
kbonds__.kb = kb;
kbonds__.kb5 = kb5;
kbonds__.kb4 = kb4;
kbonds__.kb3 = kb3;
kbonds__.kel = kel;
}

kchrge_st kchrge__;

void set_kchrge_data_ (double* chg) {
kchrge__.chg = chg;
}

kcpen_st kcpen__;

void set_kcpen_data_ (double* cpele, double* cpalp) {
kcpen__.cpele = cpele;
kcpen__.cpalp = cpalp;
}

kctrn_st kctrn__;

void set_kctrn_data_ (double* ctchg, double* ctdmp) {
kctrn__.ctchg = ctchg;
kctrn__.ctdmp = ctdmp;
}

kdipol_st kdipol__;

void set_kdipol_data_ (int* maxnd, int* maxnd5, int* maxnd4, int* maxnd3, double* dpl, double* dpl5, double* dpl4, double* dpl3, double* pos, double* pos5, double* pos4, double* pos3, char* kd, char* kd5, char* kd4, char* kd3) {
kdipol__.maxnd = *maxnd;
kdipol__.maxnd5 = *maxnd5;
kdipol__.maxnd4 = *maxnd4;
kdipol__.maxnd3 = *maxnd3;
kdipol__.dpl = dpl;
kdipol__.dpl5 = dpl5;
kdipol__.dpl4 = dpl4;
kdipol__.dpl3 = dpl3;
kdipol__.pos = pos;
kdipol__.pos5 = pos5;
kdipol__.pos4 = pos4;
kdipol__.pos3 = pos3;
kdipol__.kd = kd;
kdipol__.kd5 = kd5;
kdipol__.kd4 = kd4;
kdipol__.kd3 = kd3;
}

kdsp_st kdsp__;

void set_kdsp_data_ (double* dspsix, double* dspdmp) {
kdsp__.dspsix = dspsix;
kdsp__.dspdmp = dspdmp;
}

keys_st keys__;

void set_keys_data_ (int* maxkey, int* nkey, char* keyline) {
keys__.maxkey = *maxkey;
keys__.nkey = nkey;
keys__.keyline = keyline;
}

khbond_st khbond__;

void set_khbond_data_ (int* maxnhb, double* radhb, double* epshb, char* khb) {
khbond__.maxnhb = *maxnhb;
khbond__.radhb = radhb;
khbond__.epshb = epshb;
khbond__.khb = khb;
}

kiprop_st kiprop__;

void set_kiprop_data_ (int* maxndi, double* dcon, double* tdi, char* kdi) {
kiprop__.maxndi = *maxndi;
kiprop__.dcon = dcon;
kiprop__.tdi = tdi;
kiprop__.kdi = kdi;
}

kitors_st kitors__;

void set_kitors_data_ (int* maxnti, double* ti1, double* ti2, double* ti3, char* kti) {
kitors__.maxnti = *maxnti;
kitors__.ti1 = ti1;
kitors__.ti2 = ti2;
kitors__.ti3 = ti3;
kitors__.kti = kti;
}

kmulti_st kmulti__;

void set_kmulti_data_ (int* maxnmp, double* multip, char* mpaxis, char* kmp) {
kmulti__.maxnmp = *maxnmp;
kmulti__.multip = multip;
kmulti__.mpaxis = mpaxis;
kmulti__.kmp = kmp;
}

kopbnd_st kopbnd__;

void set_kopbnd_data_ (int* maxnopb, double* opbn, char* kopb) {
kopbnd__.maxnopb = *maxnopb;
kopbnd__.opbn = opbn;
kopbnd__.kopb = kopb;
}

kopdst_st kopdst__;

void set_kopdst_data_ (int* maxnopd, double* opds, char* kopd) {
kopdst__.maxnopd = *maxnopd;
kopdst__.opds = opds;
kopdst__.kopd = kopd;
}

korbs_st korbs__;

void set_korbs_data_ (int* maxnpi, int* maxnpi5, int* maxnpi4, double* sslope, double* sslope5, double* sslope4, double* tslope, double* tslope5, double* tslope4, double* electron, double* ionize, double* repulse, char* kpi, char* kpi5, char* kpi4) {
korbs__.maxnpi = *maxnpi;
korbs__.maxnpi5 = *maxnpi5;
korbs__.maxnpi4 = *maxnpi4;
korbs__.sslope = sslope;
korbs__.sslope5 = sslope5;
korbs__.sslope4 = sslope4;
korbs__.tslope = tslope;
korbs__.tslope5 = tslope5;
korbs__.tslope4 = tslope4;
korbs__.electron = electron;
korbs__.ionize = ionize;
korbs__.repulse = repulse;
korbs__.kpi = kpi;
korbs__.kpi5 = kpi5;
korbs__.kpi4 = kpi4;
}

kpitor_st kpitor__;

void set_kpitor_data_ (int* maxnpt, double* ptcon, char* kpt) {
kpitor__.maxnpt = *maxnpt;
kpitor__.ptcon = ptcon;
kpitor__.kpt = kpt;
}

kpolr_st kpolr__;

void set_kpolr_data_ (int* pgrp, double* polr, double* athl, double* ddir) {
kpolr__.pgrp = pgrp;
kpolr__.polr = polr;
kpolr__.athl = athl;
kpolr__.ddir = ddir;
}

krepl_st krepl__;

void set_krepl_data_ (double* prsiz, double* prdmp, double* prele) {
krepl__.prsiz = prsiz;
krepl__.prdmp = prdmp;
krepl__.prele = prele;
}

kstbnd_st kstbnd__;

void set_kstbnd_data_ (int* maxnsb, double* stbn, char* ksb) {
kstbnd__.maxnsb = *maxnsb;
kstbnd__.stbn = stbn;
kstbnd__.ksb = ksb;
}

ksttor_st ksttor__;

void set_ksttor_data_ (int* maxnbt, double* btcon, char* kbt) {
ksttor__.maxnbt = *maxnbt;
ksttor__.btcon = btcon;
ksttor__.kbt = kbt;
}

ktorsn_st ktorsn__;

void set_ktorsn_data_ (int* maxnt, int* maxnt5, int* maxnt4, double* t1, double* t2, double* t3, double* t4, double* t5, double* t6, double* t15, double* t25, double* t35, double* t45, double* t55, double* t65, double* t14, double* t24, double* t34, double* t44, double* t54, double* t64, char* kt, char* kt5, char* kt4) {
ktorsn__.maxnt = *maxnt;
ktorsn__.maxnt5 = *maxnt5;
ktorsn__.maxnt4 = *maxnt4;
ktorsn__.t1 = t1;
ktorsn__.t2 = t2;
ktorsn__.t3 = t3;
ktorsn__.t4 = t4;
ktorsn__.t5 = t5;
ktorsn__.t6 = t6;
ktorsn__.t15 = t15;
ktorsn__.t25 = t25;
ktorsn__.t35 = t35;
ktorsn__.t45 = t45;
ktorsn__.t55 = t55;
ktorsn__.t65 = t65;
ktorsn__.t14 = t14;
ktorsn__.t24 = t24;
ktorsn__.t34 = t34;
ktorsn__.t44 = t44;
ktorsn__.t54 = t54;
ktorsn__.t64 = t64;
ktorsn__.kt = kt;
ktorsn__.kt5 = kt5;
ktorsn__.kt4 = kt4;
}

ktrtor_st ktrtor__;

void set_ktrtor_data_ (int* maxntt, int* maxtgrd, int* maxtgrd2, int* tnx, int* tny, double* ttx, double* tty, double* tbf, double* tbx, double* tby, double* tbxy, char* ktt) {
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

kurybr_st kurybr__;

void set_kurybr_data_ (int* maxnu, double* ucon, double* dst13, char* ku) {
kurybr__.maxnu = *maxnu;
kurybr__.ucon = ucon;
kurybr__.dst13 = dst13;
kurybr__.ku = ku;
}

kvdwpr_st kvdwpr__;

void set_kvdwpr_data_ (int* maxnvp, double* radpr, double* epspr, char* kvpr) {
kvdwpr__.maxnvp = *maxnvp;
kvdwpr__.radpr = radpr;
kvdwpr__.epspr = epspr;
kvdwpr__.kvpr = kvpr;
}

kvdws_st kvdws__;

void set_kvdws_data_ (double* rad, double* eps, double* rad4, double* eps4, double* reduct) {
kvdws__.rad = rad;
kvdws__.eps = eps;
kvdws__.rad4 = rad4;
kvdws__.eps4 = eps4;
kvdws__.reduct = reduct;
}

light_st light__;

void set_light_data_ (int* nlight, int* kbx, int* kby, int* kbz, int* kex, int* key, int* kez, int* locx, int* locy, int* locz, int* rgx, int* rgy, int* rgz) {
light__.nlight = nlight;
light__.kbx = kbx;
light__.kby = kby;
light__.kbz = kbz;
light__.kex = kex;
light__.key = key;
light__.kez = kez;
light__.locx = locx;
light__.locy = locy;
light__.locz = locz;
light__.rgx = rgx;
light__.rgy = rgy;
light__.rgz = rgz;
}

limits_st limits__;

void set_limits_data_ (double* vdwcut, double* repcut, double* dispcut, double* chgcut, double* dplcut, double* mpolecut, double* ctrncut, double* vdwtaper, double* reptaper, double* disptaper, double* chgtaper, double* dpltaper, double* mpoletaper, double* ctrntaper, double* ewaldcut, double* dewaldcut, double* usolvcut, int* use_ewald, int* use_dewald, int* use_lights, int* use_list, int* use_vlist, int* use_dlist, int* use_clist, int* use_mlist, int* use_ulist) {
limits__.vdwcut = vdwcut;
limits__.repcut = repcut;
limits__.dispcut = dispcut;
limits__.chgcut = chgcut;
limits__.dplcut = dplcut;
limits__.mpolecut = mpolecut;
limits__.ctrncut = ctrncut;
limits__.vdwtaper = vdwtaper;
limits__.reptaper = reptaper;
limits__.disptaper = disptaper;
limits__.chgtaper = chgtaper;
limits__.dpltaper = dpltaper;
limits__.mpoletaper = mpoletaper;
limits__.ctrntaper = ctrntaper;
limits__.ewaldcut = ewaldcut;
limits__.dewaldcut = dewaldcut;
limits__.usolvcut = usolvcut;
limits__.use_ewald = use_ewald;
limits__.use_dewald = use_dewald;
limits__.use_lights = use_lights;
limits__.use_list = use_list;
limits__.use_vlist = use_vlist;
limits__.use_dlist = use_dlist;
limits__.use_clist = use_clist;
limits__.use_mlist = use_mlist;
limits__.use_ulist = use_ulist;
}

linmin_st linmin__;

void set_linmin_data_ (int* intmax, double* stpmin, double* stpmax, double* cappa, double* slpmax, double* angmax) {
linmin__.intmax = intmax;
linmin__.stpmin = stpmin;
linmin__.stpmax = stpmax;
linmin__.cappa = cappa;
linmin__.slpmax = slpmax;
linmin__.angmax = angmax;
}

math_st math__;

void set_math_data_ (double* pi, double* elog, double* radian, double* logten, double* twosix, double* sqrtpi, double* sqrttwo, double* sqrtthree) {
math__.pi = *pi;
math__.elog = *elog;
math__.radian = *radian;
math__.logten = *logten;
math__.twosix = *twosix;
math__.sqrtpi = *sqrtpi;
math__.sqrttwo = *sqrttwo;
math__.sqrtthree = *sqrtthree;
}

mdstuf_st mdstuf__;

void set_mdstuf_data_ (int* nfree, int* irest, int* bmnmix, double* arespa, int* dorest, char* integrate) {
mdstuf__.nfree = nfree;
mdstuf__.irest = irest;
mdstuf__.bmnmix = bmnmix;
mdstuf__.arespa = arespa;
mdstuf__.dorest = dorest;
mdstuf__.integrate = integrate;
}

merck_st merck__;

void set_merck_data_ (int* nlignes, int* bt_1, int* eqclass, int* crd, int* val, int* pilp, int* mltb, int* arom, int* lin, int* sbmb, int* mmffarom, int* mmffaromc, int* mmffaroma, double* rad0, double* paulel, double* r0ref, double* kbref, double* mmff_kb, double* mmff_kb1, double* mmff_b0, double* mmff_b1, double* mmff_ka, double* mmff_ka1, double* mmff_ka2, double* mmff_ka3, double* mmff_ka4, double* mmff_ka5, double* mmff_ka6, double* mmff_ka7, double* mmff_ka8, double* mmff_ang0, double* mmff_ang1, double* mmff_ang2, double* mmff_ang3, double* mmff_ang4, double* mmff_ang5, double* mmff_ang6, double* mmff_ang7, double* mmff_ang8, double* stbn_abc, double* stbn_cba, double* stbn_abc1, double* stbn_cba1, double* stbn_abc2, double* stbn_cba2, double* stbn_abc3, double* stbn_cba3, double* stbn_abc4, double* stbn_cba4, double* stbn_abc5, double* stbn_cba5, double* stbn_abc6, double* stbn_cba6, double* stbn_abc7, double* stbn_cba7, double* stbn_abc8, double* stbn_cba8, double* stbn_abc9, double* stbn_cba9, double* stbn_abc10, double* stbn_cba10, double* stbn_abc11, double* stbn_cba11, double* defstbn_abc, double* defstbn_cba, double* t1_1, double* t2_1, double* t3_1, double* t1_2, double* t2_2, double* t3_2, char* kt_1, char* kt_2, double* g, double* alph, double* nn, char* da, double* bci, double* bci_1, double* pbci, double* fcadj) {
merck__.nlignes = nlignes;
merck__.bt_1 = bt_1;
merck__.eqclass = eqclass;
merck__.crd = crd;
merck__.val = val;
merck__.pilp = pilp;
merck__.mltb = mltb;
merck__.arom = arom;
merck__.lin = lin;
merck__.sbmb = sbmb;
merck__.mmffarom = mmffarom;
merck__.mmffaromc = mmffaromc;
merck__.mmffaroma = mmffaroma;
merck__.rad0 = rad0;
merck__.paulel = paulel;
merck__.r0ref = r0ref;
merck__.kbref = kbref;
merck__.mmff_kb = mmff_kb;
merck__.mmff_kb1 = mmff_kb1;
merck__.mmff_b0 = mmff_b0;
merck__.mmff_b1 = mmff_b1;
merck__.mmff_ka = mmff_ka;
merck__.mmff_ka1 = mmff_ka1;
merck__.mmff_ka2 = mmff_ka2;
merck__.mmff_ka3 = mmff_ka3;
merck__.mmff_ka4 = mmff_ka4;
merck__.mmff_ka5 = mmff_ka5;
merck__.mmff_ka6 = mmff_ka6;
merck__.mmff_ka7 = mmff_ka7;
merck__.mmff_ka8 = mmff_ka8;
merck__.mmff_ang0 = mmff_ang0;
merck__.mmff_ang1 = mmff_ang1;
merck__.mmff_ang2 = mmff_ang2;
merck__.mmff_ang3 = mmff_ang3;
merck__.mmff_ang4 = mmff_ang4;
merck__.mmff_ang5 = mmff_ang5;
merck__.mmff_ang6 = mmff_ang6;
merck__.mmff_ang7 = mmff_ang7;
merck__.mmff_ang8 = mmff_ang8;
merck__.stbn_abc = stbn_abc;
merck__.stbn_cba = stbn_cba;
merck__.stbn_abc1 = stbn_abc1;
merck__.stbn_cba1 = stbn_cba1;
merck__.stbn_abc2 = stbn_abc2;
merck__.stbn_cba2 = stbn_cba2;
merck__.stbn_abc3 = stbn_abc3;
merck__.stbn_cba3 = stbn_cba3;
merck__.stbn_abc4 = stbn_abc4;
merck__.stbn_cba4 = stbn_cba4;
merck__.stbn_abc5 = stbn_abc5;
merck__.stbn_cba5 = stbn_cba5;
merck__.stbn_abc6 = stbn_abc6;
merck__.stbn_cba6 = stbn_cba6;
merck__.stbn_abc7 = stbn_abc7;
merck__.stbn_cba7 = stbn_cba7;
merck__.stbn_abc8 = stbn_abc8;
merck__.stbn_cba8 = stbn_cba8;
merck__.stbn_abc9 = stbn_abc9;
merck__.stbn_cba9 = stbn_cba9;
merck__.stbn_abc10 = stbn_abc10;
merck__.stbn_cba10 = stbn_cba10;
merck__.stbn_abc11 = stbn_abc11;
merck__.stbn_cba11 = stbn_cba11;
merck__.defstbn_abc = defstbn_abc;
merck__.defstbn_cba = defstbn_cba;
merck__.t1_1 = t1_1;
merck__.t2_1 = t2_1;
merck__.t3_1 = t3_1;
merck__.t1_2 = t1_2;
merck__.t2_2 = t2_2;
merck__.t3_2 = t3_2;
merck__.kt_1 = kt_1;
merck__.kt_2 = kt_2;
merck__.g = g;
merck__.alph = alph;
merck__.nn = nn;
merck__.da = da;
merck__.bci = bci;
merck__.bci_1 = bci_1;
merck__.pbci = pbci;
merck__.fcadj = fcadj;
}

minima_st minima__;

void set_minima_data_ (int* maxiter, int* nextiter, double* fctmin, double* hguess) {
minima__.maxiter = maxiter;
minima__.nextiter = nextiter;
minima__.fctmin = fctmin;
minima__.hguess = hguess;
}

molcul_st molcul__;

void set_molcul_data_ (int* nmol, int* imol, int* kmol, int* molcule, double* totmass, double* molmass) {
molcul__.nmol = nmol;
molcul__.imol = imol;
molcul__.kmol = kmol;
molcul__.molcule = molcule;
molcul__.totmass = totmass;
molcul__.molmass = molmass;
}

moldyn_st moldyn__;

void set_moldyn_data_ (double* v, double* a, double* aalt) {
moldyn__.v = v;
moldyn__.a = a;
moldyn__.aalt = aalt;
}

moment_st moment__;

void set_moment_data_ (double* netchg, double* netdpl, double* netqdp, double* xdpl, double* ydpl, double* zdpl, double* xxqdp, double* xyqdp, double* xzqdp, double* yxqdp, double* yyqdp, double* yzqdp, double* zxqdp, double* zyqdp, double* zzqdp) {
moment__.netchg = netchg;
moment__.netdpl = netdpl;
moment__.netqdp = netqdp;
moment__.xdpl = xdpl;
moment__.ydpl = ydpl;
moment__.zdpl = zdpl;
moment__.xxqdp = xxqdp;
moment__.xyqdp = xyqdp;
moment__.xzqdp = xzqdp;
moment__.yxqdp = yxqdp;
moment__.yyqdp = yyqdp;
moment__.yzqdp = yzqdp;
moment__.zxqdp = zxqdp;
moment__.zyqdp = zyqdp;
moment__.zzqdp = zzqdp;
}

mplpot_st mplpot__;

void set_mplpot_data_ (double* m2scale, double* m3scale, double* m4scale, double* m5scale, int* use_chgpen, char* pentyp) {
mplpot__.m2scale = m2scale;
mplpot__.m3scale = m3scale;
mplpot__.m4scale = m4scale;
mplpot__.m5scale = m5scale;
mplpot__.use_chgpen = use_chgpen;
mplpot__.pentyp = pentyp;
}

mpole_st mpole__;

void set_mpole_data_ (int* maxpole, int* npole, int* ipole, int* polsiz, int* pollist, int* zaxis, int* xaxis, int* yaxis, double* pole, double* rpole, double* spole, double* srpole, char* polaxe) {
mpole__.maxpole = *maxpole;
mpole__.npole = npole;
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

mrecip_st mrecip__;

void set_mrecip_data_ (double* vmxx, double* vmyy, double* vmzz, double* vmxy, double* vmxz, double* vmyz, double* cmp, double* fmp, double* cphi, double* fphi) {
mrecip__.vmxx = vmxx;
mrecip__.vmyy = vmyy;
mrecip__.vmzz = vmzz;
mrecip__.vmxy = vmxy;
mrecip__.vmxz = vmxz;
mrecip__.vmyz = vmyz;
mrecip__.cmp = cmp;
mrecip__.fmp = fmp;
mrecip__.cphi = cphi;
mrecip__.fphi = fphi;
}

mutant_st mutant__;

void set_mutant_data_ (int* nmut, int* vcouple, int* imut, int* type0, int* class0, int* type1, int* class1, double* lambda, double* tlambda, double* vlambda, double* elambda, double* scexp, double* scalpha, int* mut) {
mutant__.nmut = nmut;
mutant__.vcouple = vcouple;
mutant__.imut = imut;
mutant__.type0 = type0;
mutant__.class0 = class0;
mutant__.type1 = type1;
mutant__.class1 = class1;
mutant__.lambda = lambda;
mutant__.tlambda = tlambda;
mutant__.vlambda = vlambda;
mutant__.elambda = elambda;
mutant__.scexp = scexp;
mutant__.scalpha = scalpha;
mutant__.mut = mut;
}

neigh_st neigh__;

void set_neigh_data_ (int* maxvlst, int* maxelst, int* maxulst, int* nvlst, int* vlst, int* nelst, int* elst, int* nulst, int* ulst, double* lbuffer, double* pbuffer, double* lbuf2, double* pbuf2, double* vbuf2, double* vbufx, double* dbuf2, double* dbufx, double* cbuf2, double* cbufx, double* mbuf2, double* mbufx, double* ubuf2, double* ubufx, double* xvold, double* yvold, double* zvold, double* xeold, double* yeold, double* zeold, double* xuold, double* yuold, double* zuold, int* dovlst, int* dodlst, int* doclst, int* domlst, int* doulst) {
neigh__.maxvlst = maxvlst;
neigh__.maxelst = maxelst;
neigh__.maxulst = maxulst;
neigh__.nvlst = nvlst;
neigh__.vlst = vlst;
neigh__.nelst = nelst;
neigh__.elst = elst;
neigh__.nulst = nulst;
neigh__.ulst = ulst;
neigh__.lbuffer = lbuffer;
neigh__.pbuffer = pbuffer;
neigh__.lbuf2 = lbuf2;
neigh__.pbuf2 = pbuf2;
neigh__.vbuf2 = vbuf2;
neigh__.vbufx = vbufx;
neigh__.dbuf2 = dbuf2;
neigh__.dbufx = dbufx;
neigh__.cbuf2 = cbuf2;
neigh__.cbufx = cbufx;
neigh__.mbuf2 = mbuf2;
neigh__.mbufx = mbufx;
neigh__.ubuf2 = ubuf2;
neigh__.ubufx = ubufx;
neigh__.xvold = xvold;
neigh__.yvold = yvold;
neigh__.zvold = zvold;
neigh__.xeold = xeold;
neigh__.yeold = yeold;
neigh__.zeold = zeold;
neigh__.xuold = xuold;
neigh__.yuold = yuold;
neigh__.zuold = zuold;
neigh__.dovlst = dovlst;
neigh__.dodlst = dodlst;
neigh__.doclst = doclst;
neigh__.domlst = domlst;
neigh__.doulst = doulst;
}

nonpol_st nonpol__;

void set_nonpol_data_ (double* epso, double* epsh, double* rmino, double* rminh, double* awater, double* slevy, double* solvprs, double* surften, double* spcut, double* spoff, double* stcut, double* stoff, double* rcav, double* rdisp, double* cdisp) {
nonpol__.epso = *epso;
nonpol__.epsh = *epsh;
nonpol__.rmino = *rmino;
nonpol__.rminh = *rminh;
nonpol__.awater = *awater;
nonpol__.slevy = *slevy;
nonpol__.solvprs = solvprs;
nonpol__.surften = surften;
nonpol__.spcut = spcut;
nonpol__.spoff = spoff;
nonpol__.stcut = stcut;
nonpol__.stoff = stoff;
nonpol__.rcav = rcav;
nonpol__.rdisp = rdisp;
nonpol__.cdisp = cdisp;
}

nucleo_st nucleo__;

void set_nucleo_data_ (int* pucker, double* glyco, double* bkbone, int* dblhlx, int* deoxy, char* hlxform) {
nucleo__.pucker = pucker;
nucleo__.glyco = glyco;
nucleo__.bkbone = bkbone;
nucleo__.dblhlx = dblhlx;
nucleo__.deoxy = deoxy;
nucleo__.hlxform = hlxform;
}

omega_st omega__;

void set_omega_data_ (int* nomega, int* iomega, int* zline, double* dihed) {
omega__.nomega = nomega;
omega__.iomega = iomega;
omega__.zline = zline;
omega__.dihed = dihed;
}

opbend_st opbend__;

void set_opbend_data_ (int* nopbend, int* iopb, double* opbk) {
opbend__.nopbend = nopbend;
opbend__.iopb = iopb;
opbend__.opbk = opbk;
}

opdist_st opdist__;

void set_opdist_data_ (int* nopdist, int* iopd, double* opdk) {
opdist__.nopdist = nopdist;
opdist__.iopd = iopd;
opdist__.opdk = opdk;
}

openmm_st openmm__;

void set_openmm_data_ (unsigned long long* ommhandle, char* cudaprecision, char* ommplatform, char* cudadevice) {
openmm__.ommhandle = ommhandle;
openmm__.cudaprecision = cudaprecision;
openmm__.ommplatform = ommplatform;
openmm__.cudadevice = cudadevice;
}

openmp_st openmp__;

void set_openmp_data_ (int* nproc, int* nthread) {
openmp__.nproc = nproc;
openmp__.nthread = nthread;
}

orbits_st orbits__;

void set_orbits_data_ (double* qorb, double* worb, double* emorb) {
orbits__.qorb = qorb;
orbits__.worb = worb;
orbits__.emorb = emorb;
}

output_st output__;

void set_output_data_ (int* archive, int* noversion, int* overwrite, int* cyclesave, int* velsave, int* frcsave, int* uindsave, char* coordtype) {
output__.archive = archive;
output__.noversion = noversion;
output__.overwrite = overwrite;
output__.cyclesave = cyclesave;
output__.velsave = velsave;
output__.frcsave = frcsave;
output__.uindsave = uindsave;
output__.coordtype = coordtype;
}

params_st params__;

void set_params_data_ (int* maxprm, int* nprm, char* prmline) {
params__.maxprm = *maxprm;
params__.nprm = nprm;
params__.prmline = prmline;
}

paths_st paths__;

void set_paths_data_ (double* pnorm, double* acoeff, double* pc0, double* pc1, double* pvect, double* pstep, double* pzet, double* gc) {
paths__.pnorm = pnorm;
paths__.acoeff = acoeff;
paths__.pc0 = pc0;
paths__.pc1 = pc1;
paths__.pvect = pvect;
paths__.pstep = pstep;
paths__.pzet = pzet;
paths__.gc = gc;
}

pbstuf_st pbstuf__;

void set_pbstuf_data_ (int* maxion, int* ionn, int* dime, int* ionq, double* pbe, double* pdie, double* sdie, double* srad, double* swin, double* sdens, double* smin, double* grid, double* gcent, double* cgrid, double* cgcent, double* fgrid, double* fgcent, double* ionr, double* ionc, double* apbe, double* pbr, double* pbep, double* pbfp, double* pbtp, double* pbeuind, double* pbeuinp, char* pbtyp, char* pbsoln, char* bcfl, char* chgm, char* srfm) {
pbstuf__.maxion = *maxion;
pbstuf__.ionn = ionn;
pbstuf__.dime = dime;
pbstuf__.ionq = ionq;
pbstuf__.pbe = pbe;
pbstuf__.pdie = pdie;
pbstuf__.sdie = sdie;
pbstuf__.srad = srad;
pbstuf__.swin = swin;
pbstuf__.sdens = sdens;
pbstuf__.smin = smin;
pbstuf__.grid = grid;
pbstuf__.gcent = gcent;
pbstuf__.cgrid = cgrid;
pbstuf__.cgcent = cgcent;
pbstuf__.fgrid = fgrid;
pbstuf__.fgcent = fgcent;
pbstuf__.ionr = ionr;
pbstuf__.ionc = ionc;
pbstuf__.apbe = apbe;
pbstuf__.pbr = pbr;
pbstuf__.pbep = pbep;
pbstuf__.pbfp = pbfp;
pbstuf__.pbtp = pbtp;
pbstuf__.pbeuind = pbeuind;
pbstuf__.pbeuinp = pbeuinp;
pbstuf__.pbtyp = pbtyp;
pbstuf__.pbsoln = pbsoln;
pbstuf__.bcfl = bcfl;
pbstuf__.chgm = chgm;
pbstuf__.srfm = srfm;
}

pdb_st pdb__;

void set_pdb_data_ (int* npdb, int* nres, int* resnum, int* resatm, int* npdb12, int* ipdb12, int* pdblist, double* xpdb, double* ypdb, double* zpdb, char* altsym, char* pdbres, char* pdbatm, char* pdbtyp, char* chnsym, char* instyp) {
pdb__.npdb = npdb;
pdb__.nres = nres;
pdb__.resnum = resnum;
pdb__.resatm = resatm;
pdb__.npdb12 = npdb12;
pdb__.ipdb12 = ipdb12;
pdb__.pdblist = pdblist;
pdb__.xpdb = xpdb;
pdb__.ypdb = ypdb;
pdb__.zpdb = zpdb;
pdb__.altsym = altsym;
pdb__.pdbres = pdbres;
pdb__.pdbatm = pdbatm;
pdb__.pdbtyp = pdbtyp;
pdb__.chnsym = chnsym;
pdb__.instyp = instyp;
}

phipsi_st phipsi__;

void set_phipsi_data_ (int* chiral, int* disulf, double* phi, double* psi, double* omg, double* chi) {
phipsi__.chiral = chiral;
phipsi__.disulf = disulf;
phipsi__.phi = phi;
phipsi__.psi = psi;
phipsi__.omg = omg;
phipsi__.chi = chi;
}

piorbs_st piorbs__;

void set_piorbs_data_ (int* norbit, int* nconj, int* reorbit, int* nbpi, int* ntpi, int* iorbit, int* iconj, int* kconj, int* piperp, int* ibpi, int* itpi, double* pbpl, double* pnpl, int* listpi) {
piorbs__.norbit = norbit;
piorbs__.nconj = nconj;
piorbs__.reorbit = reorbit;
piorbs__.nbpi = nbpi;
piorbs__.ntpi = ntpi;
piorbs__.iorbit = iorbit;
piorbs__.iconj = iconj;
piorbs__.kconj = kconj;
piorbs__.piperp = piperp;
piorbs__.ibpi = ibpi;
piorbs__.itpi = itpi;
piorbs__.pbpl = pbpl;
piorbs__.pnpl = pnpl;
piorbs__.listpi = listpi;
}

pistuf_st pistuf__;

void set_pistuf_data_ (double* bkpi, double* blpi, double* kslope, double* lslope, double* torsp2) {
pistuf__.bkpi = bkpi;
pistuf__.blpi = blpi;
pistuf__.kslope = kslope;
pistuf__.lslope = lslope;
pistuf__.torsp2 = torsp2;
}

pitors_st pitors__;

void set_pitors_data_ (int* npitors, int* ipit, double* kpit) {
pitors__.npitors = npitors;
pitors__.ipit = ipit;
pitors__.kpit = kpit;
}

pme_st pme__;

void set_pme_data_ (int* nfft1, int* nfft2, int* nfft3, int* nefft1, int* nefft2, int* nefft3, int* ndfft1, int* ndfft2, int* ndfft3, int* bsorder, int* bseorder, int* bsporder, int* bsdorder, int* igrid, double* bsmod1, double* bsmod2, double* bsmod3, double* bsbuild, double* thetai1, double* thetai2, double* thetai3, double* qgrid, double* qfac) {
pme__.nfft1 = nfft1;
pme__.nfft2 = nfft2;
pme__.nfft3 = nfft3;
pme__.nefft1 = nefft1;
pme__.nefft2 = nefft2;
pme__.nefft3 = nefft3;
pme__.ndfft1 = ndfft1;
pme__.ndfft2 = ndfft2;
pme__.ndfft3 = ndfft3;
pme__.bsorder = bsorder;
pme__.bseorder = bseorder;
pme__.bsporder = bsporder;
pme__.bsdorder = bsdorder;
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

polar_st polar__;

void set_polar_data_ (int* npolar, int* ipolar, double* polarity, double* thole, double* dirdamp, double* pdamp, double* udir, double* udirp, double* udirs, double* udirps, double* uind, double* uinp, double* uinds, double* uinps, double* uexact, int* douind) {
polar__.npolar = npolar;
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

polgrp_st polgrp__;

void set_polgrp_data_ (int* maxp11, int* maxp12, int* maxp13, int* maxp14, int* np11, int* np12, int* np13, int* np14, int* ip11, int* ip12, int* ip13, int* ip14) {
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

polopt_st polopt__;

void set_polopt_data_ (int* maxopt, int* optorder, int* optlevel, double* copt, double* copm, double* uopt, double* uoptp, double* uopts, double* uoptps, double* fopt, double* foptp) {
polopt__.maxopt = *maxopt;
polopt__.optorder = optorder;
polopt__.optlevel = optlevel;
polopt__.copt = copt;
polopt__.copm = copm;
polopt__.uopt = uopt;
polopt__.uoptp = uoptp;
polopt__.uopts = uopts;
polopt__.uoptps = uoptps;
polopt__.fopt = fopt;
polopt__.foptp = foptp;
}

polpcg_st polpcg__;

void set_polpcg_data_ (int* mindex, double* pcgpeek, double* minv, int* pcgprec, int* pcgguess) {
polpcg__.mindex = mindex;
polpcg__.pcgpeek = pcgpeek;
polpcg__.minv = minv;
polpcg__.pcgprec = pcgprec;
polpcg__.pcgguess = pcgguess;
}

polpot_st polpot__;

void set_polpot_data_ (int* politer, double* poleps, double* p2scale, double* p3scale, double* p4scale, double* p5scale, double* p2iscale, double* p3iscale, double* p4iscale, double* p5iscale, double* d1scale, double* d2scale, double* d3scale, double* d4scale, double* u1scale, double* u2scale, double* u3scale, double* u4scale, double* w2scale, double* w3scale, double* w4scale, double* w5scale, double* udiag, int* dpequal, int* use_thole, int* use_dirdamp, char* poltyp) {
polpot__.politer = politer;
polpot__.poleps = poleps;
polpot__.p2scale = p2scale;
polpot__.p3scale = p3scale;
polpot__.p4scale = p4scale;
polpot__.p5scale = p5scale;
polpot__.p2iscale = p2iscale;
polpot__.p3iscale = p3iscale;
polpot__.p4iscale = p4iscale;
polpot__.p5iscale = p5iscale;
polpot__.d1scale = d1scale;
polpot__.d2scale = d2scale;
polpot__.d3scale = d3scale;
polpot__.d4scale = d4scale;
polpot__.u1scale = u1scale;
polpot__.u2scale = u2scale;
polpot__.u3scale = u3scale;
polpot__.u4scale = u4scale;
polpot__.w2scale = w2scale;
polpot__.w3scale = w3scale;
polpot__.w4scale = w4scale;
polpot__.w5scale = w5scale;
polpot__.udiag = udiag;
polpot__.dpequal = dpequal;
polpot__.use_thole = use_thole;
polpot__.use_dirdamp = use_dirdamp;
polpot__.poltyp = poltyp;
}

poltcg_st poltcg__;

void set_poltcg_data_ (int* tcgorder, int* tcgnab, double* tcgpeek, double* uad, double* uap, double* ubd, double* ubp, int* tcgguess) {
poltcg__.tcgorder = tcgorder;
poltcg__.tcgnab = tcgnab;
poltcg__.tcgpeek = tcgpeek;
poltcg__.uad = uad;
poltcg__.uap = uap;
poltcg__.ubd = ubd;
poltcg__.ubp = ubp;
poltcg__.tcgguess = tcgguess;
}

potent_st potent__;

void set_potent_data_ (int* use_bond, int* use_angle, int* use_strbnd, int* use_urey, int* use_angang, int* use_opbend, int* use_opdist, int* use_improp, int* use_imptor, int* use_tors, int* use_pitors, int* use_strtor, int* use_angtor, int* use_tortor, int* use_vdw, int* use_repuls, int* use_disp, int* use_charge, int* use_chgdpl, int* use_dipole, int* use_mpole, int* use_polar, int* use_chgtrn, int* use_rxnfld, int* use_solv, int* use_metal, int* use_geom, int* use_extra, int* use_born, int* use_orbit) {
potent__.use_bond = use_bond;
potent__.use_angle = use_angle;
potent__.use_strbnd = use_strbnd;
potent__.use_urey = use_urey;
potent__.use_angang = use_angang;
potent__.use_opbend = use_opbend;
potent__.use_opdist = use_opdist;
potent__.use_improp = use_improp;
potent__.use_imptor = use_imptor;
potent__.use_tors = use_tors;
potent__.use_pitors = use_pitors;
potent__.use_strtor = use_strtor;
potent__.use_angtor = use_angtor;
potent__.use_tortor = use_tortor;
potent__.use_vdw = use_vdw;
potent__.use_repuls = use_repuls;
potent__.use_disp = use_disp;
potent__.use_charge = use_charge;
potent__.use_chgdpl = use_chgdpl;
potent__.use_dipole = use_dipole;
potent__.use_mpole = use_mpole;
potent__.use_polar = use_polar;
potent__.use_chgtrn = use_chgtrn;
potent__.use_rxnfld = use_rxnfld;
potent__.use_solv = use_solv;
potent__.use_metal = use_metal;
potent__.use_geom = use_geom;
potent__.use_extra = use_extra;
potent__.use_born = use_born;
potent__.use_orbit = use_orbit;
}

potfit_st potfit__;

void set_potfit_data_ (int* nconf, int* namax, int* ngatm, int* nfatm, int* npgrid, int* ipgrid, double* resp, double* xdpl0, double* ydpl0, double* zdpl0, double* xxqdp0, double* xyqdp0, double* xzqdp0, double* yyqdp0, double* yzqdp0, double* zzqdp0, double* fit0, double* fchg, double* fpol, double* pgrid, double* epot, int* use_dpl, int* use_qdp, int* fit_mpl, int* fit_dpl, int* fit_qdp, int* fitchg, int* fitpol, int* gatm, int* fatm) {
potfit__.nconf = nconf;
potfit__.namax = namax;
potfit__.ngatm = ngatm;
potfit__.nfatm = nfatm;
potfit__.npgrid = npgrid;
potfit__.ipgrid = ipgrid;
potfit__.resp = resp;
potfit__.xdpl0 = xdpl0;
potfit__.ydpl0 = ydpl0;
potfit__.zdpl0 = zdpl0;
potfit__.xxqdp0 = xxqdp0;
potfit__.xyqdp0 = xyqdp0;
potfit__.xzqdp0 = xzqdp0;
potfit__.yyqdp0 = yyqdp0;
potfit__.yzqdp0 = yzqdp0;
potfit__.zzqdp0 = zzqdp0;
potfit__.fit0 = fit0;
potfit__.fchg = fchg;
potfit__.fpol = fpol;
potfit__.pgrid = pgrid;
potfit__.epot = epot;
potfit__.use_dpl = use_dpl;
potfit__.use_qdp = use_qdp;
potfit__.fit_mpl = fit_mpl;
potfit__.fit_dpl = fit_dpl;
potfit__.fit_qdp = fit_qdp;
potfit__.fitchg = fitchg;
potfit__.fitpol = fitpol;
potfit__.gatm = gatm;
potfit__.fatm = fatm;
}

ptable_st ptable__;

void set_ptable_data_ (int* maxele, double* atmass, double* vdwrad, double* covrad, char* elemnt) {
ptable__.maxele = *maxele;
ptable__.atmass = atmass;
ptable__.vdwrad = vdwrad;
ptable__.covrad = covrad;
ptable__.elemnt = elemnt;
}

qmstuf_st qmstuf__;

void set_qmstuf_data_ (int* ngatom, double* egau, double* gx, double* gy, double* gz, double* gfreq, double* gforce, double* gh) {
qmstuf__.ngatom = ngatom;
qmstuf__.egau = egau;
qmstuf__.gx = gx;
qmstuf__.gy = gy;
qmstuf__.gz = gz;
qmstuf__.gfreq = gfreq;
qmstuf__.gforce = gforce;
qmstuf__.gh = gh;
}

refer_st refer__;

void set_refer_data_ (int* nref, int* refltitle, int* refleng, int* reftyp, int* n12ref, int* i12ref, double* xboxref, double* yboxref, double* zboxref, double* alpharef, double* betaref, double* gammaref, double* xref, double* yref, double* zref, char* refnam, char* reffile, char* reftitle) {
refer__.nref = nref;
refer__.refltitle = refltitle;
refer__.refleng = refleng;
refer__.reftyp = reftyp;
refer__.n12ref = n12ref;
refer__.i12ref = i12ref;
refer__.xboxref = xboxref;
refer__.yboxref = yboxref;
refer__.zboxref = zboxref;
refer__.alpharef = alpharef;
refer__.betaref = betaref;
refer__.gammaref = gammaref;
refer__.xref = xref;
refer__.yref = yref;
refer__.zref = zref;
refer__.refnam = refnam;
refer__.reffile = reffile;
refer__.reftitle = reftitle;
}

repel_st repel__;

void set_repel_data_ (int* nrep, double* sizpr, double* dmppr, double* elepr) {
repel__.nrep = nrep;
repel__.sizpr = sizpr;
repel__.dmppr = dmppr;
repel__.elepr = elepr;
}

reppot_st reppot__;

void set_reppot_data_ (double* r2scale, double* r3scale, double* r4scale, double* r5scale) {
reppot__.r2scale = r2scale;
reppot__.r3scale = r3scale;
reppot__.r4scale = r4scale;
reppot__.r5scale = r5scale;
}

resdue_st resdue__;

void set_resdue_data_ (int* maxamino, int* maxnuc, int* ntyp, int* catyp, int* ctyp, int* hntyp, int* otyp, int* hatyp, int* cbtyp, int* nntyp, int* cantyp, int* cntyp, int* hnntyp, int* ontyp, int* hantyp, int* nctyp, int* cactyp, int* cctyp, int* hnctyp, int* octyp, int* hactyp, int* o5typ, int* c5typ, int* h51typ, int* h52typ, int* c4typ, int* h4typ, int* o4typ, int* c1typ, int* h1typ, int* c3typ, int* h3typ, int* c2typ, int* h21typ, int* o2typ, int* h22typ, int* o3typ, int* ptyp, int* optyp, int* h5ttyp, int* h3ttyp, char* amino1, char* nuclz1, char* amino, char* nuclz) {
resdue__.maxamino = *maxamino;
resdue__.maxnuc = *maxnuc;
resdue__.ntyp = ntyp;
resdue__.catyp = catyp;
resdue__.ctyp = ctyp;
resdue__.hntyp = hntyp;
resdue__.otyp = otyp;
resdue__.hatyp = hatyp;
resdue__.cbtyp = cbtyp;
resdue__.nntyp = nntyp;
resdue__.cantyp = cantyp;
resdue__.cntyp = cntyp;
resdue__.hnntyp = hnntyp;
resdue__.ontyp = ontyp;
resdue__.hantyp = hantyp;
resdue__.nctyp = nctyp;
resdue__.cactyp = cactyp;
resdue__.cctyp = cctyp;
resdue__.hnctyp = hnctyp;
resdue__.octyp = octyp;
resdue__.hactyp = hactyp;
resdue__.o5typ = o5typ;
resdue__.c5typ = c5typ;
resdue__.h51typ = h51typ;
resdue__.h52typ = h52typ;
resdue__.c4typ = c4typ;
resdue__.h4typ = h4typ;
resdue__.o4typ = o4typ;
resdue__.c1typ = c1typ;
resdue__.h1typ = h1typ;
resdue__.c3typ = c3typ;
resdue__.h3typ = h3typ;
resdue__.c2typ = c2typ;
resdue__.h21typ = h21typ;
resdue__.o2typ = o2typ;
resdue__.h22typ = h22typ;
resdue__.o3typ = o3typ;
resdue__.ptyp = ptyp;
resdue__.optyp = optyp;
resdue__.h5ttyp = h5ttyp;
resdue__.h3ttyp = h3ttyp;
resdue__.amino1 = amino1;
resdue__.nuclz1 = nuclz1;
resdue__.amino = amino;
resdue__.nuclz = nuclz;
}

restrn_st restrn__;

void set_restrn_data_ (int* npfix, int* ndfix, int* nafix, int* ntfix, int* ngfix, int* nchir, int* ipfix, int* kpfix, int* idfix, int* iafix, int* itfix, int* igfix, int* ichir, double* depth, double* width, double* rwall, double* xpfix, double* ypfix, double* zpfix, double* pfix, double* dfix, double* afix, double* tfix, double* gfix, double* chir, int* use_basin, int* use_wall) {
restrn__.npfix = npfix;
restrn__.ndfix = ndfix;
restrn__.nafix = nafix;
restrn__.ntfix = ntfix;
restrn__.ngfix = ngfix;
restrn__.nchir = nchir;
restrn__.ipfix = ipfix;
restrn__.kpfix = kpfix;
restrn__.idfix = idfix;
restrn__.iafix = iafix;
restrn__.itfix = itfix;
restrn__.igfix = igfix;
restrn__.ichir = ichir;
restrn__.depth = depth;
restrn__.width = width;
restrn__.rwall = rwall;
restrn__.xpfix = xpfix;
restrn__.ypfix = ypfix;
restrn__.zpfix = zpfix;
restrn__.pfix = pfix;
restrn__.dfix = dfix;
restrn__.afix = afix;
restrn__.tfix = tfix;
restrn__.gfix = gfix;
restrn__.chir = chir;
restrn__.use_basin = use_basin;
restrn__.use_wall = use_wall;
}

rgddyn_st rgddyn__;

void set_rgddyn_data_ (double* xcmo, double* ycmo, double* zcmo, double* vcm, double* wcm, double* lm, double* vc, double* wc, int* linear) {
rgddyn__.xcmo = xcmo;
rgddyn__.ycmo = ycmo;
rgddyn__.zcmo = zcmo;
rgddyn__.vcm = vcm;
rgddyn__.wcm = wcm;
rgddyn__.lm = lm;
rgddyn__.vc = vc;
rgddyn__.wc = wc;
rgddyn__.linear = linear;
}

rigid_st rigid__;

void set_rigid_data_ (double* xrb, double* yrb, double* zrb, double* rbc, int* use_rigid) {
rigid__.xrb = xrb;
rigid__.yrb = yrb;
rigid__.zrb = zrb;
rigid__.rbc = rbc;
rigid__.use_rigid = use_rigid;
}

ring_st ring__;

void set_ring_data_ (int* nring3, int* nring4, int* nring5, int* nring6, int* nring7, int* iring3, int* iring4, int* iring5, int* iring6, int* iring7) {
ring__.nring3 = nring3;
ring__.nring4 = nring4;
ring__.nring5 = nring5;
ring__.nring6 = nring6;
ring__.nring7 = nring7;
ring__.iring3 = iring3;
ring__.iring4 = iring4;
ring__.iring5 = iring5;
ring__.iring6 = iring6;
ring__.iring7 = iring7;
}

rotbnd_st rotbnd__;

void set_rotbnd_data_ (int* nrot, int* rot, int* use_short) {
rotbnd__.nrot = nrot;
rotbnd__.rot = rot;
rotbnd__.use_short = use_short;
}

rxnfld_st rxnfld__;

void set_rxnfld_data_ (int* ijk, double* b1, double* b2) {
rxnfld__.ijk = ijk;
rxnfld__.b1 = b1;
rxnfld__.b2 = b2;
}

rxnpot_st rxnpot__;

void set_rxnpot_data_ (int* rfterms, double* rfsize, double* rfbulkd) {
rxnpot__.rfterms = rfterms;
rxnpot__.rfsize = rfsize;
rxnpot__.rfbulkd = rfbulkd;
}

scales_st scales__;

void set_scales_data_ (double* scale, int* set_scale) {
scales__.scale = scale;
scales__.set_scale = set_scale;
}

sequen_st sequen__;

void set_sequen_data_ (int* nseq, int* nchain, int* ichain, int* seqtyp, char* chnnam, char* seq, char* chntyp) {
sequen__.nseq = nseq;
sequen__.nchain = nchain;
sequen__.ichain = ichain;
sequen__.seqtyp = seqtyp;
sequen__.chnnam = chnnam;
sequen__.seq = seq;
sequen__.chntyp = chntyp;
}

shunt_st shunt__;

void set_shunt_data_ (double* off, double* off2, double* cut, double* cut2, double* c0, double* c1, double* c2, double* c3, double* c4, double* c5, double* f0, double* f1, double* f2, double* f3, double* f4, double* f5, double* f6, double* f7) {
shunt__.off = off;
shunt__.off2 = off2;
shunt__.cut = cut;
shunt__.cut2 = cut2;
shunt__.c0 = c0;
shunt__.c1 = c1;
shunt__.c2 = c2;
shunt__.c3 = c3;
shunt__.c4 = c4;
shunt__.c5 = c5;
shunt__.f0 = f0;
shunt__.f1 = f1;
shunt__.f2 = f2;
shunt__.f3 = f3;
shunt__.f4 = f4;
shunt__.f5 = f5;
shunt__.f6 = f6;
shunt__.f7 = f7;
}

sizes_st sizes__;

void set_sizes_data_ (int* maxatm, int* maxtyp, int* maxclass, int* maxval, int* maxref, int* maxgrp, int* maxres, int* maxfix) {
sizes__.maxatm = *maxatm;
sizes__.maxtyp = *maxtyp;
sizes__.maxclass = *maxclass;
sizes__.maxval = *maxval;
sizes__.maxref = *maxref;
sizes__.maxgrp = *maxgrp;
sizes__.maxres = *maxres;
sizes__.maxfix = *maxfix;
}

socket_st socket__;

void set_socket_data_ (int* skttyp, int* cstep, double* cdt, double* cenergy, int* sktstart, int* sktstop, int* use_socket) {
socket__.skttyp = skttyp;
socket__.cstep = cstep;
socket__.cdt = cdt;
socket__.cenergy = cenergy;
socket__.sktstart = sktstart;
socket__.sktstop = sktstop;
socket__.use_socket = use_socket;
}

solute_st solute__;

void set_solute_data_ (double* doffset, double* p1, double* p2, double* p3, double* p4, double* p5, double* rsolv, double* asolv, double* rborn, double* drb, double* drbp, double* drobc, double* gpol, double* shct, double* aobc, double* bobc, double* gobc, double* vsolv, double* wace, double* s2ace, double* uace, char* solvtyp, char* borntyp) {
solute__.doffset = doffset;
solute__.p1 = p1;
solute__.p2 = p2;
solute__.p3 = p3;
solute__.p4 = p4;
solute__.p5 = p5;
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
solute__.solvtyp = solvtyp;
solute__.borntyp = borntyp;
}

stodyn_st stodyn__;

void set_stodyn_data_ (double* friction, double* fgamma, int* use_sdarea) {
stodyn__.friction = friction;
stodyn__.fgamma = fgamma;
stodyn__.use_sdarea = use_sdarea;
}

strbnd_st strbnd__;

void set_strbnd_data_ (int* nstrbnd, int* isb, double* sbk) {
strbnd__.nstrbnd = nstrbnd;
strbnd__.isb = isb;
strbnd__.sbk = sbk;
}

strtor_st strtor__;

void set_strtor_data_ (int* nstrtor, int* ist, double* kst) {
strtor__.nstrtor = nstrtor;
strtor__.ist = ist;
strtor__.kst = kst;
}

syntrn_st syntrn__;

void set_syntrn_data_ (double* tpath, double* ppath, double* xmin1, double* xmin2, double* xm) {
syntrn__.tpath = tpath;
syntrn__.ppath = ppath;
syntrn__.xmin1 = xmin1;
syntrn__.xmin2 = xmin2;
syntrn__.xm = xm;
}

tarray_st tarray__;

void set_tarray_data_ (int* ntpair, int* tindex, double* tdipdip) {
tarray__.ntpair = ntpair;
tarray__.tindex = tindex;
tarray__.tdipdip = tdipdip;
}

titles_st titles__;

void set_titles_data_ (int* ltitle, char* title) {
titles__.ltitle = ltitle;
titles__.title = title;
}

torpot_st torpot__;

void set_torpot_data_ (double* idihunit, double* itorunit, double* torsunit, double* ptorunit, double* storunit, double* atorunit, double* ttorunit) {
torpot__.idihunit = idihunit;
torpot__.itorunit = itorunit;
torpot__.torsunit = torsunit;
torpot__.ptorunit = ptorunit;
torpot__.storunit = storunit;
torpot__.atorunit = atorunit;
torpot__.ttorunit = ttorunit;
}

tors_st tors__;

void set_tors_data_ (int* ntors, int* itors, double* tors1, double* tors2, double* tors3, double* tors4, double* tors5, double* tors6) {
tors__.ntors = ntors;
tors__.itors = itors;
tors__.tors1 = tors1;
tors__.tors2 = tors2;
tors__.tors3 = tors3;
tors__.tors4 = tors4;
tors__.tors5 = tors5;
tors__.tors6 = tors6;
}

tortor_st tortor__;

void set_tortor_data_ (int* ntortor, int* itt) {
tortor__.ntortor = ntortor;
tortor__.itt = itt;
}

tree_st tree__;

void set_tree_data_ (int* maxpss, int* nlevel, double* etree, double* ilevel) {
tree__.maxpss = *maxpss;
tree__.nlevel = nlevel;
tree__.etree = etree;
tree__.ilevel = ilevel;
}

units_st units__;

void set_units_data_ (double* avogadro, double* lightspd, double* boltzmann, double* gasconst, double* elemchg, double* vacperm, double* emass, double* planck, double* joule, double* ekcal, double* bohr, double* hartree, double* evolt, double* efreq, double* coulomb, double* debye, double* prescon) {
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

uprior_st uprior__;

void set_uprior_data_ (int* maxpred, int* nualt, int* maxualt, double* gear, double* aspc, double* bpred, double* bpredp, double* bpreds, double* bpredps, double* udalt, double* upalt, double* usalt, double* upsalt, int* use_pred, char* polpred) {
uprior__.maxpred = *maxpred;
uprior__.nualt = nualt;
uprior__.maxualt = maxualt;
uprior__.gear = gear;
uprior__.aspc = aspc;
uprior__.bpred = bpred;
uprior__.bpredp = bpredp;
uprior__.bpreds = bpreds;
uprior__.bpredps = bpredps;
uprior__.udalt = udalt;
uprior__.upalt = upalt;
uprior__.usalt = usalt;
uprior__.upsalt = upsalt;
uprior__.use_pred = use_pred;
uprior__.polpred = polpred;
}

urey_st urey__;

void set_urey_data_ (int* nurey, int* iury, double* uk, double* ul) {
urey__.nurey = nurey;
urey__.iury = iury;
urey__.uk = uk;
urey__.ul = ul;
}

urypot_st urypot__;

void set_urypot_data_ (double* cury, double* qury, double* ureyunit) {
urypot__.cury = cury;
urypot__.qury = qury;
urypot__.ureyunit = ureyunit;
}

usage_st usage__;

void set_usage_data_ (int* nuse, int* iuse, int* use) {
usage__.nuse = nuse;
usage__.iuse = iuse;
usage__.use = use;
}

valfit_st valfit__;

void set_valfit_data_ (int* fit_bond, int* fit_angle, int* fit_strbnd, int* fit_urey, int* fit_opbend, int* fit_tors, int* fit_struct, int* fit_force) {
valfit__.fit_bond = fit_bond;
valfit__.fit_angle = fit_angle;
valfit__.fit_strbnd = fit_strbnd;
valfit__.fit_urey = fit_urey;
valfit__.fit_opbend = fit_opbend;
valfit__.fit_tors = fit_tors;
valfit__.fit_struct = fit_struct;
valfit__.fit_force = fit_force;
}

vdw_st vdw__;

void set_vdw_data_ (int* nvdw, int* ivdw, int* jvdw, int* ired, double* kred, double* xred, double* yred, double* zred, double* radmin, double* epsilon, double* radmin4, double* epsilon4, double* radhbnd, double* epshbnd) {
vdw__.nvdw = nvdw;
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

vdwpot_st vdwpot__;

void set_vdwpot_data_ (int* maxgauss, int* ngauss, double* igauss, double* abuck, double* bbuck, double* cbuck, double* ghal, double* dhal, double* v2scale, double* v3scale, double* v4scale, double* v5scale, int* use_vcorr, char* vdwindex, char* radtyp, char* radsiz, char* gausstyp, char* radrule, char* epsrule, char* vdwtyp) {
vdwpot__.maxgauss = *maxgauss;
vdwpot__.ngauss = ngauss;
vdwpot__.igauss = igauss;
vdwpot__.abuck = abuck;
vdwpot__.bbuck = bbuck;
vdwpot__.cbuck = cbuck;
vdwpot__.ghal = ghal;
vdwpot__.dhal = dhal;
vdwpot__.v2scale = v2scale;
vdwpot__.v3scale = v3scale;
vdwpot__.v4scale = v4scale;
vdwpot__.v5scale = v5scale;
vdwpot__.use_vcorr = use_vcorr;
vdwpot__.vdwindex = vdwindex;
vdwpot__.radtyp = radtyp;
vdwpot__.radsiz = radsiz;
vdwpot__.gausstyp = gausstyp;
vdwpot__.radrule = radrule;
vdwpot__.epsrule = epsrule;
vdwpot__.vdwtyp = vdwtyp;
}

vibs_st vibs__;

void set_vibs_data_ (double* rho, double* rhok, double* rwork) {
vibs__.rho = rho;
vibs__.rhok = rhok;
vibs__.rwork = rwork;
}

virial_st virial__;

void set_virial_data_ (double* vir, int* use_virial) {
virial__.vir = vir;
virial__.use_virial = use_virial;
}

warp_st warp__;

void set_warp_data_ (double* deform, double* difft, double* diffv, double* diffc, double* m2, int* use_smooth, int* use_dem, int* use_gda, int* use_tophat, int* use_stophat) {
warp__.deform = deform;
warp__.difft = difft;
warp__.diffv = diffv;
warp__.diffc = diffc;
warp__.m2 = m2;
warp__.use_smooth = use_smooth;
warp__.use_dem = use_dem;
warp__.use_gda = use_gda;
warp__.use_tophat = use_tophat;
warp__.use_stophat = use_stophat;
}

xtals_st xtals__;

void set_xtals_data_ (int* maxlsq, int* maxrsd, int* nxtal, int* nvary, int* ivary, int* iresid, int* vary, double* e0_lattice, char* vartyp, char* rsdtyp) {
xtals__.maxlsq = *maxlsq;
xtals__.maxrsd = *maxrsd;
xtals__.nxtal = nxtal;
xtals__.nvary = nvary;
xtals__.ivary = ivary;
xtals__.iresid = iresid;
xtals__.vary = vary;
xtals__.e0_lattice = e0_lattice;
xtals__.vartyp = vartyp;
xtals__.rsdtyp = rsdtyp;
}

zclose_st zclose__;

void set_zclose_data_ (int* nadd, int* ndel, int* iadd, int* idel) {
zclose__.nadd = nadd;
zclose__.ndel = ndel;
zclose__.iadd = iadd;
zclose__.idel = idel;
}

zcoord_st zcoord__;

void set_zcoord_data_ (int* iz, double* zbond, double* zang, double* ztors) {
zcoord__.iz = iz;
zcoord__.zbond = zbond;
zcoord__.zang = zang;
zcoord__.ztors = ztors;
}

