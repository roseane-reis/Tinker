      subroutine ommdata ()
      use action
      use align
      use analyz
      use angang
      use angbnd
      use angpot
      use angtor
      use argue
      use ascii
      use atmlst
      use atomid
      use atoms
      use bath
      use bitor
      use bndpot
      use bndstr
      use bound
      use boxes
      use cell
      use charge
      use chgpen
      use chgpot
      use chgtrn
      use chrono
      use chunks
      use couple
      use ctrpot
      use deriv
      use dipole
      use disgeo
      use disp
      use dma
      use domega
      use dsppot
      use energi
      use ewald
      use faces
      use fft
      use fields
      use files
      use fracs
      use freeze
      use gkstuf
      use group
      use hescut
      use hessn
      use hpmf
      use ielscf
      use improp
      use imptor
      use inform
      use inter
      use iounit
      use kanang
      use kangs
      use kantor
      use katoms
      use kbonds
      use kchrge
      use kcpen
      use kctrn
      use kdipol
      use kdsp
      use keys
      use khbond
      use kiprop
      use kitors
      use kmulti
      use kopbnd
      use kopdst
      use korbs
      use kpitor
      use kpolr
      use krepl
      use kstbnd
      use ksttor
      use ktorsn
      use ktrtor
      use kurybr
      use kvdwpr
      use kvdws
      use light
      use limits
      use linmin
      use math
      use mdstuf
      use merck
      use minima
      use molcul
      use moldyn
      use moment
      use mplpot
      use mpole
      use mrecip
      use mutant
      use neigh
      use nonpol
      use nucleo
      use omega
      use opbend
      use opdist
      use openmm
      use openmp
      use orbits
      use output
      use params
      use paths
      use pbstuf
      use pdb
      use phipsi
      use piorbs
      use pistuf
      use pitors
      use pme
      use polar
      use polgrp
      use polopt
      use polpcg
      use polpot
      use poltcg
      use potent
      use potfit
      use ptable
      use qmstuf
      use refer
      use repel
      use reppot
      use resdue
      use restrn
      use rgddyn
      use rigid
      use ring
      use rotbnd
      use rxnfld
      use rxnpot
      use scales
      use sequen
      use shunt
      use sizes
      use socket
      use solute
      use stodyn
      use strbnd
      use strtor
      use syntrn
      use tarray
      use titles
      use torpot
      use tors
      use tortor
      use tree
      use units
      use uprior
      use urey
      use urypot
      use usage
      use valfit
      use vdw
      use vdwpot
      use vibs
      use virial
      use warp
      use xtals
      use zclose
      use zcoord
      implicit none
      call set_action_data (neb,nea,neba,neub,neaa,neopb,neopd,neid,
     &                      neit,net,nept,nebt,neat,nett,nev,ner,nedsp,
     &                      nec,necd,ned,nem,nep,nect,new,nerxf,nes,
     &                      nelf,neg,nex)
      call set_align_data (nfit,ifit,wfit)
      call set_analyz_data (aesum,aeb,aea,aeba,aeub,aeaa,aeopb,aeopd,
     &                      aeid,aeit,aet,aept,aebt,aeat,aett,aev,aer,
     &                      aedsp,aec,aecd,aed,aem,aep,aect,aerxf,aes,
     &                      aelf,aeg,aex)
      call set_angang_data (nangang,iaa,kaa)
      call set_angbnd_data (nangle,iang,ak,anat,afld)
      call set_angpot_data (angunit,stbnunit,aaunit,opbunit,opdunit,
     &                      cang,qang,pang,sang,copb,qopb,popb,sopb,
     &                      copd,qopd,popd,sopd,opbtyp,angtyp)
      call set_angtor_data (nangtor,iat,kant)
      call set_argue_data (maxarg,narg,listarg,arg)
      call set_ascii_data (null,tab,linefeed,formfeed,carriage,escape,
     &                     space,exclamation,quote,pound,dollar,percent,
     &                     ampersand,apostrophe,asterisk,plus,comma,
     &                     minus,period,frontslash,colon,semicolon,
     &                     equal,question,atsign,backslash,caret,
     &                     underbar,vertical,tilde)
      call set_atmlst_data (bndlist,anglist)
      call set_atomid_data (tag,class,atomic,valence,mass,name,story)
      call set_atoms_data (n,type,x,y,z)
      call set_bath_data (maxnose,voltrial,kelvin,atmsph,tautemp,
     &                    taupres,compress,collide,eta,volmove,vbar,
     &                    qbar,gbar,vnh,qnh,gnh,isothermal,isobaric,
     &                    anisotrop,volscale,barostat,thermostat)
      call set_bitor_data (nbitor,ibitor)
      call set_bndpot_data (cbnd,qbnd,bndunit,bndtyp)
      call set_bndstr_data (nbond,ibnd,bk,bl)
      call set_bound_data (polycut,polycut2,use_bounds,use_replica,
     &                     use_polymer)
      call set_boxes_data (xbox,ybox,zbox,alpha,beta,gamma,xbox2,ybox2,
     &                     zbox2,box34,volbox,beta_sin,beta_cos,
     &                     gamma_sin,gamma_cos,beta_term,gamma_term,
     &                     lvec,recip,orthogonal,monoclinic,triclinic,
     &                     octahedron,spacegrp)
      call set_cell_data (ncell,icell,xcell,ycell,zcell,xcell2,ycell2,
     &                    zcell2)
      call set_charge_data (nion,iion,jion,kion,pchg)
      call set_chgpen_data (ncp,pcore,pval,palpha)
      call set_chgpot_data (electric,dielec,ebuffer,c2scale,c3scale,
     &                      c4scale,c5scale,neutnbr,neutcut)
      call set_chgtrn_data (nct,chgct,dmpct)
      call set_chrono_data (twall,tcpu)
      call set_chunks_data (nchunk,nchk1,nchk2,nchk3,ngrd1,ngrd2,ngrd3,
     &                      nlpts,nrpts,grdoff,pmetable)
      call set_couple_data (n12,n13,n14,n15,i12,i13,i14,i15)
      call set_ctrpot_data (ctrntyp)
      call set_deriv_data (desum,deb,dea,deba,deub,deaa,deopb,deopd,
     &                     deid,deit,det,dept,debt,deat,dett,dev,der,
     &                     dedsp,dec,decd,ded,dem,dep,dect,derxf,des,
     &                     delf,deg,dex)
      call set_dipole_data (ndipole,idpl,bdpl,sdpl)
      call set_disgeo_data (vdwmax,compact,pathmax,dbnd,georad,
     &                      use_invert,use_anneal)
      call set_disp_data (ndisp,idisp,csixpr,csix,adisp)
      call set_dma_data (mp,dpx,dpy,dpz,q20,q21c,q21s,q22c,q22s)
      call set_domega_data (tesum,teb,tea,teba,teub,teaa,teopb,teopd,
     &                      teid,teit,tet,tept,tebt,teat,tett,tev,ter,
     &                      tedsp,tec,tecd,ted,tem,tep,tect,terxf,tes,
     &                      telf,teg,tex)
      call set_dsppot_data (dsp2scale,dsp3scale,dsp4scale,dsp5scale,
     &                      use_dcorr)
      call set_energi_data (esum,eb,ea,eba,eub,eaa,eopb,eopd,eid,eit,et,
     &                      ept,ebt,eat,ett,ev,er,edsp,ec,ecd,ed,em,ep,
     &                      ect,erxf,es,elf,eg,ex)
      call set_ewald_data (aewald,aeewald,apewald,adewald,boundary)
      call set_faces_data (maxcls,maxtt,maxt,maxp,maxv,maxen,maxfn,maxc,
     &                     maxep,maxfs,maxcy,mxcyep,maxfp,mxfpcy,na,pr,
     &                     ar,axyz,skip,nosurf,afree,abur,cls,clst,acls,
     &                     ntt,ttfe,ttle,enext,tta,ttbur,ttfree,nt,tfe,
     &                     ta,tr,t,tax,tfree,np,pa,p,nv,va,vp,vxyz,nen,
     &                     nfn,env,fnen,nc,ca,ct,cr,c,nep,epc,epv,afe,
     &                     ale,epnext,nfs,fsen,fsep,ncy,cynep,cyep,nfp,
     &                     fpa,fpncy,fpcy)
      call set_fft_data (maxprime,iprime,planf,planb,ffttable,ffttyp)
      call set_fields_data (maxbio,biotyp,forcefield)
      call set_files_data (nprior,ldir,leng,filename,outfile)
      call set_fracs_data (xfrac,yfrac,zfrac)
      call set_freeze_data (nrat,nratx,iratx,kratx,irat,rateps,krat,
     &                      use_rattle,ratimage)
      call set_gkstuf_data (gkc,gkr)
      call set_group_data (ngrp,kgrp,grplist,igrp,grpmass,wgrp,
     &                     use_group,use_intra,use_inter)
      call set_hescut_data (hesscut)
      call set_hessn_data (hessx,hessy,hessz)
      call set_hpmf_data (rcarbon,rwater,acsurf,safact,tslope,toffset,
     &                    hpmfcut,h1,h2,h3,c1,c2,c3,w1,w2,w3,npmf,ipmf,
     &                    rpmf,acsa)
      call set_ielscf_data (nfree_aux,tautemp_aux,kelvin_aux,uaux,upaux,
     &                      vaux,vpaux,aaux,apaux,use_ielscf)
      call set_improp_data (niprop,iiprop,kprop,vprop)
      call set_imptor_data (nitors,iitors,itors1,itors2,itors3)
      call set_inform_data (maxask,digits,iprint,iwrite,isend,silent,
     &                      verbose,debug,holdup,abort)
      call set_inter_data (einter)
      call set_iounit_data (input,iout)
      call set_kanang_data (anan)
      call set_kangs_data (maxna,maxna5,maxna4,maxna3,maxnap,maxnaf,
     &                     acon,acon5,acon4,acon3,aconp,aconf,ang,ang5,
     &                     ang4,ang3,angp,angf,ka,ka5,ka4,ka3,kap,kaf)
      call set_kantor_data (maxnat,atcon,kat)
      call set_katoms_data (atmcls,atmnum,ligand,weight,symbol,describe)
      call set_kbonds_data (maxnb,maxnb5,maxnb4,maxnb3,maxnel,bcon,
     &                      bcon5,bcon4,bcon3,blen,blen5,blen4,blen3,
     &                      dlen,kb,kb5,kb4,kb3,kel)
      call set_kchrge_data (chg)
      call set_kcpen_data (cpele,cpalp)
      call set_kctrn_data (ctchg,ctdmp)
      call set_kdipol_data (maxnd,maxnd5,maxnd4,maxnd3,dpl,dpl5,dpl4,
     &                      dpl3,pos,pos5,pos4,pos3,kd,kd5,kd4,kd3)
      call set_kdsp_data (dspsix,dspdmp)
      call set_keys_data (maxkey,nkey,keyline)
      call set_khbond_data (maxnhb,radhb,epshb,khb)
      call set_kiprop_data (maxndi,dcon,tdi,kdi)
      call set_kitors_data (maxnti,ti1,ti2,ti3,kti)
      call set_kmulti_data (maxnmp,multip,mpaxis,kmp)
      call set_kopbnd_data (maxnopb,opbn,kopb)
      call set_kopdst_data (maxnopd,opds,kopd)
      call set_korbs_data (maxnpi,maxnpi5,maxnpi4,sslope,sslope5,
     &                     sslope4,tslope,tslope5,tslope4,electron,
     &                     ionize,repulse,kpi,kpi5,kpi4)
      call set_kpitor_data (maxnpt,ptcon,kpt)
      call set_kpolr_data (pgrp,polr,athl,ddir)
      call set_krepl_data (prsiz,prdmp,prele)
      call set_kstbnd_data (maxnsb,stbn,ksb)
      call set_ksttor_data (maxnbt,btcon,kbt)
      call set_ktorsn_data (maxnt,maxnt5,maxnt4,t1,t2,t3,t4,t5,t6,t15,
     &                      t25,t35,t45,t55,t65,t14,t24,t34,t44,t54,t64,
     &                      kt,kt5,kt4)
      call set_ktrtor_data (maxntt,maxtgrd,maxtgrd2,tnx,tny,ttx,tty,tbf,
     &                      tbx,tby,tbxy,ktt)
      call set_kurybr_data (maxnu,ucon,dst13,ku)
      call set_kvdwpr_data (maxnvp,radpr,epspr,kvpr)
      call set_kvdws_data (rad,eps,rad4,eps4,reduct)
      call set_light_data (nlight,kbx,kby,kbz,kex,key,kez,locx,locy,
     &                     locz,rgx,rgy,rgz)
      call set_limits_data (vdwcut,repcut,dispcut,chgcut,dplcut,
     &                      mpolecut,ctrncut,vdwtaper,reptaper,
     &                      disptaper,chgtaper,dpltaper,mpoletaper,
     &                      ctrntaper,ewaldcut,dewaldcut,usolvcut,
     &                      use_ewald,use_dewald,use_lights,use_list,
     &                      use_vlist,use_dlist,use_clist,use_mlist,
     &                      use_ulist)
      call set_linmin_data (intmax,stpmin,stpmax,cappa,slpmax,angmax)
      call set_math_data (pi,elog,radian,logten,twosix,sqrtpi,sqrttwo,
     &                    sqrtthree)
      call set_mdstuf_data (nfree,irest,bmnmix,arespa,dorest,integrate)
      call set_merck_data (nlignes,bt_1,eqclass,crd,val,pilp,mltb,arom,
     &                     lin,sbmb,mmffarom,mmffaromc,mmffaroma,rad0,
     &                     paulel,r0ref,kbref,mmff_kb,mmff_kb1,mmff_b0,
     &                     mmff_b1,mmff_ka,mmff_ka1,mmff_ka2,mmff_ka3,
     &                     mmff_ka4,mmff_ka5,mmff_ka6,mmff_ka7,mmff_ka8,
     &                     mmff_ang0,mmff_ang1,mmff_ang2,mmff_ang3,
     &                     mmff_ang4,mmff_ang5,mmff_ang6,mmff_ang7,
     &                     mmff_ang8,stbn_abc,stbn_cba,stbn_abc1,
     &                     stbn_cba1,stbn_abc2,stbn_cba2,stbn_abc3,
     &                     stbn_cba3,stbn_abc4,stbn_cba4,stbn_abc5,
     &                     stbn_cba5,stbn_abc6,stbn_cba6,stbn_abc7,
     &                     stbn_cba7,stbn_abc8,stbn_cba8,stbn_abc9,
     &                     stbn_cba9,stbn_abc10,stbn_cba10,stbn_abc11,
     &                     stbn_cba11,defstbn_abc,defstbn_cba,t1_1,t2_1,
     &                     t3_1,t1_2,t2_2,t3_2,kt_1,kt_2,g,alph,nn,da,
     &                     bci,bci_1,pbci,fcadj)
      call set_minima_data (maxiter,nextiter,fctmin,hguess)
      call set_molcul_data (nmol,imol,kmol,molcule,totmass,molmass)
      call set_moldyn_data (v,a,aalt)
      call set_moment_data (netchg,netdpl,netqdp,xdpl,ydpl,zdpl,xxqdp,
     &                      xyqdp,xzqdp,yxqdp,yyqdp,yzqdp,zxqdp,zyqdp,
     &                      zzqdp)
      call set_mplpot_data (m2scale,m3scale,m4scale,m5scale,use_chgpen,
     &                      pentyp)
      call set_mpole_data (maxpole,npole,ipole,polsiz,pollist,zaxis,
     &                     xaxis,yaxis,pole,rpole,spole,srpole,polaxe)
      call set_mrecip_data (vmxx,vmyy,vmzz,vmxy,vmxz,vmyz,cmp,fmp,cphi,
     &                      fphi)
      call set_mutant_data (nmut,vcouple,imut,type0,class0,type1,class1,
     &                      lambda,tlambda,vlambda,elambda,scexp,
     &                      scalpha,mut)
      call set_neigh_data (maxvlst,maxelst,maxulst,nvlst,vlst,nelst,
     &                     elst,nulst,ulst,lbuffer,pbuffer,lbuf2,pbuf2,
     &                     vbuf2,vbufx,dbuf2,dbufx,cbuf2,cbufx,mbuf2,
     &                     mbufx,ubuf2,ubufx,xvold,yvold,zvold,xeold,
     &                     yeold,zeold,xuold,yuold,zuold,dovlst,dodlst,
     &                     doclst,domlst,doulst)
      call set_nonpol_data (epso,epsh,rmino,rminh,awater,slevy,solvprs,
     &                      surften,spcut,spoff,stcut,stoff,rcav,rdisp,
     &                      cdisp)
      call set_nucleo_data (pucker,glyco,bkbone,dblhlx,deoxy,hlxform)
      call set_omega_data (nomega,iomega,zline,dihed)
      call set_opbend_data (nopbend,iopb,opbk)
      call set_opdist_data (nopdist,iopd,opdk)
      call set_openmm_data (ommhandle,cudaprecision,ommplatform,
     &                      cudadevice)
      call set_openmp_data (nproc,nthread)
      call set_orbits_data (qorb,worb,emorb)
      call set_output_data (archive,noversion,overwrite,cyclesave,
     &                      velsave,frcsave,uindsave,coordtype)
      call set_params_data (maxprm,nprm,prmline)
      call set_paths_data (pnorm,acoeff,pc0,pc1,pvect,pstep,pzet,gc)
      call set_pbstuf_data (maxion,ionn,dime,ionq,pbe,pdie,sdie,srad,
     &                      swin,sdens,smin,grid,gcent,cgrid,cgcent,
     &                      fgrid,fgcent,ionr,ionc,apbe,pbr,pbep,pbfp,
     &                      pbtp,pbeuind,pbeuinp,pbtyp,pbsoln,bcfl,chgm,
     &                      srfm)
      call set_pdb_data (npdb,nres,resnum,resatm,npdb12,ipdb12,pdblist,
     &                   xpdb,ypdb,zpdb,altsym,pdbres,pdbatm,pdbtyp,
     &                   chnsym,instyp)
      call set_phipsi_data (chiral,disulf,phi,psi,omega,chi)
      call set_piorbs_data (norbit,nconj,reorbit,nbpi,ntpi,iorbit,iconj,
     &                      kconj,piperp,ibpi,itpi,pbpl,pnpl,listpi)
      call set_pistuf_data (bkpi,blpi,kslope,lslope,torsp2)
      call set_pitors_data (npitors,ipit,kpit)
      call set_pme_data (nfft1,nfft2,nfft3,nefft1,nefft2,nefft3,ndfft1,
     &                   ndfft2,ndfft3,bsorder,bseorder,bsporder,
     &                   bsdorder,igrid,bsmod1,bsmod2,bsmod3,bsbuild,
     &                   thetai1,thetai2,thetai3,qgrid,qfac)
      call set_polar_data (npolar,ipolar,polarity,thole,dirdamp,pdamp,
     &                     udir,udirp,udirs,udirps,uind,uinp,uinds,
     &                     uinps,uexact,douind)
      call set_polgrp_data (maxp11,maxp12,maxp13,maxp14,np11,np12,np13,
     &                      np14,ip11,ip12,ip13,ip14)
      call set_polopt_data (maxopt,optorder,optlevel,copt,copm,uopt,
     &                      uoptp,uopts,uoptps,fopt,foptp)
      call set_polpcg_data (mindex,pcgpeek,minv,pcgprec,pcgguess)
      call set_polpot_data (politer,poleps,p2scale,p3scale,p4scale,
     &                      p5scale,p2iscale,p3iscale,p4iscale,p5iscale,
     &                      d1scale,d2scale,d3scale,d4scale,u1scale,
     &                      u2scale,u3scale,u4scale,w2scale,w3scale,
     &                      w4scale,w5scale,udiag,dpequal,use_thole,
     &                      use_dirdamp,poltyp)
      call set_poltcg_data (tcgorder,tcgnab,tcgpeek,uad,uap,ubd,ubp,
     &                      tcgguess)
      call set_potent_data (use_bond,use_angle,use_strbnd,use_urey,
     &                      use_angang,use_opbend,use_opdist,use_improp,
     &                      use_imptor,use_tors,use_pitors,use_strtor,
     &                      use_angtor,use_tortor,use_vdw,use_repuls,
     &                      use_disp,use_charge,use_chgdpl,use_dipole,
     &                      use_mpole,use_polar,use_chgtrn,use_rxnfld,
     &                      use_solv,use_metal,use_geom,use_extra,
     &                      use_born,use_orbit)
      call set_potfit_data (nconf,namax,ngatm,nfatm,npgrid,ipgrid,resp,
     &                      xdpl0,ydpl0,zdpl0,xxqdp0,xyqdp0,xzqdp0,
     &                      yyqdp0,yzqdp0,zzqdp0,fit0,fchg,fpol,pgrid,
     &                      epot,use_dpl,use_qdp,fit_mpl,fit_dpl,
     &                      fit_qdp,fitchg,fitpol,gatm,fatm)
      call set_ptable_data (maxele,atmass,vdwrad,covrad,elemnt)
      call set_qmstuf_data (ngatom,egau,gx,gy,gz,gfreq,gforce,gh)
      call set_refer_data (nref,refltitle,refleng,reftyp,n12ref,i12ref,
     &                     xboxref,yboxref,zboxref,alpharef,betaref,
     &                     gammaref,xref,yref,zref,refnam,reffile,
     &                     reftitle)
      call set_repel_data (nrep,sizpr,dmppr,elepr)
      call set_reppot_data (r2scale,r3scale,r4scale,r5scale)
      call set_resdue_data (maxamino,maxnuc,ntyp,catyp,ctyp,hntyp,otyp,
     &                      hatyp,cbtyp,nntyp,cantyp,cntyp,hnntyp,ontyp,
     &                      hantyp,nctyp,cactyp,cctyp,hnctyp,octyp,
     &                      hactyp,o5typ,c5typ,h51typ,h52typ,c4typ,
     &                      h4typ,o4typ,c1typ,h1typ,c3typ,h3typ,c2typ,
     &                      h21typ,o2typ,h22typ,o3typ,ptyp,optyp,h5ttyp,
     &                      h3ttyp,amino1,nuclz1,amino,nuclz)
      call set_restrn_data (npfix,ndfix,nafix,ntfix,ngfix,nchir,ipfix,
     &                      kpfix,idfix,iafix,itfix,igfix,ichir,depth,
     &                      width,rwall,xpfix,ypfix,zpfix,pfix,dfix,
     &                      afix,tfix,gfix,chir,use_basin,use_wall)
      call set_rgddyn_data (xcmo,ycmo,zcmo,vcm,wcm,lm,vc,wc,linear)
      call set_rigid_data (xrb,yrb,zrb,rbc,use_rigid)
      call set_ring_data (nring3,nring4,nring5,nring6,nring7,iring3,
     &                    iring4,iring5,iring6,iring7)
      call set_rotbnd_data (nrot,rot,use_short)
      call set_rxnfld_data (ijk,b1,b2)
      call set_rxnpot_data (rfterms,rfsize,rfbulkd)
      call set_scales_data (scale,set_scale)
      call set_sequen_data (nseq,nchain,ichain,seqtyp,chnnam,seq,chntyp)
      call set_shunt_data (off,off2,cut,cut2,c0,c1,c2,c3,c4,c5,f0,f1,f2,
     &                     f3,f4,f5,f6,f7)
      call set_sizes_data (maxatm,maxtyp,maxclass,maxval,maxref,maxgrp,
     &                     maxres,maxfix)
      call set_socket_data (skttyp,cstep,cdt,cenergy,sktstart,sktstop,
     &                      use_socket)
      call set_solute_data (doffset,p1,p2,p3,p4,p5,rsolv,asolv,rborn,
     &                      drb,drbp,drobc,gpol,shct,aobc,bobc,gobc,
     &                      vsolv,wace,s2ace,uace,solvtyp,borntyp)
      call set_stodyn_data (friction,fgamma,use_sdarea)
      call set_strbnd_data (nstrbnd,isb,sbk)
      call set_strtor_data (nstrtor,ist,kst)
      call set_syntrn_data (tpath,ppath,xmin1,xmin2,xm)
      call set_tarray_data (ntpair,tindex,tdipdip)
      call set_titles_data (ltitle,title)
      call set_torpot_data (idihunit,itorunit,torsunit,ptorunit,
     &                      storunit,atorunit,ttorunit)
      call set_tors_data (ntors,itors,tors1,tors2,tors3,tors4,tors5,
     &                    tors6)
      call set_tortor_data (ntortor,itt)
      call set_tree_data (maxpss,nlevel,etree,ilevel)
      call set_units_data (avogadro,lightspd,boltzmann,gasconst,elemchg,
     &                     vacperm,emass,planck,joule,ekcal,bohr,
     &                     hartree,evolt,efreq,coulomb,debye,prescon)
      call set_uprior_data (maxpred,nualt,maxualt,gear,aspc,bpred,
     &                      bpredp,bpreds,bpredps,udalt,upalt,usalt,
     &                      upsalt,use_pred,polpred)
      call set_urey_data (nurey,iury,uk,ul)
      call set_urypot_data (cury,qury,ureyunit)
      call set_usage_data (nuse,iuse,use)
      call set_valfit_data (fit_bond,fit_angle,fit_strbnd,fit_urey,
     &                      fit_opbend,fit_tors,fit_struct,fit_force)
      call set_vdw_data (nvdw,ivdw,jvdw,ired,kred,xred,yred,zred,radmin,
     &                   epsilon,radmin4,epsilon4,radhbnd,epshbnd)
      call set_vdwpot_data (maxgauss,ngauss,igauss,abuck,bbuck,cbuck,
     &                      ghal,dhal,v2scale,v3scale,v4scale,v5scale,
     &                      use_vcorr,vdwindex,radtyp,radsiz,gausstyp,
     &                      radrule,epsrule,vdwtyp)
      call set_vibs_data (phi,phik,pwork)
      call set_virial_data (vir,use_virial)
      call set_warp_data (deform,difft,diffv,diffc,m2,use_smooth,
     &                    use_dem,use_gda,use_tophat,use_stophat)
      call set_xtals_data (maxlsq,maxrsd,nxtal,nvary,ivary,iresid,vary,
     &                     e0_lattice,vartyp,rsdtyp)
      call set_zclose_data (nadd,ndel,iadd,idel)
      call set_zcoord_data (iz,zbond,zang,ztors)
      return
      end
