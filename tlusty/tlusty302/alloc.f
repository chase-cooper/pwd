

!*********************************************************************


      subroutine alloc_ini

      use accura
      use basics
      use atomic
      use hydodf
      use modelq

      INTEGER :: MIO0   =  170  ! initial max.num. of explicit ions
      INTEGER :: MLEVE0 = 1134  ! initial max.num. of explicit levels
      INTEGER :: MME0   =   12  ! initial max.num. of merged levels

      IF(IOPTAB.LT.0) THEN
         MIO0=1
         MLEVEL0=1
         MME0=1
      END IF

      write(*,*) 'alloc_ini MATOM,MDEPTH,MIO0,MLEVE0,MME0    ',
     *                      MATOM,MDEPTH,MIO0,MLEVE0,MME0
      write(66,*) 'alloc_ini MATOM,MDEPTH,MIO0,MLEVE0,MME0    ',
     *                      MATOM,MDEPTH,MIO0,MLEVE0,MME0
!     from common/atopar
      ALLOCATE (AMASS(MATOM))
      ALLOCATE (ABUND(MATOM,MDEPTH))
      ALLOCATE (ABNDD(MATOM,MDEPTH))
      ALLOCATE (NUMAT(MATOM))
      ALLOCATE (N0A(MATOM))
      ALLOCATE (NKA(MATOM))
      ALLOCATE (NREF(MATOM))
      ALLOCATE (IATEX(MATOM))
      ALLOCATE (NREFS(MATOM,MDEPTH))
      ALLOCATE (IADOP(MATOM))
      ALLOCATE (IIFIX(MATOM))
      ALLOCATE (ENEV(MATOM,26))
      ALLOCATE (LGR(MATOM),LRM(MATOM))
      ALLOCATE (VTURB(MDEPTH),VTURBS(MDEPTH))

!     The rest unnecessary for a pure-table mode
 
      IF(IOPTAB.LT.0) RETURN

!     from common/ionpar
      ALLOCATE (FF(MIO0))
      ALLOCATE (CHARG2(MIO0))
      ALLOCATE (NFIRST(MIO0))
      ALLOCATE (NLAST(MIO0))
      ALLOCATE (NNEXT(MIO0))
      ALLOCATE (MODFF(MIO0))
      ALLOCATE (IZ(MIO0))
      ALLOCATE (IUPSUM(MIO0))
      ALLOCATE (ICUP(MIO0))
      ALLOCATE (ILTE(MIO0))
      ALLOCATE (ILTION(MIO0))
      ALLOCATE (IATI(MIO0))
      ALLOCATE (IZI(MIO0))
      ALLOCATE (NLEVS(MIO0))
      ALLOCATE (NLLIM(MIO0))
      ALLOCATE (IKOBS(MIO0))
      ALLOCATE (IFMETA(MIO0))

      ALLOCATE (INODF1(MIO0))
      ALLOCATE (INODF2(MIO0))
      ALLOCATE (INBFCS(MIO0))

!     from common/levpar
      ALLOCATE (ENION(MLEVE0))
      ALLOCATE (G(MLEVE0))
      ALLOCATE (NQUANT(MLEVE0))
      ALLOCATE (IMODL(MLEVE0),IMODL0(MLEVE0))
      ALLOCATE (IFWOP(MLEVE0))
      ALLOCATE (ILTLEV(MLEVE0))
      ALLOCATE (IMRG(MLEVE0))
      ALLOCATE (FRODF(MLEVE0))
      ALLOCATE (IATM(MLEVE0))
      ALLOCATE (IEL(MLEVE0))
      ALLOCATE (ILK(MLEVE0))
      ALLOCATE (IGUIDE(MLEVE0))
      ALLOCATE (IIFOR(MLEVE0))
      ALLOCATE (IIEXP(MLEVE0))
      ALLOCATE (INDLEV(MLEVE0))

      ALLOCATE (IIMER(MME0))

      end subroutine alloc_ini

!********************************************************


      subroutine alloc_tra

      use accura
      use basics
      use atomic
      use hydodf
      use modelq

      IF(IOPTAB.LT.0) RETURN
  
      write(*,*) 'alloc_tra MTRAN0,MLEVEL,MBF,MFIT,MCROSS    ',
     *                      MTRAN0,MLEVEL,MBF,MFIT,MCROSS
      write(*,*) 'alloc_tra MION,MXTCOL,MCFIT,MVOIG0         ',
     *                      MION,MXTCOL,MCFIT,MVOIG0
      write(66,*) 'alloc_tra MTRAN0,MLEVEL,MBF,MFIT,MCROSS    ',
     *                      MTRAN0,MLEVEL,MBF,MFIT,MCROSS
      write(66,*) 'alloc_tra MION,MXTCOL,MCFIT,MVOIG0         ',
     *                      MION,MXTCOL,MCFIT,MVOIG0
!
      ALLOCATE (FR0(MTRAN0))
      ALLOCATE (OSC0(MTRAN0))
      ALLOCATE (CPAR(MTRAN0))
      ALLOCATE (FRQMX(MTRAN0))
      ALLOCATE (FR0PC(MTRAN0))
      ALLOCATE (OMECOL(MLEVEL,MLEVEL))
      ALLOCATE (ILOW(MTRAN0))
      ALLOCATE (IUP(MTRAN0))
      ALLOCATE (INDEXP(MTRAN0))
      ALLOCATE (KFR0(MTRAN0))
      ALLOCATE (KFR1(MTRAN0))
      ALLOCATE (ILUCTR(MTRAN0))
      ALLOCATE (IFC0(MTRAN0))
      ALLOCATE (IFC1(MTRAN0))
      ALLOCATE (IFR0(MTRAN0))
      ALLOCATE (IFR1(MTRAN0))
      ALLOCATE (ITRA(MLEVEL,MLEVEL))
      ALLOCATE (IPROF(MTRAN0))
      ALLOCATE (IPROF0(MTRAN0))
      ALLOCATE (ICOL(MTRAN0))
      ALLOCATE (INTMOD(MTRAN0))
      ALLOCATE (ITRCON(MTRAN0))
      ALLOCATE (IDIEL(MTRAN0))
      ALLOCATE (IJTF(MTRAN0))
      ALLOCATE (JNDODF(MTRAN0))
      ALLOCATE (MCDW(MTRAN0))
      ALLOCATE (LCOMP(MTRAN0))
      ALLOCATE (LINE(MTRAN0))
      ALLOCATE (LALI(MTRAN0))
      ALLOCATE (LEXP(MTRAN0))

      ALLOCATE (S0CS(MBF))
      ALLOCATE (ALFCS(MBF))
      ALLOCATE (BETCS(MBF))
      ALLOCATE (GAMCS(MBF))
      ALLOCATE (IBF(MBF))
      ALLOCATE (ITRBF(MBF))
CC    ALLOCATE (ITRBF(MLEVEL))

      ALLOCATE (IPZERT(MLEVEL),IGZERT(MLEVEL),INDLGZ(MLEVEL))
      ALLOCATE (IINONZ(MLEVEL))

      ALLOCATE (CTOP(MFIT,MCROSS)) !sigma=alog10(sigma/10^-18) of fit point
      ALLOCATE (XTOP(MFIT,MCROSS)) ! x = alog10(nu/nu0) of fit point

      ALLOCATE (CTEMP(MXTCOL,MCFIT,MTRAN0)) ! temperature vs.
      ALLOCATE (CRATE(MXTCOL,MCFIT,MTRAN0)) ! collisional rates

      ALLOCATE (GAMAR(MVOIG0))
      ALLOCATE (STARK1(MVOIG0))
      ALLOCATE (STARK2(MVOIG0))
      ALLOCATE (STARK3(MVOIG0))
      ALLOCATE (VDWH(MVOIG0))

      ALLOCATE (NEVKU(MION))
      ALLOCATE (NODKU(MION))
      ALLOCATE (XEV(MLEVEL,MION))
      ALLOCATE (XOD(MLEVEL,MION))

      ALLOCATE (IJFL(MLEVEL))

      end subroutine alloc_tra

C*******************************************************************

      subroutine alloc_freq0

      use accura
      use basics
      use modelq
      use atomic

      MFREQ0=655000
      if(IOPTAB.LT.0) MFREQ0=MFREQ

      write(*,*) 'alloc_freq MFREQ0                          ',
     *                        MFREQ0
      write(66,*) 'alloc_freq MFREQ0                          ',
     *                        MFREQ0

      IF(IOPTAB.GE.0) THEN
         ALLOCATE (FREQ(MFREQ0),W(MFREQ0))
         ALLOCATE (IJALI(MFREQ0),IJX(MFREQ0),JIK(MFREQ0))
      END IF
      ALLOCATE (PROF(MFREQ0))
      ALLOCATE (WCH(MFREQ0))
      ALLOCATE (IFREQB(MFREQ0),NLINES(MFREQ0))

      end subroutine alloc_freq0

C
C     ********************************************************************
C
      subroutine alloc_modelq

      use accura
      use basics
      use atomic
      use modelq
      use hydodf 
      use array1
      use molec

      write(*,*) 'alloc_models MDEPTH,MION,MLEVEL,MLVEXP,MMER',
     *                         MDEPTH,MION,MLEVEL,MLVEXP,MMER
      write(*,*) 'alloc_modelq MTRANS,MZZ,NLMX,MMCDW         ',
     *                         MTRANS,MZZ,NLMX,MMCDW
      write(66,*) 'alloc_models MDEPTH,MION,MLEVEL,MLVEXP,MMER',
     *                         MDEPTH,MION,MLEVEL,MLVEXP,MMER
      write(66,*) 'alloc_modelq MTRANS,MZZ,NLMX,MMCDW         ',
     *                         MTRANS,MZZ,NLMX,MMCDW
      ALLOCATE (DM(MDEPTH),TEMP(MDEPTH),ELEC(MDEPTH),DENS(MDEPTH))
      ALLOCATE (TOTN(MDEPTH))
      ALLOCATE (ANTO(MDEPTH),ANMA(MDEPTH),ANH1(MDEPTH),ZD(MDEPTH))
      ALLOCATE (HKT1(MDEPTH),TK1(MDEPTH),HKT21(MDEPTH),SQT1(MDEPTH))
      ALLOCATE (TEMP1(MDEPTH),ELEC1(MDEPTH),DENS1(MDEPTH))
      ALLOCATE (DENSI(MDEPTH),DENSIM(MDEPTH))
      ALLOCATE (ELSCAT(MDEPTH),ALAB(MDEPTH))
      ALLOCATE (DELDM(MDEPTH),DELDMZ(MDEPTH))
      ALLOCATE (THETAV(MDEPTH),VISCD(MDEPTH))

      ALLOCATE (POPUL(MLEVEL,MDEPTH),POPINV(MLEVEL,MDEPTH))
      ALLOCATE (POPGRP(MLEVEL))
      ALLOCATE (POP(MLEVEL),SBF(MLEVEL),DSBF(MLEVEL),USUM(MION))
      ALLOCATE (POPUL0(MLEVEL,MDEPTH)) 

      ALLOCATE (ANH2(MDEPTH),ANHM(MDEPTH)) 
      ALLOCATE (RPOP0(MLVEXP,MDEPTH))
!     ALLOCATE (IPZERO(MLEVEL,MDEPTH),IGZERO(MLVEXP,MDEPTH))
      ALLOCATE (WNHINT(NLMX,MDEPTH),WNHEII(NLMX,MDEPTH))
      ALLOCATE (WOP(MLEVEL,MDEPTH))
!     ALLOCATE (BFAC(MLEVEL,MDEPTH))

      ALLOCATE (REINT(MDEPTH),REDIF(MDEPTH))

!     ALLOCATE (VTURB(MDEPTH),VTURBS(MDEPTH))

      ALLOCATE (SBPSI(MLEVEL,MDEPTH),SBLPSI(MLEVEL,MDEPTH))
      ALLOCATE (DSBPST(MLEVEL,MDEPTH),DSBPSN(MLEVEL,MDEPTH))
!     ALLOCATE (ILTREF(MLEVEL,MDEPTH),IGUIDE(MLEVEL))
      ALLOCATE (ILTERF(MLEVEL,MDEPTH))
      ALLOCATE (PT(MLEVEL,MDEPTH),PN(MLEVEL,MDEPTH))
      ALLOCATE (PP(MLEVEL,MDEPTH))
      ALLOCATE (USUMS(MION,MDEPTH))
      ALLOCATE (DUSMT(MION,MDEPTH))
      ALLOCATE (DUSMN(MION,MDEPTH))
      ALLOCATE (DIESIG(MION,MDEPTH))
      ALLOCATE (SGM0(MMER),FRCH(MMER),SGEXT1(MMER,MDEPTH))
      ALLOCATE (GMER(MMER,MDEPTH))
      ALLOCATE (SGMSUM(NLMX,MMER,MDEPTH))
      ALLOCATE (SGMSUD(NLMX,MMER,MDEPTH))
      ALLOCATE (SGMG(MMER,MDEPTH))
      ALLOCATE (DUSUMT(MION),DUSUMN(MION))

      ALLOCATE (ELEC23(MDEPTH))
      ALLOCATE (ACOR(MDEPTH))
      ALLOCATE (Z3(MZZ))
      ALLOCATE (DWC1(MZZ,MDEPTH))
      ALLOCATE (DWC2(MDEPTH))
      ALLOCATE (DWF1(MMCDW,MDEPTH))

      ALLOCATE (GF0(MDEPTH),GF1(MDEPTH),GF2(MDEPTH))
      ALLOCATE (GF3(MDEPTH),GF4(MDEPTH),GF5(MDEPTH),GF6(MDEPTH))
      ALLOCATE (GF0D(MDEPTH),GF1D(MDEPTH),GF2D(MDEPTH))
      ALLOCATE (GF3D(MDEPTH),GF4D(MDEPTH),GF5D(MDEPTH))
      ALLOCATE (GF6D(MDEPTH),DELTT(MDEPTH))

      ALLOCATE (SFF3(MION,MDEPTH))
      ALLOCATE (SFF2(MION,MDEPTH))
      ALLOCATE (DSFF(MION,MDEPTH))
      ALLOCATE (CFFN(MDEPTH),CFFT(MDEPTH))

      ALLOCATE (XKF(MDEPTH),XKF1(MDEPTH),XKFB(MDEPTH))
      ALLOCATE (ABSO1(MDEPTH),EMIS1(MDEPTH),SCAT1(MDEPTH))
      ALLOCATE (ABSOT(MDEPTH))
      ALLOCATE (EMEL1(MDEPTH),ABSO1L(MDEPTH),EMIS1L(MDEPTH))
      ALLOCATE (ABSOPR(MDEPTH),EMISPR(MDEPTH))
      ALLOCATE (ABSO(MFREQ),EMIS(MFREQ),SCAT(MFREQ))

      ALLOCATE (DABT1(MDEPTH),DEMT1(MDEPTH))
      ALLOCATE (DABN1(MDEPTH),DEMN1(MDEPTH))
      ALLOCATE (DABM1(MDEPTH),DEMM1(MDEPTH))
      ALLOCATE (DABX1(MDEPTH),DEMX1(MDEPTH))
      ALLOCATE (DABP1(MLEVEL,MDEPTH),DEMP1(MLEVEL,MDEPTH))
      ALLOCATE (DRCH1(MLEVEL,MDEPTH),DRET1(MLEVEL,MDEPTH))
      ALLOCATE (ABSFF(MDEPTH),DABFT(MDEPTH),DABFN(MDEPTH))
      ALLOCATE (DSFDT(MDEPTH),DSFDN(MDEPTH),DSFDM(MDEPTH))
      ALLOCATE (DSFDP(MLVEXP,MDEPTH))
      ALLOCATE (DSFDTM(MDEPTH),DSFDNM(MDEPTH))
      ALLOCATE (DSFDPM(MLVEXP,MDEPTH))
      ALLOCATE (DSFDTP(MDEPTH),DSFDNP(MDEPTH))
      ALLOCATE (DSFDPP(MLVEXP,MDEPTH))
      ALLOCATE (DSFP1D(MLVEXP))
      ALLOCATE (ALIM1(MDEPTH),ALIP1(MDEPTH))
      ALLOCATE (DSCT1(MDEPTH),DSCN1(MDEPTH),DRHODT(MDEPTH))
 
      ALLOCATE (PTOTAL(MDEPTH),PGS(MDEPTH),PRADT(MDEPTH))
      ALLOCATE (PRADA(MDEPTH),FLRD(MDEPTH))
      ALLOCATE (ABROSD(MDEPTH),SUMDPL(MDEPTH))
      ALLOCATE (ABPLAD(MDEPTH))
      ALLOCATE (QFIX(MDEPTH))
      ALLOCATE (ALBE(MFREQ))
      ALLOCATE (DDN(MLEVEL))
      ALLOCATE (FLXC(MDEPTH),DELTA(MDEPTH))
      ALLOCATE (RSAX(MLEVEL),RSBX(MLEVEL))
      ALLOCATE (TAUROS(MDEPTH),TAUFLX(MDEPTH))
      ALLOCATE (TAUTHE(MDEPTH),THETA(MDEPTH))
      ALLOCATE (TROSS(MDEPTH))
      ALLOCATE (FCOOL(MDEPTH),FPRD(MDEPTH))
      ALLOCATE (FLTOT(MDEPTH),FLFIX(MDEPTH),FLEXP(MDEPTH))
      ALLOCATE (FCOOLI(MDEPTH))
      ALLOCATE (FPRAD(MDEPTH),GRAD(MDEPTH))

      ALLOCATE (PHMOL(MDEPTH),QADD(MDEPTH))
      ALLOCATE (ANATO(100,MDEPTH),ANION(100,MDEPTH))

      ALLOCATE (RR(99,99))
!     ALLOCATE (ABNDD(99,MDEPTH),ENEV(99,30))
!     ALLOCATE (IONIZ(99))
!     ALLOCATE (LGR(99),LRM(99))

      ALLOCATE (CHANT(MDEPTH))
CC    ALLOCATE (XK0(MLINH))

      ALLOCATE (NQLODF(MLEVEL),I1ODF(MLEVEL),I2ODF(MLEVEL))

      ALLOCATE (ABTRA(MTRANS,MDEPTH),EMTRA(MTRANS,MDEPTH))
      ALLOCATE (DEMLT(MTRANS,MDEPTH))

      ALLOCATE (DT(MDEPTH),TAURS(MDEPTH))
      ALLOCATE (FLUX(MFREQ),FH(MFREQ),Q0(MFREQ),UU0(MFREQ))

!     ALLOCATE (PSY0(MTOT,MDEPTH))
c     ALLOCATE (PSY2(MTOT,MDEPTH),PSY3(MTOT,MDEPTH))

      ALLOCATE (GAMJ(MDEPTH))

      IF(IFMOL.GT.0) ALLOCATE (ANMOL(600,MDEPTH))

      end subroutine alloc_modelq
C
C     ***************************************************
c
      subroutine alloc_iniex

      use accura
      use basics
      use modelq
      use array1

      write(*,*) 'alloc_iniex MFREX,MTOT                     ',
     *                         MFREX,MTOT
      write(66,*) 'alloc_iniex MFREX,MTOT                     ',
     *                         MFREX,MTOT
      IF(IFRYB.EQ.0) THEN
         ALLOCATE (RADEX(MFREX,MDEPTH),FAKEX(MFREX,MDEPTH))
         ALLOCATE (ABSOEX(MFREX,MDEPTH))
      END IF

      ALLOCATE (PSY0(MTOT,MDEPTH))

      end subroutine alloc_iniex
C
C     ***************************************************
c
      subroutine alloc_radt

      use accura
      use basics
      use modelq

      MFREQ1=1
      IF(ISPLIN.GE.5) MFREQ1=MFREQ
      MTRAN3=MTRANS
      IF(IFALI.EQ.6) MTRAN3=MTRANS
      write(*,*) 'alloc_radt MTRANS,MDEPTH,MMCDW      ',
     *                       MTRANS,MDEPTH,MMCDW
      write(66,*) 'alloc_radt MTRANS,MDEPTH,MMCDW      ',
     *                       MTRANS,MDEPTH,MMCDW

      ALLOCATE (RRU(MTRANS,MDEPTH),RRD(MTRANS,MDEPTH))
      ALLOCATE (DRDT(MTRANS,MDEPTH))
      ALLOCATE (COLRAT(MTRANS,MDEPTH),COLTAR(MTRANS,MDEPTH))
      ALLOCATE (COL(MTRANS),CLOC(MTRANS))
      ALLOCATE (GRD(MDEPTH),PRA(MDEPTH),PGS0(MDEPTH),ANTP(MDEPTH))

      ALLOCATE (FHD(MFREQ))

      ALLOCATE (RAD1(MDEPTH),ALI1(MDEPTH),FAK1(MDEPTH))
      ALLOCATE (RADCM(MFREQ,MDEPTH))
      ALLOCATE (ALIH1(MDEPTH))
      ALLOCATE (GRADF(MDEPTH,MFREQ))

      ALLOCATE (ESEMAT(MLEVEL,MLEVEL),BESE(MLEVEL))

      IF(ISPLIN.GE.5) THEN
         ALLOCATE (RAD(MFREQ,MDEPTH))
         ALLOCATE (FAK(MFREQ,MDEPTH),RADK(MFREQ,MDEPTH))
         ALLOCATE (RDDP(MTRANS,MDEPTH),RDDM(MTRANS,MDEPTH))
      END IF

      end subroutine alloc_radt

C     
C     ***************************************************      
c
      subroutine alloc_prf

      use accura
      use basics
      use modelq
      use odfpar

      write(*,*) 'alloc_prf MFREQ,MFREQP,MDODF,MCFE,MDEPTH   ',
     *                      MFREQ,MFREQP,MDODF,MCFE,MDEPTH
      write(66,*) 'alloc_prf MFREQ,MFREQP,MDODF,MCFE,MDEPTH   ',
     *                      MFREQ,MFREQP,MDODF,MCFE,MDEPTH
      ALLOCATE (PRFLIN(MDEPTH,MFREQP),PRF(MFREQP))
CC    ALLOCATE (SIGFE(MDODF,MCFE))

      ALLOCATE (KIJ(MFREQ))

      end subroutine alloc_prf


!*********************************************************************
      
      
      subroutine alloc_tra2
   
      use accura
      use basics
      use atomic
      use modelq

      write(*,*) 'alloc_tra2 MTRANS,MMCDW                    ',
     *           MTRANS,MMCDW
      write(66,*) 'alloc_tra2 MTRANS,MMCDW                    ',
     *           MTRANS,MMCDW

      allocate (linexp(mtrans))
      ALLOCATE (ITRCDW(MMCDW))
      
! in levset
      ALLOCATE (ILTREF(MLEVEL,MDEPTH))
      ALLOCATE (IPZERO(MLEVEL,MDEPTH),IGZERO(MLEVEL,MDEPTH))
      ALLOCATE (BFAC(MLEVEL,MDEPTH))
      
      end subroutine alloc_tra2


!*********************************************************************



      subroutine alloc_freq

      use accura
      use basics
      use modelq

      write(*,*) 'alloc_freq2 MFREQ,MFREQC,MCROSS,MDEPTH     ',
     *               MFREQ,MFREQC,MCROSS,MDEPTH
      write(66,*) 'alloc_freq2 MFREQ,MFREQC,MCROSS,MDEPTH     ',
     *               MFREQ,MFREQC,MCROSS,MDEPTH

      ALLOCATE (AIJBF(MFREQ))
      ALLOCATE (BFCS(MCROSS,MFREQC))
      ALLOCATE (IJBF(MFREQ),IJEX(MFREQ),IJFR(MFREQ))
      ALLOCATE (LSKIP(MDEPTH,MFREQ))
      ALLOCATE (W0E(MFREQ),BNUE(MFREQ),WC(MFREQ))

      ALLOCATE (IJLIN(MFREQ))
      allocate (sigec(mfreq))


      end subroutine alloc_freq




!*********************************************************************

      subroutine alloc_iter

      use accura
      use basics
      use iterat

      MITER=NITER+1
      MLAMBD=NLAMBD+1

      ALLOCATE (NITLAM(MITER))

      ALLOCATE (IFFIX(MITER),NETEXP(MITER),NETFIX(MITER))
      ALLOCATE (IETEXP(MITER,MLAMBD),IETFIX(MITER,MLAMBD))
      ALLOCATE (INHE0(MITER),INRE0(MITER),INPC0(MITER),
     *         INDL0(MITER),INSE0(MITER),INMP0(MITER),
     *         NN00(MITER),NDRE0(MITER),KANT(MITER))

      end subroutine alloc_iter


!*********************************************************************


      subroutine alloc_accel

      use accura
      use basics
      use iterat
      use modelq
      use accel

      write(*,*) 'alloc_accel iacpp,mtot',iacpp,mtot

      write(*,*) 'alloc_accel MLEVEL,MDEPTH,MTOT             ',
     *                         MLEVEL,MDEPTH,MTOT
      write(66,*) 'alloc_accel MLEVEL,MDEPTH,MTOT             ',
     *                         MLEVEL,MDEPTH,MTOT
      if(iacpp.gt.0) then
         allocate (popul1(mlevel,mdepth),popul2(mlevel,mdepth),
     *             popul3(mlevel,mdepth))
      end if

      write(*,*) 'alloc_accel iacpp,mtot',iacpp,mtot

      if(iacc.gt.0) then
         ALLOCATE (PSY1(MTOT,MDEPTH))
         ALLOCATE (PSY2(MTOT,MDEPTH),PSY3(MTOT,MDEPTH))
      end if

!     ALLOCATE (ISNG(MTOT))

      end subroutine alloc_accel


!*********************************************************************


      subroutine alloc_therm
  
      use accura
      use thermo

      allocate (SL(330,100),PL(330,100))

      end subroutine alloc_therm


!*********************************************************************

      
      subroutine alloc_molec

      use accura
      use molec

      ALLOCATE (C(600,5),PPMOL(600),APMLOG(600),
     *          XIP(100),XIP2(100),CCOMP(100), UIIDUI(100),
     *          P(100),FP(100),XKP(100),XK2(100))

      ALLOCATE (NELEM(5,600),NATO(5,600),MMAX(600),NELEMX(100))

      ALLOCATE (anion2(30,mdepth))
      ALLOCATE (entato(100),ention(100),entmol(600))
      ALLOCATE (uelem(100),ull(100),anden(800),
     *          aelem(100),ammol(600),
     *          anat0(100),anio0(100),anmo0(600),pfmol(600),
     *          denso(mdepth),eleco(mdepth),wmmo(mdepth))


      end subroutine alloc_molec




!*********************************************************************


      subroutine alloc_alipar

      use accura
      use basics
      use modelq
      use alipar

      write(*,*) 'alloc_alipar MDEPTH,MLVEXP                  ',
     *                         MDEPTH,MLVEXP

      ALLOCATE (REIT(MDEPTH),REIN(MDEPTH),REIM(MDEPTH))
      ALLOCATE (AREIT(MDEPTH),AREIN(MDEPTH),AREIM(MDEPTH))
      ALLOCATE (CREIT(MDEPTH),CREIN(MDEPTH),CREIM(MDEPTH))
      ALLOCATE (REIX(MDEPTH),CREIX(MDEPTH))
      ALLOCATE (REDX(MDEPTH),REDT(MDEPTH),REDN(MDEPTH))
      ALLOCATE (REDM(MDEPTH))
      ALLOCATE (REDXM(MDEPTH),REDTM(MDEPTH),REDNM(MDEPTH))
      ALLOCATE (REDMM(MDEPTH))
      ALLOCATE (REDTP(MDEPTH),REDNP(MDEPTH),REDXP(MDEPTH))
      ALLOCATE (REDMP(MDEPTH))


CC    IF(IFRYB.EQ.0) THEN
         ALLOCATE (HEIT(MDEPTH),HEIN(MDEPTH),HEIM(MDEPTH))
         ALLOCATE (HEIP(MLVEXP,MDEPTH))
         ALLOCATE (HEITM(MDEPTH),HEINM(MDEPTH),HEIMM(MDEPTH))
         ALLOCATE (HEIPM(MLVEXP,MDEPTH))
         ALLOCATE (HEITP(MDEPTH),HEINP(MDEPTH),HEIMP(MDEPTH))
         ALLOCATE (HEIPP(MLVEXP,MDEPTH))
         ALLOCATE (REIP(MLVEXP,MDEPTH),AREIP(MLVEXP,MDEPTH),
     *             CREIP(MLVEXP,MDEPTH),REDP(MLVEXP,MDEPTH))
         ALLOCATE (REDPM(MLVEXP,MDEPTH),REDPP(MLVEXP,MDEPTH))
         ALLOCATE (EHET(MDEPTH),EHEN(MDEPTH),ERET(MDEPTH),EREN(MDEPTH))
         ALLOCATE (APT(MLVEXP,MDEPTH),APN(MLVEXP,MDEPTH))
         ALLOCATE (APM(MLVEXP,MDEPTH))
         ALLOCATE (AAPT(MLVEXP,MDEPTH),AAPN(MLVEXP,MDEPTH))
         ALLOCATE (AAPM(MLVEXP,MDEPTH))
         ALLOCATE (APP(MLVEXP,MLVEXP,MDEPTH))
         ALLOCATE (AAPP(MLVEXP,MLVEXP,MDEPTH))
CC    END IF

CC    IF(IFALI.GT.5) THEN
         ALLOCATE (EHEP(MLVEXP,MDEPTH),EREP(MLVEXP,MDEPTH))
         ALLOCATE (CAPT(MLVEXP,MDEPTH),CAPN(MLVEXP,MDEPTH))
CC       ALLOCATE (CAPP(MLVEXP,MLVEXP,MDEPTH))
CC    END IF

      end subroutine alloc_alipar


!*********************************************************************

      subroutine alloc_array1

      use accura
      use basics
      use array1
      use modelq

      write(66,*) 'alloc_array1 MTOT,MFREX,MLEVEL,MDEPTH      ',
     *                          MTOT,MFREX,MLEVEL,MDEPTH

      ALLOCATE (A(MTOT,MTOT),  B(MTOT,MTOT),  C(MTOT,MTOT))
      ALLOCATE (E(MTOT,MTOT))
      ALLOCATE (VECL(MTOT),    Y1(MTOT),     Y2(MTOT))
      ALLOCATE (ALF(MTOT,MTOT),BET(MTOT,MDEPTH),DPSI(MTOT))
      ALLOCATE (PSI0(MTOT),    PSIM(MTOT),   PSIP(MTOT))
      ALLOCATE (RAD0(MTOT),    RADM(MTOT),   RADP(MTOT))
      ALLOCATE (FKM(MFREX),    FK0(MFREX),   FKP(MFREX))
      ALLOCATE (ABSOM(MFREX),  ABSO0(MFREX), ABSOP(MFREX))
      ALLOCATE (EMISM(MFREX),  EMIS0(MFREX), EMISP(MFREX))
      ALLOCATE (SCATM(MFREX),  SCAT0(MFREX), SCATP(MFREX))
      ALLOCATE (DABTM(MFREX),  DABT0(MFREX), DABTP(MFREX))
      ALLOCATE (DEMTM(MFREX),  DEMT0(MFREX), DEMTP(MFREX))
      ALLOCATE (DABNM(MFREX),  DABN0(MFREX), DABNP(MFREX))
      ALLOCATE (DEMNM(MFREX),  DEMN0(MFREX), DEMNP(MFREX))
      ALLOCATE (DABMM(MFREX),  DABM0(MFREX), DABMP(MFREX))
      ALLOCATE (DEMMM(MFREX),  DEMM0(MFREX), DEMMP(MFREX))
      ALLOCATE (WDEPM(MFREX),  WDEP0(MFREX), WDEPP(MFREX))
      ALLOCATE (SBFM(MLEVEL),  SBF0(MLEVEL), SBFP(MLEVEL))
      ALLOCATE (HEX(MLEVEL),   REX(MLEVEL),  REXA(MLEVEL))
      ALLOCATE (DSBFM(MLEVEL),DSBF0(MLEVEL),DSBFP(MLEVEL))
      ALLOCATE (SUMDCH(MLEVEL))
      ALLOCATE (DRCHM(MLEVEL,MFREX),   DRETM(MLEVEL,MFREX))
      ALLOCATE (DRCH0(MLEVEL,MFREX),   DRET0(MLEVEL,MFREX))
      ALLOCATE (DRCHP(MLEVEL,MFREX),   DRETP(MLEVEL,MFREX))

!     ALLOCATE (ABSOEX(MFREX,MDEPTH))
      ALLOCATE (EMISEX(MFREX,MDEPTH),SCATEX(MFREX,MDEPTH))
      ALLOCATE (DABTEX(MFREX,MDEPTH),DEMTEX(MFREX,MDEPTH))
      ALLOCATE (DABNEX(MFREX,MDEPTH),DEMNEX(MFREX,MDEPTH))
      ALLOCATE (DABMEX(MFREX,MDEPTH),DEMMEX(MFREX,MDEPTH))
      ALLOCATE (DRCHEX(MLVEXP,MFREX,MDEPTH))
      ALLOCATE (DRETEX(MLVEXP,MFREX,MDEPTH))

      end subroutine alloc_array1


!*********************************************************************


      subroutine alloc_rybpar

      use accura
      use basics
      use rybpar

!     write(66,*) 'alloc_rybpar MLEVEL           ',
!    *                          MLEVEL

      ALLOCATE (RA(MDEPTH),RB(MDEPTH),RC(MDEPTH),VR(MDEPTH),
     *          UA(MDEPTH),UB(MDEPTH),UC(MDEPTH),
     *          VA(MDEPTH),VB(MDEPTH),VC(MDEPTH),WR(MDEPTH),
     *          WM(MDEPTH,MDEPTH))

      ALLOCATE (ABSOPP(MFREQ,MDEPTH),SCATPP(MFREQ,MDEPTH),
     *          EMISPP(MFREQ,MDEPTH),BFABS(MLEVEL,MDEPTH))

      ALLOCATE (CSND(MDEPTH),PRAD2D(MDEPTH))

      END SUBROUTINE ALLOC_RYBPAR


!*********************************************************************


      subroutine alloc_topbase

      use accura
      use topbase

      write(66,*) ' alloc_topbase MOP,MMAXOP          ',
     &                            MOP,MMAXOP        

      ALLOCATE (XFIT(MOP) ,! local array containing x     for OP data
     +          SFIT(MOP)) ! local array containing sigma for OP data

      ALLOCATE (SOP(MOP,MMAXOP) ,! sigma = alog10(sigma/10^-18) of fit point
     +          XOP(MOP,MMAXOP)) ! x = alog10(nu/nu0) of fit point
      ALLOCATE (NOP(MMAXOP))     ! number of fit points for current level

      END SUBROUTINE ALLOC_TOPBASE

!-----------------------------------------------------------------------



! in iroset

!     ALLOCATE (EMKU(MLEVEL,2),YMKU(MLEVEL,2),EU(2*MLEVEL))
!     ALLOCATE (EEV(MKULEV),AEV(MKULEV),SEV(MKULEV),WEV(MKULEV))
!     ALLOCATE (EOD(MKULEV),AOD(MKULEV),SOD(MKULEV),WOD(MKULEV))
!     ALLOCATE (KSEV(MKULEV),KSOD(MKULEV),JEN(2*MLEVEL))
!     ALLOCATE (XJID(MDEPTH),JIDI(MDEPTH),JIDR(MDODF))
!     ALLOCATE (OMES(100,100),EKU(2*MKULEV),GKU(2*MKULEV))
!     ALLOCATE (KKU(2*MKULEV))

! in irost0

!     ALLOCATE (SIGFE(MDODF,MCFE))

! in srtfrq

!     ALLOCATE (ITRLIN(MITJ,MFREQ))

! in  rdatax

!     allocate (iex(mtrx),itrind(mtrx),izx0(mtrx),izx1(mtrx),
!    *          nmaxx(mtrx),izx(mtrx),nshx(mtrx),nax(mtrx),icx(mtrx))
!     allocate (etx(mtrx),ssx(mtrx),dx(mtrx))
!     allocate (aphx(11,5,mtrx),bphx(5,mtrx))


! in dkini  

CC    allocate (tabden(mtab,mden),tabtem(mtab,mtem,mden),
!    *          tablam(mtab,mlam,mden),prftab(mtab,mlam,mtem,mden))
!     allocate (nden(mtab),numtem(mtab,mden),numlam(mtab,mden))

! in sigave

!     ALLOCATE (FRINSG(MFREQ),CRIN(MFREQ))
!     ALLOCATE (JKF(MFREQ))

! in topbas

!     ALLOCATE (XFIT(MOP) ,! local array containing x     for OP data
!    +          SFIT(MOP)) ! local array containing sigma for OP data

! in opdata

!     ALLOCATE (SOP(MOP,MMAXOP) ,! sigma = alog10(sigma/10^-18) of fit point
!    +          XOP(MOP,MMAXOP)) ! x = alog10(nu/nu0) of fit point
!     ALLOCATE (NOP(MMAXOP))     ! number of fit points for current level

! in opahst

!     ALLOCATE (M1FILE(NLMX,MHOD),M2FILE(NLMX,MHOD))

! in iterini

!     ALLOCATE (LSNG(MTOT))

! in irrad

!     ALLOCATE (EXTRAD(MFREQ),EXTINT(MFREQ,MMU),HEXTRD(MFREQ))
!     ALLOCATE (FST(MFREQ),EXTRD0(MFREQ))
!     DEALLOCATE (FST,EXTRD0)

C
C     ***************************************************
c

!     in prdini

!     ALLOCATE (DOPTR(MTRPRD,MDEPTH),COHER(MTRPRD,MDEPTH),
!    *         PJBAR(MTRPRD,MDEPTH),RJBAR(MTRPRD,MDEPTH))
!     ALLOCATE (IPRD(MTRANS),ITRTOT(MTRPRD))

