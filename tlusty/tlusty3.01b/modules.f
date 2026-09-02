      module accura
!!!   integer,parameter::dp=selected_real_kind(15,300)
      integer,parameter::dp=selected_real_kind(15,307)
      end module accura

!     *******************************************************
!
      module basics
!     module replacement for BASICS.FOR

      use accura
CC    implicit none
!
!     Parameters that specify dimensions of arrays
!
      INTEGER, PARAMETER :: MATOM  =   99  ! max.num. of explicit atoms
      INTEGER            :: MION
      INTEGER            :: MLEVEL
      INTEGER            :: MLVEXP
      INTEGER            :: MTRAN0   
      INTEGER, PARAMETER :: MDEPT0 =  100  ! max.num. of depth points
      INTEGER            :: MTRANS         ! max.num. of transitions
      INTEGER            :: MDEPTH         ! max num. of depth points
      INTEGER            :: MFREQ0         ! initial .num. of frequency points
      INTEGER            :: MFREQ          ! max. num. of frequency points
      INTEGER            :: MFREQP         ! working arrays of frequency
      INTEGER            :: MFREQC         ! max num. of cont. frequency pts.
      INTEGER            :: MFREX          ! max num. of explicit freqs.
      INTEGER            :: MFREQL         ! max.num. of frequencies per line
      INTEGER            :: MTOT           ! max.num. of linearized parameters
      INTEGER, PARAMETER :: MMU    =    6  ! max.num. of angle points
      INTEGER, PARAMETER :: MFIT   =  357  ! max.num. of fit points (OP b-f c.s)
      INTEGER            :: MITJ            ! max.num. of overlapping transitions
      INTEGER, PARAMETER :: MMCDW  =   46  ! max.num. of levels with pseudocont.
      INTEGER, PARAMETER :: MMER   =   12  ! max.num. of merged levels
      INTEGER, PARAMETER :: MVOIGT = 8080  ! max.num. of lines with Voigt profile
      INTEGER            :: MVOIG0
      INTEGER, PARAMETER :: MZZ    =   30  ! maximum charge for occup.prob. ions
      INTEGER, PARAMETER :: NLMX   =   80  ! highest hydrogenic level considered
      INTEGER, PARAMETER :: MSMX   =    1  ! size of matrix kept in memory in SOLVE
      INTEGER            :: MFREQ1
      INTEGER, PARAMETER :: MDEPTC =    2  ! max.num. of depth points (Compton)
      INTEGER, PARAMETER :: MMUC   =    2  ! max.num. of angle points (Compton)
      INTEGER            :: MLEVE3         ! =1 for diag.prec; =MLEVEL for trid.prec.
      INTEGER            :: MLEVEX3        ! =1 for diag.oper.; =MLVEXP for tridiag.
      INTEGER            :: MTRAN3         ! =1 for diag.oper; =MTRANS for tridiag.
      INTEGER            :: MCROSS         ! max.num. of b-f cross.secs.
      INTEGER            :: MBF
      INTEGER, PARAMETER :: MCFIT  =   10  ! max.num. of collision fit points
      INTEGER, PARAMETER :: MXTCOL =    3  ! max.num. of collision types (CE, CP, CH)
      INTEGER            :: MCORAT         ! max.num. of col excitation transitions
!     
!     Basic physical constants
!
      REAL(DP), PARAMETER :: H     = 6.6256e-27_dp      ! Planck constant     h 
      REAL(DP), PARAMETER :: BOLK  = 1.38054e-16_dp     ! Boltzmann constant  k
      REAL(DP), PARAMETER :: HK    = 4.79928144e-11_dp  ! h/k
      REAL(DP), PARAMETER :: CAS   = 2.997925e18_dp     ! light speed c (A/s)
      REAL(DP), PARAMETER :: EH    = 2.17853041e-11_dp  ! ionizaton energy of hydrogen
      REAL(DP), PARAMETER :: BN    = 1.4743e-2_dp       ! 2*h/c**3, c -light speed
      REAL(DP), PARAMETER :: SIGE  = 6.6516e-25_dp      ! Thomson scattering c-s
      REAL(DP), PARAMETER :: SIG4P = 4.5114062e-6_dp    ! Stefan-Boltzmann const/4pi
      REAL(DP), PARAMETER :: PI4H  = 1.8966e27_dp       ! 4pi/h
      REAL(DP), PARAMETER :: PCK   = 4.19168946e-10_dp  ! 4pi/c
      REAL(DP), PARAMETER :: HMASS = 1.67333e-24_dp     ! mass of hydrogen atom
!
!     Basic mathematical constants
!
      REAL(DP), PARAMETER :: UN    = 1.0_dp
      REAL(DP), PARAMETER :: HALF  = 0.5_dp
      REAL(DP), PARAMETER :: TWO   = 2.0_dp
      REAL(DP), PARAMETER :: ZERO  = 0.0_dp
!
!     Unit number
!
      INTEGER, PARAMETER :: IBUFF = 95
!
!     from common/basnum
      INTEGER ::    NATOM,NION,NLEVEL,NTRANS,ND,NFREQ,NFREQC,NFREQE
      INTEGER ::    IOPTAB,IDISK,IZSCAL,IDMFIX,IHESO6,IFMOL,IFENTR
      INTEGER ::    NFREQL,NLEV0,ICOLHN,IOSCOR,ILGDER,IFRYB,IFRSET
      INTEGER ::    NFREAD,NELSC,NTRANC,IOVER,JALI,IBC,IUBC,INTENS
      INTEGER ::    IRDER,ILMCOR,IFDIEL,IFALI,IFTENE,ITNDRE,NTRX
      INTEGER ::    ILPSCT,ILASCT,IRTE,IDLTE,IBFINT,INTRPL,ICHANG
      INTEGER ::    NATOMS,IPSLTE,ISPODF,ITLUCY,NRETC,IFRAYL,IFPRAD
      INTEGER ::    IOPFR,IFRAL,ICHCI,IFPZEV,IOPFRC,IFALIH
      INTEGER ::    NLVEXP,NLVFOR,NLVEXZ,NVOIGT,NMER,NFHOD
!     from common/inppar
      REAL(DP) ::     TEFF,GRAV
      REAL(DP) ::     YTOT(MDEPT0),WMM(MDEPT0),WMY(MDEPT0)
      REAL(DP) ::     TMOLIM,EE0
      REAL(DP) ::     xmstar,xmdot,rstar,alpha0,reynum,adist
      REAL(DP) ::     QGRAV,EDISC,DZETA,RELDST
      REAL(DP) ::     visc,zeta0,zeta1,dmvisc,fractv
      REAL(DP) ::     omeg32,wbarm,wbar,alphav,pgas0
      REAL(DP) ::     bergfc,cutlym,cutbal,aneut
      INTEGER  ::    ISPLIN,IRSPLT,ivisc,ibche
      LOGICAL  ::    LTE,LTGREY,LCHC,LRESC,LASV
!     from common/matkey
      INTEGER ::    NN,NN0,INHE,INRE,INPC,INSE,INZD,INMP,NDRE,NDS,IFIXDE
      REAL(DP)::    XI2(NLMX),XI3(NLMX)
      REAL(DP)::    CHMAX
      INTEGER ::    ITER,NITER,NITZER,INIT,NLAMBD
      LOGICAL ::    LAC2,LFIN
!     from common/conkey
      REAL(DP)::     HMIX0,crflim
CC    INTEGER ::    NCONIT,ICONV,INDL,IPRESS,ITEMP,ICBEG,ICEND
      INTEGER ::    NCONIT,ICONV,INDL,IPRESS,ITEMP
      INTEGER ::    itmcor,iconre,ideepc,ndcgap,IDCONZ
!     from common/opckey
      INTEGER ::    NCON,IOPHL1,IOPHL2,IPHE2C,IFMOFF
!     from common/prints
      INTEGER ::    IPRINT,IPRING,IPRIND,IPRINP,ICOOLP,ICHCKP
      INTEGER ::    IPOPAC,IPRINI
!     from common/psilim
      REAL(DP)::     DPSILG,DPSILT,DPSILN,DPSILD
!     from common/centrl
      REAL(DP)::     ZND,IFZ0
!     from common/opadd   
      INTEGER ::    IOPADD,IOPHMI,IOPH2P,IOPHEM,IOPCH,IOPOH,IOPH2M
      INTEGER ::    IOH2H2,IOH2HE,IOH2H,IOHHE
      INTEGER ::    IOPHLI,IOPLYM
      INTEGER ::    IRSCT,IRSCHE,IRSCH2
      INTEGER ::    KEEPOP,IOPOLD
!     from comomon/angles
      REAL(DP)::    AMU(MMU),WTMU(MMU),FMU(MMU)
      INTEGER ::    NMU
!     from common/comptn, compti,comite
      REAL(DP)::     amuc(mmuc),wtmuc(mmuc)
      REAL(DP)::     amuc1(mmuc),amuc2(mmuc),amuc3(mmuc)
      REAL(DP)::     amuj(mmuc),amuk(mmuc),amuh(mmuc),amun(mmuc)
      REAL(DP)::     calph(mmuc,mmuc),cbeta(mmuc,mmuc)
      REAL(DP)::     cgamm(mmuc,mmuc),RADZER,FRLCOM
      INTEGER ::    nmuc
      INTEGER ::    nedd,nsti,islab,ilbc,icompt,icomst,icomde
      INTEGER ::    icombc,icmdra,knish,itcomp,icomve,icomrt
      INTEGER ::    ichcoo,icomgr
      INTEGER ::    ncfor1,ncfor2,nccoup,ncitot,ncfull
!     from common/mlcons
      REAL(DP)::     aconml,bconml,cconml
!     from common iprkey and others,ifratp,irbefa,iprat
      INTEGER ::    iprkey,ifratp,irbefa,iprat
      INTEGER ::    iprybh,ipelch,ipeldo,ipconf
!     from common/iunit
      INTEGER ::    IUNIT

!from commons in nstpar 
      INTEGER   ::    IIRWIN,MOLTAB,IRWTAB,IOVERR
      INTEGER   ::    IHYDDK,IGCOMB,IQUASI,NUNALP,NUNBET,NUNGAM,NUNBAL
      INTEGER   ::    ICONRS,NCONRF,IMUCON,ICHANM
      INTEGER   ::    IPRCRS,NPRCRS,IPTRAN
      INTEGER   ::    iflskp,intrho,nodiss,ipzev
      INTEGER   ::    NFRECL,ITGMAX,ITGMX0
      INTEGER   ::    NLASTE,NHOD
      INTEGER   ::    ifcaut,ifcrea
      INTEGER   ::    irgrad,iprybc,itgrad
      REAL(DP)  ::    ARH,BRH,CRH,DRH

!     REAL(DP)  ::    T4,PGAS,PRAD,PGM,PRADM,PGSGRA,Q0MAX
      REAL(DP)  ::    T4,PGSGRA,Q0MAX

      REAL(DP)  ::    frmin,frmax
      REAL(DP)  ::    tqmprf
      REAL(DP)  ::    hcmass,radstr
      REAL(DP)  ::    tcoaut
      REAL(DP)  ::    TSTD, VTB, TFLOOR
      REAL(DP)  ::    GRDADB,GRDAD0,DIFT,DIFP,DERT
!    allocatable arrays
      REAL(DP), ALLOCATABLE  ::    SIGEC(:)
      INTEGER, ALLOCATABLE ::    IJORIG(:)
      end module basics

!*************************************************************************

      module alipar
!     module replacement for ALIPAR.FOR
      use accura
c     usa basics

      REAL(DP), ALLOCATABLE :: REIT(:),REIN(:),REIM(:)
      REAL(DP), ALLOCATABLE :: REIP(:,:)
      REAL(DP), ALLOCATABLE :: AREIT(:),AREIN(:),AREIM(:)
      REAL(DP), ALLOCATABLE :: AREIP(:,:)
      REAL(DP), ALLOCATABLE :: CREIT(:),CREIN(:),CREIM(:)
      REAL(DP), ALLOCATABLE :: CREIP(:,:)
      REAL(DP), ALLOCATABLE :: REIX(:),CREIX(:)
      REAL(DP), ALLOCATABLE :: REDX(:),REDT(:),REDN(:)
      REAL(DP), ALLOCATABLE :: REDM(:),REDP(:,:)
      REAL(DP), ALLOCATABLE :: REDXM(:),REDTM(:),REDNM(:)
      REAL(DP), ALLOCATABLE :: REDMM(:),REDPM(:,:)
      REAL(DP), ALLOCATABLE :: REDTP(:),REDNP(:),REDXP(:)
      REAL(DP), ALLOCATABLE :: REDMP(:),REDPP(:,:)
      REAL(DP), ALLOCATABLE :: HEIT(:),HEIN(:),HEIM(:)
      REAL(DP), ALLOCATABLE :: HEIP(:,:)
      REAL(DP), ALLOCATABLE :: HEITM(:),HEINM(:),HEIMM(:)
      REAL(DP), ALLOCATABLE :: HEIPM(:,:)
      REAL(DP), ALLOCATABLE :: HEITP(:),HEINP(:),HEIMP(:)
      REAL(DP), ALLOCATABLE :: HEIPP(:,:)
      REAL(DP), ALLOCATABLE :: EHET(:),EHEN(:),ERET(:),EREN(:)
      REAL(DP), ALLOCATABLE :: EHEP(:,:),EREP(:,:)
      REAL(DP), ALLOCATABLE :: APT(:,:),APN(:,:)
      REAL(DP), ALLOCATABLE :: APM(:,:)
      REAL(DP), ALLOCATABLE :: AAPT(:,:),AAPN(:,:)
      REAL(DP), ALLOCATABLE :: AAPM(:,:)
      REAL(DP), ALLOCATABLE :: CAPT(:,:),CAPN(:,:)
      REAL(DP), ALLOCATABLE :: APP(:,:,:)
      REAL(DP), ALLOCATABLE :: AAPP(:,:,:)
      REAL(DP), ALLOCATABLE :: CAPP(:,:,:)
      REAL(DP)              :: QTLAS
      INTEGER             :: IFPOPR,irprec,ifprec,itold1,itold2
      INTEGER             :: itlas
      end module alipar
 
!*************************************************************************

      module array1
!     module replacement for ARRAY1.FOR

      use accura
CCC   implicit none

      REAL(DP), ALLOCATABLE :: A(:,:),    B(:,:),    C(:,:) 
      REAL(DP), ALLOCATABLE :: E(:,:) 
      REAL(DP), ALLOCATABLE :: VECL(:),   Y1(:),     Y2(:) 
      REAL(DP), ALLOCATABLE :: ALF(:,:),  BET(:,:),  DPSI(:)
      REAL(DP), ALLOCATABLE :: PSI0(:),   PSIM(:),   PSIP(:) 
      REAL(DP), ALLOCATABLE :: RAD0(:),   RADM(:),   RADP(:) 
      REAL(DP), ALLOCATABLE :: FKM(:),    FK0(:),    FKP(:) 
      REAL(DP), ALLOCATABLE :: ABSOM(:),  ABSO0(:),  ABSOP(:) 
      REAL(DP), ALLOCATABLE :: EMISM(:),  EMIS0(:),  EMISP(:) 
      REAL(DP), ALLOCATABLE :: SCATM(:),  SCAT0(:),  SCATP(:) 
      REAL(DP), ALLOCATABLE :: DABTM(:),  DABT0(:),  DABTP(:) 
      REAL(DP), ALLOCATABLE :: DEMTM(:),  DEMT0(:),  DEMTP(:) 
      REAL(DP), ALLOCATABLE :: DABNM(:),  DABN0(:),  DABNP(:) 
      REAL(DP), ALLOCATABLE :: DEMNM(:),  DEMN0(:),  DEMNP(:) 
      REAL(DP), ALLOCATABLE :: DABMM(:),  DABM0(:),  DABMP(:) 
      REAL(DP), ALLOCATABLE :: DEMMM(:),  DEMM0(:),  DEMMP(:) 
      REAL(DP), ALLOCATABLE :: WDEPM(:),  WDEP0(:),  WDEPP(:) 
      REAL(DP), ALLOCATABLE :: SBFM(:),  SBF0(:),  SBFP(:) 
      REAL(DP), ALLOCATABLE :: HEX(:),   REX(:),   REXA(:) 
      REAL(DP), ALLOCATABLE :: DSBFM(:), DSBF0(:), DSBFP(:) 
      REAL(DP), ALLOCATABLE :: SUMDCH(:) 
      REAL(DP), ALLOCATABLE :: DRCHM(:,:),   DRETM(:,:) 
      REAL(DP), ALLOCATABLE :: DRCH0(:,:),   DRET0(:,:) 
      REAL(DP), ALLOCATABLE :: DRCHP(:,:),   DRETP(:,:)
!     from common/exprad
CC    REAL(DP), ALLOCATABLE :: ABSOEX(:,:),EMISEX(:,:) 
CC    REAL(DP), ALLOCATABLE :: SCATEX(:,:) 
      REAL(DP), ALLOCATABLE :: DABTEX(:,:),DEMTEX(:,:) 
      REAL(DP), ALLOCATABLE :: DABNEX(:,:),DEMNEX(:,:) 
      REAL(DP), ALLOCATABLE :: DABMEX(:,:),DEMMEX(:,:) 
      REAL(DP), ALLOCATABLE :: DRCHEX(:,:,:) 
      REAL(DP), ALLOCATABLE :: DRETEX(:,:,:)

      REAL(DP)               :: CZZ,CZN,CZE,CZM
!     from common/bpocom
CC    REAL(DP)              :: BESE(:),ATT(:),ANN(:)

      end module array1

!*************************************************************************

      module atomic
!     module replacement for ATOMIC.FOR

      use accura
      use basics, only: MATOM, MION, MLEVEL, MDEPTH
CCC   implicit none

      INTEGER, PARAMETER :: MIOI = 170
      INTEGER, PARAMETER :: MLEVEI  = 1200
      CHARACTER(LEN=40)  :: FIDATA(MIOI),FIODF1(MIOI),FIODF2(MIOI)
      CHARACTER(LEN=40)  :: FIBFCS(MIOI)
      CHARACTER(LEN=10)  :: TYPLEV(MLEVEI)
      CHARACTER(LEN=4)   :: TYPION(MIOI)

!     from common/atopar 
      REAL(DP), ALLOCATABLE :: AMASS(:),ABUND(:,:),ABNDD(:,:),ENEV(:,:)
      REAL(DP)              :: AMAS(100)
      INTEGER,ALLOCATABLE   :: NUMAT(:)
      INTEGER,ALLOCATABLE   :: N0A(:),NKA(:),nref(:),iatex(:)
      INTEGER,ALLOCATABLE   :: nrefs(:,:),iadop(:)
      INTEGER,ALLOCATABLE   :: iifix(:)
      INTEGER               :: iatref,modref
      INTEGER               :: IONIZ(MATOM),MODPF(MATOM)
      LOGICAL,ALLOCATABLE   :: LRM(:),LGR(:)

!     from common/ionpar
      REAL(DP),ALLOCATABLE  :: FF(:),CHARG2(:)
      INTEGER,ALLOCATABLE   :: NFIRST(:),NLAST(:),NNEXT(:)
      INTEGER,ALLOCATABLE   :: IZ(:),IUPSUM(:),ICUP(:),ILTE(:)
      INTEGER,ALLOCATABLE   :: ILTION(:),IKOBS(:),MODFF(:)
      INTEGER,ALLOCATABLE   :: IFMETA(:)
      INTEGER,ALLOCATABLE   :: INODF1(:),INODF2(:),INBFCS(:)
      INTEGER,ALLOCATABLE   :: IMRG(:),IIMER(:)
      INTEGER,ALLOCATABLE   :: NEVKU(:),NODKU(:)
      REAL(DP), ALLOCATABLE :: XEV(:,:),XOD(:,:),EU(:)
      INTEGER               :: NLEVKU,NLINKU,KEVE,KODD
      INTEGER               :: N0HN

!     from common/levpar
      REAL(DP), ALLOCATABLE :: ENION(:),G(:)
      INTEGER,ALLOCATABLE   :: NQUANT(:)
      INTEGER,ALLOCATABLE   :: IATM(:),IEL(:),ILK(:),ilin(:)
      INTEGER,ALLOCATABLE   :: iltlev(:),indlev(:)
      INTEGER,ALLOCATABLE   :: imodl(:),iiexp(:),iifor(:),imodl0(:)
      INTEGER,ALLOCATABLE   :: ipzert(:)
      INTEGER,ALLOCATABLE   :: igzert(:),indlgz(:),iinonz(:)
      INTEGER,ALLOCATABLE   :: IFWOP(:)
      LOGICAL               :: LBPFX

!     from common/trapar
      REAL(DP), ALLOCATABLE :: FR0(:),OSC0(:),CPAR(:)
      REAL(DP), ALLOCATABLE :: FRQMX(:),FR0PC(:),OMECOL(:,:)
      REAL(DP)              :: XGRAD,STRL1,STRL2,STRLX
      INTEGER,ALLOCATABLE :: ILOW(:),IUP(:),INDEXP(:)
      INTEGER,ALLOCATABLE :: KFR0(:),KFR1(:),ILUCTR(:)
      INTEGER,ALLOCATABLE :: IFC0(:),IFC1(:)
      INTEGER,ALLOCATABLE :: IFR0(:),IFR1(:),ITRA(:,:)
      INTEGER,ALLOCATABLE :: IPROF(:),ICOL(:),INTMOD(:)
      INTEGER,ALLOCATABLE :: IPROF0(:),ITRCON(:),IDIEL(:)
      INTEGER,ALLOCATABLE :: IJTF(:)
      INTEGER,ALLOCATABLE :: ITRBF(:)
      INTEGER,ALLOCATABLE :: ILTREF(:,:),IGUIDE(:)
      LOGICAL,ALLOCATABLE :: LINEXP(:)
      LOGICAL,ALLOCATABLE :: LCOMP(:),LINE(:)

!     from common/phoset
      REAL(DP), ALLOCATABLE :: S0CS(:),ALFCS(:),BETCS(:)
      REAL(DP), ALLOCATABLE :: GAMCS(:)  
      INTEGER,ALLOCATABLE :: IBF(:)

      REAL(DP), ALLOCATABLE :: CTOP(:,:),XTOP(:,:)

      REAL(DP), ALLOCATABLE :: CTEMP(:,:,:),CRATE(:,:,:)

      REAL(DP), ALLOCATABLE :: GAMAR(:),STARK1(:),STARK2(:),STARK3(:)
      REAL(DP), ALLOCATABLE :: VDWH(:)

      REAL(DP)              :: COLHE1(19,19)

      LOGICAL,ALLOCATABLE :: LEXP(:),LALI(:)

      INTEGER             :: NFFIX,IFSUB,IFLEV
      INTEGER             :: IATH,IATHE,IELH,IELHM,IELHE1,IELHE2
      INTEGER,ALLOCATABLE :: IATI(:),IZI(:),NLEVS(:),NLLIM(:)
    
      REAL(DP)            :: OSH(20,20)
      INTEGER             :: IOSH,JOSH
      DATA ((OSH(IOSH,JOSH),IOSH=1,20),JOSH=1,16)/20*0.,
     * 0.4162,19*0.,7.910E-2,0.6407,18*0.,2.899E-2,0.1193,
     * 0.8421,17*0.,1.394E-2,4.467E-2,0.1506,1.038,16*0.,7.799E-3,
     * 2.209E-2,5.584E-2,0.1793,1.231,15*0.,4.814E-3,1.270E-2,2.768E-2,
     * 6.549E-2,0.2069,1.424,14*0.,3.183E-3,8.036E-3,1.604E-2,3.23E-2,
     * 7.448E-2,0.234,1.616,13*0.,2.216E-3,5.429E-3,1.023E-2,1.87E-2,
     * 3.645E-2,8.315E-2,0.2609,1.807,12*0.,1.605E-3,3.851E-3,6.98E-3,
     * 1.196E-2,2.104E-2,4.038E-2,9.163E-2,0.2876,1.999,11*0.,1.201E-3,
     * 2.835E-3,4.996E-3,8.187E-3,1.344E-2,2.32E-2,4.416E-2,0.1,0.3143,
     * 2.19,10*0.,9.214E-4,2.151E-3,3.711E-3,5.886E-3,9.209E-3,1.479E-2,
     * 2.525E-2,4.787E-2,0.1083,0.3408,2.381,9*0.,7.227E-4,1.672E-3,
     * 2.839E-3,4.393E-3,6.631E-3,1.012E-2,1.605E-2,2.724E-2,5.152E-2,
     * 0.1166,0.3673,2.572,8*0.,5.744E-4,1.326E-3,2.224E-3,3.375E-3,
     * 4.959E-3,7.289E-3,1.097E-2,1.726E-2,2.918E-2,5.513E-2,0.1248,
     * 0.3938,2.763,7*0.,4.686E-4,1.07E-3,1.776E-3,2.656E-3,3.821E-3,
     * 5.455E-3,7.891E-3,1.177E-2,1.843E-2,3.109E-2,5.872E-2,0.133,
     * 0.4202,2.954,6*0.,3.856E-4,8.764E-4,1.443E-3,2.131E-3,3.014E-3,
     * 4.207E-3,5.905E-3,8.456E-3,1.254E-2,1.958E-2,3.298E-2,6.228E-2,
     * 0.1412,0.4467,3.145,5*0./
      DATA ((OSH(IOSH,JOSH),IOSH=1,20),JOSH=17,20)/3.211E-4,
     * 7.270E-4,1.188E-3,1.739E-3,
     * 2.425E-3,3.324E-3,4.556E-3,6.323E-3,8.995E-3,.01328,.0207,.03486,
     * .06584,.1494,0.4731,3.336,4*0.,2.702E-4,6.099E-4,9.916E-4,
     * 1.439E-3,1.984E-3,2.679E-3,3.602E-3,4.877E-3,6.719E-3,9.515E-3,
     * 0.01402,.02182,.03672,.06938,.1575,.4995,3.527,3*0.,2.296E-4,
     * 5.167E-4,8.361E-4,1.204E-3,1.646E-3,2.196E-3,2.905E-3,3.856E-3,
     * 5.180E-3,7.099E-3,.01002,.01474,.02292,.03858,.07292,.1657,.5259,
     * 3.718,2*0.,1.967E-4,4.416E-4,7.118E-4,1.019E-3,1.382E-3,1.825E-3,
     * 2.383E-3,3.112E-3,4.094E-3,5.468E-3,7.468E-3,.01052,.01545,
     * .02402,.04043,.07644,0.1738,.5523,3.909,0./


      REAL(DP)              :: FRTABM

!     from common/prdpar
      INTEGER,PARAMETER     :: MTRPRD=5
      REAL(DP), ALLOCATABLE :: DOPTR(:,:),COHER(:,:)
      REAL(DP), ALLOCATABLE :: PJBAR(:,:),RJBAR(:,:)
      REAL(DP)              :: XPDIV
      INTEGER,ALLOCATABLE   :: IPRD(:),ITRTOT(:)
      INTEGER               :: NTRPRD,IFPRD
 
      INTEGER,ALLOCATABLE   :: IJFL(:),IJALI(:)
 
      end module atomic

!*************************************************************************

      module odfpar
!     module replacementfor ODFPAR.FOR
      
      use accura
CCC   implicit none

      INTEGER, PARAMETER    :: MDODF =        3
      INTEGER, PARAMETER    :: MFODF =      180
      INTEGER, PARAMETER    :: MKULEV=     7000
      INTEGER, PARAMETER    :: MLINE =  1140000
c     INTEGER, PARAMETER    :: MCFE  = 34824000
      INTEGER               :: MCFE
 
      REAL, ALLOCATABLE     :: SIGFE(:,:)
      REAL, ALLOCATABLE     :: AGAF(:,:),VDOP(:,:),SIG0(:,:),WAVE(:)

      REAL(DP)              :: FRS1,FRS2,FRS3,DXNU,GST
      REAL(DP), ALLOCATABLE :: XJID(:)
      INTEGER,ALLOCATABLE   :: JIDI(:)
      INTEGER,ALLOCATABLE   :: JIDR(:)
      INTEGER               :: JIDS,JIDN,NFRS1,NFTT
      INTEGER,ALLOCATABLE   :: JEN(:),JTR(:,:)
      REAL(DP), ALLOCATABLE :: EMKU(:,:),YMKU(:,:)
      REAL(DP),ALLOCATABLE  :: EEV(:),AEV(:),SEV(:),WEV(:)
      REAL(DP),ALLOCATABLE  :: EOD(:),AOD(:),SOD(:),WOD(:)
      INTEGER ,ALLOCATABLE  :: KSEV(:),KSOD(:)
      REAL(DP),ALLOCATABLE  :: OMES(:,:)
      REAL(DP), ALLOCATABLE :: EKU(:),GKU(:)
      INTEGER,  ALLOCATABLE :: KKU(:)

      REAL(DP), ALLOCATABLE :: OFR(:),OW(:),OWSUB(:)
      REAL(DP), ALLOCATABLE :: ODFL0(:,:),ODF2(:)
      INTEGER, ALLOCATABLE  :: IFTRA(:),IDODF(:)
      INTEGER               :: NDODF       

      end module odfpar


!*************************************************************************


      module hydodf

      use accura
      use basics, only: NLMX
CCC   implicit none

      INTEGER, PARAMETER :: MHOD  =       3
      INTEGER, PARAMETER :: MFRO  =     5000    

      REAL(DP), ALLOCATABLE :: FRODF(:)
      INTEGER               :: NFRODF(MHOD)
      INTEGER,ALLOCATABLE   :: INDODF(:),JNDODF(:)
      REAL(DP), ALLOCATABLE :: FROS(:,:),WNUS(:,:),FFRO(:)
      REAL(DP)              :: XDO(3,MHOD)
      INTEGER               :: KDO(4,MHOD)
      INTEGER,ALLOCATABLE   :: I1ODF(:),I2ODF(:),NQLODF(:)
      REAL(DP)              :: XKIJ(MHOD,NLMX),WL0(MHOD,NLMX)
      REAL(DP)              :: FIJ(MHOD,NLMX)
      INTEGER,ALLOCATABLE   :: M1FILE(:,:),M2FILE(:,:)
      INTEGER               :: IMERG
      REAL(DP)              :: ALLIM1,ABLIM1,ABLIM2,ABLIM3

      end module hydodf


!*************************************************************************


      module iterat
!     module replacement for ITERAT.FOR
      use accura
!     use basics, only: MTOT

!     implicit none

!     INTEGER, PARAMETER :: MITER  = 100
!     INTEGER, PARAMETER :: MLAMBD =  50

      INTEGER, ALLOCATABLE :: NITLAM(:)
      INTEGER, ALLOCATABLE :: IFFIX(:),NETEXP(:),NETFIX(:)
      INTEGER, ALLOCATABLE :: IETEXP(:,:),IETFIX(:,:)
      INTEGER, ALLOCATABLE :: INHE0(:),INRE0(:),INPC0(:)
      INTEGER, ALLOCATABLE :: INDL0(:),INSE0(:),INMP0(:)
      INTEGER, ALLOCATABLE :: NN00(:),NDRE0(:),KANT(:)
      INTEGER              :: ITEK,IACC,IACC0,IACD,KSNG
      INTEGER              :: ILAM,IACPP,IACC0P,IACDP
      INTEGER              :: IACLT,IACLDT
      INTEGER              :: NLAMT,ILDER,IBPOPE,IELCOR
      REAL(DP)             :: CHMAXT,ORELAX   
CCC   LOGICAL              :: LSNG(MTOT)
      INTEGER,ALLOCATABLE  :: ISNG(:)
      LOGICAL ::              LASO,LRES2,LCHMAT,LIROST,LAC2P
     
      end module iterat

!*************************************************************************

      module hyglin
      use accura

      integer, parameter :: mlinn     =  12
      integer, parameter :: mtabtn    =   5
      integer, parameter :: mtaben    =  13
      integer, parameter :: mfhtab    =1500
      integer, parameter :: mlilor    =  30
      integer, parameter :: mtelor    =  15

      INTEGER, PARAMETER :: MLINH  = 78,
     *                       MHT    = 7,
     *                       MHE    = 20,
     *                       MHWL   = 90


      REAL(DP)               :: hglin
      INTEGER                :: ihgon
      REAL(DP)               :: ENTAB1(mlinn), ENTAB2(mlinn)
      REAL(DP)               :: TNTAB1(mlinn),TNTAB2(mlinn)
      INTEGER                :: nunfreq(mlinn)
      INTEGER                :: nunele(mlinn),nuntem(mlinn)
      REAL(DP)               :: temn(mtabtn,mlinn),elen(mtaben,mlinn)
      REAL(DP), ALLOCATABLE  :: hygcrs(:,:,:,:)
      REAL(DP), ALLOCATABLE  :: hygcr0(:,:,:,:)
      REAL(DP), ALLOCATABLE  :: ynint(:,:)
      INTEGER, ALLOCATABLE   :: jnint(:,:)

!     original Gomez hydrogen line profiles

      REAL(DP)              :: HGLIM
      INTEGER               :: IHGOM
      REAL(DP)              :: FRGTB1,FRGTB2,EGTAB1,EGTAB2,TGTAB1,TGTAB2
      INTEGER               :: nugfreq,nugele,nugtemp
      REAL(DP), ALLOCATABLE :: temvec(:),elevec(:)
      REAL(DP), ALLOCATABLE :: frgtab(:),hydcrs(:,:,:)
      REAL(DP), ALLOCATABLE :: YGINT(:)
      INTEGER, ALLOCATABLE  :: JGINT(:)

      REAL(DP), ALLOCATABLE  :: temm(:),gamm(:,:),shif(:,:)
      INTEGER                :: NTEMM

!     Lemke or Tremblay tables

      REAL(DP), ALLOCATABLE :: PRFHYD(:,:,:,:)
      REAL(DP), ALLOCATABLE :: WLHYD(:,:)
      REAL(DP), ALLOCATABLE :: WLH(:,:)
      REAL(DP), ALLOCATABLE :: XTLEM(:,:)
      REAL(DP), ALLOCATABLE :: XNELEM(:,:)
      REAL(DP), ALLOCATABLE :: XK0(:)
      INTEGER,ALLOCATABLE   :: NWLHYD(:)
      INTEGER,ALLOCATABLE   :: NWLH(:)
      INTEGER,ALLOCATABLE   :: NTH(:)
      INTEGER,ALLOCATABLE   :: NEH(:)
      INTEGER               :: ILINH(4,22)
      INTEGER               :: IHYDPR

!     Gomze Xenomorph tables

      REAL(DP), ALLOCATABLE :: PRFXB(:,:,:,:)
      REAL(DP), ALLOCATABLE :: PRFXR(:,:,:,:)
      REAL(DP), ALLOCATABLE :: ALXEN(:,:)
      REAL(DP), ALLOCATABLE :: XTXEN(:,:)
      REAL(DP), ALLOCATABLE :: XNEXEN(:,:)
      REAL(DP)              :: XNEMIN
      INTEGER,ALLOCATABLE   :: NWLXEN(:)
      INTEGER,ALLOCATABLE   :: NTHXEN(:)
      INTEGER,ALLOCATABLE   :: NEHXEN(:)
      INTEGER               :: ILXEN(4,22)
      INTEGER               :: IHXENB



      end module hyglin

!*************************************************************************
      
      module dktabl

      use accura

      integer,parameter     :: mtab=5,
     *                         mden=10,mtem=10,mlam=1000
      real(dp), allocatable :: tabden(:,:),tabtem(:,:,:)
      real(dp), allocatable :: tablam(:,:,:),prftab(:,:,:,:)
      integer, allocatable  :: nden(:),numtem(:,:),numlam(:,:)

      end module dktabl


!*************************************************************************

      module modelq

      use accura

      REAL(DP), ALLOCATABLE :: DM(:),TEMP(:),ELEC(:),DENS(:)
      REAL(DP), ALLOCATABLE :: TOTN(:)
      REAL(DP), ALLOCATABLE :: ANTO(:),ANMA(:),ANH1(:),ZD(:)
      REAL(DP), ALLOCATABLE :: HKT1(:),TK1(:),HKT21(:),SQT1(:)
      REAL(DP), ALLOCATABLE :: TEMP1(:),ELEC1(:),DENS1(:)
      REAL(DP), ALLOCATABLE :: DENSI(:),DENSIM(:)
      REAL(DP), ALLOCATABLE :: ELSCAT(:),ALAB(:)
      REAL(DP), ALLOCATABLE :: DELDM(:),DELDMZ(:)
      REAL(DP), ALLOCATABLE :: THETAV(:),VISCD(:)
      REAL(DP)              :: DEDM1,DALPMX,XHYD
      REAL(DP)              :: DMTOT,RRDIL,TEMPBD,ALPTAV,ALPGAV
      INTEGER               :: NALP,IBETA

      REAL(DP), ALLOCATABLE :: POPUL(:,:),BFAC(:,:)
      REAL(DP), ALLOCATABLE :: POPINV(:,:),POPUL0(:,:)
      REAL(DP), ALLOCATABLE :: POPGRP(:)
      REAL(DP), ALLOCATABLE :: POP(:),SBF(:),DSBF(:),USUM(:)
      REAL(DP), ALLOCATABLE :: ANH2(:),ANHM(:)

      REAL(DP), ALLOCATABLE :: RPOP0(:,:)
      REAL(DP)              :: POPZER,POPZR2,POPZCH
      INTEGER,ALLOCATABLE   :: IPZERO(:,:),IGZERO(:,:)
      REAL(DP), ALLOCATABLE :: WNHINT(:,:),WNHEII(:,:),WOP(:,:)

      REAL(DP), ALLOCATABLE :: REINT(:),REDIF(:)

      REAL(DP)              :: TAUDIV
      INTEGER               :: IDLST

      REAL(DP), ALLOCATABLE :: VTURB(:),VTURBS(:)
      INTEGER               :: IPTURB

      REAL(DP), ALLOCATABLE :: SBPSI(:,:),SBLPSI(:,:)
      REAL(DP), ALLOCATABLE :: DSBPST(:,:),DSBPSN(:,:)
      INTEGER,ALLOCATABLE   :: ILTERF(:,:)
      REAL(DP), ALLOCATABLE :: PT(:,:),PN(:,:)
      REAL(DP), ALLOCATABLE :: PP(:,:)
      REAL(DP), ALLOCATABLE :: USUMS(:,:)
      REAL(DP), ALLOCATABLE :: DUSMT(:,:)
      REAL(DP), ALLOCATABLE :: DUSMN(:,:)
      REAL(DP), ALLOCATABLE :: DIESIG(:,:)
      REAL(DP), ALLOCATABLE :: SGM0(:),FRCH(:),SGEXT1(:,:)
      REAL(DP), ALLOCATABLE :: GMER(:,:)
      REAL(DP), ALLOCATABLE :: SGMSUM(:,:,:)
      REAL(DP), ALLOCATABLE :: SGMSUD(:,:,:)
      REAL(DP), ALLOCATABLE :: SGMG(:,:)
      REAL(DP), ALLOCATABLE :: DUSUMT(:),DUSUMN(:)

      REAL(DP)              :: XK,DBETA,BETAD,ADH,DIVH

      REAL(DP), ALLOCATABLE :: ELEC23(:)
      REAL(DP), ALLOCATABLE :: ACOR(:)
      REAL(DP), ALLOCATABLE :: Z3(:)
      REAL(DP), ALLOCATABLE :: DWC1(:,:)
      REAL(DP), ALLOCATABLE :: DWC2(:)
      REAL(DP), ALLOCATABLE :: DWF1(:,:)

      REAL(DP), ALLOCATABLE :: GF0(:),GF1(:),GF2(:)
      REAL(DP), ALLOCATABLE :: GF3(:),GF4(:),GF5(:),GF6(:)
      REAL(DP), ALLOCATABLE :: GF0D(:),GF1D(:),GF2D(:)
      REAL(DP), ALLOCATABLE :: GF3D(:),GF4D(:),GF5D(:)
      REAL(DP), ALLOCATABLE :: GF6D(:),DELTT(:)

      REAL(DP), ALLOCATABLE :: SFF3(:,:)
      REAL(DP), ALLOCATABLE :: SFF2(:,:)
      REAL(DP), ALLOCATABLE :: DSFF(:,:)
      REAL(DP), ALLOCATABLE :: CFFN(:),CFFT(:)

      REAL(DP), ALLOCATABLE :: XKF(:),XKF1(:),XKFB(:)
      REAL(DP), ALLOCATABLE :: ABSO1(:),EMIS1(:),SCAT1(:)
      REAL(DP), ALLOCATABLE :: ABSOT(:),ABSOE1(:)
      REAL(DP), ALLOCATABLE :: EMEL1(:),ABSO1L(:),EMIS1L(:)
      REAL(DP), ALLOCATABLE :: ABSOPR(:),EMISPR(:)
      REAL(DP), ALLOCATABLE :: ABSO(:),EMIS(:),SCAT(:)

      REAL(DP), ALLOCATABLE :: DABT1(:),DEMT1(:)
      REAL(DP), ALLOCATABLE :: DABN1(:),DEMN1(:)
      REAL(DP), ALLOCATABLE :: DABM1(:),DEMM1(:)
      REAL(DP), ALLOCATABLE :: DABX1(:),DEMX1(:)
      REAL(DP), ALLOCATABLE :: DABP1(:,:),DEMP1(:,:)
      REAL(DP), ALLOCATABLE :: DRCH1(:,:),DRET1(:,:)
      REAL(DP), ALLOCATABLE :: ABSFF(:),DABFT(:),DABFN(:)
      REAL(DP), ALLOCATABLE :: DSFDT(:),DSFDN(:),DSFDM(:)
      REAL(DP), ALLOCATABLE :: DSFDP(:,:)
      REAL(DP), ALLOCATABLE :: DSFDTM(:),DSFDNM(:)
      REAL(DP), ALLOCATABLE :: DSFDPM(:,:)
      REAL(DP), ALLOCATABLE :: DSFDTP(:),DSFDNP(:)
      REAL(DP), ALLOCATABLE :: DSFDPP(:,:),DSFP1D(:)
      REAL(DP)              :: DSFT1D,DSFN1D,DSFM1D

      REAL(DP), ALLOCATABLE :: ALIM1(:),ALIP1(:)
      REAL(DP), ALLOCATABLE :: DSCT1(:),DSCN1(:),DRHODT(:)

      REAL(DP), ALLOCATABLE :: FLTOT(:),FLFIX(:),FLEXP(:),FCOOLI(:) 
      REAL(DP), ALLOCATABLE :: FCOOL(:),FPRD(:),FLRD(:)
      REAL(DP), ALLOCATABLE :: FPRAD(:),GRAD(:)

      REAL(DP), ALLOCATABLE :: PTOTAL(:),PGS(:),PRADT(:)
      REAL(DP), ALLOCATABLE :: PRADA(:)
      REAL(DP)              :: PRD0
      INTEGER               :: IHECOR
      REAL(DP), ALLOCATABLE :: ABROSD(:),SUMDPL(:)
      REAL(DP), ALLOCATABLE :: ABPLAD(:)
      REAL(DP)              :: ABPMIN
      REAL(DP), ALLOCATABLE :: QFIX(:)
      REAL(DP), ALLOCATABLE :: ALBE(:)
      INTEGER               :: IWINBL
      REAL(DP)              :: DJMAX
      INTEGER               :: NTRALI
      REAL(DP)              :: ABAD,EMAD,SCAD,DAT,DAN,DET,DEN,DST,DSN
      REAL(DP), ALLOCATABLE :: DDN(:)
      REAL(DP), ALLOCATABLE :: FLXC(:),DELTA(:)
      REAL(DP)              :: RSAT,RSBT,RSAN,RSBN
      REAL(DP), ALLOCATABLE :: RSAX(:),RSBX(:)
      REAL(DP), ALLOCATABLE :: TAUROS(:),TAUFLX(:)
      REAL(DP), ALLOCATABLE :: TAUTHE(:),THETA(:)
      REAL(DP), ALLOCATABLE :: TROSS(:)
      REAL(DP), ALLOCATABLE :: TAURS(:),DT(:)
      REAL(DP)              :: TDISK
      INTEGER               :: ITCONS

      REAL(DP), ALLOCATABLE :: PHMOL(:),QADD(:)
      REAL(DP)              :: ANEREL
      INTEGER               :: IHM,IH2,IH2P
      REAL(DP)              :: ANP,AHTOT,AHMOL

      REAL(DP), ALLOCATABLE :: RR(:,:)
      INTEGER               :: IREF,IREFA
!     LOGICAL,ALLOCATABLE   :: LGR(:),LRM(:)
      REAL(DP)              :: Q,QM,DQT,DQN,DQM,ENER,ENTR,QREF,DQTR,DQNR
      REAL(DP)              :: PFHYD

      REAL(DP), ALLOCATABLE :: CHANT(:)
      REAL(DP)              :: ELSTD
      INTEGER               :: IDSTD

      REAL(DP)              :: DELMDE


      REAL(DP), ALLOCATABLE :: GRD(:),PRA(:),PGS0(:),ANTP(:)
      REAL(DP), ALLOCATABLE :: GRADF(:,:)

! psy

      REAL(DP), ALLOCATABLE :: PSY0(:,:),PSY1(:,:)
      REAL(DP), ALLOCATABLE :: PSY2(:,:),PSY3(:,:)



! mtrans
  
      INTEGER,ALLOCATABLE   :: MCDW(:),ITRCDW(:)
      REAL(DP), ALLOCATABLE :: ABTRA(:,:),EMTRA(:,:)
      REAL(DP), ALLOCATABLE :: DEMLT(:,:)

      REAL(DP), ALLOCATABLE :: RRU(:,:),RRD(:,:)
      REAL(DP), ALLOCATABLE :: DRDT(:,:)
      REAL(DP), ALLOCATABLE :: RDDP(:,:),RDDM(:,:)
      REAL(DP), ALLOCATABLE :: COLRAT(:,:),COLTAR(:,:)
      REAL(DP), ALLOCATABLE :: COL(:),CLOC(:)

      REAL(DP), ALLOCATABLE :: ESEMAT(:,:),BESE(:) 


! mfreq,mfreql

      REAL(DP), ALLOCATABLE :: RAD(:,:),FHD(:)
      REAL(DP), ALLOCATABLE :: FAK(:,:),RADK(:,:)
      REAL(DP), ALLOCATABLE :: EXTRAD(:),EXTINT(:,:),HEXTRD(:)
      REAL(DP), ALLOCATABLE :: EXTJ(:),EXTH(:)
      REAL(DP)              :: TRAD,WDIL,EXTOT,TSTAR,WANGLE
      REAL(DP), ALLOCATABLE :: RAD1(:),ALI1(:),FAK1(:)
      REAL(DP), ALLOCATABLE :: radcm(:,:)
      REAL(DP), ALLOCATABLE :: alih1(:)
      REAL(DP), ALLOCATABLE :: FREQ(:),W(:),PROF(:),WCH(:)
      REAL(DP),ALLOCATABLE  :: FREQC(:)
      INTEGER,ALLOCATABLE   :: JIK(:),IJX(:),IJBF(:)
      INTEGER               :: IFS0
      INTEGER,ALLOCATABLE   :: KIJ(:)
      LOGICAL,ALLOCATABLE   :: LSKIP(:,:)
      REAL(DP)              :: FRCMAX,FRCMIN,FRLMAX,FRLMIN,CFRMAX,DFTAIL
      REAL(DP)              :: TSNU,VTNU,DDNU,DDNU2,CNU1,CNU2,CNU3
      INTEGER               :: IELNU,NFTAIL
      INTEGER,ALLOCATABLE   :: NITJ(:),IJLIN(:),ITRLIN(:,:)
      INTEGER,ALLOCATABLE   :: NLINES(:)
      REAL(DP), ALLOCATABLE :: AIJBF(:),BFCS(:,:)
      INTEGER,ALLOCATABLE   :: IFREQB(:)
      REAL(DP), ALLOCATABLE :: W0E(:),BNUE(:),WC(:)

      INTEGER,ALLOCATABLE   :: IJTC(:),IJEX(:),IJFR(:)

      REAL(DP), ALLOCATABLE :: FLUX(:),FH(:),Q0(:),UU0(:)
      REAL(DP), ALLOCATABLE :: RADEX(:,:),FAKEX(:,:)
      REAL(DP), ALLOCATABLE :: ABSOEX(:,:),EMISEX(:,:)
      REAL(DP), ALLOCATABLE :: SCATEX(:,:)

      REAL(DP), ALLOCATABLE :: PRFLIN(:,:),PRF(:)

! compton

      REAL(DP), ALLOCATABLE :: DLNFR(:),BNUS(:)
      REAL(DP), ALLOCATABLE :: CDER10(:),CDER1P(:),CDER1M(:)
      REAL(DP), ALLOCATABLE :: CDER20(:),CDER2P(:),CDER2M(:)
      REAL(DP), ALLOCATABLE :: DELJ(:,:)


! ltegr

      REAL(DP)              :: TAUFIR,TAULAS,ABROS0,TSURF,ALBAVE,DION0
      REAL(DP)              :: DM1,ABPLA0
      INTEGER               :: NNEWD,NDGREY,IDGREY


! disk

      REAL(DP), ALLOCATABLE :: TVISC(:),DTVIST(:),DTVISN(:)
      REAL(DP), ALLOCATABLE :: DTVISR(:)
      REAL(DP), ALLOCATABLE :: GAMJ(:)
      REAL(DP)              :: GAMH,FAK0

! molec

      REAL(DP), ALLOCATABLE :: ANATO(:,:),ANION(:,:),ANMOL(:,:)

      end module modelq

!*************************************************************************


      module optabl
!     opacity tables

      use accura
      use basics
      use atomic

      INTEGER               :: MTABT, MTABR
      INTEGER               :: MFRTAB
      REAL(DP)              :: FRTB1,FRTB2,RTAB1,RTAB2,TTAB1,TTAB2
      INTEGER               :: numfreq,numrho,numtemp
      INTEGER,  ALLOCATABLE :: numrh(:)
      REAL(DP), ALLOCATABLE :: tempvec(:),rhovec(:,:)
      REAL(DP), ALLOCATABLE :: rhomat(:,:)
      REAL(DP), ALLOCATABLE :: frtab(:),frlt(:)
      REAL(DP)              :: frtlim
      REAL, ALLOCATABLE     :: absopac(:,:,:)
      REAL(DP), ALLOCATABLE :: raytab(:,:)
      REAL(DP), ALLOCATABLE :: raysc(:)
      INTEGER               :: ibinop
      INTEGER, ALLOCATABLE  :: JINT(:)
      REAL(DP), ALLOCATABLE :: YINR(:)
      REAL(DP)              :: abunt(matom),abuno(matom),tmolit
      INTEGER               :: iophmt,ioph2t,iophet,iopcht,iopoht,
     *                         ioh2mt,ih2h2t,ih2het,ioh2ht,iohhet,
     *                         ifmolt
      REAL(DP), ALLOCATABLE :: elecgr(:,:)
      REAL(DP), ALLOCATABLE :: RCS(:),RCHE(:),RCH2(:)

      end module optabl

!*************************************************************************



      module molec

      use accura
      use basics

      REAL(DP), ALLOCATABLE :: C(:,:),PPMOL(:),APMLOG(:),
     *                         XIP(:),XIP2(:),CCOMP(:),UIIDUI(:),
     *                         P(:),FP(:),XKP(:),XK2(:)

      INTEGER,  ALLOCATABLE :: NELEM(:,:),NATO(:,:),MMAX(:),
     *                         NELEMX(:)
      INTEGER               :: NMOLEC,NIMAX
      INTEGER, PARAMETER    :: NMETAL=92

      REAL(DP), ALLOCATABLE :: anion2(:,:)
      REAL(DP), ALLOCATABLE :: entato(:),ention(:),entmol(:)
      REAL(DP), ALLOCATABLE :: uelem(:),ull(:),anden(:),
     *                         aelem(:),ammol(:),
     *                         anat0(:),anio0(:),anmo0(:),pfmol(:),
     *                         denso(:),eleco(:),wmmo(:)

      CHARACTER(LEN=8)      :: CMOL(600)


      REAL(DP)              :: RHOTER,ANTA,ENTRP

      end module molec


!*************************************************************************

      module topbase
  
      use accura
      
      INTEGER, PARAMETER    :: MMAXOP=200, ! maximum number of levels in OP data
     *                         MOP=15      ! maximum number of fit points per level
      REAL(DP), ALLOCATABLE :: SOP(:,:),XOP(:,:),XFIT(:),SFIT(:)
      INTEGER, ALLOCATABLE  :: NOP(:)

      INTEGER               :: NTOTOP  ! total number of levels in OP data
      CHARACTER(LEN=10)     :: IDLVOP(MMAXOP) ! level identifyer Opacity-Project data
      LOGICAL               :: LOPREA  ! .T. OP data read in; .F. OP data not yer read in

      REAL(DP), ALLOCATABLE :: FRINSG(:),CRIN(:)
      INTEGER, ALLOCATABLE  :: JKF(:) 

      end module topbase


!*************************************************************************

      module allarp

      use accura

      integer, parameter :: NXMAX=1400,NNMAX=5,NTAMAX=6
      real(dp),parameter ::  xnorma=8.8528e-29*1215.6*1215.6*0.41618,
     *           xnormb=8.8528e-29*1025.73*1025.7*0.0791,
     *           xnormg=8.8528e-29*972.53*972.53*0.0290,
     *           xnormc=8.8528e-29*6562.*6562.*0.6407
      REAL(DP), ALLOCATABLE :: XLALP(:),PLALP(:,:)
      REAL(DP), ALLOCATABLE :: XLBET(:),PLBET(:,:)
      REAL(DP), ALLOCATABLE :: XLGAM(:),PLGAM(:,:)
      REAL(DP), ALLOCATABLE :: XLBAL(:),PLBAL(:,:)
      REAL(DP)              :: STNNEA,STNCHA,VNEUA,VCHAA
      REAL(DP)              :: STNNEB,STNCHB,VNEUB,VCHAB
      REAL(DP)              :: STNNEG,STNCHG,VNEUG,VCHAG
      REAL(DP)              :: STNNEC,STNCHC,VNEUC,VCHAC
      INTEGER               :: NXALP,IWARNA
      INTEGER               :: NXBET,IWARNB
      INTEGER               :: NXGAM,IWARNG
      INTEGER               :: NXBAL,IWARNC
      REAL(DP), ALLOCATABLE :: XLALPD(:,:),PLALPD(:,:,:)
      REAL(DP), ALLOCATABLE :: STNEAD(:),STNCHD(:),TALPD(:)
      REAL(DP), ALLOCATABLE :: VNEUAD(:),VCHAAD(:)
      INTEGER, ALLOCATABLE  :: NXALPD(:)

      INTEGER               :: NTALPD
!     INTEGER               :: tqmprf,iquasi,nunalp,nunbet,nungam,nunbal

      end module allarp

      
      
!*************************************************************************
      
      module accel

      use accura

      REAL(DP), ALLOCATABLE :: POPUL1(:,:),POPUL2(:,:),POPUL3(:,:)
!     REAL(DP), ALLOCATABLE :: PSY1(:,:),PSY2(:,:),PSY3(:,:)

      end module accel


!*************************************************************************


      module thermo

      use accura

      REAL(DP),ALLOCATABLE :: SL(:,:),PL(:,:)
      REAL(DP) :: redge,pedge(100),sedge(100),cvedge(100),
     *              cpedge(100),gammaedge(100),tedge(100)
       REAL(DP):: R1,R2,T1,T2,T12,T22
       INTEGER :: INDEX,JON

       end module thermo


!*************************************************************************

      module conveq

      use accura

      real(dp) :: heatcp,dlrdlt,hscale,flco,vco,rhod,flxtot,gravd
      real(dp) :: acnv,bcnv,ddel
      INTEGER  :: ICBEG,ICEND,ICBEGP,ICENDP,ICBEG0,ICBEGD

      end module conveq


!*************************************************************************

      module xdata

      use accura

      INTEGER, parameter   :: mtrx=1000

      integer, allocatable :: iex(:),itrind(:),izx0(:),izx1(:),
     *                        nmaxx(:),izx(:),nshx(:),nax(:),icx(:)
      real(dp),allocatable :: etx(:),ssx(:),dx(:),aphx(:,:,:),bphx(:,:)

      end module xdata



!*************************************************************************


      module rybpar

      use accura

      REAL(DP), ALLOCATABLE :: RA(:),RB(:),RC(:),VR(:),
     *                         UA(:),UB(:),UC(:),
     *                         VA(:),VB(:),VC(:),WR(:),
     *                         WM(:,:)

      REAL(DP), ALLOCATABLE :: ABSOPP(:,:),SCATPP(:,:),
     *                         EMISPP(:,:),BFABS(:,:)

      REAL(DP), ALLOCATABLE :: CSND(:),PRAD2D(:)
      REAL(DP)              :: F1HE
      
      end module rybpar
 


!*************************************************************************

      module disaux

      use accura
      use basics

      REAL(DP),ALLOCATABLE :: DDM(:),DDP(:),DD0(:),
     *              DDMIN(:),DDPLU(:),DDA(:),
     *              DDC(:),DDB(:)
      REAL(DP),ALLOCATABLE :: VSND2(:)

      end module disaux

      
!*************************************************************************
      


      module ctdata
!
!      real CTIon
!     second dimension is ionization stage,
!     1=+0 for parent, etc
!     third dimension is atomic number of atom
!      real CTRecomb
!     second dimension is ionization stage,
!     1=+1 for parent, etc
!     third dimension is atomic number of atom
!
!     digital form of the fits to the charge transfer
!     ionization rate coefficients
!
!     Note: First parameter is in units of 1e-9!
!     Note: Seventh parameter is in units of 1e4 K
!     ionization
      
      use accura
      
      implicit none

      real(dp) :: CTIon(7,4,30),CTRecomb(6,4,30)
      integer, private  :: i

      data (CTIon(i,1,3),i=1,7)/2.84e-3,1.99,375.54,-54.07,1e2,1e4,0.0/
      data (CTIon(i,2,3),i=1,7)/7*0./
      data (CTIon(i,3,3),i=1,7)/7*0./
      data (CTIon(i,1,4),i=1,7)/7*0./
      data (CTIon(i,2,4),i=1,7)/7*0./
      data (CTIon(i,3,4),i=1,7)/7*0./
      data (CTIon(i,1,5),i=1,7)/7*0./
      data (CTIon(i,2,5),i=1,7)/7*0./
      data (CTIon(i,3,5),i=1,7)/7*0./
      data (CTIon(i,1,6),i=1,7)/1.07e-6,3.15,176.43,-4.29,1e3,1e5,0.0/
      data (CTIon(i,2,6),i=1,7)/7*0./
      data (CTIon(i,3,6),i=1,7)/7*0./
      data (CTIon(i,1,7),i=1,7)/4.55e-3,-0.29,-0.92,-8.38,1e2,5e4,1.086/
      data (CTIon(i,2,7),i=1,7)/7*0./
      data (CTIon(i,3,7),i=1,7)/7*0./
      data (CTIon(i,1,8),i=1,7)/7.40e-2,0.47,24.37,-0.74,1e1,1e4,0.023/
      data (CTIon(i,2,8),i=1,7)/7*0./
      data (CTIon(i,3,8),i=1,7)/7*0./
      data (CTIon(i,1,9),i=1,7)/7*0./
      data (CTIon(i,2,9),i=1,7)/7*0./
      data (CTIon(i,3,9),i=1,7)/7*0./
      data (CTIon(i,1,10),i=1,7)/7*0./
      data (CTIon(i,2,10),i=1,7)/7*0./
      data (CTIon(i,3,10),i=1,7)/7*0./
      data (CTIon(i,1,11),i=1,7)/3.34e-6,9.31,2632.31,-3.04,1e3,2e4,0.0/
      data (CTIon(i,2,11),i=1,7)/7*0./
      data (CTIon(i,3,11),i=1,7)/7*0./
      data (CTIon(i,1,12),i=1,7)/9.76e-3,3.14,55.54,-1.12,5e3,3e4,0.0/
      data (CTIon(i,2,12),i=1,7)/7.60e-5,0.00,-1.97,-4.32,1e4,3e5,1.670/
      data (CTIon(i,3,12),i=1,7)/7*0./
      data (CTIon(i,1,13),i=1,7)/7*0./
      data (CTIon(i,2,13),i=1,7)/7*0./
      data (CTIon(i,3,13),i=1,7)/7*0./
      data (CTIon(i,1,14),i=1,7)/0.92,1.15,0.80,-0.24,1e3,2e5,0.0/
      data (CTIon(i,2,14),i=1,7)/2.26,7.36e-2,-0.43,-0.11,2e3,1e5,
     & 3.031/
      data (CTIon(i,3,14),i=1,7)/7*0./
      data (CTIon(i,1,15),i=1,7)/7*0./
      data (CTIon(i,2,15),i=1,7)/7*0./
      data (CTIon(i,3,15),i=1,7)/7*0./
      data (CTIon(i,1,16),i=1,7)/1.00e-5,0.00,0.00,0.00,1e3,1e4,0.0/
      data (CTIon(i,2,16),i=1,7)/7*0./
      data (CTIon(i,3,16),i=1,7)/7*0./
      data (CTIon(i,1,17),i=1,7)/7*0./
      data (CTIon(i,2,17),i=1,7)/7*0./
      data (CTIon(i,3,17),i=1,7)/7*0./
      data (CTIon(i,1,18),i=1,7)/7*0./
      data (CTIon(i,2,18),i=1,7)/7*0./
      data (CTIon(i,3,18),i=1,7)/7*0./
      data (CTIon(i,1,19),i=1,7)/7*0./
      data (CTIon(i,2,19),i=1,7)/7*0./
      data (CTIon(i,3,19),i=1,7)/7*0./
      data (CTIon(i,1,20),i=1,7)/7*0./
      data (CTIon(i,2,20),i=1,7)/7*0./
      data (CTIon(i,3,20),i=1,7)/7*0./
      data (CTIon(i,1,21),i=1,7)/7*0./
      data (CTIon(i,2,21),i=1,7)/7*0./
      data (CTIon(i,3,21),i=1,7)/7*0./
      data (CTIon(i,1,22),i=1,7)/7*0./
      data (CTIon(i,2,22),i=1,7)/7*0./
      data (CTIon(i,3,22),i=1,7)/7*0./
      data (CTIon(i,1,23),i=1,7)/7*0./
      data (CTIon(i,2,23),i=1,7)/7*0./
      data (CTIon(i,3,23),i=1,7)/7*0./
      data (CTIon(i,1,24),i=1,7)/7*0./
      data (CTIon(i,2,24),i=1,7)/4.39,0.61,-0.89,-3.56,1e3,3e4,3.349/
      data (CTIon(i,3,24),i=1,7)/7*0./
      data (CTIon(i,1,25),i=1,7)/7*0./
      data (CTIon(i,2,25),i=1,7)/2.83e-1,6.80e-3,6.44e-2,-9.70,1e3,3e4,
     & 2.368/
      data (CTIon(i,3,25),i=1,7)/7*0./
      data (CTIon(i,1,26),i=1,7)/7*0./
      data (CTIon(i,2,26),i=1,7)/2.10,7.72e-2,-0.41,-7.31,1e4,1e5,3.005/
      data (CTIon(i,3,26),i=1,7)/7*0./
      data (CTIon(i,1,27),i=1,7)/7*0./
      data (CTIon(i,2,27),i=1,7)/1.20e-2,3.49,24.41,-1.26,1e3,3e4,4.044/
      data (CTIon(i,3,27),i=1,7)/7*0./
      data (CTIon(i,1,28),i=1,7)/7*0./
      data (CTIon(i,2,28),i=1,7)/7*0./
      data (CTIon(i,3,28),i=1,7)/7*0./
      data (CTIon(i,1,29),i=1,7)/7*0./
      data (CTIon(i,2,29),i=1,7)/7*0./
      data (CTIon(i,3,29),i=1,7)/7*0./
      data (CTIon(i,1,30),i=1,7)/7*0./
      data (CTIon(i,2,30),i=1,7)/7*0./
      data (CTIon(i,3,30),i=1,7)/7*0./
!
!     digital form of the fits to the charge transfer
!     recombination rate coefficients (total)
!
!     Note: First parameter is in units of 1e-9!
!     recombination
      data (CTRecomb(i,1,2),i=1,6)/7.47e-6,2.06,9.93,-3.89,6e3,1e5/
      data (CTRecomb(i,2,2),i=1,6)/1.00e-5,0.,0.,0.,1e3,1e7/
      data (CTRecomb(i,1,3),i=1,6)/6*0./
      data (CTRecomb(i,2,3),i=1,6)/1.26,0.96,3.02,-0.65,1e3,3e4/
      data (CTRecomb(i,3,3),i=1,6)/1.00e-5,0.,0.,0.,2e3,5e4/
      data (CTRecomb(i,1,4),i=1,6)/6*0./
      data (CTRecomb(i,2,4),i=1,6)/1.00e-5,0.,0.,0.,2e3,5e4/
      data (CTRecomb(i,3,4),i=1,6)/1.00e-5,0.,0.,0.,2e3,5e4/
      data (CTRecomb(i,4,4),i=1,6)/5.17,0.82,-0.69,-1.12,2e3,5e4/
      data (CTRecomb(i,1,5),i=1,6)/6*0./
      data (CTRecomb(i,2,5),i=1,6)/2.00e-2,0.,0.,0.,1e3,1e9/
      data (CTRecomb(i,3,5),i=1,6)/1.00e-5,0.,0.,0.,2e3,5e4/
      data (CTRecomb(i,4,5),i=1,6)/2.74,0.93,-0.61,-1.13,2e3,5e4/
      data (CTRecomb(i,1,6),i=1,6)/4.88e-7,3.25,-1.12,-0.21,5.5e3,1e5/
      data (CTRecomb(i,2,6),i=1,6)/1.67e-4,2.79,304.72,-4.07,5e3,5e4/
      data (CTRecomb(i,3,6),i=1,6)/3.25,0.21,0.19,-3.29,1e3,1e5/
      data (CTRecomb(i,4,6),i=1,6)/332.46,-0.11,-9.95e-1,-1.58e-3,1e1,
     & 1e5/
      data (CTRecomb(i,1,7),i=1,6)/1.01e-3,-0.29,-0.92,-8.38,1e2,5e4/
      data (CTRecomb(i,2,7),i=1,6)/3.05e-1,0.60,2.65,-0.93,1e3,1e5/
      data (CTRecomb(i,3,7),i=1,6)/4.54,0.57,-0.65,-0.89,1e1,1e5/
      data (CTRecomb(i,4,7),i=1,6)/2.95,0.55,-0.39,-1.07,1e3,1e6/
      data (CTRecomb(i,1,8),i=1,6)/1.04,3.15e-2,-0.61,-9.73,1e1,1e4/
      data (CTRecomb(i,2,8),i=1,6)/1.04,0.27,2.02,-5.92,1e2,1e5/
      data (CTRecomb(i,3,8),i=1,6)/3.98,0.26,0.56,-2.62,1e3,5e4/
      data (CTRecomb(i,4,8),i=1,6)/2.52e-1,0.63,2.08,-4.16,1e3,3e4/
      data (CTRecomb(i,1,9),i=1,6)/6*0./
      data (CTRecomb(i,2,9),i=1,6)/1.00e-5,0.,0.,0.,2e3,5e4/
      data (CTRecomb(i,3,9),i=1,6)/9.86,0.29,-0.21,-1.15,2e3,5e4/
      data (CTRecomb(i,4,9),i=1,6)/7.15e-1,1.21,-0.70,-0.85,2e3,5e4/
      data (CTRecomb(i,1,10),i=1,6)/6*0./
      data (CTRecomb(i,2,10),i=1,6)/1.00e-5,0.,0.,0.,5e3,5e4/
      data (CTRecomb(i,3,10),i=1,6)/14.73,4.52e-2,-0.84,-0.31,5e3,5e4/
      data (CTRecomb(i,4,10),i=1,6)/6.47,0.54,3.59,-5.22,1e3,3e4/
      data (CTRecomb(i,1,11),i=1,6)/6*0./
      data (CTRecomb(i,2,11),i=1,6)/1.00e-5,0.,0.,0.,2e3,5e4/
      data (CTRecomb(i,3,11),i=1,6)/1.33,1.15,1.20,-0.32,2e3,5e4/
      data (CTRecomb(i,4,11),i=1,6)/1.01e-1,1.34,10.05,-6.41,2e3,5e4/
      data (CTRecomb(i,1,12),i=1,6)/6*0./
      data (CTRecomb(i,2,12),i=1,6)/8.58e-5,2.49e-3,2.93e-2,-4.33,1e3,
     & 3e4/
      data (CTRecomb(i,3,12),i=1,6)/6.49,0.53,2.82,-7.63,1e3,3e4/
      data (CTRecomb(i,4,12),i=1,6)/6.36,0.55,3.86,-5.19,1e3,3e4/
      data (CTRecomb(i,1,13),i=1,6)/6*0./
      data (CTRecomb(i,2,13),i=1,6)/1.00e-5,0.,0.,0.,1e3,3e4/
      data (CTRecomb(i,3,13),i=1,6)/7.11e-5,4.12,1.72e4,-22.24,1e3,3e4/
      data (CTRecomb(i,4,13),i=1,6)/7.52e-1,0.77,6.24,-5.67,1e3,3e4/
      data (CTRecomb(i,1,14),i=1,6)/6*0./
      data (CTRecomb(i,2,14),i=1,6)/6.77,7.36e-2,-0.43,-0.11,5e2,1e5/
      data (CTRecomb(i,3,14),i=1,6)/4.90e-1,-8.74e-2,-0.36,-0.79,1e3,
     & 3e4/
      data (CTRecomb(i,4,14),i=1,6)/7.58,0.37,1.06,-4.09,1e3,5e4/
      data (CTRecomb(i,1,15),i=1,6)/6*0./
      data (CTRecomb(i,2,15),i=1,6)/1.74e-4,3.84,36.06,-0.97,1e3,3e4/
      data (CTRecomb(i,3,15),i=1,6)/9.46e-2,-5.58e-2,0.77,-6.43,1e3,3e4/
      data (CTRecomb(i,4,15),i=1,6)/5.37,0.47,2.21,-8.52,1e3,3e4/
      data (CTRecomb(i,1,16),i=1,6)/3.82e-7,11.10,2.57e4,-8.22,1e3,1e4/
      data (CTRecomb(i,2,16),i=1,6)/1.00e-5,0.,0.,0.,1e3,3e4/
      data (CTRecomb(i,3,16),i=1,6)/2.29,4.02e-2,1.59,-6.06,1e3,3e4/
      data (CTRecomb(i,4,16),i=1,6)/6.44,0.13,2.69,-5.69,1e3,3e4/
      data (CTRecomb(i,1,17),i=1,6)/6*0./
      data (CTRecomb(i,2,17),i=1,6)/1.00e-5,0.,0.,0.,1e3,3e4/
      data (CTRecomb(i,3,17),i=1,6)/1.88,0.32,1.77,-5.70,1e3,3e4/
      data (CTRecomb(i,4,17),i=1,6)/7.27,0.29,1.04,-10.14,1e3,3e4/
      data (CTRecomb(i,1,18),i=1,6)/6*0./
      data (CTRecomb(i,2,18),i=1,6)/1.00e-5,0.,0.,0.,1e3,3e4/
      data (CTRecomb(i,3,18),i=1,6)/4.57,0.27,-0.18,-1.57,1e3,3e4/
      data (CTRecomb(i,4,18),i=1,6)/6.37,0.85,10.21,-6.22,1e3,3e4/
      data (CTRecomb(i,1,19),i=1,6)/6*0./
      data (CTRecomb(i,2,19),i=1,6)/1.00e-5,0.,0.,0.,1e3,3e4/
      data (CTRecomb(i,3,19),i=1,6)/4.76,0.44,-0.56,-0.88,1e3,3e4/
      data (CTRecomb(i,4,19),i=1,6)/1.00e-5,0.,0.,0.,1e3,3e4/
      data (CTRecomb(i,1,20),i=1,6)/6*0./
      data (CTRecomb(i,2,20),i=1,6)/0.,0.,0.,0.,1e1,1e9/
      data (CTRecomb(i,3,20),i=1,6)/3.17e-2,2.12,12.06,-0.40,1e3,3e4/
      data (CTRecomb(i,4,20),i=1,6)/2.68,0.69,-0.68,-4.47,1e3,3e4/
      data (CTRecomb(i,1,21),i=1,6)/6*0./
      data (CTRecomb(i,2,21),i=1,6)/0.,0.,0.,0.,1e1,1e9/
      data (CTRecomb(i,3,21),i=1,6)/7.22e-3,2.34,411.50,-13.24,1e3,3e4/
      data (CTRecomb(i,4,21),i=1,6)/1.20e-1,1.48,4.00,-9.33,1e3,3e4/
      data (CTRecomb(i,1,22),i=1,6)/6*0./
      data (CTRecomb(i,2,22),i=1,6)/0.,0.,0.,0.,1e1,1e9/
      data (CTRecomb(i,3,22),i=1,6)/6.34e-1,6.87e-3,0.18,-8.04,1e3,3e4/
      data (CTRecomb(i,4,22),i=1,6)/4.37e-3,1.25,40.02,-8.05,1e3,3e4/
      data (CTRecomb(i,1,23),i=1,6)/6*0./
      data (CTRecomb(i,2,23),i=1,6)/1.00e-5,0.,0.,0.,1e3,3e4/
      data (CTRecomb(i,3,23),i=1,6)/5.12,-2.18e-2,-0.24,-0.83,1e3,3e4/
      data (CTRecomb(i,4,23),i=1,6)/1.96e-1,-8.53e-3,0.28,-6.46,1e3,3e4/
      data (CTRecomb(i,1,24),i=1,6)/6*0./
      data (CTRecomb(i,2,24),i=1,6)/5.27e-1,0.61,-0.89,-3.56,1e3,3e4/
      data (CTRecomb(i,3,24),i=1,6)/10.90,0.24,0.26,-11.94,1e3,3e4/
      data (CTRecomb(i,4,24),i=1,6)/1.18,0.20,0.77,-7.09,1e3,3e4/
      data (CTRecomb(i,1,25),i=1,6)/6*0./
      data (CTRecomb(i,2,25),i=1,6)/1.65e-1,6.80e-3,6.44e-2,-9.70,1e3,
     & 3e4/
      data (CTRecomb(i,3,25),i=1,6)/14.20,0.34,-0.41,-1.19,1e3,3e4/
      data (CTRecomb(i,4,25),i=1,6)/4.43e-1,0.91,10.76,-7.49,1e3,3e4/
      data (CTRecomb(i,1,26),i=1,6)/6*0./
      data (CTRecomb(i,2,26),i=1,6)/1.26,7.72e-2,-0.41,-7.31,1e3,1e5/
      data (CTRecomb(i,3,26),i=1,6)/3.42,0.51,-2.06,-8.99,1e3,1e5/
      data (CTRecomb(i,4,26),i=1,6)/14.60,3.57e-2,-0.92,-0.37,1e3,3e4/
      data (CTRecomb(i,1,27),i=1,6)/6*0./
      data (CTRecomb(i,2,27),i=1,6)/5.30,0.24,-0.91,-0.47,1e3,3e4/
      data (CTRecomb(i,3,27),i=1,6)/3.26,0.87,2.85,-9.23,1e3,3e4/
      data (CTRecomb(i,4,27),i=1,6)/1.03,0.58,-0.89,-0.66,1e3,3e4/
      data (CTRecomb(i,1,28),i=1,6)/6*0./
      data (CTRecomb(i,2,28),i=1,6)/1.05,1.28,6.54,-1.81,1e3,1e5/
      data (CTRecomb(i,3,28),i=1,6)/9.73,0.35,0.90,-5.33,1e3,3e4/
      data (CTRecomb(i,4,28),i=1,6)/6.14,0.25,-0.91,-0.42,1e3,3e4/
      data (CTRecomb(i,1,29),i=1,6)/6*0./
      data (CTRecomb(i,2,29),i=1,6)/1.47e-3,3.51,23.91,-0.93,1e3,3e4/
      data (CTRecomb(i,3,29),i=1,6)/9.26,0.37,0.40,-10.73,1e3,3e4/
      data (CTRecomb(i,4,29),i=1,6)/11.59,0.20,0.80,-6.62,1e3,3e4/
      data (CTRecomb(i,1,30),i=1,6)/6*0./
      data (CTRecomb(i,2,30),i=1,6)/1.00e-5,0.,0.,0.,1e3,3e4/
      data (CTRecomb(i,3,30),i=1,6)/6.96e-4,4.24,26.06,-1.24,1e3,3e4/
      data (CTRecomb(i,4,30),i=1,6)/1.33e-2,1.56,-0.92,-1.20,1e3,3e4/
c
      end module ctdata
