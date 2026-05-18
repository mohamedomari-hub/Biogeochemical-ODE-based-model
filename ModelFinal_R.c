#include <R.h>
#include<math.h> 

static double parms[69];

#define Kn     parms[0]
#define Kp       parms[1]
#define Kc       parms[2]
#define Iop         parms[3]
#define Gmax         parms[4]
#define bm       parms[5]
#define ArrhBM         parms[6]
#define Vcph            parms[7]
#define ArrhGr   parms[8]
#define NCPhy   parms[9]
#define PCPhy   parms[10]


#define    ArrhFluxDOC parms[11]
#define    ArrhHyd parms[12]
#define    ArrhMinDOC parms[13]

#define    k0  parms[14]
#define    k1  parms[15]
#define    k2  parms[16]
#define    k3  parms[17]
#define    k4  parms[18]

#define    Isat  parms[19]
#define    alpha  parms[20]
#define    beta parms[21]
#define    lossweed  parms[22]

#define    Arrhloss  parms[23]
#define    Tmin  parms[24]
#define    KnDWeed  parms[25]
#define    KpDWeed  parms[26]
#define    KBiomass  parms[27]
#define    GmaxDWeed  parms[28]    
#define    resDW  parms[29]
#define    ArrhresDW  parms[30]

#define    hydpoc  parms[31]
#define    Vpoc parms[32]



#define    mindoc  parms[33]

#define    DOsatdoc  parms[34]
#define    SedfluxDOC parms[35]
#define    Krefrespdoc  parms[36]
#define    KHoxresp  parms[37]

#define    K parms[38]
#define    SedfluxDIC  parms[39]


#define  Vpon parms[40]
#define  mindon parms[41]
#define  DOsatdon  parms[42]
#define  SedfluxNH4  parms[43]

#define    ArrhNit parms[44]
#define    Nitrmax parms[45]
#define    KhnitDo  parms[46]


#define    denit  parms[47]
#define    Arrhdent parms[48]
#define    SedfluxNO3 parms[49]

#define    Vpop  parms[50]

#define    mindop  parms[51]
#define    DOsatdop  parms[52]
#define    SedfluxP  parms[53]
#define    ArrhPflux  parms[54]
#define    Abso parms[55]
#define    Vset parms[56]

#define    SOD  parms[57]
#define    ArrhDO  parms[58]
#define    Rea parms[59]
#define    DOsatsod parms[60]    


#define    co parms[61]
#define    ArrhFluxDIC parms[62]
#define    ArrhFluxNH4 parms[63]
#define    ArrhFluxNO3 parms[64]
#define    ArrhFluxPO4 parms[65]

#define    CCDW parms[66]
#define    NCDW parms[67]
#define    PCDW parms[68]

#define    ArrhMinON parms[69]
#define    ArrhMinOP parms[70]




// ArrhDo and XXX and YYY are repeated


static double forc[9];

#define temp            forc[0]

//static double forcs[1];
#define H            forc[1]

//static double forcss[1];
#define I0            forc[2]

//static double forcs1[1];
#define PH            forc[3]

//static double forcs2[1];
#define Ca            forc[4]

//static double forcs3[1];
#define Fe            forc[5]

#define Mn            forc[6]

#define Wind            forc[7]

#define Vol            forc[8]


#define CPhytoplankton         y[0]
#define DWeed                  y[1]
#define DO                     y[2]

#define CDetritus              y[3]
#define DOC                    y[4]
#define DIC                    y[5]
#define TC                     y[6]

#define NDetritus              y[7]
#define NH4                    y[8]  
#define NO3                    y[9]
#define TN                     y[10]

#define PDetritus              y[11]
#define PO4                    y[12]
#define TP                     y[13]

#define POC                    y[14]
#define ON                     y[15]
#define OP                     y[16]

//#define NPhytoplankton                     y[17]
  //#define PPhytoplankton                     y[18]
  //#define PON                    y[16]
  //#define POP                    y[17]
  
  //#define DOCSediment           y[19]
  
  
  //#define POCSediment           y[22]
  //#define PONSediment           y[23]
  //#define POPSediment           y[24]
  
  
  
  
  #define dCPhytoplanktondt         ydot[0]
  #define dDWeeddt                  ydot[1]
  #define dDOdt                     ydot[2]
  
  #define dCDetritusdt              ydot[3]
  #define dDOCdt                    ydot[4]
  #define dDICdt                    ydot[5]
#define dTCdt                     ydot[6]

#define dNDetritusdt              ydot[7]
#define dNH4dt                    ydot[8]
#define dNO3dt                    ydot[9]
#define dTNdt                     ydot[10]        

#define dPDetritusdt              ydot[11]
#define dPO4dt                    ydot[12]
#define dTPdt                     ydot[13]

#define dPOCdt                    ydot[14]
#define dONdt                     ydot[15]
#define dOPdt                     ydot[16]

//#define dNPhytoplanktondt                     ydot[17]
  //#define dPPhytoplanktondt                     ydot[18]
  /*=======================================================================
  c=======================================================================
  c                      Model initialisation
c=======================================================================
  c=======================================================================
  
  c=======================================================================
  c Initialise parameter common block
c=======================================================================
  */
  void initmod(void (* odeparms)(int *, double *))
{
  int N =  71;
  odeparms(&N, parms);
  }

/*
  c=======================================================================
  c Initialise forcing function common block
c=======================================================================
  */
  
  void forcc(void (* odeforcs)(int *, double *))
{
  int N=9;
  odeforcs(&N, forc);
  }




void derivs (int *neq, double *t, double *y, double *ydot)

{     
  
  double k, UTem,Ulight,GN,GNutN,GP,GNutP,Gr,BM , IN, UTemDWeed, GDWeedN, GDWeedP, LimEffBiomass, GrDWeed, loss, 
  HYD, Setpoc, Min , DOCsf, H0, CO2PH, CO2water, DICsf,  Setpon, Mindon, Nitr, NH4sf, 
  r, NO3sf , Setpop, Mindop, Psf, DOSed, DOsat, Ligth, TemNit, RESDW; 
  
  
  double  BMPOC= 0.25,BMDOC=0.5, BMDIC= 0.25,BMDON= 0.5,BMNH4= 0.5 , BMDOP= 0.5, BMPO4= 0.5, FDOM=0.5;
  //#-----------------------------------------      PHYTOPLANKTON     ---------------------------------------------------------------
  
  //    #Temeprature term:
    //    #----------------
  
  
  //if (temp <= 20)   {UTem  = exp(-minTphy*pow((temp - 20),2));}
  
  //         else                    {UTem  = exp(-maxTphy*pow((temp - 20),2));}   
  
  //if (temp <= 10)   {UTem =exp(-0.004*pow((temp - 20),2));}
  //else if (10 < temp && temp < 18)   UTem =0.0;
  //else                    {UTem =exp(-0.004*pow((temp - 20),2));}
  
  
  //if (temp <= 20)   {Tem  = exp(-minTeff*pow((temp - 20),2));}
  
  //       else                    {Tem  = exp(-maxTeff*pow((temp - 20),2));}     
  
  
  double TemFxDOC, TemFxDIC, TemFxNH4, TemFxNO3, TemFxPO4, TemHy, TemMinDOC, TemMinON, TemMinOP;
  
  //UTem = exp(bmTeff*(temp - 20));
  UTem =   pow(ArrhGr,(temp - 20));
  
  TemFxDOC  = pow(ArrhFluxDOC,(temp - 20));
  TemFxDIC  = pow(ArrhFluxDIC,(temp - 20));
  TemFxNH4  = pow(ArrhFluxNH4,(temp - 20));
  TemFxNO3  = pow(ArrhFluxNO3,(temp - 20));
  TemFxPO4  = pow(ArrhFluxPO4,(temp - 20));
  
  TemHy  = pow(ArrhHyd,(temp - 20));       
  
  TemMinDOC  = pow(ArrhMinDOC,(temp - 20));  
  TemMinON   = pow(ArrhMinON,(temp - 20));
  TemMinOP   = pow(ArrhMinOP,(temp - 20));
  
  TemNit = pow(ArrhNit,(temp - 20));
  
  
  //    #Light term:
    //    #-----------
  
  k = k0 + k1*CPhytoplankton + k4*DWeed + k2*CDetritus + k3*DOC;/* + k4*DWeed ;*/  //# this is Extinction coefficient
    
    //double par = 0.43;
  
  //double J = par*I0*exp(-k*H);
  
  //double Light = (J/Iop)*exp(1 - (J/(Iop)));
  
  double Light = (2.718/(k*H))*(exp(-(I0/Iop)*exp(-k*H)) - exp(-I0/Iop));
  
  Ulight = fmax(0,fmin(1, Light));
  
  //   #Nutrient term:
    //   #-------------
  
  IN = NH4 + NO3;
  
  GN   =  (IN/(IN+Kn)) ;
  //GN   =  ((NPhytoplankton - minNC)/(maxNC - minNC)) ;
  // GN   =  ((maxNC)/(maxNC - minNC))*(1 - (minNC/NPhytoplankton)) ;
  
  GNutN = fmax(0,fmin(1, GN));
  
  
  
  GP   =  (PO4/(PO4+Kp)) ;
  //GP   =  ((PPhytoplankton - minPC)/(maxPC - minPC));
  // GP   =  ((maxPC)/(maxPC - minPC))*(1 - (minPC/PPhytoplankton)) ;
  
  GNutP   =  fmax(0,fmin(1, GP));
  
  double GC, GNutC;
  
  GC   =  (DIC/(DIC+Kc)) ;
  GNutC   =  fmax(0,fmin(1, GC));
  
  //  double kDIC = 0.020;
  //   #Growth function:
    //   #---------------
  double Nut;
  
  Nut = fmin(GNutC,fmin(GNutP,GNutN));
  
  Gr =  Gmax*UTem*Ulight*Nut;//
    
    //   #Mortality function:
    //   #-------------------
  
  
  
  BM  = bm*pow(ArrhBM,(temp - 20)); //bm*exp(bmTeff*(temp - 20)) ;     
  
  //   #Settling function:
    //   #-----------------
  
  double Set;
  
  Set   =  (Vcph/H) ;                                                                                 
  
  //   #Respiration function:
    
    //MOR = mor;
  
  
  
  //   # Phytoplankton Stoichiometries:
    //   # ------------------------------
  
  //   #Uptake:
    //   #-------            
  
  //IN = NH4 + NO3;
  
  //UN       =  (IN/(Kn + IN))*((maxNC - NPhytoplankton)/(maxNC - minNC));
  //UP       =  (PO4/(Kp + PO4))*((maxPC - PPhytoplankton)/(maxPC - minPC));
  
  //UpN      =  upmaxN*UN;
  //UpP      =  upmaxP*UP;
  
  //#---------------------------------------------------------------------------------------------------------
  
  
  //# ODE for Phytoplankton:             
    //# ----------------------
  
  
  dCPhytoplanktondt = (Gr - BM - Set)*CPhytoplankton;                                            
  
  
  //  # Phytoplankton Stoichiometries:
    //  # ------------------------------    
  
  
  //dNPhytoplanktondt = UpN  - Gr*NPhytoplankton ;
  //dNPhytoplanktondt = UpN*CPhytoplankton  - MOR*NPhytoplankton - BM*NPhytoplankton  ; 
  
  //dPPhytoplanktondt = UpP  - Gr*PPhytoplankton ;
  //dPPhytoplanktondt = UpP*CPhytoplankton  - MOR*PPhytoplankton - BM*PPhytoplankton  ; 
  
  //   #-----------------------------------------        DuckWeed        ---------------------------------------------------------------
  
  //  # Temeprature term:
    //  #-----------------
  
  
  
  /* if (tempair <= 20) (UTemDWeed  = exp(-log(10)*(pow((tempair - 20),2)/pow((Tmin - 20),2))));
  else (UTemDWeed  = exp(-log(10)*(pow((tempair - 20),2)/pow((Tmax - 20),2))));*/
    
    UTemDWeed  = fmin(fmax(0,(temp - Tmin)/(26 - Tmin)),1);               
  
  //  # Temeprature term:
    //  #-----------------
  
  
  if (I0 <= Isat) (Ligth = alpha*I0 + beta);
  else(Ligth = 1);
  
  
  //# Nutrient term:
    //#-----------------
  
  GDWeedN = IN/(IN + KnDWeed);
  
  GDWeedP = PO4/(PO4 + KpDWeed);
  
  //# limiting effect of biomass term:
    //#----------------------------
  
  LimEffBiomass =  KBiomass/(KBiomass + DWeed);
  
  //# Growth and loss functions:
    //#---------------------------
  
  GrDWeed     =  GmaxDWeed*UTemDWeed*GDWeedN*GDWeedP*LimEffBiomass*Ligth;   
  
  loss         =  lossweed*pow(Arrhloss,((temp - 20)));
  
  RESDW = resDW*pow(ArrhresDW,(temp - 20));    
  
  
  
  //# ODE for DuckWeed:             
    //# -----------------
  
  dDWeeddt = GrDWeed*DWeed  - 0*RESDW*DWeed  - loss*DWeed; 
  
  
  // #-----------------------------------------         Water plants     ------------------------------------------------------------
  
  //#----------------------------------------     POC     -------------------------------------------------------------           
  
  //#Hydrolysis and Settling:
    //#------------------------  
  
  HYD      =  hydpoc*TemHy;                                                     
  
  Setpoc  =  (Vpoc/H);    
  
  //double resus= 0.01/48, granulation=0.4;
  
  //double Resus = resus*granulation*(Hs/H) ;
  
  //# ODE for POC:
    //# ------------
  
  double Resus, m = 1.2 ; 
  double Winds = Wind*0.27;
  
  if (Winds > 2.)   {Resus = co*pow(Winds,m)/H;}
  
  else              {Resus  = 0;}    
  
  double POCinflow;
  
  if ( 2500 <*t && *t< 3000)   POCinflow = 500000/48;  // 500 gram
  
  else POCinflow = 0;
  
  double  mort=0.25; 
  
  dCDetritusdt = BMPOC*BM*CPhytoplankton  - HYD*CDetritus - Setpoc*CDetritus 
  + 0*Resus*1e-3 + 0*POCinflow/Vol + (1/H)*CCDW*mort*loss*DWeed;/*+ Resus*POCSediment*/ ; //+ POCinflow;    
  
  
  dPOCdt = dCDetritusdt + dCPhytoplanktondt ;
  //#----------------------------------------     DOC     -------------------------------------------------------------     
  
  //# fraction exuded as DOC:
    //#------------------------
  
  //bmDOC    =  BMDOC + (1-BMDOC)*(KHR/(KHR + DO));
  
  
  //# Mineralisation, Denitrification and sediment flux respectively:
    //#---------------------------------------------------------------
  
  Min =  mindoc*TemMinDOC*(DO/(DOsatdoc + DO));
  
  //Kres = Krefrespdoc*TemMin;
  
  //Khr = Kres*(DO/(DO+KHoxresp));
  
  
  
  DOCsf = SedfluxDOC*TemFxDOC*(1/H);   
  //DOCsf = (SedfluxDOC/(Hs/2))*(poro/H) ;
  
  //# ODE for DOC:
    //# ------------
  
  dDOCdt = BMDOC*BM*CPhytoplankton + HYD*CDetritus - Min*DOC  + DOCsf ; //+ DOCinflow; 
  
  
  //#----------------------------------------     DIC     -------------------------------------------------------------     
  
  //# Air exchange:
    //#-------------
  
  
  //      K0 = 3.7*12*(10^(-5)),  Kc1 = 4.3*1e-7,  Kc2 = 5.61*1e-11
  double  K0 = 0.000444,  Kc1 = 4.3e-07,  Kc2 = 5.61e-11, pCO2atm = 360;
  
  H0        =  pow(10,-PH);                                                            
  
  CO2PH    =  pow(H0,2)/(pow(H0,2) + (H0*Kc1) + (Kc1*Kc2));                                        
  
  CO2water  =  (DIC*CO2PH)/K0;                                                        
  
  
  
  //#Diffusion fluxes:
    //#----------------
  
  DICsf = SedfluxDIC*TemFxDIC*(1/H);      
  //DICsf = (SedfluxDIC/(Hs/2))*(poro/H) ;
  
  //# ODE for DIC:
    //# ------------
  
  // dDICdt = (BM - Gr)*CPhytoplankton + (1/H)*(resDICvgs*RESVgs - (1 - Frprp)*GrVgs)*MacrophyteShoot 
  //           + (1/H)*(DICDW*RESDW - GrDWeed)*DWeed + fdenitDIC*(Denitr/DOC)*NO3  + Min*DOC 
  //           + DICsf*(DICSediment - DIC) - K*(CO2water - pCO2atm)   - (Abso/H)*Vset*Ca*DIC ; // + DICinflow;
  
  dDICdt =  (BMDIC*BM  - Gr)*CPhytoplankton 
  + Min*DOC + K*(pCO2atm - CO2water) + DICsf  ; // + DICinflow;
  
  
  dTCdt = dDICdt + dDOCdt + dPOCdt  ;       
  
  //#----------------------------------------     PON     -------------------------------------------------------------     
  
  //# Hydrolysis and Settling:
    //#------------------------
  
  //HYDpon  =  hydpon*TemHy;                                             
  
  
  
  Setpon  =  (1 - FDOM)*(Vpon/H);                                                        
  
  
  //# ODE for PON:
    //# ------------
  
  
  
  //dPONdt  = - HYDpon*PON + BMPON*NCPhy*BM*CPhytoplankton  - 0*Setpon*PON /*+ Resus*PONSediment*/  ; //+ PONinflow;
  
  //#----------------------------------------     DON     -------------------------------------------------------------     
  
  //# Mineralisation:
    //#---------------  
  
  Mindon = mindon*TemMinON*(DO/(DOsatdon+DO));
  
  
  //# ODE for DON:
    //# ------------
  
  double ONinflow;
  
  if ( 2500 <*t && *t< 3000)   ONinflow = 50000/48;  // 50 gram
  
  else ONinflow = 0;
  
  //double  PONDW= 0.25;
  
  dNDetritusdt = BMDON*BM*NCPhy*CPhytoplankton - Mindon*NDetritus - Setpon*NDetritus 
  + 0*Resus*1e-3 + 0*ONinflow/Vol + (1/H)*NCDW*mort*loss*DWeed /*0.05*loss*DWeed*/; // + DONinflow;
  
  
  dONdt = dNDetritusdt + NCPhy*dCPhytoplanktondt ;
  
  //#----------------------------------------     NH4     -------------------------------------------------------------     
  
  //# Nitrificaton, sediment fluxes and prerence term for Phytoplankton and macrophyte respectively:
    //#------------------------------------------------------------------------------
  
  // if (temp <= 26)   {TemNit  = exp(-minTNit*pow((temp - 26),2));}
  
  //           else    {TemNit  = exp(-maxTNit*pow((temp - 26),2));}
  
  
  Nitr    =   Nitrmax*(DO/(KhnitDo + DO))*TemNit;
  
  NH4sf  =   SedfluxNH4*TemFxNH4*(1/H);  
  //NH4sf = (SedfluxNH4/(Hs/2))*(poro/H)*Tem ;
  
  
  r       =   (NH4*NO3)/((NH4+Kn)*(NO3+Kn)) +  (NH4*Kn)/((NH4+NO3)*(NO3+Kn));
  
  
  //double psi = 3000;
  
  // r        =   1 - exp(- psi*NH4);
  
  
  
  //# ODE for NH4:
    //# ------------
  double NH4inflow;
  
  if(1000 <*t && *t<  1790)   NH4inflow = 3/48;  // 96 gram
  else NH4inflow = 0;
  
  dNH4dt  =  (BMNH4*BM - r*Gr)*NCPhy*CPhytoplankton + Mindon*NDetritus - Nitr*NH4 + NH4sf  + 0*NH4inflow/Vol;
  
  //#-----------------------------------------      NO3     ---------------------------------------------------------------
  
  //# Sediment fluxes:
    //#----------------
  
  NO3sf = SedfluxNO3*TemFxNO3*(1/H)   ;
  //NO3sf = (SedfluxNO3/(Hs/2))*(poro/H) ;
  
  //# ODE for NO3:
    //# ------------
  
  
  double Denit = denit*pow(Arrhdent,(temp - 20))*(KHoxresp/(KHoxresp+DO));
  
  double NO3inflow;
  
  if(1000 <*t && *t<  1790)   NO3inflow = 96000/48;  // 96 gram  from 20/03/2013  -- 06/04/2013
  else NO3inflow = 0;
  
  dNO3dt = - (1-r)*Gr*NCPhy*CPhytoplankton + Nitr*NH4 - Denit*NO3 + NO3sf + NO3inflow/Vol; //;
  
  
  
  dTNdt =  dNO3dt + dNH4dt + dONdt  ;
  //#-----------------------------------------      POP     ---------------------------------------------------------------
  
  //# Hydrolysis and Settling:
    //#------------------------
  
  
  //HYDpop = hydpop*TemHy;                                                   
  
  Setpop = (1- FDOM)*(Vpop/H);         
  
  // 31*1/106*12
  
  //# ODE for POP:
    //# ------------
  
  // dPOPdt =  BMPOP*BM*PCPhy*CPhytoplankton  - HYDpop*POP - 0*Setpop*POP  /*+ Resus*POPSediment*/ ; //+ POPinflow;
  
  //#-----------------------------------------      DOP     ---------------------------------------------------------------
  
  //# Mineralisation:
    //#---------------
  
  Mindop = mindop*TemMinOP*(DO/(DOsatdop+DO)) ;                                      
  
  //# ODE for DOP:
    //# ------------
  
  //double  POPDW = 0.25;
  
  dPDetritusdt = BMDOP*BM*PCPhy*CPhytoplankton - Mindop*PDetritus 
  - Setpop*PDetritus + 0*Resus*1e-3  + (1/H)*PCDW*mort*loss*DWeed /*0.001*loss*DWeed*/; //+ DOPinflow;
  
  dOPdt = dPDetritusdt +  PCPhy*dCPhytoplanktondt ;
  
  //#-----------------------------------------      PO4     ---------------------------------------------------------------
  
  //# Sediment fluxes:
    //#-----------------
  
  TemFxPO4  = pow(ArrhFluxPO4,(temp - 20));
  Psf    = SedfluxP*TemFxPO4*(1/H)   ;
  //Psf = (SedfluxP/(Hs/2))*(poro/H);
  
  
  
  //#ODE for PO4:
    //#------------
  //if(1000 <*t && *t<  1790)   PO4inflow = 5e-4;
  //else PO4inflow = 0;
  
  dPO4dt =  (BMPO4*BM - Gr)*PCPhy*CPhytoplankton + Mindop*PDetritus + Psf - (Abso/H)*Vset*(Fe + Mn)*PO4 ; //+ PO4inflow;
  
  
  dTPdt =  dPO4dt + dOPdt ;
  
  //#-----------------------------------------      DO     ----------------------------------------------------------------
  
  //# DO consumed in sediment:
    //#-------------------------
  
  double ResDOC = 2.67;
  
  DOSed  =  SOD*pow(ArrhDO,(temp - 20))/H;                                           
  
  //#DO saturation:
    //#-------------
  
  
  //DOsat     = do1 - do2*temp + do3*pow(temp,2) - do4*pow(temp,3);
  DOsat    = - 0.000063*pow(temp,3) + 0.0054258*pow(temp,2) - 0.68217*temp + 10.5532;                                                                                                
  
  //10:291 0:2809T ˛ 0:006009T2 :000063T3
  
  //# ODE for DO:
    //# -----------
  
  /* dDOdt = (Rea/Hsurf)*(DOsat - DO) + 0*((1.3 - 0.3*r)*Gr - ((1-BMDOC)*(DO/(KHR + DO)))*BM)*ResDOC*CPhytoplankton +
    (1/H)*(GrVgs - RESVgs)*ResDOC*MacrophyteShoot - 
    Khr*ResDOC*DOC - DOSed - 4.57*Nitr*NH4 ;*/
    
    
    //if ( 2500 <*t && *t< 3000)   Rea = 0;  // 500 gram
  
  //else Rea = 0.02/48;
  
  
  dDOdt =  (Gr - BMDIC*BM)*ResDOC*CPhytoplankton + Rea*(DOsat - DO) - Min*ResDOC*DOC - DOSed - 4.57*Nitr*NH4 ;
  
  
  
  
}


