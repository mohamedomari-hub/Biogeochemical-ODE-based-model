

parms <- c(
  
  #Phytoplankton
  #--------------
  
  
  Kn = 0.065, Kp = 0.006, Kc = 0.7, Iop = 100., Gmax = 1.2/48,  bm = 0.06/48, ArrhBM = 1.02, 
  
  Vcph=0.02/48, ArrhGr= 1.001, NCPhy =0.17, PCPhy = 0.024,     ArrhFluxDOC=1.05, 
  
  ArrhHyd= 1.08, ArrhMinDOC=1.001, k0=0.01, k1=0.02, k2=0.02, k3=0.01, k4= 0.16, 
  
  
  
  #Dweed Parameters:
  #-----------------
  
  #alpha=0.00432, beta=0.03
  
  Isat=173.6, alpha=0.15, beta=0.03, lossweed = 0.15/48, Arrhloss = 1.01,  Tmin = 5, 
  
  KnDWeed = 0.005, KpDWeed = 0.0043, KBiomass =  26, GmaxDWeed = 0.41/48, resDW= 0.12/48, ArrhresDW = 1.02,
  
  
  #POC Parameters:
  #--------------------------
  
  hydpoc = 0.02/48,  Vpoc = 0.6/48, 
  
  #DOC Parameters:
  #--------------------------
  
  mindoc = 0.0057/48,  DOsatdoc = 1.5,   SedfluxDOC = 0.025/48 , 
  
  Krefrespdoc = 0.0024/48, KHoxresp=0.1,
  
  #DIC:
  #----
  
  K = (12*0.015*1e-3)/48, SedfluxDIC = 0.02/48, 
  
  #PON Parameters:
  #--------------------------
  
  Vpon=0.3/48,  
  
  #DON Parameters:
  #--------------------------
  
  mindon=0.01/48,  DOsatdon = 0.5,  
  
  #NH4 Parameters:
  #--------------------------
  
  SedfluxNH4 = 0.01/48 , ArrhNit= 1.03, Nitrmax= 0.1/48, KhnitDo = 2.0, 
  
  #NO3 Parameters:
  #--------------------------
  
  denit = 0.06, Arrhdent=1.08, SedfluxNO3 = 0.005/48,
  
  
  
  #POP Parameters:
  #--------------------------
  
  Vpop=0.2/48,  
  
  #DOP Parameters:
  #--------------------------
  
  mindop= 0.02/48, DOsatdop = 0.5, 
  
  
  #PO4 Parameters:
  #--------------------------
  
  SedfluxP = 0.04/48 , ArrhPflux=1.05 , Abso  = 0.2 , Vset = 0.3/48,
  
  #DO Parameters:
  #--------------------------
  
  SOD = 0.04/48, ArrhDO = 1.05,  Rea = 0.1/48, DOsatsod = 1.5,
  
  
  #Additional Parameters:
  #----------------------
  
  co = 8/48, ArrhFluxDIC=1.05, ArrhFluxNH4=1.05, ArrhFluxNO3=1.05, ArrhFluxPO4=1.15,
  
  CCDW=0.3, NCDW =0.03, PCDW = 0.0003, ArrhMinON=1.01, ArrhMinOP=1.01)

