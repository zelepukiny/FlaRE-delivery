{ Tumour cord model using a Krogh cylinder. Bolus or short-time infusion DOX.
Solve for 1D time-varying drug release source
Treats as time-domain problem
AZ April 2021}

TITLE 'TumourCord'     
COORDINATES cylinder1("r")  
VARIABLES       
Cf   (threshold = 1e-3)  ! Concentration of free DOX in tissue, measured in ug/L
Cb  (threshold = 1e-3)  ! Concentration of bound DOX in tissue, measured in ug/L
Ci	  (threshold = 1e-3)  ! Concentration of free DOX in cells, measured in ng/10^5 cells

DEFINITIONS    
 

 A = 7.46e-2 		! Compartment 1 parameter, L-3
 Dose = 1.5e5		! Total dose injected, ug
 alpha = 0.161  	! Compartment 1 clearance rate, min-1
 Tinf = 180				! Period of infusion, min (bolus)
  
 Rc = 10 				! vessel radius, um
 Rt = 120 			! Tissue radius, um
 Df = 7000 			! diffusivity free DOX, um2/min
 Db = 530			! diffusivity albumin-bound DOX, um2/min
 Pf = 160				! diffusive permeability free DOX, um/min
 Pb = 0.5				! diffusive permeability bound DOX, um/min
 delta = 0.75		! Fraction of plasma DOX bound
 theta = 0.6			! fraction of blood that in plasma 
 ka = 50				! Free DOX-albumin binding rate, min-1
 kd = 17				! DOX-albumin dissociation rate, min-1
 
 KE = 219			! Michaelis constant, extracellular, ug/L
 KI = 1.37				! Michaelis constant, intracellular, ng/10^5 cells
 Vmax = 0.28		! Rate for transmembrane transport, ng/10^5 cells/min
 rhoc = 1e12		! Density of tumour cells, cells/L
 zeta = 1e-8			! scaling factor, ug (10^5 cells)/(ng cell)
 phi = 0.4				! Tumour fraction extracellular space
 ksi = 0.1				! fractional blood volume in the lungs
 iota = 0.8			! fractional MOF-DOX delivery to the lungs
 
 
 !kappa = 8e-3 ! drug release rate, s-1, [6.8e-5 for VP], [8e-3 for MOF]
 
 Cv =  IF ( t < Tinf) 
          THEN   Dose*A/alpha/Tinf* (1 - exp(-alpha* t))
          ELSE   Dose*A/alpha/Tinf* (exp(alpha*Tinf) - 1)* exp(-alpha*t) 
 
mu = Vmax*(Cf + Cb)/(Cf + Cb + KE*phi)
nu = Vmax*Ci/(Ci + KI)

! INITIAL VALUES

EQUATIONS        { PDE's, one for each variable }
  
  Cf : 	-div(Df*grad(Cf)) + dt(Cf) + rhoc*zeta*(mu - nu) + ka*Cf - kd*Cb = 0
  Cb:  	-div(Db*grad(Cb)) + dt(Cb) + kd*Cb - ka*Cf = 0
  Ci:		dt(Ci) = mu - nu
  
  
BOUNDARIES       { The domain definition }
 Region 1 "Tumour Cord"
 START(Rc) 		
 POINT LOAD(Cf) = -Pf*(Cv*theta*(1 - delta) - Cf)*iota/ksi
 POINT LOAD(Cb) = -Pb*(Cv*theta*delta - Cb)*iota/ksi
 LINE TO (Rt) 	
  
 TIME 0 TO 300    !  min

MONITORS   
FOR t = 0 BY 10 TO endtime
elevation(Cf) from (Rc) to (Rt) 
elevation(Cb) from (Rc) to (Rt)
elevation(Ci) from (Rc) to (Rt)

PLOTS
FOR t = 1 BY 20 TO endtime
 elevation(Cb) from (Rc) to (Rt) AS "Concentration Cb, ug/L" !export format "#x#b#1" file="C_Cb_3min_MOF.txt" 
 elevation(Ci) from (Rc) to (Rt) AS "Concentration Ci, ng/(10^5 cells)"   
 export format "#x#b#1" file="C_Ci_180min_MOF.txt"
  

HISTORIES
	! 	HISTORY(Cf) AT (Rt/2) AS "Concentration Cf Tumour, ug/L" 
    HISTORY(Cb) AT (Rt/2) AS "Concentration Cb Tumour, 50 um from capillary, ug/L" 
    HISTORY(Ci) AT (Rt/2) AS "Concentration Ci Tumour, 50 um from capillary, ng/10^5 cells"  !export format "#t#r,#i" file="Ci_t_3min_MOF.txt"
    HISTORY(mu/Vmax) AT (Rt/2) AS "mu/Vmax" 

END
