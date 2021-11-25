#=============================================================================
# Start code + packages
#=============================================================================

require(deSolve)

#=============================================================================
# State and parameter formulation
#=============================================================================
x <- 0

#compartments/state variables go here, units are mmolC/m3
state <- c(BACTERIA = x, FUNGI = x, A.BISPORUS = x, SUGARS = x, COMPOST = , CO2 = x)

#parameters that comprise rate laws of processes go here
parms <- c(
  k20.bac.mort * BACTERIA
  k21.fun.mort * FUNGI
  k22.bisp.mort * A.BISPORUS
  k11.bac.uptake * SUGARS/(SUGARS+kSUGARS.bac) * BACTERIA * kBISPORUS/(A.BISPORUS+kBISPORUS)
  kSUGARS.bac = 
  kBISPORUS = 
  k12.fun.uptake * SUGARS/(SUGARS+kSUGARS.fun) * FUNGI * kBISPORUS/(A.BISPORUS+kBISPORUS)
  kSUGARS.fun = 
  k13.bisp.uptake * SUGARS/(SUGARS+kSUGARS.bisp) * A.BISPORUS * (1-A.BISPORUS/MAX.A.BISPORUS)
  kSUGARS.bisp = 
  MAX.A.BISPORUS = 
  k14.fun.killing.bac * BACTERIA * FUNGI
  k15.bisp.killing.bac * BACTERIA * A.BISPORUS
  k16.bisp.killing.fun * FUNGI * A.BISPORUS
  k1.deg.bac * BACTERIA
  k2.deg.fun * FUNGI
  k3.deg.bisp * A.BISPORUS
  Ef.fun = 
  Ef.bac = 
  Ef.bisp = 
)

#=============================================================================
# Model formulation
#=============================================================================

WHITEBUTTON <- function(t, state, parameters) { 
  with(as.list(c(state, parameters)),{
 
  # Rate expressions - all in units of [mmolC/m3/day] go here
      R20.bac.mort <- k20.bac.mort * BACTERIA
      R21.fun.mort <- k21.fun.mort * FUNGI
      R22.bisp.mort <- k22.bisp.mort * A.BISPORUS
      R11.bac.uptake <- k11.bac.uptake * SUGARS/(SUGARS+kSUGARS.bac) * BACTERIA * kBISPORUS/(A.BISPORUS+kBISPORUS)
      R12.fun.uptake <- k12.fun.uptake * SUGARS/(SUGARS+kSUGARS.fun) * FUNGI * kBISPORUS/(A.BISPORUS+kBISPORUS)
      R13.bisp.uptake <- k13.bisp.uptake * SUGARS/(SUGARS+kSUGARS.bisp) * A.BISPORUS * (1-A.BISPORUS/MAX.A.BISPORUS)
      R14.fun.killing.bac <- k14.fun.killing.bac * BACTERIA * FUNGI
      R15.bisp.killing.bac <- k15.bisp.killing.bac * BACTERIA * A.BISPORUS
      R16.bisp.killing.fun <- k16.bisp.killing.fun * FUNGI * A.BISPORUS
      R1.deg.bac <- k1.deg.bac * BACTERIA
      R2.deg.fun <- k2.deg.fun * FUNGI
      R3.deg.bisp <- k3.deg.bisp * A.BISPORUS
      
  # Mass balances [molN/m3/day] go here
      dSUGARS     <- R20.bac.mort + 
                     R21.fun.mort + 
                     R22.bisp.mort - 
                     R11.bac.uptake - 
                     R12.fun.uptake - 
                     R13.bisp.uptake + 
                     R14.fun.killing.bac + 
                     R15.bisp.killing.bac + 
                     R16.bisp.killing.fun +
                     R1.deg.bac + 
                     R2.deg.fun + 
                     R3.deg.bisp
      
      dBACTERIA   <- Ef.bac*R11.bac.uptake - 
                     R15.bisp.killing.bac - 
                     R20.bac.mort - 
                     R14.fun.killing.bac
      
      dFUNGI      <- Ef.fun*R12.fun.uptake - 
                     R16.bisp.killing.fun - 
                     R21.fun.mort
      
      dA.BISPORUS <- Ef.bisp*R13.bisp.uptake - 
                     R22.bisp.mort
      
      dCO2        <- (1-Ef.bac)*R11.bac.uptake + 
                     (1-Ef.fun)*R12.fun.uptake + 
                     (1-Ef.bisp)*R13.bisp.uptake
      
      dCOMPOST    <- -R1.deg.bac - R2.deg.fun - R3.deg.bisp
      
      TOTAL_C     <- COMPOST + BACTERIA + SUGARS + FUNGI + A.BISPORUS + CO2
      
  # Model output goes here
      return (list(c(dBACTERIA, dFUNGI, dA.BISPORUS, dSUGARS, dCOMPOST, dCO2),
                     TOTAL_C = TOTAL_C)
              )
      })
}

#=============================================================================
# Model run
#=============================================================================

outtimes     <- seq(from = 0, to = 365, length.out = 1000)

model_output <- ode(y = state, parms = parms, func = WHITEBUTTON, times = outtimes)

plot(model_output, mfrow=c(2,3))

#=============================================================================
# End code
#=============================================================================
