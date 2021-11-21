#=============================================================================
# Start code + packages
#=============================================================================

require(deSolve)

#=============================================================================
# State and parameter formulation
#=============================================================================

#compartments/state variables go here, units are mmolC/m3
state <- c(BACTERIA = x, FUNGI = x, A.BISPORUS = x, NEMATODES = x, COMPOST = x)

#parameters that comprise rate laws of processes go here
parms <- c(
  fluxTRACER    = x,          # [mmolC/m3]
  pResp_2       = x,          # [-]
  pResp_3       = x,          # [-]
  pResp_4       = x,          # [-]
  pResp_5       = x,          # [-]
  pResp_6       = x,          # [-]
  pResp_7       = x,          # [-]
  pResp_8       = x,          # [-]
  pResp_9       = x,          # [-]
  pResp_10      = x,          # [-]
  rupmax.B      = x,          # [/day]
  rupmax.F      = x,          # [/day]
  rupmax.AB     = x,          # [/day]
  kCOMPOST      = x,          # [mmolC/m3]
  kA.BISPORUS   = x,          # [mmolC/m3]
  kBACTERIA     = x,          # [mmolC/m3]
  maxPred.B_F   = x,          # [/day]
  maxPred.B_AB  = x,          # [/day]
  maxPred.F_AB  = x,          # [/day]
  maxPred.B_N   = x,          # [/day]
  maxPred.F_N   = x,          # [/day]
  maxPred.AB_N  = x,          # [/day]
  r.Mort.B      = x,          # [/(mmolC/m3)/day]
  r.Mort.F      = x,          # [/(mmolC/m3)/day]
  r.Mort.AB     = x,          # [/(mmolC/m3)/day]
  r.Mort.N      = x           # [/(mmolC/m3)/day]
)

#=============================================================================
# Model formulation
#=============================================================================

WHITEBUTTON <- function(t, state, parameters) { 
  with(as.list(c(state, parameters)),{
 
  # Rate expressions - all in units of [mmolC/m3/day] go here
      R1   <-   fluxTRACER
      R2   <-   pResp_2      * R11 
      R3   <-   pResp_3      * R12
      R4   <-   pResp_4      * R13 
      R5   <-   pResp_5      * R14 
      R6   <-   pResp_6      * R15 
      R7   <-   pResp_7      * R16 
      R8   <-   pResp_8      * R17 
      R9   <-   pResp_9      * R18 
      R10  <-   pResp_10     * R19 
      R11  <-   rupmax.B     * COMPOST/(COMPOST+kCOMPOST) * BACTERIA
      R12  <-   rupmax.F     * COMPOST/(COMPOST+kCOMPOST) * kA.BISPORUS/(A.BISPORUS+kA.BISPORUS) * FUNGI
      R13  <-   rupmax.AB    * COMPOST/(COMPOST+kCOMPOST) * kFUNGI/(FUNGI+kFUNGI) * A.BISPORUS
      R14  <-   maxPred.B_F  * BACTERIA/(BACTERIA+kBACTERIA) * kA.BISPORUS/(A.BISPORUS+kA.BISPORUS) * FUNGI
      R15  <-   maxPred.B_AB * BACTERIA/(BACTERIA+kBACTERIA) * kFUNGI/(FUNGI+kFUNGI) * A.BISPORUS
      R16  <-   maxPred.F_AB * FUNGI/(FUNGI+kFUNGI) * A.BISPORUS
      R17  <-   maxPred.B_N  * BACTERIA/(BACTERIA+kBACTERIA) * NEMATODES
      R18  <-   maxPred.F_N  * FUNGI/(FUNGI+kFUNGI) * NEMATODES
      R19  <-   maxPred.AB_N * A.BISPORUS/(A.BISPORUS+kA.BISPORUS) * NEMATODES
      R20  <-   r.Mort.B     * BACTERIA * BACTERIA
      R21  <-   r.Mort.F     * FUNGI * FUNGI
      R22  <-   r.Mort.AB    * A.BISPORUS * A.BISPORUS
      R23  <-   r.Mort.N     * NEMATODES * NEMATODES
    
  # Mass balances [molN/m3/day] go here
      dCOMPOST    = R1 + R20 + R21 + R22 + R23 - R11 - R12 - R13
      dBACTERIA   = R11 - R2 - R14 - R15 - R17 - R20
      dFUNGI      = R12 + R14 - R3 - R5 - R16 - R18 - R21
      dA.BISPORUS = R13 + R15 + R16 - R4 - R6 - R7 - R19 - R22
      dNEMATODES  = R17 + R18 + R19 - R8 - R9 - R10 - R23
      TOTAL_C     = COMPOST + BACTERIA + FUNGI + A.BISPORUS + NEMATODES            # [mmolC/m3]
      
  # Model output goes here
      return (list(c(dCOMPOST, dBACTERIA, dFUNGI, dA.BISPORUS, dNEMATODES),
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
