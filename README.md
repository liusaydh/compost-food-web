# compost-food-web
Mushroom compost food web modelling in R, supplemented by nematode analysis.

<h2> Natural abundance (carbon flow based) food web model </h2>

This is the most basic extant model for an "artificial" mushroom compost model. Artificial, due to the fact that it is a food web not commonly found in nature, but rather an achievement of dozens and hundreds of years of selective breeding and composting to achieve perfect conditions in which to raise _Agaricus bisporus_ to a commercially viable extent. 

As an artificial environment, specifically created to maximally propagate a single living organism, this food web is likewise truncated in comparison to any natural counterpart. The microbiota (pre-extant bacteria and fungi) are altered through processes of casing and aerating, changing the humidity and temperature and affecting entire generations and types of microbial fauna present. Likewise, some larger critters are not present (centipedes, mites, beetles, snails), as well as entire populations of common pests (aphids, flies). 

<p align="center">
  <img 
    width="800"
    src="flowcharts/diagramMar22.png"
  >
</p>

As such, the food web comes down to several factors: fungal populations, bacterial populations, _A. bisporus_ populations, nematode populations, recalcitrant sources of carbohydrates and lignin, readily available sources of broken down sugars, and carbon dioxide respired.

In its first functional form, the food web is modeled to exclude nematode populations, as they add an additional layer of higher trophic predatory complexity as both first and second level consumers. This documentation serves to both explain the model and serve as a repository of parameters. Some of the parameters are relied upon from (scarce) literature on the subject, while others are heavily inferred from a combination of literature sources or inferred from the outputs of the model itself. The table below serves to denote these sources and any comments on parameters found and their overall reliability.

<h2> <sup>13</sup>C model addition </h2>

An extension of this basic model comes from labelling experiments, where labelled, or <sup>13</sup>C carbon is added, as a one-off injection in the form of sugars, to the pool of what is otherwise predominantly naturally abundant <sup>12</sup>C monosaccharide carbon. The surplus <sup>13</sup>C starts to be uptaken by the living biomass in the compost.

<p align="center">
  <img 
    width="800"
    src="flowcharts/injection.png"
  >
</p>

As the time of the experiment progresses, surplus <sup>13</sup>C accumulates in the biomass of the living organisms that uptake it. For the second trophic step, living organisms that eat other living organisms incorporate their biomass labelling, as well.

<p align="center">
  <img 
    width="800"
    src="flowcharts/incorporation.png"
  >
</p>

Finally, surplus <sup>13</sup>C gets respired, incorporated and accumulates further in the biomass.

<p align="center">
  <img 
    width="800"
    src="flowcharts/spread.png"
  >
</p>

<h2> Natural abundance parameters and state variables </h2>

Parameter | Description | Unit | Value | Source(s) | Comments |
--- | --- | --- | --- |--- |--- |
Ef.bac | fraction of C assimilated by bacteria into their biomass as a result of metabolisation of sugars | (-) | 0.3 | Sinsabaugh, 2016. + Aanderud, 2018. | max. 30% is assimilated in biomass for soil dwellers, and is based on their C:N biomass ratio (around 8)
Ef.fun | fraction of C assimilated by fungi into their biomass as a result of metabolisation of sugars | (-) | 0.3 | Sinsabaugh, 2016. + Aanderud, 2018. | should be less than bacteria though, up to 10 times less 
Ef.bisp | fraction of C assimilated by _A. bisporus_ into their biomass as a result of metabolisation of sugars | (-) | 0.016 - 0.3 | Krakowska, 2021., Sinsabaugh, 2016. + Aanderud, 2018. |
kBISPORUS | Monod rate constant for _A. bisporus_ biomass | mmol C per cubic meter | UNKNOWN | (to be) calibrated |  
MAX.A.BISPORUS | maximum possible unlimited _A. bisporus_ biomass growth | mmol C per cubic meter | UNKNOWN | (to be) calibrated |
kSUGARS.bac | bacterial growth limitation due to sugars availability (C.INI * %) | mmol C per cubic meter | 0.1 - 0.2 | Vîtă (van Dam), 2022. | check star conditions in Femke's thesis
kSUGARS.fun | fungal growth limitation due to sugars availability (C.INI * %) | mmol C per cubic meter | 0.1 - 0.2 | Vîtă (van Dam), 2022. | check star conditions in Femke's thesis    
kSUGARS.bisp | _A. bisporus_ growth limitation due to sugars availability (C.INI * %) | mmol C per cubic meter | UNKNOWN | (to be) calibrated | expected similar to other fungi 
k1.deg.bac | degradation of recalcitrant materials from compost into sugars by bacteria | /d | 0.05 - 0.15 | Soares & Rousk, 2019. + Vîtă (van Dam), 2022. | 
k2.deg.fun | degradation of recalcitrant materials from compost into sugars by fungi | /d | 0.05 - 0.15 | Soares & Rousk, 2019. + Vîtă (van Dam), 2022. | fungi are expected to be better degraders than bacteria, in this context 
k3.deg.bisp | degradation of recalcitrant materials from compost into sugars by _A. bisporus_ | /d | 0.075 - 0.13 | Andlar et al. 2018. | values found correspond to fungi whose active enzymes are incredibly similar in activity to _A. bisporus_
k4.bac.uptake | maximum uptake rate of carbohydrates by bacteria | /d | 0.04 - 0.1 | Bore et al. 2017. |
k5.fun.uptake | maximum uptake rate of carbohydrates by fungi | /d | 0.01 - 0.05 | Bore et al. 2017. |
k6.bisp.uptake | maximum uptake rate of carbohydrates by _A. bisporus_ | /d | UNKNOWN | (to be) calibrated | expected similar to other fungi
k7.fun.killing.bac | maximum predation rate constant for fungi predating on bacteria | /d | UNKNOWN | (to be) calibrated |
k8.bisp.killing.bac | maximum predation rate constant for _A. bisporus_ predating on bacteria | /d | UNKNOWN | (to be) calibrated |
k9.bisp.killing.fun | maximum predation rate constant for _A. bisporus_ predating on fungi | /d | UNKNOWN | (to be) calibrated |
k10.bac.mort | bacterial linear mortality rate constant | (mmol C per cubic meter)/d | 0.24 - 0.72 | Servais et al., 1985. | in water environments |
k11.fun.mort | fungal linear mortality rate constant | (mmol C per cubic meter)/d | 0.01 - 0.1 | Lamour, 2002. | generally slower than bacteria |
k12.bisp.mort | _A. bisporus_ linear mortality rate constant | (mmol C per cubic meter)/d | 0 | Koch, 1958. | death under suitable growth conditions is uncommon |

In addition, the system has starting state variables which are conditions that are in place the moment the model is ran. These most often refer to the amounts of biomass at the beginning of the run and are viewable in the table below.

State variable | Description | Unit | Value | Source(s) | Comments |
--- | --- | --- | --- |--- |--- |
BACTERIA | biomass of bacteria present at Phase III start | (-) | 1.4% - 2.2% of compost | Vos, 2017., Vîtă, 2022. |
FUNGI | biomass of fungi (excluding _A. bisporus_) present at Phase III start | (-) | 1.77% of compost | Vos, 2017. |
A.BISPORUS | biomass of _A. bisporus_ present at Phase III start | (-) | 0.01% of compost | Vos, 2017., Vîtă, 2022. | technically, _A. bisporus_ biomass is zero because only he inoculate on rye grain is present, not the mycelium
SUGARS | amount of sugars present at Phase III start as available monosaccharides | (-) | 26% of compost (w/w) | Jurak, 2013. |
COMPOST | amount of total carbon present at Phase III start, mostly as undegraded polysaccharides | grams | 12 - 15 | Vîtă, 2022. | 30% carbon in 45 grams compost
CO2 | amount of carbon dioxide present at Phase III start as a result of respiration | (-) | 0 | Vîtă (van Dam), 2022. |
C.INI | amount of total carbon present in the system at Phase III start | grams | 12 - 15 | Vîtă, 2022. | same value as compost, used for easier calculation only




