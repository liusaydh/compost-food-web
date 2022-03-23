# compost-food-web
Mushroom compost food web modelling in R, supplemented by nematode analysis.

This is the most basic extant model for an "artificial" mushroom compost model. Artificial, due to the fact that it is a food web not commonly found in nature, but rather an achievement of dozens and hundreds of years of selective breeding and composting to achieve perfect conditions in wbich to raise _Agaricus bisporus_ to a commercially viable extent. 

As an artificial environment, specifically created to maximally propagate a single living organism, this food web is likewise truncated in comparison to any natural counterpart. The microbiota (pre-extant bacteria and fungi) are altered through processes of casing and aerating, changing the humidity and temperature and affecting entire generations and types of microbial fauna present. Likewise, some larger critters are not present (centipedes, mites, beetles, snails), as well as entire populations of common pests (aphids, flies). 

As such, the food web comes down to several factors: fungal populations, bacterial populations, _A. bisporus_ populations, nematode populations, recalcitrant sources of carbohydrates and lignin, readily available sources of broken down sugars, and carbon dioxide respired.

In its first functional form, the food web is modelled to exclude nematode populations, as they add an additional layer of higher trophic predatory complexity as both first and second level consumers. This documentation serves to both explain the model and serve as a repository of parameters. Some of the parameters are relied upon from (scarce) literature on the subject, while others are heavily inferred from a combination of literature sources or inferred from the outputs of the model itself. The table below serves to denote these sources and any comments on parameters found and their overall reliability.

Parameter | Description | Unit | Value | Source(s) | Comments |
--- | --- | --- | --- |--- |--- |
k10.bac.mort | bacterial mortality rate | /d | 0.24 - 0.72 | Servais et al., 1985. | in water environments |
k11.fun.mort | fungal mortality rate | 283 | /d | 0 | Koch, 1958. | death under suitable growth conditions is uncommon |
k12.bisp.mort | _A. bisporus_ mortality rate | /d | 0.01 | Lamour, 2002. | generally slower than bacteria, also due to predation |
