# Owners
Antoine van Kampen, Elena Merino Tejero

# Project
## MSM_PCdifferentiation

This repository includes code ONLY of the multiscale model of plasma cell diferentiation in Germinal Centers. The Agent-based
model is based on Mafalda, which is based on Hyphasma (e.g., Michael-Meyer Hermann, 2012).  
The GRN is based on Martinez et al., 2012. It is posible to simulate 2 scenarios: Plasma cell diferentiation based on 
Ag inheritance or on BLIMP1 level.

CD40 signal can be modeled 2 ways: Fixed CD40 (=50, 10) or CD40 dependent on affinity (=affinity*50). 
The switch from one to the other has to be done inside the code (for now).

## Software
All software is written in C++

## References
* Martínez, M. R., Corradin, A., Klein, U., Álvarez, M. J., Toffolo, G. M., di Camillo, B., … Stolovitzky, G. A. (2012). Quantitative modeling of the terminal differentiation of B * cells and mechanisms of lymphomagenesis. Proceedings of the National Academy of Sciences, 109(7), 2672–2677. 
Meyer-Hermann, M., Mohr, E., Pelletier, N., Zhang, Y., Victora, G. D., & Toellner, K. M. (2012). A theory of germinal center b cell selection, division, and exit. Cell Reports, 2(1), 162–174. 

