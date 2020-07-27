# Publication
_In preparation_

This software accompanies the publication

**Multiscale modelling of germinal center plasma cell differentiation Germinal center multiscale model shows temporal switch from memory B cells to plasma cells with affinity-based CD40 signaling**. 

Elena Merino1,&, Danial Lashgari1,#, Rodrigo Garcia Valiente1, Xuefeng Gao2, Fabien Crauste2,3, Olivier Grandillon2, Philip A Robert4, Michael Meyer-Hermann4, Maria Rodriguez Martinez5, Marieke van Ham6, Jeroen Guikema7, Huub Hoefsloot8, Antoine H.C. van Kampen1,8,@

1 Bioinformatics Laboratory, Epidemiology and Data Science, Amsterdam Public Health research institute, Amsterdam Institute for Infection and Immunity Amsterdam, the Netherlands.

2 Laboratory of Biology and Modeling of the Cell, Ecole normale supérieure de Lyon, Lyon France.

3 Institute for Mathematics, University of Bordeaux, France

4 Department for Systems Immunology, Helmholtz Centre for Infection Research, Inhoffenstrasse 7, 38124 Braunschweig, Germany.

5 Zurich Research Laboratory, IBM, Zurich, Switzerland

6 Department of Immunopathology, Division Research, Sanquin, Amsterdam

7 Department of Pathology, Amsterdam Institute for Infection and Immunity Amsterdam, the Netherlands

8 Biosystems Data Analysis, Swammerdam Institute for Life Sciences, University of Amsterdam, Amsterdam, the Netherlands.

&co-first author

@Corresponding author: 
Antoine H.C. van Kampen
Bioinformatics Laboratory
Epidemiology & Data Science
Amsterdam University Medical Centers
Meibergdreef 9, 1105 AZ Amsterdam, the Netherlands
a.h.vankampen@amsterdamumc.nl
tel. +31-20-5667096


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

