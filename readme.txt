
Authors: Anna Melnykova, Irene Tubikanec
Date:   2023-08-18

Description: This code has been developed jointly by Anna Melnykova and Irene Tubikanec. It provides an implementation of the algorithm MMLH, proposed in the paper:
         
Granger Causal Inference in Multivariate Hawkes Processes by Minimum Message Length, 
by K. Hlavackova-Schindler, A. Melnykova and I. Tubikanec

In particular, it reproduces the estimation results of MMLH-u and MMLH-e in the "Sparse Setting 1 - Cascade structure" for p=7 and T=200 (see Table 1 in the paper) and can be easily adapted to other settings.

How does the code work?

1. Install relevant packages (in particular, see the top of the R-file main_MML_cascade).

2. Adjust the section "Choose setting" (e.g., set the number of trials N, the dimension p and the time horizon T). 

3. Run the R-file main_MML_cascade (Algorithm MMLH).
After termination of the algorithm, the obtained F1-scores (averaged over the N trials) are shown.

Code description:
For a detailed description of the code, please consider the respective files.

Licence information:
Please consider the txt-files LICENCE and COPYING.GPL.v3.0.




