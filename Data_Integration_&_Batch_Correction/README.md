This folder contains the code demonstrating how to integrate all the 17 spatial transcriptomics (ST) datasets mentioned in the current study.
The ST datasets belong to the following experimental groups namely:
1. Cold Ischemia Injury (CIS)
2. Warm Ischemia-Reperfusion Injury (AKI) (female mice)
3. Native kidneys (CTRL)
4. 24-hours Warm Ischemia-Reperfusion Injury (AKI24) (male mice)
   Note: AKI24 is same as IRL group (as mentioned in the codes)

The code makes use of Harmony package to remove batch effects among the ST datasets and perform subsequent clustering which facilitates idenitification
of shared compartments compartments across all the 17 ST datasets namely: inner medulla, outer medulla and cortex.

Note: 1. inner medulla is referred to as medulla (within the code)
      2. outer medulla is referred to as interface (within the code)
