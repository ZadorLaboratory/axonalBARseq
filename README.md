# axonalBARseq
These codes are for mouse 65A


## TBS
Necessary, contains functions for all scripts

## (Optional) flatmap
Construct flatmap using Allen CCF.
Output: ctxAP/MP/DepthPrctile; can be directly used for 65A, no need to rerun.

## register2Ref
Conatains two parts:
  1. register 3D volume to Allen CCF (Optional)
       Output: ref2AllenTform.mat; can be directly used for 65A, no need to rerun.
       Note: angle adjustment procedure is not good, will have big modification for future experiments.
  2. Visualize data in registered coordinates/flatmap (contains code for Fig 1E, SupFig 4E-K)
       
## analysis_65A
Script for analysis of 65A data
