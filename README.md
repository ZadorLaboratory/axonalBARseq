# axonalBARseq
These codes are for mouse 65A


## TBS
Necessary, contains functions for all scripts

## imProcess_5d
Image processing (i.e. max projection) for both sequencing and antibody data.

Output image available, no need to rerun.

Note: not recommend, updated script available.

## alignmentTest27
Script alignment and basecalling for every batch of experiments. 

Output variable available, no need to refun.

Output data need to be cleaned for the downstream analysis.

## combineBC_8C

## abStitching_2
Stithing DAPI/antibody stanining whole brain images.

## (Optional) flatmap_2
Construct flatmap using Allen CCF.

Output: ctxAP/MP/DepthPrctile; can be directly used for 65A, no need to rerun.

## register2Ref_3
Conatains two parts:
  1. register 3D volume to Allen CCF (Optional)
       Output: ref2AllenTform.mat; can be directly used for 65A, no need to rerun.
       Note: angle adjustment procedure is not good, will have big modification for future experiments.
  2. Visualize data in registered coordinates/flatmap (contains code for Fig 1E, SupFig 4E-K)
       
## analysis_65A_2
Script for analysis of 65A data
