# axonalBARseq
These codes are for mouse 65A

# Structure of the codes & analysis
<img src='https://user-images.githubusercontent.com/60980561/167521108-25124b75-708b-4e21-b69b-089ad1110f1b.png' width='600'>


## TBS
Necessary, contains functions for all scripts


## imProcess_5d
Image processing (i.e. max projection) for both sequencing and antibody data.

Output image available, no need to rerun.

Note: not recommend, updated script available.


## abStitching_2
Stithing DAPI/antibody stanining whole brain images.

Output image available, no need to rerun.

Note: not recommend, updated script available.


## alignmentTest27
Alignment and basecalling for every batch of experiments. 

Output variable available, no need to refun.

Output data need to be cleaned for the downstream analysis.


## register2Volume

Registration experiments to a 3-D volume.
Contains two parts:
  1. Register individual areas into whole brain coronal sections using DAPI channel.
  2. Register individual brain section into a 3-D volume.


## combineBC_8C
Construct codebook and lookup table, exclude nonspecific barcodes.


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
