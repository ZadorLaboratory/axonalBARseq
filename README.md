# axonalBARseq

biorxiv: https://doi.org/10.1101/2023.02.18.528865

Note:   
  1. The mouse used in the dataset is called EF65A or 65A.
  2. The numbers in the file name are version IDs. 
  3. The processing and analysis pipelines are still under active development. Major changes are likely for the future axonal BARseq experiments.
  4. Questions and comments are welcome! 
  
  (Ongoing) deposit raw and processed data to public database


# Structure of the codes & analysis
<img src='https://user-images.githubusercontent.com/60980561/167521108-25124b75-708b-4e21-b69b-089ad1110f1b.png' width='600'>


## TBS
Necessary, contains functions for all scripts


## imProcess_5d
Image processing (i.e. max projection) for both sequencing and antibody data.

Processed image available, no need to rerun.


## abStitching_2
Stithing DAPI/antibody stanining whole brain images.

Processed image available, no need to rerun.


## alignmentTest27
Alignment and basecalling for every batch of experiments. 

Output variables available, no need to refun.

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
  2. Visualize data in registered coordinates/flatmap (contains code for Fig 1E, SupFig 4G-I)
       
## analysis_65A_2
Script for analysis of 65A data
