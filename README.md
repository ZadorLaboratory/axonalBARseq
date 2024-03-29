# Codes for axonal BARseq

biorxiv: https://doi.org/10.1101/2023.02.18.528865

Note:   
  1. The mouse used in the dataset is called EF65A or 65A.
  2. All the scripts were writen in MATLAB (i.e. 2021a), additioanl toolboxes are needed (parallel computing, image processing, curve fitting, statistics and machine learning, computer vision, bioinformatics toolbox).
  3. FIJI and MIJ (http://bigwww.epfl.ch/sage/soft/mij/) were used for data visualization.
  4. The numbers in the file names are version IDs. 
  5. The processing and analysis pipelines are still under active development. Major changes are likely for the future axonal BARseq experiments.
  6. Questions and comments are welcome! 
  
  (Ongoing) deposit raw and processed data to public database    
  
<br>

# Variables for quick testing
We included four variables (xyz*/mlapd*.mat) for those who want to take a quick look at the data:

* xyzDot and xyzSoma are the xyz-coordinates for rolonies (cell array) and soma (matrix) in CCFv3, with a voxel size of 25 µm.
* mlapdDot and mlapdSoma are the ML-AP-depth coordinates in the cortex for rolonies (cell array) and soma (matrix). The units for ML and AP (columns 1 and 2) are microns, and the unit for depth (column 3) is the percentage of depth within the cortical column.

Each row in these variables represents a barcoded cell, and you can locate the corresponding rolonies and soma locations using the row number. Please note that not all somas were identified due to technical issues (further details can be found in the paper), and these unidentified somas are represented as [0 0 0] in the soma variables. Additionally, only cortical rolonies were included in mlapdDot.

Example: the three signle-cell examples in Fig. 5, barcode (row) 505, 1, 431.

<br>

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
Construct flatmap using Allen CCF (25 um/voxel).

Output: ctxAP/MP/DepthPrctile; can be directly used for 65A, no need to rerun.

## register2Ref_3
Conatains two parts:
  1. register 3D volume to Allen CCF (Optional)
       Output: ref2AllenTform.mat; can be directly used for 65A, no need to rerun.
  2. Visualize data in registered coordinates/flatmap (contains code for Fig 1E, SupFig 4G-I)
       
## analysis_65A_2
Script for analysis of 65A data

## evaluateCCFalignment
Script for CCF-registration, by manually selecting edges of brain regions.

Output: roiCCFacc. roi selected for 65A-evaluation
