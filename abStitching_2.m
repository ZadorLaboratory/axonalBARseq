% abStitching
% Stitiching for antibody/flourescent images, using combination of
% different channels
% Updated 03152021 LY

clear
clc

% Regular setting==========================================================
% TBS path
addpath('C:\Users\Li Yuan\MATLAB Drive\Script65A');

% Directory
directory.main = 'F:\Li Yuan\65A_Ab_corrected';
directory.imFolder = directory.main;%fullfile(directory.main,'Ab01');
directory.stitched = 'F:\Li Yuan\65A_AbImage';

% System setting ----------------------------------------------------------
sysSetting = TBS.getSysSetting;
% 1-sectionElement; 2-seqElement; 3-regionElement
sysSetting.nameElement = 1:3;

% Image setting -----------------------------------------------------------
imageSetting = TBS.getImageSetting(sysSetting,[]);

% Aligment setting --------------------------------------------------------
alignmentSetting = [];
% Channel sequence for alignment
% Channel combination: allow use multipe channel per alignment
% ie. [{1:3},{4}]: align 1 to 3 at the same time, the tile doesnt works use
% channel 4
alignmentSetting.ch = [{1},{2},{3},{4}];

% Each channel combination can have its own methods
% Method (cell) for images (ie. @(X) max(X,[],3), @(X) MIN(X,[],3))
alignmentSetting.method = {};

% Prepend of image for stitching
alignmentSetting.imAppend = '_pixelCorrected';

% Whether only do translaiton in imregcorr
alignmentSetting.translaitonOnly = false;

% Whether skip aligned image or redo everything
alignmentSetting.skipAlignedIm = true;

% Output name on disk
alignmentSetting.tformVarName = 'tformStitchTable';

% Padding for output image (for stitching)
% (i.e for tileSize*2, padding will be one tile size for each side)
alignmentSetting.paddingSize = imageSetting.tileSize.*2;

% Output table ------------------------------------------------------------

if ~exist(directory.stitched)
    mkdir(directory.stitched);
end

% Load the output table if it exist
cd(directory.stitched);
if ~exist([alignmentSetting.tformVarName,'.mat'])
    tformStitchTable = table();
else
    load([alignmentSetting.tformVarName,'.mat']); 
end

%% Get tform for stitched image

tformStitchTable = TBS.stitchImInfo(tformStitchTable,alignmentSetting,...
    imageSetting,sysSetting,directory);

%% Fuse image

tformStitchTable = TBS.fuseImage(tformStitchTable,alignmentSetting,...
    imageSetting,sysSetting,directory);



