% Folder for the Allen nrrd data
directoryAllen = 'I:\AllenTest\AllenConnectivity';

% Setting for reference map
refSetting = TBS.getRefSetting(directory.main);
refScale = refSetting.refScale;

zFlatmap = 0:2.5:95;

% Minimum intensity to be considered as signal per ML-AP pixel
minInten = 0.5;

% Get cortex pixels on flatmap --------------------------------------------
ctxTF = any(layerFlat,3);

% Exclude edges
SE = strel('disk',9);
ctxTF = imerode(ctxTF,SE);

%% Plot CT/PT/WT bundle location ==========================================
% SupFig.7A-B

MIJ.run('Close All');

cd(directoryAllen);
im = [];
for i = 1:numel(fileName)
    iFile = fileName(i);
    [iIm,iFile] = TBS.getNrrdAllen(iFile);
    
    im(:,:,:,i) = TBS.correctInjSide(iIm);
    
    disp(['Got: ',num2str(iFile)]);
end

% Combine all the CT and PT registered images
% imPT = max(im,[],4);
% imCT = max(im,[],4);

% Manually draw a line on the same plate for profiling intensity along the
% line in ImageJ
return

MIJ.createImage(cat(3,imPT,imCT));
MIJ.run("Stack to Hyperstack...", "order=xyzct channels=2 slices=528 frames=1 display=Composite");

% % Bundle localization of WT mice
% MIJ.createImage(cat(3,imPT,imCT,im));
% MIJ.run("Stack to Hyperstack...", "order=xyzct channels=3 slices=528 frames=1 display=Grayscale");

%% (Stat) Visualization of PT/CT bundle location different ----------------
% stat, cell, with 2-3 column of PT/CT intensity profiling

% Input channel number
z = size(stat{1},2);

n = cellfun(@(X) size(X,1),stat);
n = max(n);
stat2 = [];
for i = 1:numel(stat)
    istat = [];
    for j = 1:z
        % Normalize to the same height
        istat(:,j) = imresize(stat{i}(:,j),[n,1]);
    end
    % Normalize max intensity to 1
    istat = max(istat,0);
    istat = istat./max(istat);
    
    stat2(:,i,:) = istat;    
end

stat2 = imresize(stat2,[100 280],'Method','nearest');

% Output in ImageJ
MIJ.createImage(stat2);
MIJ.run("Stack to Hyperstack...", strcat("order=xyzct channels=",...
    num2str(z)," slices=1 frames=1 display=Composite"));

%% SupFig. 10 Mannually select ROI for Lat projections ====================

MIJ.run('Close All');

cd(directoryAllen);
iFile = fileName(i+1);
[iIm,iFile] = TBS.getNrrdAllen(iFile);
iIm = TBS.correctInjSide(iIm);

i = i+1; disp(['Current image: ',iFile]);

% Get flapmap
flatIm = TBS.im2flatmapAllen(iIm,ctxML,ctxAP,ctxDepthPrctile,refScale);

% Correct for edge effect
flatIm = flatIm.*ctxTF;

% Depth: 2.5% per voxel, total 40
flatIm = flatIm(:,:,1:zFlatmap(end)/2.5,:);

% Get mirror flatmap for ROI selection
test = sum(flatIm,3);
test = cat(3,test,fliplr(test));
MIJ.createImage(test);
MIJ.run("Make Composite","Display Mode = Composite");
% Mannually choose ROI of homotopic projection

% Projection symmetry across hemisphere -----------------------------------
% Pause here: crop using roi
v = MIJ.getCurrentImage(); 
v = reshape(v,[],2);
% Exclude background pixel pairs 
TF = v >= minInten;
TF = any(TF,2);
v = v(TF,:);
% Correlation uses rank
[rho,pval] = corr(v(:,1),v(:,2),'type','Spearman')
tableOut.corr{iFile} = [rho,pval];

% Projection intensity across in deep layer --------------------------------------
MIJ.run('Close All');
l = layerFlat(:,:,1:zFlatmap(end)/2.5,:);
MIJ.createImage(l);
% Pause here: crop using roi
v = MIJ.getCurrentImage();
% threshold defined in analysis_65A
l = v < threshold & v > 0;
l = double(l);
l(v >= threshold) = 2;
lScale = [sum(l == 1,'all'),sum(l == 2,'all')];
lScale = lScale./sum(lScale);

% Projection intensity across in deep layer --------------------------------------
MIJ.run('Close All');
MIJ.createImage(flatIm);
% Pause here: crop using roi
v = MIJ.getCurrentImage();
v2 = [sum(v(l == 1),'all'),sum(v(l == 2),'all')];
% Normalization to voxel counts & %
v2 = v2./lScale;
v2 = v2./sum(v2).*100;
tableOut.deepProj{iFile} = v2

% save(fullfile(directory.main,'roiC_WT_Allen.mat'),'tableOut');

%% Mannual quantification of fiber bundles

MIJ.run('Close All');

cd(directoryAllen);
iFile = fileName(i+1);
[iIm,iFile] = TBS.getNrrdAllen(iFile);
iIm = TBS.correctInjSide(iIm);

i = i+1; disp(['Current image: ',num2str(iFile)]);

% ROI selection for fiber intensity measurement
iIm = iIm(85:285,80:236,:);
MIJ.createImage(iIm); MIJ.run('Fire');

% Mannually measure upper lower fiber ratio -------------------------------
% Upper fiber: enter thalamus; lower fiber: enter internal capsul
% Pause here: Mean * length
userIn = userIn(:,3).*userIn(:,7);
userIn = reshape(userIn,3,2);

userIn = mean(userIn,1);
tableOut.ulInten{iFile} = userIn
% save(fullfile(directory.main,'roiC_WT_Allen.mat'),'tableOut');

%% Stat -------------------------------------------------------------------
% File name of the region
if ~iscell(fileName)
    fileName = cellstr(num2str(fileName));
end
TF = tableOut.Properties.RowNames;
TF = cellfun(@(X) contains(TF,X),fileName,'UniformOutput',false);
TF = horzcat(TF{:});
TF = any(TF,2);

% Have correlation rho >= 0.5
TF2 = tableOut.corr;
TF2 = cellfun(@(X) X(1),TF2);
TF2 = TF2 > 0.5;
TF = TF & TF2;

% Measured fiber intensity
TF2 = tableOut.ulInten;
TF2 = ~cellfun(@isempty,TF2);
TF = TF & TF2;

tableOut2 = tableOut(TF,:);

% Correlation -------------------------------------------------------------
% Upper/all fiber intensity
x = cellfun(@(X) X(1)/sum(X).*100,tableOut2.ulInten);
% > 60% in contralateral side
y = cellfun(@(X) X(2),tableOut2.deepProj);   
[x,y]

%% SupFig. 8. Med/Lat bulk projection in Allen
% Description: Quantify the Med/Lat intensity along depth/layer

% Settings
zL = 1:0.25:6;  %'MedLatAUD_WT_Allen.mat'

% Binned layerFlat, with bin number ---------------------------------------
layerFlatBin = zeros(size(layerFlat));
for i = 1:numel(zL)-1
    TF = layerFlat >= zL(i) & layerFlat < zL(i+1);
    layerFlatBin(TF) = i;
end
layerFlatBin = layerFlatBin(:,:,1:rolonyDepthEdges(end)/2.5,:);

% ML-Boundaries (edges) pixel locations -----------------------------------
mlBoundary2 = [-mlBoundary 0 mlBoundary]';
mlBoundary2(:,2:3) = repmat([nan nan],3,1);

% mlapd2flatmapXYZ(mlapd,refSetting,flatmap)
mlBoundary2 = TBS.mlapd2flatmapXYZ(mlBoundary2,refSetting,ctxTF);
mlBoundary2 = mlBoundary2(:,1)'; % 195   335   475

% AP-Boundaries (edges) pixel locations -----------------------------------
flatmapYLim2 = nan(2,3);
flatmapYLim2(:,2) = flatmapYLim';

% mlapd2flatmapXYZ(mlapd,refSetting,flatmap)
flatmapYLim2 = TBS.mlapd2flatmapXYZ(flatmapYLim2,refSetting,ctxTF);
flatmapYLim2 = round(flatmapYLim2(:,2))';
flatmapYLim2 = flatmapYLim2(1):flatmapYLim2(2);  % 255   414

%% % ----------------------------------------------------------------------
% Need to mannually copy the file name
cd(directoryAllen);
iFile = fileName(i+1);
[iIm,iFile] = TBS.getNrrdAllen(iFile);
% make the injection side on the left, basing on intensity
iIm = TBS.correctInjSide(iIm);

i = i+1; disp(['Current image: ',iFile]);

% Get flapmap
flatIm = TBS.im2flatmapAllen(iIm,ctxML,ctxAP,ctxDepthPrctile,refScale);
% Correct for edge effect
flatIm = flatIm.*ctxTF;

% Depth: 2.5% per voxel, total 40
flatIm = flatIm(:,:,1:rolonyDepthEdges(end)/2.5,:);

% Exclude pixels near saturated pixels, 8.*25 um-200 um range
injSite = any(flatIm == 1,3);
injSite(:,mlBoundary2(1):end,:) = false;
flatIm = injSiteExclusion(flatIm,injSite);

% Try to eliminate blood vessles/noise using imreconstruction -------------
TF = sum(flatIm,3);
% Thresholding on sum
TF = TF >= minInten;
TF = flatIm.*TF;
flatIm0 = flatIm;
flatIm = imreconstruct(TF,flatIm);
% % (Check point)
% MIJ.createImage(cat(3,sum(flatIm0,3),sum(flatIm,3)));
MIJ.createImage(sum(flatIm(flatmapYLim2,:,:),3)); ij.IJ.setMinAndMax(0,3);

% Sum of intensity along depth in LatI/MedI/MedC/LatC ---------------------
[intenDepth,intenLayer] = profProjAUD(flatIm(flatmapYLim2,:,:),...
    mlBoundary2,layerFlatBin(flatmapYLim2,:,:));

tableOut.intenDepth{iFile} = intenDepth;
tableOut.intenLayer{iFile} = intenLayer;

% save(fullfile(directory.main,'MedLatAUD_WT_Allen.mat'),'tableOut');
% save(fullfile(directory.main,'MedLatAUD_Cre_Allen.mat'),'tableOut');

%% Function:    injSiteExclusion
% Description:  injection site and nearby region exclusion on flatmap
function flatIm = injSiteExclusion(flatIm,injSite)
% Input & output: flatIm, mat, flat map stack
%           injectSide, mat/logical, with injection site as non-zero

% ~200 um radius for 25 um/voxel
SE = strel('disk',8);

injSite = imdilate(injSite,SE);
flatIm = flatIm.*(~injSite);

end

%% Function:    profProjAUD
% Description:  profile projection intensity along the depth & layer for
% AUD projection
function [intenDepth,intenLayer] = profProjAUD(flatIm,mlBoundary2,layerFlat2)
% Input:    flatIm, mat, image stack
%           mlBoundary2, vector, left and right boundaries for Lat/Med in
%           each hemisphere, in pixel space
%           layerFlatBin, mat, layer position on flatmap, binned
% Output:   intenDepth, mat, sum intensity along depth, in Lat and Med (from left to right)
%           intenLayer, mat, sum intensity along layer

mlBoundary2 = [0 mlBoundary2 size(flatIm,2)];

intenDepth = {}; intenLayer = {};
for i = 1:4
    col = [mlBoundary2(i)+1,mlBoundary2(i+1)];
    col = col(1):col(2);
    
    iIm = flatIm(:,col,:);
    intenDepth{1,i} = sum(iIm,1:2);
    
    iLayer = layerFlat2(:,col,:);    
    TF = iLayer > 0;
    intenLayer{1,i} = accumarray(iLayer(TF),iIm(TF));   
end

intenDepth = cellfun(@squeeze,intenDepth,'UniformOutput',false);

intenDepth = cell2mat(intenDepth);
intenLayer = cell2mat(intenLayer);

end
