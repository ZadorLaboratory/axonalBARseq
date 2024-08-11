% Code for figures and analysis for mouse 65A
% Need 
%   bioinformative toolbox
%   nrrdread.m (https://www.mathworks.com/matlabcentral/fileexchange/34653-nrrd-format-file-reader)
%   MIJ (http://bigwww.epfl.ch/sage/soft/mij/)
%
% Other scripts:
%   rotating brain for fig. 2A and supVideo, 'rotatingBrain_12072023.m'
%   evaluation of alignment to CCF, 'evaluateCCFalignment.m'
%   Get data from Allen Mouse Brain Connectivity Atlas,
%   'fiberBundle_Allen_09282023.m'
% 
% Last update: 12182023 LY

clc

% Settings

% Directory with all files
directory.main = 'D:\Li Yuan';

sysSetting = TBS.getSysSetting;
imageSetting = TBS.getImageSetting(sysSetting,[]);

refSetting = TBS.getRefSetting(directory.main);

scaleFactor = 0.055; % (change way to set this)

mouseID = imageSetting.mouseID;
msBrainName = [mouseID,'_',num2str(scaleFactor),'.tif'];

% Depth/ML/AP lookup table
cd(directory.main);
load('ctxML.mat'); load('ctxAP.mat'); load('ctxDepthPrctile.mat');

% Edges for histogram along cortical depth: soma & rolony
binSize = 5;
somaDepthEdges = -binSize:binSize:100;

binSize = 1;
rolonyDepthEdges = -binSize:binSize:95;

% Boundaries for median & lateral region
mlBoundary = 3500;

% Colormaps ===============================================================
% Colormap for CT, PT
cCorn = [1 0.82 0; 0.25 0.75 0.1];

% Colormap for iTc, iTi, CT, PT
cmap4 = [0.9 0.32 0; 0  0.45 0.74];
cmap4(3:4,:) = cCorn;
cmap4 = cmap4([4 3 1 2],:);

% meganta & green
cmapMG = [235 0 139; 13 177 75]./255;

% Red, orange & blue
cmapROB = [0.9 0 0; 1 0.45 0; 0 0 0.9];
% for black background in ImageJ: cmapROB = [0.9 0 0; 1 0.45 0; 0 0.22 1];

% Annotation Info =========================================================

% Annotation map
annoMap = TBS.getRefMap('anno',refSetting);

% Annotation structure
annoStruct = refSetting.annoStruct;

% Default for 65A: CtxI, CtxC, Thal, Str, Mb
regVoxel = TBS.getRegVoxel(annoMap,annoStruct);

% Column number
regCol = [];
regCol.ctxI = 1;
regCol.ctxC = 2;
regCol.thal = 3;
regCol.str = 4;
regCol.mb = 5;

% Get annoMap with combined layers (level 11)------------------------------
% Annotated region in 3D
regionOutline = TBS.getAnnoRegion(refSetting);

% Get outline of annotated regions
regionOutline = TBS.annoRegionOutline(regionOutline);

% Layers ==================================================================
regName = {'layer 1';'layer 2';'layer 3';'layer 4'; 'layer 5'; 'layer 6'};
id = cellfun(@(X) TBS.findAnnoStrID(annoStruct,X),regName,'UniformOutput',false);
% Combine layer 2 & 3
id{2} = [id{2};id{3}]; id(3) = [];

layer = zeros(size(annoMap),'uint8');
for i = 1:numel(id)
    TF = ismember(annoMap,id{i});
    layer(TF) = i;
end

% Get layer on the flatmap ------------------------------------------------
layerFlat = TBS.getLayerFlatmap(layer,ctxML,ctxAP,ctxDepthPrctile,refSetting.refScale);
% Get gradience within each layer
layerFlat = TBS.layerFlatmapGradience(layerFlat);
% % (Check point)
% MIJ.createImage(layerFlat);

% Function: mlapd2flatmapV(mlapd,refSetting,flatmap,edgeCorrect);

% Current brain, DAPI stack ===============================================
moving = TBS.getStack(msBrainName,[]);

sz = size(moving);
moving = reshape(moving,sz(1),sz(2),4,[]);
% Get DAPI: channel 1
moving = moving(:,:,1,:);
moving = squeeze(moving);

% Transform to registrated voxels
load('reg2AllenTform.mat');
% Transformation matrix and displacement field
tform = reg2AllenTform.tform;
D = reg2AllenTform.D;

% Output size as annotation map
sz = size(annoMap);
R = imref3d(sz);

tfMoving = imwarp(moving,tform,'OutputView',R);
tfMoving = imwarp(tfMoving,D);

% Rolony & soma location coordinates ======================================
cd(directory.main);
load('codeBookref.mat');

% Rolony location ---------------------------------------------------------
load('axonBCref.mat'); 

xyzDot = TBS.BCtable2cell(axonBC,'xyzRef');

n = cellfun(@(X) size(X,1),xyzDot);
disp(['Total rolony number: ',num2str(sum(n))]);
disp(['Median number: ',num2str(median(n))]);

% Cortical ML/AP/Depth
% getMLAPD(xyz,ctxML,ctxAP,ctxDepthPrctile)
mlapdDot = TBS.getMLAPD(xyzDot,ctxML,ctxAP,ctxDepthPrctile);

% Only include rolony within the edges
fh = @(X) X(:,3)>= min(rolonyDepthEdges) & X(:,3)<= max(rolonyDepthEdges);
TF = cellfun(@(X) fh(X),mlapdDot,'UniformOutput',false);
mlapdDot = cellfun(@(X,Y) X(Y,:),mlapdDot,TF,'UniformOutput',false);

% Soma location -----------------------------------------------------------
load('somaLocationRef.mat');
xyzSoma = somaLocation;

TF = any(xyzSoma,2);
mlapdSoma = zeros(size(xyzSoma));
mlapdSoma(TF,:) = TBS.getMLAPD(xyzSoma(TF,:),ctxML,ctxAP,ctxDepthPrctile);

TF = any(mlapdSoma,2);
disp(['Soma included in analysis: ',num2str(sum(TF))]);
disp('08302022, two soma were excluded during flatmap registration.')

% (02142023, median distance (um) to the injection center)
TF = any(mlapdSoma,2);
injCenter = median(xyzSoma(TF,:));
D = pdist2(xyzSoma(TF,:),injCenter);
D = median(D);
% Convert from ccf to micron
D = D/refSetting.refScale;
disp(['Soma median distance to injection center (um): ',num2str(D)]);

% mlapd map settings ======================================================
% Area boundaries only include the middle layer

% Area boundaries on flatmap
[y,x,z,~] = TBS.find3(regionOutline);

regionOutlineFlat = TBS.getMLAPD([x,y,z],ctxML,ctxAP,ctxDepthPrctile);
% Only get the middle layers
TF = regionOutlineFlat(:,3) > 45 & regionOutlineFlat(:,3) <= 55;
regionOutlineFlat = regionOutlineFlat(TF,:);

% (Check point) 
figure; TBS.plotRegionRef(mlapdSoma,regionOutlineFlat);

% Other settngs
% AP limit for figure
flatmapYLim = vertcat(mlapdDot{:});
flatmapYLim = nonzeros(flatmapYLim(:,2));
flatmapYLim = [min(flatmapYLim), max(flatmapYLim)];

%% % Cell grouping ===========================================================

% Range from the injection center (median soma ML/AP location)
localRng = 1000;

% Min count to be defined as with projection
regionMinCount = [3 5];

% CT/PT -------------------------------------------------------------------
% Note: cannot use annotation to get the ROI, registration is not exactly
% like CCF. Therefore, ROI need to be defined mannually
% (Already thresholding for minimum rolony in thalamus, i.e t)

% Thal+ neurons
TF = cellfun(@(X) TBS.xyz2v(X,regVoxel{3}),xyzDot,'UniformOutput',false);
TF = cellfun(@sum,TF) >= max(regionMinCount);

% Group CT and PT basing on Thal Projection
CF = nan(size(TF));
% groupCTPT(xyzDot,thalReg,strReg,mbReg)
CF(TF) = TBS.groupCTPT(xyzDot(TF),regVoxel{3},regVoxel{4},regVoxel{5});

% thal-mb+ cells
TF = cellfun(@(X) TBS.xyz2v(X,regVoxel{5}),xyzDot,'UniformOutput',false);
TF = cellfun(@sum,TF) >= min(regionMinCount);
TF = isnan(CF) & TF;

disp(['Thal-MB+ cells: ', num2str(sum(TF))]);

CF(TF) = 1;
CT = CF == 0;
PT = CF == 1;
disp(['CF-cells: ', num2str(sum(~isnan(CF)))]);
disp(['PT-cells: ', num2str(sum(PT))]);
disp(['CT-cells: ', num2str(sum(CT))]);

% IT-cells ----------------------------------------------------------------
IT = isnan(CF);

disp(['IT-cells: ', num2str(sum(IT))]);

% ITc & ITi ---------------------------------------------------------------
% Whether has dots in the contra hemisphere
TF = cellfun(@(X) TBS.xyz2v(X,regVoxel{2}),xyzDot,'UniformOutput',false);
% With rolony more than minimum counts
TF = cellfun(@sum,TF);
TF = TF >= max(regionMinCount);
ITc = IT & TF;

disp(['ITc-cells: ', num2str(sum(ITc))]);

ITi = IT & ~ITc;
disp(['ITi-cells: ', num2str(sum(ITi))]);

% Cell type ---------------------------------------------------------------
% 1-4: PT, CT, ITc, ITi
cellType = zeros(size(xyzDot));
cellType(PT) = 1;
cellType(CT) = 2;
cellType(ITc) = 3;
cellType(ITi) = 4;

% (Report) CF with IT projection ------------------------------------------
% CFIT: CF cells with non-local IT projection

% Cells has long-range cortical rolony (>= regionMinCount)
TF = cellfun(@(X) ~TBS.isLocalProj(X,mlapdSoma,localRng) & X(:,1) ~= 0,...
    mlapdDot,'UniformOutput',false);
TF = cellfun(@(X) sum(X),TF) >= max(regionMinCount);

CFIT = TF & (CT | PT);

disp(['CF-cells project out of local area (CFIT): ', num2str(sum(CFIT))]);
disp(['CFIT CT counts: ', num2str(sum(TF & CT))]);
disp(['CFIT PT counts: ', num2str(sum(TF & PT))]);

%% Identify IT cells innvervate the same cortical cluster =================
% Output: cellITcluster, cell, BC number for each cortical cluster
%         sorted using projection center to local region ratio

% Mannually identify cluster of projection center
center = [6400 8300;6900 7600; 4850 8500];

figOutTF = true;

% Range from the projection center for within/surrounding area
rng = 300;

% Minimum barcode within & in the surrounding range
minCount = 10;

% Min ratio (cluster/local) to categorize a cell project to a cluster
minRatio = 75;

cellITcluster = {};
for i = 1:size(center,1)
    
    iCenter = center(i,:);
    
    % Whether rolony is within cluster/local area
    % Local: same size as the cluster
    % inCtxCluster(mlapdDot,center,rng,figOutTF,mlapdSoma,regionOutlineFlat)
    [clusterRng, localRng] = TBS.inCtxCluster(mlapdDot,iCenter,rng,...
        figOutTF,mlapdSoma,regionOutlineFlat);
    
    % Number of rolony within cluster or local area
    nCluster = cellfun(@sum,clusterRng);
    nLocal = cellfun(@sum,localRng);
    
    % within/local cluster ratio
    ratio = nCluster./nLocal.*100;
    
    % IT cell only
    TF = ismember(cellType,3:4);
    
    % Cell with more than min count in local area
    TF = TF & nLocal >= minCount;
    
    % (Check Point) Stat for histogram
    stat = ratio(TF);
    
    % Cells project to the cortical cluster
    TF = TF & ratio >= minRatio;
    
    disp(['Cortical cluster center: ', num2str(iCenter),...
        ', positive cell: ',num2str(sum(TF))]);
    
    % Sort the cells using ratio
    row = find(TF);
    ratio = ratio(TF);
    [~,I] = sort(ratio,'descend');
    cellITcluster{i} = row(I);
end

return

%% Fig 1. Connect dots ====================================================
% To show how the reconstruction were made

% Example neuron from ITc with the highest counts (or pick mannually)
row = cellType == 3 & any(mlapdSoma,2);
row = find(row);

n = cellfun(@(X) size(X,1),xyzDot(row));
[~,I] = max(n);
I = row(I);

xyzSoma2 = xyzSoma(I,:); xyzDot2 = xyzDot(I);
% Reconstruction
[im, outline] = TBS.BCmodel(xyzDot2,xyzSoma2,[],1,annoMap,refSetting);

% Rolony location
im2 = TBS.xyzv2im(size(im),xyzDot2{:},[]);
SE = strel('sphere',1);
im2 = imdilate(im2,SE);

im = single(im); im2 = single(im2);

% Show the whole brain here
MIJ.createImage(cat(3,im,im2,outline));
MIJ.run("Stack to Hyperstack...", "order=xyzct channels=3 slices=528 frames=1 display=Composite");

%% Fig 2. Plot BC neuron model ============================================
% Discription: Plot 100 most rolony & soma+ neurons from 4 groups

nCell = 25;

% BC id for plotting
id = [];
for i = 1:4
    row = cellType == i & any(mlapdSoma,2);
    row = find(row);
    
    % Sort using rolony count
    n = cellfun(@(X) size(X,1),xyzDot(row));
    [~,I] = sort(n,'descend');
    row = row(I,:);
    
    id(:,i) = row(1:nCell);
end

% Plot --------------------------------------------------------------------
% Shuffle order for random colorcode
id = reshape(id,[],1);
id = TBS.shuffleRows(id);

xyzSoma2 = xyzSoma(id,:); 
xyzDot2 = xyzDot(id);
% Use barcode id as color
% Make it bigger (2X ref) for finner details in image output
% [im, outline] = TBS.BCmodel(xyzDot,xyzSoma,c,scaleFactor,annoMap,refSetting)
[im, outline] = TBS.BCmodel(xyzDot2,xyzSoma2,[],2,annoMap,refSetting);

% Output In imageJ
MIJ.createImage(cat(3,im,outline));
MIJ.run("Stack to Hyperstack...", "order=xyzct channels=2 slices=1056 frames=1 display=Composite");

%% Fig 2. Simulated retrograd injection ===================================
% Get cells project to each cluster with minimum rolony number

% Minimum count with cluster range
minCount = 10;

% Maximum cells to plot
nCell = 25;

im = {};

% Cortical cluster to plot (i.e. 2: [7050, 6950])
for iCluster = 1:size(center,1)
    
    % Get rolony within the cluster
    iCenter = center(iCluster,:);
        
    [clusterRng, ~] = TBS.inCtxCluster(mlapdDot,iCenter,rng,...
        figOutTF,mlapdSoma,regionOutlineFlat);
    
    % Number of rolony within cluster or local area
    nCluster = cellfun(@sum,clusterRng);
    
    % IT cell only
    TF = ismember(cellType,3:4);
    
    % Cell with more than min count in local area
    TF = TF & nCluster >= minCount;

    row = find(TF);
    
    % Soma+ cells
    TF = mlapdSoma(row,:);
    TF = any(TF,2);
    row = row(TF);
    
    disp(['Cells meet the criteria: ', num2str(size(row,1))]);
    % continue
    
    % Randomly select N cells to plot
    if numel(row) > nCell
        row = TBS.shuffleRows(row);
        row = row(1:nCell);
    end
        
    % Plot BC model -----------------------------------------------------------
    % Shuffle order for random colorcode
    row = TBS.shuffleRows(row);
    xyzSoma2 = xyzSoma(row,:); xyzDot2 = xyzDot(row);
    [im{iCluster}, outline] = TBS.BCmodel(xyzDot2,xyzSoma2,[],2,annoMap,refSetting);
end

% cell2stack(cellArray,dimension4stackCell)
im = TBS.cell2stack(im,3);
im = uint8(im);

MIJ.createImage(cat(3,im,outline));
MIJ.run("Stack to Hyperstack...", "order=xyzct channels=4 slices=1056 frames=1 display=Composite");

%% Fig 2. Plot cells solely inervate single cortical cluster ==============
% 06112022
% Uses a subset of cells for axon model
% Plot cells projects to a cluster, and most of the projection in the
% hemisphere is within the cluster; also soma-positive

im = {};

% Cortical cluster to plot (i.e. 2: [7050, 6950])
for iCluster = 1:size(center,1)
    
    % Get rolony within the cluster
    iCenter = center(iCluster,:);
    [clusterRng, ~] = TBS.inCtxCluster(mlapdDot,iCenter,rng,...
        figOutTF,mlapdSoma,regionOutlineFlat);
    
    % Dot in the ctxC
    if iCenter(1) > 0
        row = cellfun(@(X) sum(X(:,1)> 0),mlapdDot);
    else
        row = cellfun(@(X) sum(X(:,1)< 0),mlapdDot);
    end
    nCluster = cellfun(@sum,clusterRng);
    % More than half of the ctxC projection is in the cluster
    row = (nCluster./row).*100 >= 75;
        
    % Get cells belongs to the cluster also most of the projection in the
    % hemisphere is within the cluster
    row = find(row);
    row = intersect(row, cellITcluster{iCluster});
    
    % Soma+ cells
    TF = mlapdSoma(row,:);
    TF = any(TF,2);
    row = row(TF);
    
    disp(['Cells meet the criteria: ', num2str(size(row,1))]);
%     continue
        
    % Plot BC model -----------------------------------------------------------
    % Shuffle order for random colorcode
    row = TBS.shuffleRows(row);
    xyzSoma2 = xyzSoma(row,:); xyzDot2 = xyzDot(row);
    [im{iCluster}, outline] = TBS.BCmodel(xyzDot2,xyzSoma2,[],2,annoMap,refSetting);
end

% cell2stack(cellArray,dimension4stackCell)
im = TBS.cell2stack(im,3);
im = uint8(im);

MIJ.createImage(cat(3,im,outline));
MIJ.run("Stack to Hyperstack...", "order=xyzct channels=4 slices=1056 frames=1 display=Composite");

%% SupFig 6, Rolony doesnt belongs to 5 brain regions
% Note, some of the rolony doesnt belongs the 5 region so not included in
% the heatmap

% Get brain region voxels, colorcoded
% Default for 65A: CtxI, CtxC, Thal, Str, Mb
regVoxel2 = regVoxel{1};
for i = 2:numel(regVoxel)
    regVoxel2 = regVoxel2 + regVoxel{i}.*i;
end

% Rolony in the 5 region
TF = vertcat(xyzDot{:});
TF = TBS.xyz2v(TF,regVoxel2>0);

disp(['Rolony within the 5 regions: ', num2str(sum(TF)),...
    ' in total: ', num2str(numel(TF))]);

% Visualization -----------------------------------------------------------
TF = any(xyzSoma,2);

% Rolony/soma location
im = [vertcat(xyzDot{:});xyzSoma(TF,:)];
im = TBS.xyzv2im(size(regVoxel2),im,[]);

% 5regions/regionOutline/BCs
test = cat(3,regVoxel2,regionOutline.*255,im.*255);
test = uint8(test);

MIJ.createImage(test);
MIJ.run("Stack to Hyperstack...",...
    "order=xyzct channels=3 slices=528 frames=1 display=Composite");

%% SupFig.6 Projection strength heatmap

% Projection strength in 5 regions: CtxI, CtxC, Thal, Str, Mb
stat = [];
for i = 1:numel(regVoxel)
    TF = cellfun(@(X) TBS.xyz2v(X,regVoxel{i}),xyzDot,...
        'UniformOutput',false);
    stat(:,i) = cellfun(@sum,TF);
end

% Sorting

% 1-4: PT, CT, ITc, ITi
idx = cellType;

% 1. Sorting using groups
[idx,I] = sort(idx,'ascend');
stat = stat(I,:);

% 2. Sorting using Thal counts
for i = 1:max(idx)
    row = idx == i;
    iStat = stat(row,:);
    
    iStat = sortrows(iStat,regCol.thal,'descend');
    stat(row,:) = iStat;
end

% Fig.2 Heatmap -----------------------------------------------------------
figure; 
h = heatmap(stat,'GridVisible','off','ColorScaling','log',...
    'FontSize',12,'FontName','Myriad Pro'); 

% xy labels
h.XData = {'CtxI','CtxC','Thal','Str','Mb'};
h.YDisplayLabels = nan(size(stat,1),1);

% Group lines
n = accumarray(idx,1);
n = cumsum(n);
S = struct(h); yline(S.Axes,n,'k');

colormap(flipud(gray)); caxis([0 4.7]);
set(gcf,'Position',[100 100 180 500]);

%% SupFig 6. Soma layer distribution
% Distribution of all barcoded somas (BC with soma)

TF = any(mlapdSoma,2);
stat = mlapdSoma(TF,3);
disp(['Number of barcode soma: ',num2str(size(stat,1))]);

% Disctribution along cortical depth (%)
N = histcounts(stat,somaDepthEdges,'Normalization','probability');
N = N.*100;

y = TBS.getMedEdges(somaDepthEdges);
stat = [N',y'];

% Plot in Prism

%% Fig 3. Soma layers

% Cell type 1-4: PT, CT, ITc, ITi 
mlapdSoma2 = {};
for i = 1:max(cellType)
    row = cellType == i & any(mlapdSoma,2);    
    mlapdSoma2{i,1} = mlapdSoma(row,:);
end

% Distribution of each group (%)
N = cellfun(@(X) histcounts(X,somaDepthEdges),mlapdSoma2,'UniformOutput',false);
N = vertcat(N{:});
N = N./sum(N,1).*100;

y = TBS.getMedEdges(somaDepthEdges);
stat = [N',y'];

% Visulization ------------------------------------------------------------
% Plot a portion of the soma

% Number of random cells for each group
randN = 200;

% Shuffle within group
mlapdSoma2 = cellfun(@(X) TBS.shuffleRows(X),mlapdSoma2,'UniformOutput',false);
mlapdSoma2 = cellfun(@(X) X(1:randN,:),mlapdSoma2,'UniformOutput',false);
c = [1:numel(mlapdSoma2)]';
c = TBS.repmat2cell(c,mlapdSoma2);

mlapdSoma2 = vertcat(mlapdSoma2{:});
c = vertcat(c{:});

TBS.plotSoma(mlapdSoma2,c); 
yline(0,'k:'); % Pia
colormap(cmap4); pbaspect([1.4 1 1]);
set(gcf,'Position',[100 100 400 250]);

% Layer location
% mlapd2flatmapV(mlapd,refSetting,flatmap,edgeCorrect);
c = TBS.mlapd2flatmapV(mlapdSoma2,refSetting,layerFlat,true);
% Need to floor due the the layerFlat contain precentage
c = floor(c);
TBS.plotSoma(mlapdSoma2,c); 
yline(0,'k:'); % Pia
pbaspect([1.4 1 1]); set(gcf,'Position',[100 100 400 250]);

%% Fig. 3 Plot PT & CT axon tract
% Plot axon tract of 1-PT & 2-CT neuron, with soma

% Mannually pick the cells
nCell = 60;

xyzDot2 = {}; xyzSoma2 = []; c = [];
for i = 1:2 
    row = cellType == i;
    row = find(row);
    
    row = TBS.shuffleRows(row);
    row = row(1:nCell);
    
    xyzDot2 = [xyzDot2; xyzDot(row)];
    xyzSoma2 = [xyzSoma2; xyzSoma(row,:)];
    c = [c; repmat(i,size(row))];
end

scaleFactor2 = 2;
% BCmodel(xyzDot,xyzSoma,c,scaleFactor,annoMap,refSetting)
[im, outline] = TBS.BCmodel(xyzDot2,xyzSoma2,c,scaleFactor2,annoMap,refSetting);

% %% Visualization -----------------------------------------------------------
% Use FIJI: 3D viewer; LUT: PTCT5.lut
% Need to cut part of the outline for better visualization

% Coronal view
zLim = vertcat(xyzDot{:});
zLim = zLim(:,3);
zLim = zLim.*scaleFactor2;
zLim = floor(min(zLim)):ceil(max(zLim));
imOut = cat(3,im(:,:,zLim),outline(:,:,zLim));

MIJ.createImage(imOut);
MIJ.run("Stack to Hyperstack...", "order=xyzct channels=2 slices=191 frames=1 display=Composite");
% LUT: CTPT3.lut

return

% Horizontal view (to cover the whole brain)
yLim = vertcat(xyzDot{:});
yLim = yLim(:,2);
yLim = yLim.*scaleFactor2;
yLim = floor(min(yLim)):ceil(max(yLim));
imOut = cat(3,im,outline);
imOut = imOut(yLim,:,:);

MIJ.createImage(imOut);
MIJ.run("Stack to Hyperstack...", "order=xyzct channels=2 slices=1056 frames=1 display=Composite");

%% SupFig 7. CT/PT rolony with region boundaries
% 04062022, resize xy to 2 for visulization

% Resize factor for xy, not z
resizeFactor = 2;

% Rolony location and cell type: 1-PT & 2-CT
xyzDot2 = {}; c = [];
for i = 1:2 
    row = cellType == i;
    row = find(row);
    
    xyzDot2 = [xyzDot2; xyzDot(row)];
    c = [c; repmat(i,size(row))];
end

c = TBS.repmat2cell(c,xyzDot2);

xyzDot2 = vertcat(xyzDot2{:});
c = vertcat(c{:});

% Annotated region in 3D
outline = TBS.getAnnoRegion(refSetting);
% Resize
outline = imresize(outline,2,'nearest');
% Get outline of annotated regions
outline = TBS.annoRegionOutline(outline);

% resize coordinates
xyzDot2(:,1:2) = xyzDot2(:,1:2).*resizeFactor;

% Plot rolony location
im = TBS.xyzv2im(size(outline),xyzDot2,c);

im = uint8(im);
outline = uint8(outline);
im = cat(3,im,outline);

MIJ.createImage(im);
MIJ.run("Stack to Hyperstack...", "order=xyzct channels=2 slices=528 frames=1 display=Composite");
% Plot in ImageJ

%% SupFig 7. CT/PT with vGlut2
% Using raw vGlut2 data to preserve the resolution

% Mannually select image
imName = 'EF65ASlide6L_Ab01_Whole.tif';     % 11
imName = 'EF65ASlide19L_Ab01_Whole.tif';    % 39
imName = 'EF65ASlide27R_Ab01_Whole.tif';    % 56
imName = 'EF65ASlide38R_Ab01_Whole.tif';    % 78
imName = 'EF65ASlide50L_Ab01_Whole.tif';    % 101

% Rolony and whole brain image
% getRolonyOnWholeIm(imName,varIn,cellType,ch,sysSetting,directory)
im = TBS.getRolonyOnWholeIm(imName,axonBC,cellType,4,sysSetting,directory);

% Make rolony bigger
SE = strel('disk',8);
im(:,:,1) = imdilate(im(:,:,1),SE);

% Convert to 1 micro/pixel
im = imresize(im,imageSetting.resolution,'nearest');
MIJ.createImage(im);

%% Fig 4-5. Prepare data for cortical analysis =======================================================================

% Min count per region
minCount = 5;% Fig. 4I Ctrl: 15

% Whether is a LatC control
isCtrl = false;

% Exclude rolony within 95% radius from the injection center
% (The current setting will leave some layer 1 rolony around injection
% site)
% (nearSoma function includs delete all 0 rows)
TF = TBS.nearSoma(mlapdDot,mlapdSoma,95,isCtrl);
mlapdDot2 = cellfun(@(X,Y) X(Y,:),mlapdDot,TF,'UniformOutput',false);

% Only include IT cells, with cortical projections
TF = cellfun(@(X) size(X,1),mlapdDot2);
TF = TF >= minCount & ismember(cellType,3:4);

% ML/AP/Depth coordinates
mlapdDot2 = mlapdDot2(TF);
mlapdSoma2 = mlapdSoma(TF,:);

% mlapd2flatmapV(mlapd,refSetting,flatmap,edgeCorrect);
layerDot2 = TBS.mlapd2flatmapV(mlapdDot2,refSetting,layerFlat,true);
layerSoma2 = TBS.mlapd2flatmapV(mlapdSoma2,refSetting,layerFlat,true);

% Reference xyz coordinates
xyzDot2 = xyzDot(TF);
xyzSoma2 = xyzSoma(TF,:);
codeBook2 = codeBook(TF,:);

disp(['Cells used for analysis: ',num2str(size(mlapdDot2,1))]);

% Cortical ROI ------------------------------------------------------------
% LatI, MedI, MedC, LatC
ctxROI = cellfun(@(X) X(:,1) <= -mlBoundary, mlapdDot2,'UniformOutput',false);
ctxROI(:,2) = cellfun(@(X) X(:,1) > -mlBoundary & X(:,1) < 0, mlapdDot2,'UniformOutput',false);
ctxROI(:,3) = cellfun(@(X) X(:,1) < mlBoundary & X(:,1) > 0, mlapdDot2,'UniformOutput',false);
ctxROI(:,4) = cellfun(@(X) X(:,1) >= mlBoundary, mlapdDot2,'UniformOutput',false);

% (Report) 
n = cellfun(@sum,ctxROI);
n = n >= minCount;
disp('Cells with LatI/MedI/MedC/LatC projection');
disp(sum(n));

TF = any(n(:,[1 4]),2) & ~any(n(:,2:3),2);
disp(['Cells with Lat+Med- projection: ',num2str(sum(TF)),';%: ',num2str(sum(TF)/sum(any(n,2)).*100)]);
TF = ~any(n(:,[1 4]),2) & any(n(:,2:3),2);
disp(['Cells with Lat-Med+ projection: ',num2str(sum(TF)),';%: ',num2str(sum(TF)/sum(any(n,2)).*100)]);
TF = any(n(:,[1 4]),2) & any(n(:,2:3),2);
disp(['Cells with Lat+Med+ projection: ',num2str(sum(TF)),';%: ',num2str(sum(TF)/sum(any(n,2)).*100)]);

return

%% Fig 4. Graph of corticl targets ========================================
% Seperation of Med & Lat region, outlie of region boundaries

figure; TBS.plotRegionRef(mlapdSoma,regionOutlineFlat); 
% figure; TBS.plotRegionRef(mlapdSoma,regionOutlineFlat,true);  % LatC-Ctrl
TBS.flatmapSetting(flatmapYLim);
xline([-mlBoundary 0 mlBoundary],'k:','LineWidth',2);
set(gca,'Visible','off');

%% (08082023) Fig 4. Batch projection depth/layer
% Discription: quantification rolony group distribution along the depth
% 5% depther per bin

zD = -5:5:rolonyDepthEdges(end);
zL = 1:0.5:6;
layerVoxelScale = getLayerVoxelScale(zL,mlapdDot2,ctxROI,layerFlat,refSetting);

d = cellfun(@(X,Y) X(Y,3),repmat(mlapdDot2,1,size(ctxROI,2)),...
    ctxROI,'UniformOutput',false);
l = cellfun(@(X,Y) X(Y),repmat(layerDot2,1,size(ctxROI,2)),...
    ctxROI,'UniformOutput',false);

n = cellfun(@numel,d);
TF = n >= minCount;
d(~TF) = {[]};
l(~TF) = {[]};

% Depth -------------------------------------------------------------------
stat = arrayfun(@(i) vertcat(d{:,i}),1:size(d,2),'UniformOutput',false);
stat = cellfun(@(X) histcounts(X,zD,'Normalization','probability'),...
    stat,'UniformOutput',false);
stat = vertcat(stat{:});
stat = stat.*100;

zD = TBS.getMedEdges(zD);
stat = [stat;zD]';
openvar stat

% Layer -------------------------------------------------------------------
stat = arrayfun(@(i) vertcat(l{:,i}),1:size(l,2),'UniformOutput',false);
stat = cellfun(@(X) histcounts(X,zL),stat,'UniformOutput',false);
stat = vertcat(stat{:});
% Normalized projection density per layer
stat = stat./layerVoxelScale;
stat = stat./sum(stat,2).*100;
stat = stat';
stat = imresize(stat,size(stat).*[1 2.5].*20,'nearest');
MIJ.createImage(stat); MIJ.run('Invert LUT'); ij.IJ.setMinAndMax(0,30);

%% Bulk projection coorelation with Allen AUD WT mice
% (12112023)

zD = -5:5:rolonyDepthEdges(end);

% Allen AUD WT data -------------------------------------------------------
load('MedLatAUD_WT_Allen.mat');
% Col1-4: LatI/MedI/MedC/LatC

% Projection along the depth ----------------------------------------------
d = tableOut{:,1};
d = cellfun(@(X) mat2cell(X,size(X,1),ones(size(X,2),1)),d,'UniformOutput',false);
d = vertcat(d{:});
% Change from 2.5% to 5% per bin
d = cellfun(@(X) X./sum(X).*100,d,'UniformOutput',false);
d = cellfun(@(X) reshape(X,2,[]),d,'UniformOutput',false);
d = cellfun(@(X) sum(X,1),d,'UniformOutput',false);

stat = {};
for i = 1:size(d,2)
    stat{i} = vertcat(d{:,i})';
end

% Mannually excluded MedI/C for low intensity
% WT
stat{2}(:,[5 10]) = [];
stat{3}(:,[5 8 10 11]) = [];

% all-to-all ksdistance in the area
D = {};
for i = 1:numel(stat)
    iStat = stat{i};
    iStat = cumsum(iStat,1);
    
    iD = iStat - permute(iStat,[1 3 2]);
    iD = abs(iD);
    iD = max(iD,[],1);
    iD = squeeze(iD);

    % Exclude self-self distance
    TF = eye(size(iD)) > 0;
    iD(TF) = nan;
    % mean for ploting, 95% CI, one tail
    D{i} = [mean(iD,1,'omitnan'); prctile(iD,95,1)];
end
D = horzcat(D{:})';

% Bulk projection of rolony -----------------------------------------------
d = cellfun(@(X,Y) X(Y,3),repmat(mlapdDot2,1,size(ctxROI,2)),...
    ctxROI,'UniformOutput',false);

n = cellfun(@numel,d);
TF = n >= minCount;
d(~TF) = {[]};

stat2 = arrayfun(@(i) vertcat(d{:,i}),1:size(d,2),'UniformOutput',false);
stat2 = cellfun(@(X) histcounts(X,zD,'Normalization','probability'),...
    stat2,'UniformOutput',false);
stat2 = vertcat(stat2{:});
stat2 = stat2.*100;
stat2 = stat2';
stat2 = stat2(2:end,:);

% ksdistance
D2 = {};
for i = 1:numel(stat)
    iStat = stat{i};
    iStat = cumsum(iStat,1);

    iStat2 = cumsum(stat2(:,i));

    iD = iStat2 - iStat;
    iD = abs(iD);
    iD = max(iD,[],1);
    D2{i} = iD;
end
D2 = horzcat(D2{:})';

% Plot in prism

%% Fig 4. Single cell projection depth ==================================== 
% Single cell laminar pattern in each area

zL = 1:0.5:6;
layerVoxelScale = getLayerVoxelScale(zL,mlapdDot2,ctxROI,layerFlat,refSetting);

k = 4;
% LatI, MedI, MedC, LatC
d = {}; l = {}; idxK = {}; somaL = {};
for i = 1:size(ctxROI,2)
    % Rolony depth & layer
    iD = cellfun(@(X,Y) X(Y,3),mlapdDot2,ctxROI(:,i),'UniformOutput',false);
    iL = cellfun(@(X,Y) X(Y),layerDot2,ctxROI(:,i),'UniformOutput',false);
    
    % Exclude cells without enough count
    n = cellfun(@numel,iD);
    TF = n >= minCount;
    
    % Frequency distribution (%)
    iD = cellfun(@(X) histcounts(X,rolonyDepthEdges,...
        'Normalization','probability'),iD,'UniformOutput',false);
    iD = vertcat(iD{:});
    iD = iD.*100;
    iD = iD(TF,:);
    
    iL = cellfun(@(X) histcounts(X,zL),iL,'UniformOutput',false);
    iL = vertcat(iL{:});
    % Normalize to layer voxel fraction
    iL = iL./layerVoxelScale(i,:);
    % Precentage
    iL = iL./sum(iL,2).*100;
    iL = iL(TF,:);  

    iSomaL = layerSoma2(TF);
    
    % Sup% (09112023)
    TF = sum(iL(:,1:3),2) >= 50;

    idx = zeros(size(TF));   
    % Clustering using kmeans, sort index using depth
    % kmeansDepthHist(X,k,C)
    idx(TF) = TBS.kmeansDepthHist(iD(TF,:),k);
    idx(~TF) = TBS.kmeansDepthHist(iD(~TF,:),k)+k;
    
    [~,I] = sort(idx);
    
    d{i} = iD(I,:);
    l{i} = iL(I,:);
    idxK{i} = idx(I);
    somaL{i} = iSomaL(I);
end

% For fig5 quantification
l2 = l([1 4]);
idxK = idxK([1 4]); somaL = somaL([1 4]);

% Fig 4. Plot -------------------------------------------------------------
% Mark the number for each group to draw a line
n = cellfun(@(X) size(X,1),d);
disp(['(Report) Single cell projection plotted: ', num2str(n)]); 
n = cumsum(n);

% Row, cortical depth; one column per barcode per region
d = vertcat(d{:})';
l = vertcat(l{:})';

% Exclude the 1st bin
d = d(2:end,:);

% Plot depth
figure; h = heatmap(d);
TBS.depthHeatmapSetting(h,d,rolonyDepthEdges(2:end));
S = struct(h); xline(S.Axes,n,'k-');
set(gcf,'Position',[100 100 1800 200]);

% Plot layer density
figure; h = heatmap(l);
TBS.depthHeatmapSetting(h,l,zL(1:end));
S = struct(h); 
yline(S.Axes,2.5:2:8.5,'w-','LineWidth',2);
xline(S.Axes,n,'k-'); caxis([0 50]);
ylabel('Layer'); S.Colorbar.Label.String = 'Normalized density (%)';   
set(gcf,'Position',[100 100 1800 200]);

%% Fig 5. Lat plot --------------------------------------------------------
% 12182023
% (internal reference) Mean projection intensity for each group in LatI & C
meanL = {};
for i = 1:numel(l2)
    iL = l2{i};
    idx = idxK{i};
    for j = 1:max(idx)
        TF = idx == j;
        meanL{i}(j,:) = mean(iL(TF,:),1);
    end
end

meanL = vertcat(meanL{:})';

figure; h = heatmap(meanL);
TBS.depthHeatmapSetting(h,meanL,zL(1:end));
S = struct(h); 
yline(S.Axes,2.5:2:8.5,'w-','LineWidth',2);
xline(S.Axes,8.5,'k-'); caxis([0 50]);
ylabel('Layer'); S.Colorbar.Label.String = 'Normalized density (%)';   
set(gcf,'Position',[100 100 800 200]);

%% SupFig. Report cluster-group heatmat ------------------------------------

r = 2; % 1-LatI; 2-LatC

% Cell group using soma layer (fig5threshold, defined in Fig. 5 below)
group2 = isnan(somaL{r});
for i = 1:numel(fig5threshold)-1
    TF = somaL{r} >= fig5threshold(i) & somaL{r} < fig5threshold(i+1);
    group2(:,end+1) = TF & ~group2(:,1);
end

% projection cluster precentage in each group of cells
stat = [];
for i = 1:max(idxK{r})
    TF = idxK{r} == i & group2(:,2:end);
    stat(:,i) = sum(TF,1);
end
n = sum(stat,1);
stat = stat./n.*100;

% Image output, adjust the cluster length to cell counts per group
% Output width-800, height-200
n = round(n./sum(n).*800);
for i = 1:numel(n)
    iStat = imresize(stat(:,i),[200, n(i)],'nearest');
    if i == 1
        stat2 = iStat;
    else
        stat2 = [stat2, iStat];
    end
end
MIJ.createImage(stat2);

%% Fig4. Single cell model ================================================
% Discription: coronal view of two example neurons of IT cells

% Scaling for the image output
scaleFactor = 2;

% (for manual check)
% With biLat & Med & soma 
I = group(:,2) & any(group(:,3:4),2) & any(xyzSoma2,2);
% With biMed & soma
% I = group(:,1) & any(xyzSoma2,2);

I = find(I);

% Sort using rolony number
n = cellfun(@(X) size(X,1),xyzDot2);
n = n(I);
[~,I2] = sort(n,'descend');
I = I(I2);

% Manually checked: 
% 370: 1   1   2   3   2   1   1   3   3   2   1   1   1   3   3   3   3
I = 370;

% Plot
xyzSoma3 = xyzSoma2(I,:); xyzDot3 = xyzDot2(I);
% Reconstruction
[im, outline] = TBS.BCmodel(xyzDot3,xyzSoma3,[],scaleFactor,annoMap,refSetting);
% Dot plot (updated: 12192022)
im2 = TBS.xyzv2im(size(im),xyzDot3{:}.*scaleFactor,[]);
SE = strel('sphere',2);
im2 = imdilate(im2,SE);

% Get the limits of z-axis
zLim = vertcat(xyzDot{:});
zLim = zLim(:,3);
zLim = zLim.*scaleFactor;
zLim = floor(min(zLim)): ceil(max(zLim));

% Output in ImageJ
MIJ.createImage(cat(3,im(:,:,zLim),im2(:,:,zLim),outline(:,:,zLim)));
MIJ.run("Stack to Hyperstack...", "order=xyzct channels=3 slices=191 frames=1 display=Composite");

%% Fig4. Sup projeciton comparison ======================================== 
% (09052023)
% Discription: compare superifical projection within and across areas

zL = 1:0.5:6;
layerVoxelScale = getLayerVoxelScale(zL,mlapdDot2,ctxROI,layerFlat,refSetting);

l = cellfun(@(X,Y) X(Y),repmat(layerDot2,1,size(ctxROI,2)),...
    ctxROI,'UniformOutput',false);

n = cellfun(@numel,l);

l = cellfun(@(X) histcounts(X,zL),l,'UniformOutput',false);
% Normalized density
for i = 1:size(l,2)
    iL = cellfun(@(X) X./layerVoxelScale(i,:),l(:,i),'UniformOutput',false);        
end
l = cellfun(@(X) X./sum(X).*100,l,'UniformOutput',false);

% SUP%
l = cellfun(@(X) sum(X(1:3)),l);

% Exclude ROI below min BC counts
TF = n >= minCount;
l(~TF) = nan;

% 1-4: Med, Sup and Inf; Mean of Lat, Sup and Inf
l = [mean(l(:,2:3),2,'omitnan'),mean(l(:,[1 4]),2,'omitnan')];

l = l(:,[1 1 2 2]);
l(:,[2 4]) = 100-l(:,[2 4]);

l(:,[1 3]) = l(:,[1 3]) > 50;
l(:,[2 4]) = l(:,[2 4]) >= 50;

disp('Counts: Med, Sup/Inf; Lat, Sup/Inf');
disp(sum(l));

disp('Counts Lat+Med-: Med, Sup/Inf; Lat, Sup/Inf');
disp(sum(l & ~any(l(:,1:2),2)));

disp('Counts Lat+Med+: Med, Sup/Inf; Lat, Sup/Inf');
disp(sum(l & any(l(:,1:2),2) & any(l(:,3:4),2)));

% Get combinations --------------------------------------------------------
% 12022023, delete co-projection
C = num2cell(1:4);
C = [C, {[1 3]}, {[1 4]}, {[2 3]}, {[2 4]}];

% Output to show the pattern ----------------------------------------------
C2 = false(numel(C),4);
for i = 1:numel(C)
   C2(i,C{i}) = true;
end
MIJ.createImage(imresize(C2.*255,30,'nearest')'); MIJ.run('Invert LUT');

stat = [];
for i = 1:size(C2,1)
    iL = l == C2(i,:);
    iL = all(iL,2);
    stat(i,1) = sum(iL);
end
stat = stat./sum(stat);

%% Fig5. Lat 3 groups, thresholding soma ==================================
% (09052023) Changed into grouping using layer
% (02212024) Change the name from group1/2/3 to upper/middle/deep

% soma threshold for diving into three groups
% Group 1: L1,2/3; Group 2: L4 & upper L5. Group3, lower L5 & 6
fig5threshold = [1, 3, 4.5, 6];

% Lat+ cells: either Lat pass threshold
TF = cellfun(@sum,ctxROI(:,[1 4]));
TF = TF >= minCount;
TF = any(TF,2);

% 1: No Lat-projection or no soma; 2: group1; 3: group2; 4: group3
group = ~TF | isnan(layerSoma2);
for i = 1:numel(fig5threshold)-1
    TF = layerSoma2 >= fig5threshold(i) & layerSoma2 < fig5threshold(i+1);
    group(:,i+1) = TF & ~group(:,1);
end

fig5SomaGroup = group;
% Report cell number
disp('Grouping using soma depth. cell count: excluded, group1,2,3: ')
sum(group)

% Fig. 5 Soma location laminar distribution -------------------------------

% Scatter plot
[~,c] = max(group,[],2);
TF = c > 1;

TBS.plotSoma(mlapdSoma2(TF,:),c(TF,:)); alpha 0.5;
colormap(cmapROB); pbaspect([1.4 1 1]);
set(gcf,'Position',[100 100 400 250]);

%% Fig5. single cell projection model and flatmap 
% Discription: example cortical projection on flatmap
% (09062023 Checked)

% 505:   1   1   3   2   2   2   3   1   3   2   1   3   3   1   4   2   1
% 1:  1   1   1   1   1   1   2   1   3   2   1   4   1   4   1   2   3
% 431:   1   1   2   4   4   1   1   1   3   2   2   1   3   1   1   3   1

% Manually pick BC ID
I = [505 1 431];

% Flat map
for i = 1:numel(I)    
    
    xy = mlapdDot2{I(i)};
    
    figure; scatter(xy(:,1),xy(:,2),10,cmapROB(i,:),'filled'); alpha 0.8
    hold on; TBS.plotRegionRef(mlapdSoma,regionOutlineFlat);
    TBS.flatmapSetting(flatmapYLim);
    xline([-mlBoundary 0 mlBoundary],'k:','LineWidth',2);   
end

% Single cell model -------------------------------------------------------

scaleFactor = 2;

xyzSoma3 = xyzSoma2(I,:); xyzDot3 = xyzDot2(I);
[im, outline] = TBS.BCmodel(xyzDot3,xyzSoma3,[],scaleFactor,annoMap,refSetting);

% Get the limits of z-axis
zLim = vertcat(xyzDot{:});
zLim = zLim(:,3);
zLim = zLim.*scaleFactor;
zLim = floor(min(zLim)): ceil(max(zLim));

MIJ.createImage(cat(3,im(:,:,zLim),outline(:,:,zLim)));
MIJ.run("Stack to Hyperstack...", "order=xyzct channels=2 slices=191 frames=1 display=Composite");
% Output in ImageJ

%% Fig 5. Rolony location visualization -----------------------------------
% Discription: plot rolony xyz in reference brain outline
% (updated: 12192022) changed to heatmap

% Group to visualizing; group1:3
i = 3;

% 04042022,Define using the minimum group count
n = min(sum(group(:,2:4)));
disp(['Plot number of cells: ',num2str(n)]);

% Color code
[~,c] = max(group,[],2);
c = c -1;

I = find(c == i);
I = TBS.shuffleRows(I);
I = I(1:n);

xyz = vertcat(xyzDot2{I});

% Image of rolony location
im = TBS.xyzv2im(size(annoMap),xyz,[]);

% Only plot the z with data for better visualization
% Get the limits of z-axis
zLim = vertcat(xyzDot{:});
zLim = zLim(:,3);
zLim = floor(min(zLim)): ceil(max(zLim));

% Brain outline
outline = TBS.getBrainOutline(annoMap,5);
outline = outline(:,:,zLim);
outline = single(outline);

im = im(:,:,zLim);

% Sum projection
im = sum(im,3); im = fliplr(im);
outline = sum(outline,3); outline = fliplr(outline);

% Fixed min & max for lookup table for comparison
disp(['Lookup table max: ', num2str(n./40)]);

% Visualized in ImageJ
MIJ.createImage(cat(3,im,outline));
% Output in ImageJ

%% SupFig 9. Rolony depth quantification -----------------------------------
% Discription: quantification rolony group distribution along the depth
% 5% depther per bin

r = 4; % 1-LatI; 4-LatC

% Bin edges for depth and layer
zD = -5:5:rolonyDepthEdges(end);
zL = 1:0.25:6; 

% layer voxel normalization for current layer binning
layerVoxelScale = getLayerVoxelScale(zL,mlapdDot2,ctxROI,layerFlat,refSetting);

d = cellfun(@(X,Y) X(Y,3), mlapdDot2,ctxROI(:,r),'UniformOutput',false);
l = cellfun(@(X,Y) X(Y), layerDot2,ctxROI(:,r),'UniformOutput',false);

stat = []; stat2 = [];
for i = 2:size(group,2)
    TF = group(:,i);
    d2 = vertcat(d{TF});
    l2 = vertcat(l{TF});
    
    stat(end+1,:) = histcounts(d2,zD);
    stat2(end+1,:) = histcounts(l2,zL);
end

% depth (normalized to 100%)
stat = stat./sum(stat,2).*100;

z = TBS.getMedEdges(zD);
stat = [stat',z'];
% Plot in Prism

% layer (normalized to 100%)
stat2 = stat2./layerVoxelScale(r,:);
% Exclude the last bin due to not enough data
stat2 = stat2(:,1:end-1);
stat2 = stat2./sum(stat2,2).*100;
stat2 = stat2';
stat2 = imresize(stat2,size(stat2).*[1 2.5].*20,'nearest');
MIJ.createImage(stat2); MIJ.run('Invert LUT'); ij.IJ.setMinAndMax(0,30);

%% Fig 5. Relationship between soma and projection depth
% Plot soma depth along with projection pattern for soma+ cells
% Cannot directly save it in MATLAB, need to save it in ImageJ
% (Modified 09102023)

r = 3; % LatI,MedI,MedC,LatC

% Threshold for lower projections
threshold = 4.5;

xLabel = [3 4.5];

% Projection binning
% 02282022: seems binning size = 5 is better
binSize = 5;
z = -binSize:binSize:rolonyDepthEdges(end);

d = cellfun(@(X,Y) X(Y,3), mlapdDot2,ctxROI(:,r),'UniformOutput',false);
% Lower projection (%)
l = cellfun(@(X,Y) X(Y), layerDot2,ctxROI(:,r),'UniformOutput',false);
p = cellfun(@(X) sum(X >= threshold)./numel(X),l);
p = p.*100;

n = cellfun(@sum,ctxROI(:,r));
TF = n >= minCount & ~group(:,1);
d = d(TF);
p = p(TF);
somaL = layerSoma2(TF);

% Sort using soma layer
[somaL,I] = sort(somaL,'ascend');
d = d(I);
p = p(I);

% Projection percentage
stat = cellfun(@(X) histcounts(X,z,'Normalization','probability'),...
    d,'UniformOutput',false);
stat = vertcat(stat{:});
stat = stat.*100;

disp(['Cell count: ', num2str(length(stat))]);

% Plot --------------------------------------------------------------------
stat = stat';
% Delete the first bin
stat = stat(2:end,:);

% Bin soma, get median p per bin
somaL2 = round(somaL./0.2).*0.2;
[somaL2,~,ic] = unique(somaL2);
p = accumarray(ic,p,[],@median);

n = accumarray(ic,1);
% Not plotting lower layer % with < 5 per bin
TF = n >= 5;
n = cumsum(n); 
n = n(TF);
p = p(TF);

% Normalize p-value to the same scale
p = p.*0.95./binSize+0.5;

% Soma layer, for x-axis
x = [];
for i = xLabel
    x(end+1) = find(somaL >= i,1,'first');
end

% Heatmap for projection frequency
figure; h = heatmap(stat); TBS.depthHeatmapSetting(h,stat,z(2:end));

% Plot soma location
S = struct(h); hold(S.Axes,'on'); 
plot(S.Axes,n,p,'r','LineWidth',2);

if r == 1 || r == 4
    set(gcf,'Position',[100 100 1200 200]);
else
    set(gcf,'Position',[100 100 550 200]);
end
colorbar off

% Right-axis
yyaxis(S.Axes,'right');
ylabel(S.Axes,'Deep-layer rolonies (%)');
S.Axes.YTick = [0 1];
S.Axes.YTickLabel =[{'0'}; {'100'}];
S.Axes.YDir = 'reverse';
S.Axes.YColor = [1 0 0];

% Xtick for soma depth
S.Axes.XTickLabel(x) = num2cell(xLabel);
S.Axes.XTickLabelRotation = 0;

%% SupFig 9. Relationship between soma and lower layer projection (%)

r = 4; % LatI,LatC

% Lower projection (%)
l = cellfun(@(X,Y) X(Y), layerDot2,ctxROI(:,r),'UniformOutput',false);
p = cellfun(@(X) sum(X >= threshold)./numel(X),l);
p = p.*100;

n = cellfun(@sum,ctxROI(:,r));
TF = n >= 10 & ~group(:,1);

p = p(TF);
somaD = mlapdSoma2(TF,3);

disp(['Cell plotted: ', num2str(sum(TF))]);

% Scatter plot
figure; scatter(somaD,p,5,'k','filled'); alpha 0.3;
g = gca; g.YDir = 'reverse'; g.YLim = [0 100]; g.XLim = [0 100];
xlabel('Soma depth (%)'); ylabel('Deep-layer rolonies (%)');
TBS.axLabelSettings('Myriad Pro',12);
set(gcf,'Position',[100 100 500 200]);

% Correlation -------------------------------------------------------------
p = cellfun(@(X) sum(X >= threshold)./numel(X),l);
p = p.*100;

n = cellfun(@sum,ctxROI(:,r));
TF = n >= 10 & group(:,4);

p = p(TF);
somaD = mlapdSoma2(TF,3);

disp(['Cell count: ', num2str(sum(TF))]);

stat = [somaD, p];
openvar stat
% Correlation in Prism

%% SupFig 9. Test minimum rolony number for focal projection

i = 4; % Use LatC for testing

% Minmum count to be tested
x = minCount:1:80;

% Number of resampling per cell
nIteration = 100;

% Propotion of rolony to compute focal projection distance
p = 0.33;

% ML-AP of rolony
xy = cellfun(@(X,Y) X(Y,1:2),mlapdDot2,ctxROI(:,i),'UniformOutput',false);

% BC with abv min rolony count in the region
n = cellfun(@sum,ctxROI(:,i));
TF = n >= max(x);
xy = xy(TF);

disp(['(Report) Cell number used for down sampling test: ', num2str(sum(TF))]);

% Compute focal distance using different minimum rolony number
d = {};
for k = 1:nIteration
    for j = x       
        % Random sample N rolony
        jxy = cellfun(@(X) TBS.shuffleRows(X),xy,'UniformOutput',false);
        jxy = cellfun(@(X) X(1:j,:),jxy,'UniformOutput',false);
        
        % focalProjPercentage(xy,dia)
        D = cellfun(@(X) TBS.focalProjPct(X,p),jxy,'UniformOutput',false);
        D = cellfun(@mean,D);
        
        d{end+1} = D;
    end
end

% Delete cells with not enough cell to finish the test
d = horzcat(d{:});

% Use all the value as ground true ----------------------------------------
groundTruth = cellfun(@(X) TBS.focalProjPct(X,p),xy,'UniformOutput',false);
groundTruth = cellfun(@mean,groundTruth);

% Difference between random sample and groundturth
d = (d - groundTruth)./groundTruth;
d = abs(d);

% Get the median for every cell
d = reshape(d,size(d,1),numel(x),[]);
d = median(d,3);

% Plot test result --------------------------------------------------------
figure; plot(x,d','Color',[0 0 0 0.1]);
hold on; plot(x,median(d,1),'r','LineWidth',2);   % Median
TBS.axLabelSettings('Myriad Pro',12);
xlabel('Subsample rolony number');
% ylabel('Absolute difference from ground truth');
ylabel('Error (%)');
set(gcf,'Position',[100 100 250 250]);

[x; median(d,1)]

%% Fig 5. Quantificaiton: focal projection distance for three group
% Minimum 55 rolony per cell

r = 4; % 1-LatI, 4-LatC

p = 0.33;

% Group number
[~,c] = max(group,[],2);

% ML-AP of rolony
xy = cellfun(@(X,Y) X(Y,1:2),mlapdDot2,ctxROI(:,r),'UniformOutput',false);

% BC with abv min rolony count in the region
n = cellfun(@sum,ctxROI(:,r));
TF = n >= 55 & any(mlapdSoma2,2);

xy = xy(TF);
c = c(TF);

% focalProjPct(mlap,p)
d = cellfun(@(X) TBS.focalProjPct(X,p),xy,'UniformOutput',false);
d = cellfun(@mean,d);

stat = {};
for i = 2:4
    stat{end+1} = d(c == i);
end
% Plot in Prism

disp('Focal projection cell counts: ');
cellfun(@numel,stat)

%% Fig5. Projection targets of three groups of cell

nIteration = 2000;

% Med+ cells --------------------------------------------------------------
roi = cellfun(@(X,Y) X | Y, ctxROI(:,2),ctxROI(:,3),'UniformOutput',false);

% Rolony count pass threshold
n = cellfun(@sum,roi);
TF = n >= minCount;

% (Report) Med+ percentage, plot in Prism
[stat,data] = cmptCI(TF,group(:,2:end),nIteration);

% Str+ cells --------------------------------------------------------------
% Rolony within str
roi = regVoxel{regCol.str};
roi = cellfun(@(X) TBS.xyz2v(X,roi),xyzDot2,'UniformOutput',false);

n = cellfun(@sum,roi);
TF = n >= minCount;

% (Report) Str+ percentage, plot in Prism
[stat,data] = cmptCI(TF,group(:,2:end),nIteration);

% LatC+ cells -------------------------------------------------------------
roi = ctxROI(:,4);

n = cellfun(@sum,roi);
TF = n >= minCount;

% (Report) LatC+ percentage, plot in Prism
[stat,data] = cmptCI(TF,group(:,2:end),nIteration);

%% SupFig 9. Split Lat+ no-soma cells into two groups
% Upper & lower projecting

group2 = [];
for r = [1 4]
    % Median projection depth in the region
    d = cellfun(@(X,Y) X(Y), layerDot2,ctxROI(:,r),'UniformOutput',false);
    d = cellfun(@median,d);
    d = d <= threshold;
    % Upper - 2; lower-1
    d = d + 1;
    
    % Cells with projection pass the threshold
    TF = cellfun(@sum,ctxROI(:,r));
    TF = TF >= minCount;
    
    group2(TF,end+1) = d(TF);    
end

% IT cells without soma
TF = any(mlapdSoma2,2);
group2(TF,:) = 0; 

% LatI & C, no belongs to the same group
TF = all(group2,2) & group2(:,1)~= group2(:,2);
disp(['LatI & C not belongs to the same group: ',num2str(sum(TF))]);
disp(['biLat counts: ', num2str(sum(all(group2,2)))]);

% Upper & lower projecting cells
group2 = max(group2,[],2) == [2 1];

% Exclude cells with different group
group2(TF,:) = 0;

disp('cell counts (upper/lower): ');
disp(sum(group2));

%% SupFig 9, Med/Str/LatC in cells w/o soma

nIteration = 2000;

% Med+ cells --------------------------------------------------------------
roi = cellfun(@(X,Y) X | Y, ctxROI(:,2),ctxROI(:,3),'UniformOutput',false);

% Rolony count pass threshold
n = cellfun(@sum,roi);
TF = n >= minCount;

% (Report) Med+ percentage, plot in Prism
[stat,data] = cmptCI(TF,group2,nIteration);

% Str+ cells --------------------------------------------------------------
% Rolony within str
roi = regVoxel{regCol.str};
roi = cellfun(@(X) TBS.xyz2v(X,roi),xyzDot2,'UniformOutput',false);

n = cellfun(@sum,roi);
TF = n >= minCount;

% (Report) Str+ percentage, plot in Prism
[stat,data] = cmptCI(TF,group2,nIteration);

% LatC+ cells -------------------------------------------------------------
roi = ctxROI(:,4);

n = cellfun(@sum,roi);
TF = n >= minCount;

% (Report) LatC+ percentage, plot in Prism
[stat,data] = cmptCI(TF,group2,nIteration);

%% SupFig 9. Focal projection in cells w/o soma location

r = 4; % 1-LatI; 4-LatC

n = cellfun(@sum,ctxROI(:,r));
TF = n >= 55;

xyz = cellfun(@(X,Y) X(Y,:), mlapdDot2, ctxROI(:,r),'UniformOutput',false);
xyz = xyz(TF);

d = cellfun(@(X) TBS.focalProjPct(X,p),xyz,'UniformOutput',false);
d = cellfun(@mean,d);

% Upper and lower
stat = [{d(group2(TF,1))},{d(group2(TF,2))}];
% Plot in Prism

%% SupFig. DAPI and vGlut2 image ===================================================
% (09062023 updated)

% resize for 25 um to 12.5 um/voxel
resizeFactor = 2;

scaleTform = eye(3).*resizeFactor;
scaleTform(4,4) = 1;
scaleTform = affine3d(scaleTform);

% Output size
sz = size(annoMap).*resizeFactor;
R = imref3d(sz);

% Transformation matrix and displacement field ----------------------------
load('reg2AllenTform.mat');

tform = reg2AllenTform.tform;
D = reg2AllenTform.D;

% resize transformation matrix and displacement field
tform.T = tform.T*scaleTform.T;
D2 = [];
for i = 1:size(D,4)
    D2(:,:,:,i) = imresize3(D(:,:,:,i),resizeFactor);
end
D = D2.*resizeFactor;

im = TBS.getStack(msBrainName,[]);

sz = size(im);
im = reshape(im,sz(1),sz(2),4,[]);

% Get DAPI: channel 1; vglut2: channel 3
marker = {};
for i = [1 3]
     iIm = im(:,:,i,:);
     iIm = squeeze(iIm);
     
     iIm = imwarp(iIm,tform,'OutputView',R);
     iIm = imwarp(iIm,D);
         
     marker{end+1} = iIm;
end

% Resize layer 
marker{3} = imresize3(layer,resizeFactor,'Method','nearest');

%% SupFig. 6 layer & cell location
% (checked 09092023)
% Slide selected for the output
I = 606:630;

SE = strel('disk',1);

% Get median value for the reference
im = cellfun(@(X) median(X(:,:,I),3),marker,'UniformOutput',false);

% Layer outline 
im{3} = imdilate(im{3},SE)~= im{3};
im{3} = uint8(im{3});

% PT and CT soma

sz = size(marker{1});

% Size for soma dilation
SE = strel('disk',1);

im2 = zeros(sz(1:2),'uint8');
for i = 1:2
    TF = cellType == i;
    ixyz = xyzSoma(TF,:);
    ixyz = transformPointsForward(scaleTform,ixyz);
    iIm = TBS.xyzv2im(sz,ixyz,[]);    
    
    iIm = iIm(:,:,I);
    
    % Dilatetion for visilizaiton
    iIm = imdilate(iIm,SE);
    
    % Sum projection
    im2(:,:,i) = sum(iIm,3);
end

im = cat(3,im{:},im2);

MIJ.createImage(im);

%% SupFig.9 deep layer projection =========================================
% (checked 09092023)
% Group to visualizing; group1:3
i = 3;

SE = strel('disk',1);

% Slide selected for the output
I = 606:630;

TF = group(:,i+1);
xyz = vertcat(xyzDot2{TF});
xyz = transformPointsForward(scaleTform,xyz);

[xyz,~,ic] = unique(round(xyz),'rows');
v = accumarray(ic,1);

im = {};
iIm = TBS.xyzv2im(sz,xyz,v);
% 0-10
im{1} = sum(iIm(:,:,I),3);
% DAPI: 25-220
im{2} = median(marker{1}(:,:,I),3);
% Layer outline: 0-4
im{3} = median(marker{3}(:,:,I),3);
im{3} = imdilate(im{3},SE)~= im{3};

MIJ.createImage(cat(3,im{:}));

%% SupFig11. additional flatmaps for 3 groups
% (05272023 for revision, checked 09112023)

% Mannually picked BC for each group
supFig10 = {};
supFig10{1} = [1271,3844,1197,1110,468,981,5735,1575,5102,4224];
supFig10{2} = [6547,6074,3629,4140,1260,6564,6838,6956,5275,5501];
supFig10{3} = [4944,5994,5615,3649,6882,4639,3827,779,6317,501];

% group1
% 1271: 1  2  4  3  2  2  2  1  3  2  3  3  1  4  2  4  3
% 3844: 2  3  4  1  2  2  2  3  3  2  3  2  1  2  1  3  3
% 1197: 1  2  4  1  3  1  2  2  3  2  1  2  3  1  1  1  4
% 1110: 1  2  3  1  4  3  2  1  3  2  4  4  1  4  4  3  1
% 468: 1  1  3  1  2  3  1  4  3  2  1  3  4  4  1  2  1
% 981: 1  2  2  2  2  1  3  2  3  2  3  4  2  1  3  4  4
% 5735: 3  4  3  4  1  4  2  2  3  2  3  3  4  1  4  1  2
% 1575: 1  3  2  4  2  3  1  4  3  2  1  2  1  4  2  2  3
% 5102: 3  2  4  2  3  2  1  2  3  2  2  1  2  1  2  4  2
% 4224: 2  4  3  4  3  2  2  3  3  2  2  1  4  1  3  2  4
% group2
% 6547: 4  2  2  4  4  3  1  4  3  2  3  4  1  3  3  1  2
% 6074: 4  1  2  3  3  3  1  2  3  2  2  4  4  4  4  2  1
% 3629: 2  3  2  1  3  3  1  2  3  2  1  1  3  1  2  1  1
% 4140: 2  4  2  4  3  1  2  2  3  2  3  3  4  4  3  4  1
% 1260: 1  2  4  3  1  2  3  2  3  2  1  3  4  3  3  3  4
% 6564: 4  2  3  1  2  3  1  2  3  2  4  3  2  1  2  4  4
% 6838: 4  3  2  1  4  3  1  4  3  2  3  3  2  4  2  3  4
% 6956: 4  3  3  2  1  2  2  4  3  2  1  2  2  1  1  4  3
% 5275: 3  3  2  1  2  3  1  2  3  2  2  1  2  2  3  1  1
% 5501: 3  3  4  4  1  2  1  4  3  2  4  2  1  4  3  3  3
% group3
% 4944: 3  2  2  3  2  4  1  1  3  2  4  2  4  2  3  2  1
% 5994: 4  1  2  1  1  1  2  4  3  2  2  4  1  1  4  4  1
% 5615: 3  4  2  1  3  1  1  3  3  2  3  4  1  4  1  1  1
% 3649: 2  3  2  2  1  2  4  2  3  2  4  1  4  3  4  3  4
% 6882: 4  3  2  3  2  1  4  2  3  2  3  3  2  1  4  4  2
% 4639: 3  1  3  4  1  4  2  1  3  2  3  1  1  3  3  3  2
% 3827: 2  3  3  4  2  3  3  1  3  2  4  3  4  1  1  3  1
% 779: 1  2  1  1  4  2  2  2  3  2  4  1  1  4  4  3  2
% 6317: 4  1  4  4  4  3  3  2  3  2  2  2  3  2  4  1  2
% 501: 1  1  3  2  1  4  2  4  3  2  3  2  2  3  1  4  4

g = 1;

figure;
for i = supFig10{g}
    % double check the soma depth
    d = mlapdSoma2(i,3);
    disp(['Soma depth: ', num2str(mlapdSoma2(i,3))])

    disp([num2str(i),': ',num2str(codeBook2(i,:))]);    
        
    % double check the group number
    ig = group(i,:);
    ig = find(ig)-1;
    
    xy = mlapdDot2{i};
    
    hold off;
    scatter(xy(:,1),xy(:,2),10,cmapROB(ig,:),'filled'); alpha 0.8
    hold on; TBS.plotRegionRef(mlapdSoma,regionOutlineFlat);
    TBS.flatmapSetting(flatmapYLim);
    xline([-mlBoundary 0 mlBoundary],'k:','LineWidth',2);  
end

%% Evaluate the misassignment effect of striatal rolony
% (for revision) 06122023
% There are errors in CCFv3 registration, particularly in Str
% This section is to evaluate the effact of misassignment of rolony on LatI
% analysis

% strROI were mannually drawn using registered vGlut2 images (25 um/voxel)
load('strROI.mat');

% Rows for LatI
% --> mlapdDot
test = TBS.getMLAPD(xyzDot2,ctxML,ctxAP,ctxDepthPrctile);
fh = @(X) X(:,3)>= min(rolonyDepthEdges) & X(:,3)<= max(rolonyDepthEdges);
TF = cellfun(@(X) fh(X),test,'UniformOutput',false);
test = cellfun(@(X,Y) X(Y,:),test,TF,'UniformOutput',false);
% --> mlapdDot2
TF2 = TBS.nearSoma(test,mlapdSoma,95,isCtrl);
TF = cellfun(@(X,Y) TBS.nestTF(X,Y),TF,TF2,'UniformOutput',false);

% (Check point) check size
all(cellfun(@sum,TF) == cellfun(@(X) size(X,1),mlapdDot2))

% Rolonies belong to LatI
latITF = cellfun(@(X,Y) TBS.nestTF(X,Y),TF,ctxROI(:,1),'UniformOutput',false);

% LatI-rolony located in strROI
strTF = cellfun(@(X) TBS.xyz2v(X,strROI),xyzDot2,'UniformOutput',false);
strTF = cellfun(@(X,Y) X & Y, latITF,strTF,'UniformOutput',false);

% Calculate percentage per cell
n = cellfun(@sum,latITF);
test = cellfun(@sum, strTF);
test = test./n.*100;
% Only includ LatI-projecting cells
test = test(n >= minCount);
disp(['In total LatI-projecting cells: ', num2str(size(test,1))]);
test = test(test > 0);
disp(['Cells with rolonies in the region: ', num2str(size(test,1))]);
disp('Mean & SD: ');
mean(test)
std(test)

% % (Check point) check dots & strROI location
% test = cellfun(@(X,Y) X(Y,:),xyzDot2,strTF,'UniformOutput',false);
% test = vertcat(test{:});
% test = TBS.xyzv2im(size(strROI),test,[]);
% test = cat(3,strROI,test);
% MIJ.createImage(test.*255);
% MIJ.run("Stack to Hyperstack...", "order=xyzct channels=2 slices=528 frames=1 display=Composite");

%% SupFig 8. Allen AUD WT bulk tracing

% Divide layers into 4 parts
zL = 1:0.25:6;
layerVoxelScale = getLayerVoxelScale(zL,mlapdDot2,ctxROI,layerFlat,refSetting);
% Exclude the last bin due to lack of data
layerVoxelScale = layerVoxelScale(:,1:end-1);

load('MedLatAUD_WT_Allen.mat');
% Col1-4: LatI/MedI/MedC/LatC

% Projection along the depth ----------------------------------------------
d = tableOut{:,1};
d = cellfun(@(X) mat2cell(X,size(X,1),ones(size(X,2),1)),d,'UniformOutput',false);
d = vertcat(d{:});
% Change from 2.5% to 5% per bin
d = cellfun(@(X) X./sum(X).*100,d,'UniformOutput',false);
d = cellfun(@(X) reshape(X,2,[]),d,'UniformOutput',false);
d = cellfun(@(X) sum(X,1),d,'UniformOutput',false);

stat = {};
for i = 1:size(d,2)
    stat{i} = vertcat(d{:,i})';
end

% Mannually excluded MedI/C for low intensity
% WT
stat{2}(:,[5 10]) = [];
stat{3}(:,[5 8 10 11]) = [];
% Plot in prism

% (SupFig) Variability within allen & current dataset  -----------------
% Allen vs Allen
D = {};
for i = 1:numel(stat)
    iStat = stat{i};
    iStat = cumsum(iStat,1);

    % all-to-all ksdistance
    iStat = iStat - reshape(iStat,size(iStat,1),1,[]);
    iStat = abs(iStat);
    iStat = max(iStat,[],1);
    iStat = squeeze(iStat);

    % Delete itself
    TF = eye(size(iStat)) > 0;
    iStat(TF) = nan;

    D{i} = mean(iStat,1,'omitnan');
end
D = horzcat(D{:})';

% Current dataset (stat2) to Allen
D2 = {};
for i = 1:numel(stat)
    iStat = stat{i};
    iStat = cumsum(iStat,1);

    iStat2 = stat2(:,i);
    iStat2 = cumsum(iStat2,1);

    iStat = iStat - iStat2;
    iStat = abs(iStat);
    iStat = max(iStat,[],1);

    D2{i} = iStat;
end
D2 = horzcat(D2{:})';
return
% Plot in prism

% Layer, MedI & C ------------------------------------------------------

l = tableOut{:,2};
l = cellfun(@(X) X(:,2:3),l,'UniformOutput',false);
l = cellfun(@(X) X(1:19,:),l,'UniformOutput',false);

% Scaling to density
l = cellfun(@(X) X./layerVoxelScale(2:3,:)',l,'UniformOutput',false);
l = cellfun(@(X) [{X(:,1)},{X(:,2)}],l,'UniformOutput',false);
l = vertcat(l{:});

% Mannually excluded MedI/C for low intensity
TF = true(size(l));
% WT
TF([5 10],1) = false; TF([5 8 10 11],2) = false;

l = l(TF); 

l = horzcat(l{:});
% For 0.25% per bin, exclude the last bin, due to not enough voxel
l = l(1:19,:);
% Normalization
l = l./sum(l,1).*100;
MIJ.createImage(imresize(l,size(l).*[1 2].*5,'nearest')); 
MIJ.run('Invert LUT'); ij.IJ.setMinAndMax(0,25);

% Check where to draw the line for Sup*Inf, to maximize the difference
% Split into two group, and calculate the diff of mean of both groups
% The assumption is the upper has high density, so didnt do abs()
cost = [];
for i = 1:size(l,1)-1
    iCost = mean(l(1:i,:),1)-mean(l(i+1:end,:),1);
    cost(i,:) = [mean(iCost),std(iCost)];
end

%% SupFig. 9 Allen AUD Cre bulk tracing

load('MedLatAUD_Cre_Allen.mat');
% Col1-4: LatI/MedI/MedC/LatC

% Projection along the depth ----------------------------------------------
d = tableOut{:,1};
d = cellfun(@(X) mat2cell(X,size(X,1),ones(size(X,2),1)),d,'UniformOutput',false);
d = vertcat(d{:});
% Change from 2.5% to 5% per bin
d = cellfun(@(X) X./sum(X).*100,d,'UniformOutput',false);
d = cellfun(@(X) reshape(X,2,[]),d,'UniformOutput',false);
d = cellfun(@(X) sum(X,1),d,'UniformOutput',false);

stat = {};
for i = 1:size(d,2)
    stat{i} = vertcat(d{:,i})';
end

% Mannually excluded MedI/C for low intensity
% WT
stat{2}(:,[6]) = [];
stat{3}(:,[1 3 6]) = [];
stat{4}(:,[3]) = [];

%% Function:    cmptCI
% Discription:  compute [50 97.5 2.5] for precentage using boostrap
function [stat,data] = cmptCI(TF,groupTF,nIteration)
% Input:    TF, logical vector
%           groupTF, logical matrix, one column per group
%           nIteration, number, number of iteration
% Output:   data, mat, data from boostrap
%           stat, mat, prctile of [50 2.5 97.5] per group, one group/row

data = [];

for k = 1:nIteration
    for i = 1:size(groupTF,2)
        % Find cells belongs the group
        I = find(groupTF(:,i));
        
        % Bootstrap (with replacement)
        I2 = randsample(numel(I),numel(I),true);
        I = I(I2);
        
        data(k,i) = sum(TF(I));
    end
end

% Compute precentage
data = data./sum(groupTF,1).*100;

% Percentile of broostrap data
stat = prctile(data,[50 97.5 2.5],1);
stat = stat';

end

%% Functon:     getLayerVoxelScale
% Description: Layer voxel scale for normalization 
function layerVoxelScale = getLayerVoxelScale(zL,mlapdDot2,ctxROI,layerFlat,refSetting)
% Depth range: 0-95%, for rolony quantification
% Some layer has more voxel, this probability is for plotting density
% Input:    zL, vector, bin edges for layer
%           mlapDot2, cell, rolony location on ML/AP/Depth-coordinates
%           ctxROI, cell, logicals of whether the rolony is within the ROI
%           layerFlat, mat, layer number with gradience in flatmap
%           coordinates
%           refSetting, struct
% Output:   layerVoxelScale, mat; row, ctxROI; column, bin of zL

% xyz of rolony on flatmap
xyz = TBS.mlapd2flatmapXYZ(mlapdDot2,refSetting,layerFlat);

% Get the area (CI) for datapoints on layerFlat 
layerVoxelScale = {};
for i = 1:size(ctxROI,2)
    ixyz = cellfun(@(X,Y) X(Y,:),xyz,ctxROI(:,i),'UniformOutput',false);
    ixyz = vertcat(ixyz{:});
    
    iL = layerFlat;
    % Area for most of the data: x (col) & y (row)
    iLim = ixyz(:,1);
    iLim = prctile(iLim,[2.5 97.5]);
    iLim = round(iLim);
    iL = iL(:,iLim(1):iLim(2),:);
    
    iLim = ixyz(:,2);
    iLim = prctile(iLim,[2.5 97.5]);
    iLim = round(iLim);
    iL = iL(iLim(1):iLim(2),:,:);      
    
    layerVoxelScale{i} = iL;
end

% Depth range: 0-95%
layerVoxelScale = cellfun(@(X) X(:,:,1:38),layerVoxelScale,'UniformOutput',false);
layerVoxelScale = cellfun(@(X) reshape(X,[],1),layerVoxelScale,'UniformOutput',false);
% Calculate fraction
layerVoxelScale = cellfun(@(X) histcounts(X,zL,'Normalization','probability'),...
    layerVoxelScale,'UniformOutput',false);
layerVoxelScale = vertcat(layerVoxelScale{:});
end

