% bioinformative toolbox
%   02112022, ver 1i, change grouping for CFIT
%   03222022: flatmap2
%   04022022: no soma MLAP exclusion

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

% Cortical ML/AP/Depth
% getMLAPD(xyz,ctxML,ctxAP,ctxDepthPrctile)
mlapdDot = TBS.getMLAPD(xyzDot,ctxML,ctxAP,ctxDepthPrctile);

% Only include rolony within the edges
fh = @(X) X(:,3)>= min(rolonyDepthEdges) & X(:,3)<= max(rolonyDepthEdges);
mlapdDot = cellfun(@(X) X(fh(X),:),mlapdDot,'UniformOutput',false);

% Soma location -----------------------------------------------------------
load('somaLocationRef.mat');
xyzSoma = somaLocation;

TF = any(xyzSoma,2);
mlapdSoma = zeros(size(xyzSoma));
mlapdSoma(TF,:) = TBS.getMLAPD(xyzSoma(TF,:),ctxML,ctxAP,ctxDepthPrctile);

TF = any(xyzSoma,2);
disp(['Soma included in analysis: ',num2str(sum(TF))]);

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
center = [6400 8300];  

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

%% Fig 2. Plot cells sololy inervate single cortical cluster ==============
% Uses a subset of cells for axon model
% Plot cells projects to a cluster, and most of the projection in the
% hemisphere is within the cluster; also soma-positive

% Cortical cluster to plot (i.e. 2: [7050, 6950])
iCluster = 1;

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

% Plot BC model -----------------------------------------------------------
% Shuffle order for random colorcode
row = TBS.shuffleRows(row);
xyzSoma2 = xyzSoma(row,:); xyzDot2 = xyzDot(row);
[im, outline] = TBS.BCmodel(xyzDot2,xyzSoma2,[],2,annoMap,refSetting);

MIJ.createImage(cat(3,im,outline));
MIJ.run("Stack to Hyperstack...", "order=xyzct channels=2 slices=1056 frames=1 display=Composite");

% SupFig 6. Plot rolony of individual BC ===============================
% Need to save image mannually using 3D viwer from ImageJ

% Only plot the z with data for better visualization
% Get the limits of z-axis
zLim = vertcat(xyzDot{:});
zLim = zLim(:,3);
zLim = floor(min(zLim)): ceil(max(zLim));

% Brain outline
outline = TBS.getBrainOutline(annoMap,5);
outline = outline(:,:,zLim);
outline = single(outline);

i = 14;  % cell to choose mannually

% Plot rolony & soma 
% Rolony & soma colorcode with z-axis
im = TBS.xyzv2im(size(annoMap),xyzDot2{i},xyzDot2{i}(:,3));

SE = strel('sphere',2); 
im = imdilate(im,SE); % Make rolony bigger

if any(xyzSoma2(i,:))
    im2 = TBS.xyzv2im(size(annoMap),xyzSoma2(i,:),xyzSoma2(i,3));
    SE = strel('sphere',5); % Make soma bigger
    im2 = imdilate(im2,SE);
    im = max(im, im2);
end

im = im(:,:,zLim);

MIJ.createImage(cat(3,im,outline));
MIJ.run("Stack to Hyperstack...", "order=xyzct channels=2 slices=97 frames=1 display=Composite");
% Then run macro in ImgeJ in the folder

%% SupFig 7, Rolony doesnt belongs to 5 brain regions
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

%% SupFig.7 Projection strength heatmap

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

colormap(flipud(gray));
set(gcf,'Position',[100 100 180 500]);

%% SupFig 7. Soma layer distribution
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

% Shuffle soma between groups
[mlapdSoma2,I] = TBS.shuffleRows(mlapdSoma2); c = c(I);

figure; scatter(mlapdSoma2(:,1),mlapdSoma2(:,3),10,c,'filled');
yline(0,'k:'); % Pia
ylabel('Soma depth (%)');
TBS.axLabelSettings('Myriad Pro',12);
g = gca; g.YDir = 'reverse'; g.YLim = [0 100];g.XTick = [];
colormap(cmap4); pbaspect([1.25 1 1]);
set(gcf,'Position',[100 100 400 250]);

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

% Plot rolony location
im = TBS.xyzv2im(size(outline),xyzDot2,c);

im = uint8(im);
outline = uint8(regionOutline);
im = cat(3,im,outline);

MIJ.createImage(im);
MIJ.run("Stack to Hyperstack...", "order=xyzct channels=2 slices=528 frames=1 display=Composite");
% Plot in ImageJ

%% Fig 4-5. Prepare data for cortical analysis =======================================================================

% Min count per region
minCount = 5;

% Whether is a LatC control
isCtrl = false;

% Exclude rolony within 95% radius from the injection center
% (The current setting will leave some layer 1 rolony around injection
% site)
mlapdDot2 = TBS.nearSomaExcl(mlapdDot,mlapdSoma,95,isCtrl);

% Only include IT cells, with cortical projections
TF = cellfun(@(X) size(X,1),mlapdDot2);
TF = TF >= minCount & ismember(cellType,3:4);

% ML/AP/Depth coordinates
mlapdDot2 = mlapdDot2(TF);
mlapdSoma2 = mlapdSoma(TF,:);
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

return

%% Fig 4. Graph of corticl targets ========================================
% Seperation of Med & Lat region, outlie of region boundaries

figure; TBS.plotRegionRef(mlapdSoma,regionOutlineFlat); 
TBS.flatmapSetting(flatmapYLim);
xline([-mlBoundary 0 mlBoundary],'k:','LineWidth',2);
set(gca,'Visible','off');

%% SupFig 8. Batch cortical depth =========================================

% Binning size in ML-AP plate
binSize = 250;
% Minimum rolony count per bin
minCountBin = 50;

% Binning
xy = vertcat(mlapdDot2{:}); 
z = xy(:,3); 

xy = xy(:,1:2);
xy = round(xy./binSize).*binSize;

[xy,~,ic] = unique(xy,'rows');

% Rolony distribution along the depth per bin -----------------------------
d = [];
for i = 1:size(xy,1)
    % Get all the dpeth in the bin
    iD = z(ic == i);
    d(i,:) = histcounts(iD,rolonyDepthEdges,'Normalization','probability');    
end

% Normalized to precentage
d = d.*100;

% Exclude bins with too few counts
n = accumarray(ic,1);
TF = n >= minCountBin;
xyBin = xy(TF,:); d = d(TF,:);

disp(['(Report) Cortical bins plotted: ', num2str(sum(TF))]);

% Visualize the cortical column -------------------------------------------
% Grouping 
idx = TBS.kmeansDepthHist(d,2);

figure; scatter(xyBin(:,1),xyBin(:,2),60,idx,'filled');
hold on; TBS.plotRegionRef(mlapdSoma,regionOutlineFlat); 
TBS.flatmapSetting(flatmapYLim);
colormap(cmapMG);
xline([-mlBoundary 0 mlBoundary],'k:','LineWidth',2);

% Heatmatp 
[idx,I] = sort(idx);
d = d(I,:)';
% Exclude the 1st bin
d = d(2:end,:);
figure; h = heatmap(d);
TBS.depthHeatmapSetting(h,d,rolonyDepthEdges(2:end));
S = struct(h); xline(S.Axes,find(diff(idx)>0),'k');
set(gcf,'Position',[100 100 900 200]);

%% SupFig 8. Single cell projection depth =================================
% Single cell laminar pattern in each area

% Get laminar pattern in each cortical ROI, per cell ----------------------
% LatI, MedI, MedC, LatC
d = {};
for i = 1:size(ctxROI,2)
    % Rolony depth
    iD = cellfun(@(X,Y) X(Y,3),mlapdDot2,ctxROI(:,i),'UniformOutput',false);
    
    % Exclude cells without enough count
    n = cellfun(@numel,iD);
    TF = n >= minCount;
    
    % Frequency distribution (%)
    iD = cellfun(@(X) histcounts(X,rolonyDepthEdges,...
        'Normalization','probability'),iD,'UniformOutput',false);
    iD = vertcat(iD{:});
    iD = iD.*100;
    
    d{i,1} = iD(TF,:);
end

% Sort for each region ----------------------------------------------------
for i = 1:numel(d)
    iD = d{i};
    
    % Clustering using kmeans, sort index using depth
    % kmeansDepthHist(X,k,C)
    idx = TBS.kmeansDepthHist(iD,7);
    
    [~,I] = sort(idx);
    
    d{i} = iD(I,:);    
end

% Plot --------------------------------------------------------------------
% Mark the number for each group to draw a line
n = cellfun(@(X) size(X,1),d);
disp(['(Report) Single cell projection plotted: ', num2str(n')]); 
n = cumsum(n);

% Row, cortical depth; one column per barcode per region
d = vertcat(d{:})';

% Exclude the 1st bin
d = d(2:end,:);

% Plot
figure; h = heatmap(d);
TBS.depthHeatmapSetting(h,d,rolonyDepthEdges(2:end));
S = struct(h); xline(S.Axes,n,'k-');
set(gcf,'Position',[100 100 1800 200]);

%% Fig. 4 Depth comparison in line plot ===================================

ctxROI2 = ctxROI;
% Met as 5th col
ctxROI2(:,5) = cellfun(@(X,Y) X | Y, ctxROI2(:,2),ctxROI2(:,3),'UniformOutput',false);

% Column number for each par
groupIdx = [];
groupIdx(1,:) = [2,3];  % biMed
groupIdx(2,:) = [1,4];  % biLat
groupIdx(3,:) = [5,1];  % Med-LatI
groupIdx(4,:) = [5,4];  % Med-LatC

% Cells project to both region (>= minimum counts)
n = cellfun(@sum,ctxROI2);
group = {};
for i = 1:size(groupIdx,1)    
    group{i} = n(:,groupIdx(i,:));
end
group = cellfun(@(X) all(X >= minCount,2),group,'UniformOutput',false);
group = horzcat(group{:});

disp('Cell counts: biMed, biLat, Med-LatI, Med-LatC');
disp(sum(group));

% Fig. 4 Cumultive distribution diff between two regions ==================

% Frequency distribution of depth
d = cellfun(@(X,Y) Y(X,3),ctxROI2,...
    repmat(mlapdDot2,1,size(ctxROI2,2)),'UniformOutput',false);
d = cellfun(@(X) histcounts(X,rolonyDepthEdges,'Normalization','cdf'),...
    d,'UniformOutput',false);
d = cellfun(@(X) X.*100,d,'UniformOutput',false);

% xy axis label
axisLabel = {'MedC - MedI';'LatC - LatI';'LatI - Med';'LatC - Med'};
axisLabel = cellfun(@(X) [X,' difference (%)'],axisLabel,'UniformOutput',false);

% Line alpha 
a = [0.3, 0.05, 0.075, 0.1];

for i = 1:size(groupIdx,1)
    % Depth from the region pair
    TF = groupIdx(i,:);
    iD = d(:,TF);
    
    % Depth from the BC project to both regions, abv threshold
    TF = group(:,i);
    iD = iD(TF,:);
    
    % Depth difference
    iD = cellfun(@(X,Y) Y-X,iD(:,1),iD(:,2),'UniformOutput',false);
    iD = vertcat(iD{:});
    
    z = rolonyDepthEdges(2:end);
    
    figure; hold on;
    for j = 1:size(iD,1)        
        plot(iD(j,:),z,'Color',[0 0 0 a(i)]);
    end
    xlabel(axisLabel{i}); ylabel('Projection depth (%)');
    TBS.axLabelSettings('Myriad Pro',15);
    g = gca; g.XLim = [-100 100]; g.YLim = [0 95]; g.YDir = 'reverse';
    g.XTick = -100:50:100; g.YTick = 0:20:100;
    pbaspect([1 1 1]); set(gcf,'Position',[100 100 300 300]);   
    
    % Median and CI
    iD = prctile(iD,[50 2.5 97.5]);
    plot(iD(1,:),z,'r','LineWidth',3);
    plot(iD(2:3,:),z,'r--','LineWidth',1.5);    
end

%% Fig 4. Depth comparison in scatter plot ================================
% Frequency distribution, so both region is compariable even when rolony
% number are differnct

% Color map, 5 color
cmap = [235 0 0; 0 0 66; 0 66 255; 0 192 0; 0 33 156];
cmap = cmap./255;

% Max barcode to be displayed
maxBC = 800;

% Frequency distribution of depth
d = cellfun(@(X,Y) Y(X,3),ctxROI2,...
    repmat(mlapdDot2,1,size(ctxROI2,2)),'UniformOutput',false);
d = cellfun(@(X) histcounts(X,rolonyDepthEdges,'Normalization','probability'),...
    d,'UniformOutput',false);
d = cellfun(@(X) X.*100,d,'UniformOutput',false);

for i = 1:size(groupIdx,1)
    % Depth from the region pair
    TF = groupIdx(i,:);
    iD = d(:,TF);
    
    % Depth from the BC project to both regions, abv threshold
    TF = group(:,i);
    iD = iD(TF,:);
                
    im1 = vertcat(iD{:,1});
    im2 = vertcat(iD{:,2});
    
    % Sort using the 2nd projection area, for visualization
    idx = TBS.kmeansDepthHist(im2,7);
    [~,I] = sort(idx);    
    im1 = im1(I,:);
    im2 = im2(I,:);    
    
    % Only sample 800 max otherwise image can be too long    
    im = cat(3,im1',im2');
    
    % Need to do this after sorting
    if size(im,2) > maxBC
        I = randsample(size(im,2),maxBC);
        I = sort(I,'ascend');
        im = im(:,I,:);
    end
        
    % Exclude the 1st bin
    im = im(2:end,:,:);

    % Get scatter plot from image output ----------------------------------
    [y,x,z,v] = TBS.find3(im);
    
    % Colorcode region
    c = cmap(groupIdx(i,z),:);
    
    % Set max to 10%, rounded
    v = v./10;
    v = min(v,1);
    v = round(v,1);
        
    % (To speed up)
    [v,~,ic] = unique(v);    
    figure; hold on;
    for j = 1:numel(v)
        TF = ic == j;
        scatter(x(TF),y(TF),2,c(TF,:),'filled','MarkerFaceAlpha',v(j));
    end
    ylabel('Projection depth (%)','FontSize',12);    
    g = gca; g.YDir = 'reverse'; g.YLim = [0 rolonyDepthEdges(end)];     
    yticks([0 rolonyDepthEdges(end)]);
    g.XLim = [1 max(x)]; 
    g.XTick = [1:100:max(x)-1,max(x)]; g.XTickLabel(1:end-1) = {[]};
    % The values in units normalized relative to the longest the axes
    g.TickLength = [1000/max(x)*0.01 0.01];
    g.LineWidth = 1;
    daspect([1 1 1]); set(gcf,'Position',[100 100 1200 150]);
end

%% SupFig 8. Depth difference quantification
% Discription: compare the laminar difference between two region of the
% same barcoded neuron
% 02252022: delete FDR

nIteration = 2000;

stat = {};
for i = 1:size(groupIdx,1)
    % Depth from the region pair
    TF = groupIdx(i,:);
    d = cellfun(@(X,Y) Y(X,3),ctxROI2(:,TF),...
        repmat(mlapdDot2,1,2),'UniformOutput',false);
    
    % Depth from the BC project to both regions, abv threshold
    TF = group(:,i);
    d = d(TF,:);
    
    % Compared difference using two-sample Kolmogorov-Smirnov test
    [~,~,ks2stat] = cellfun(@(X,Y) kstest2(X,Y),d(:,1),d(:,2));
        
    % Calculate p-value ---------------------------------------------------
    sz = cellfun(@numel,d);
    sz = mat2cell(sz,ones(size(sz,1),1),2);
    
    d = cellfun(@(X,Y) [X;Y],d(:,1),d(:,2),'UniformOutput',false);
    
    shuffD = [];
    parfor j = 1:nIteration % parfor
        % Shuffle data        
        jd = cellfun(@(X) TBS.shuffleRows(X),d,'UniformOUtput',false);
        jd = cellfun(@(X,Y) mat2cell(X,Y,1)',jd,sz,'UniformOUtput',false);
        jd = vertcat(jd{:});
        
        [~,~,shuffD(:,j)] = cellfun(@(X,Y) kstest2(X,Y),jd(:,1),jd(:,2));
    end    
    
    % Calculate p-value
    p = shuffD >= ks2stat;
    p = sum(p,2)./size(shuffD,2);
        
    stat{i} = p;
    disp(i);
end
% Plot in Prism

%% Fig.4 Stat sumary of cell with different projection
% Get distribution of <= 0.05 fraction by bootstraping

% Number of bootstrap
nIteration = 2000; 

bootStat = [];
for i = 1:numel(stat)
    p = stat{i};
    % Precentage is significatn different
    p = p <= 0.05;
    
    n = numel(p);
    for j = 1:nIteration
        % Sample with replacement
        jp = randsample(n,n,true);
        jp = p(jp);
        
        bootStat(j,i) = sum(jp)/n.*100;
    end
end

prctile(bootStat,[50,97.5,2.5])'

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
[im, outline] = TBS.BCmodel(xyzDot3,xyzSoma3,[],scaleFactor,annoMap,refSetting);

% Get the limits of z-axis
zLim = vertcat(xyzDot{:});
zLim = zLim(:,3);
zLim = zLim.*scaleFactor;
zLim = floor(min(zLim)): ceil(max(zLim));

% Output in ImageJ
MIJ.createImage(cat(3,im(:,:,zLim),outline(:,:,zLim)));
MIJ.run("Stack to Hyperstack...", "order=xyzct channels=2 slices=191 frames=1 display=Composite");

%% Fig5. Lat 3 groups, thresholding soma =====================================================

% soma threshold for diving into two groups
fig5threshold = [0 35 60 100];

% Lat+ cells: either Lat pass threshold
TF = cellfun(@sum,ctxROI(:,[1 4]));
TF = TF >= minCount;
TF = any(TF,2);

% Group1, No Lat-projection or no soma; 2, 0-35; 3: 35-60; 4, > 60
group = ~TF | ~any(mlapdSoma2,2);
for i = 1:numel(fig5threshold)-1
    TF = mlapdSoma2(:,3) > fig5threshold(i) & mlapdSoma2(:,3) <= fig5threshold(i+1);
    group(:,i+1) = TF & ~group(:,1);
end

fig5SomaGroup = group;
% Report cell number
disp('Grouping using soma depth. cell count: excluded, group1,2,3: ')
sum(group)

% %% Fig. 5 Soma location laminar distribution ----------------------------

% Scatter plot
[~,c] = max(group,[],2);
TF = c > 1;

TBS.plotSoma(mlapdSoma2(TF,:),c(TF,:)); alpha 0.5;
yline(fig5threshold(2:3),'k:');
colormap(cmapROB); pbaspect([1.4 1 1]);
set(gcf,'Position',[100 100 400 250]);

%% (Report) Fig.5 Upper-projecting Layer6b cells


%% Fig5. single cell projection model and flatmap 
% Discription: example cortical projection on flatmap

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

% Group to visualizing
i = 3;

% 04042022,Define using the minimum group count
n = min(sum(group(:,2:4)));

% Color code
[~,c] = max(group,[],2);
c = c -1;

I = find(c == i);
I = TBS.shuffleRows(I);
I = I(1:n);

c = TBS.repmat2cell(c,xyzDot2);

xyz = vertcat(xyzDot2{I});
c = vertcat(c{I});

% Image of rolony location
im = TBS.xyzv2im(size(annoMap),xyz,c);

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

% Visualized in ImageJ
MIJ.createImage(cat(3,im,outline));
MIJ.run("Stack to Hyperstack...", "order=xyzct channels=2 slices=97 frames=1 display=Composite");
% Output in ImageJ

%% SupFig9. Rolony depth quantification -----------------------------------
% Discription: quantification rolony group distribution along the depth
% 5% depther per bin

r = 4; % 1-LatI; 4-LatC

z = -5:5:rolonyDepthEdges(end);

d = cellfun(@(X,Y) X(Y,3), mlapdDot2,ctxROI(:,r),'UniformOutput',false);

stat = [];
for i = 2:size(group,2)
    TF = group(:,i);
    d2 = vertcat(d{TF});
    stat(end+1,:) = histcounts(d2,z);
end

% Sum of per goup is 100
stat = stat./sum(stat,2).*100;

z = TBS.getMedEdges(z);
stat = [stat',z'];
% Plot in Prism

%% Fig5. Relationship between soma and projection depth
% Plot soma depth along with projection pattern for soma+ cells
% Cannot directly save it in MATLAB, need to save it in ImageJ

r = 4; % LatI,MedI,MedC,LatC

% Threshold for lower projections
threshold = 60;

xLabel = fig5threshold;

% Projection binning
% 02282022: seems binning size = 5 is better
binSize = 5;
z = -binSize:binSize:rolonyDepthEdges(end);

d = cellfun(@(X,Y) X(Y,3), mlapdDot2,ctxROI(:,r),'UniformOutput',false);

n = cellfun(@sum,ctxROI(:,r));
TF = n >= minCount & ~group(:,1);

d = d(TF);
somaD = mlapdSoma2(TF,3);

% Sort using soma depth
[somaD,I] = sort(somaD,'ascend');
d = d(I);

% Projection percentage
stat = cellfun(@(X) histcounts(X,z,'Normalization','probability'),...
    d,'UniformOutput',false);
stat = vertcat(stat{:});
stat = stat.*100;

% Proportion above threshold
p = cellfun(@(X) sum(X > threshold)./numel(X),d);
p = p.*100;

disp(['Cell count: ', num2str(length(stat))]);

% Plot --------------------------------------------------------------------
stat = stat';
% Delete the first bin
stat = stat(2:end,:);

% Bin soma, get median p per bin
somaD2 = round(somaD./5).*5;
[somaD2,~,ic] = unique(somaD2);
p = accumarray(ic,p,[],@median);
% Delete the bins have too few cells
n = accumarray(ic,1);
TF = n >= 5;
somaD2 = somaD2(TF);
p = p(TF);

n = cumsum(n); 
n = n(TF);

% Normalize p-value to the same scale
p = p.*0.95./binSize+0.5;

% Soma depth, for x-axis
x = 1;
for i = xLabel(2:end-1)
    x(end+1) = find(somaD >= i,1,'first');
end
x(end+1) = numel(somaD);

% Heatmap for projection frequency
figure; h = heatmap(stat); TBS.depthHeatmapSetting(h,stat,z(2:end));

% Plot soma location
S = struct(h); hold(S.Axes,'on'); 
plot(S.Axes,n,p,'r','LineWidth',2);

set(gcf,'Position',[100 100 1200 200]); 
set(gcf,'Position',[100 100 550 200]); 
colorbar off

% Right-axis
yyaxis(S.Axes,'right');
ylabel(S.Axes,'Lower layer rolony (%)');
S.Axes.YTick = [0 1];
S.Axes.YTickLabel =[{'0'}; {'100'}];
S.Axes.YDir = 'reverse';
S.Axes.YColor = [1 0 0]

% Xtick for soma depth
S.Axes.XTickLabel(x) = num2cell(xLabel);
S.Axes.XTickLabelRotation = 0;

%% SupFig9. Relationship between soma and lower layer projection (%)

r = 4; % LatI,LatC

d = cellfun(@(X,Y) X(Y,3), mlapdDot2,ctxROI(:,r),'UniformOutput',false);

n = cellfun(@sum,ctxROI(:,r));
TF = n >= 10 & ~group(:,1);

d = d(TF);
somaD = mlapdSoma2(TF,3);

% Proportion above threshold
p = cellfun(@(X) sum(X > threshold)./numel(X),d);
p = p.*100;

% Scatter plot
figure; scatter(somaD,p,5,'k','filled'); alpha 0.3;
g = gca; g.YDir = 'reverse'; g.YLim = [0 100]; g.XLim = [0 100];
g.XTick = fig5threshold;
xlabel('Soma depth (%)'); ylabel('Lower layer rolony (%)');
TBS.axLabelSettings('Myriad Pro',12);
set(gcf,'Position',[100 100 500 200]);

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

r = 1; % 1-LatI, 4-LatC

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
stat = cmptCI(TF,group(:,2:end),nIteration);

% Str+ cells --------------------------------------------------------------
% Rolony within str
roi = regVoxel{regCol.str};
roi = cellfun(@(X) TBS.xyz2v(X,roi),xyzDot2,'UniformOutput',false);

n = cellfun(@sum,roi);
TF = n >= minCount;

% (Report) Str+ percentage, plot in Prism
stat = cmptCI(TF,group(:,2:end),nIteration);

% LatC+ cells -------------------------------------------------------------
roi = ctxROI(:,4);

n = cellfun(@sum,roi);
TF = n >= minCount;

% (Report) LatC+ percentage, plot in Prism
stat = cmptCI(TF,group(:,2:end),nIteration);

%% SupFig9. Split Lat+ no-soma cells into two groups
% Upper & lower projecting

group2 = [];
for r = [1 4]
    % Median projection depth in the region
    d = cellfun(@(X,Y) X(Y,3), mlapdDot2,ctxROI(:,r),'UniformOutput',false);
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

%% SupFig9, Med/Str/LatC in cells w/o soma

nIteration = 2000;

% Med+ cells --------------------------------------------------------------
roi = cellfun(@(X,Y) X | Y, ctxROI(:,2),ctxROI(:,3),'UniformOutput',false);

% Rolony count pass threshold
n = cellfun(@sum,roi);
TF = n >= minCount;

% (Report) Med+ percentage, plot in Prism
stat = cmptCI(TF,group2,nIteration);

% Str+ cells --------------------------------------------------------------
% Rolony within str
roi = regVoxel{regCol.str};
roi = cellfun(@(X) TBS.xyz2v(X,roi),xyzDot2,'UniformOutput',false);

n = cellfun(@sum,roi);
TF = n >= minCount;

% (Report) Str+ percentage, plot in Prism
stat = cmptCI(TF,group2,nIteration);

% LatC+ cells -------------------------------------------------------------
roi = ctxROI(:,4);

n = cellfun(@sum,roi);
TF = n >= minCount;

% (Report) LatC+ percentage, plot in Prism
stat = cmptCI(TF,group2,nIteration);

%% SupFig9. Focal projection in cells w/o soma location

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

%% Test, plot soma distribution

z = -5:5:95;

% Med+ cells --------------------------------------------------------------
roi = cellfun(@(X,Y) X | Y, ctxROI(:,2),ctxROI(:,3),'UniformOutput',false);

% Rolony count pass threshold
n = cellfun(@sum,roi);
TF = n >= minCount;

TF2 = any(mlapdSoma2,2);

stat = {};
stat{1} = mlapdSoma2(TF2 & TF,3);
stat{2} = mlapdSoma2(TF2 & ~TF,3);

% Soma composition
stat = cellfun(@(X) histcounts(X,z),stat,'UniformOutput',false);
stat = vertcat(stat{:});
stat = stat./sum(stat,1).*100;

y = TBS.getMedEdges(z);

stat = [stat; y]';

% Exclude bins with too few soma
TF = histcounts(mlapdSoma2(TF2,3),z);
TF = TF >= 10;

stat = stat(TF,:);

%% Function: cmptCI
% Discription:  compute [50 97.5 2.5] for precentage using boostrap
function stat = cmptCI(TF,groupTF,nIteration)
% Input:    TF, logical vector
%           groupTF, logical matrix, one column per group
%           nIteration, number, number of iteration
% Output:   stat, mat, prctile of [50 2.5 97.5] per group, one group/row

stat = [];

for k = 1:nIteration
    for i = 1:size(groupTF,2)
        % Find cells belongs the group
        I = find(groupTF(:,i));
        
        % Bootstrap (with replacement)
        I2 = randsample(numel(I),numel(I),true);
        I = I(I2);
        
        stat(k,i) = sum(TF(I));
    end
end

% Compute precentage
stat = stat./sum(groupTF,1).*100;

% Percentile of broostrap data
stat = prctile(stat,[50 97.5 2.5],1);
stat = stat';

end
