% Get flatmap
% Cortical flatmap, use ML/AP/depth percentile coordinates system
% Current annotation map: Allen CCF3, 25 um/voxel

clc

% Settings

% TBS path
addpath('C:\Users\Li Yuan\MATLAB Drive\Script65A');

directory.main = 'D:\Li Yuan';

% Annotation Info ---------------------------------------------------------

% Setting for reference map
refSetting = TBS.getRefSetting(directory.main);

refScale = refSetting.refScale;

% Annotation map
annoMap = TBS.getRefMap('anno',refSetting);

% Annotation structure
annoStruct = refSetting.annoStruct;

sz = size(annoMap);

% Flatmap settings --------------------------------------------------------
% Depth %) for ML/AP-reference plate
refDepth = 50;

% Midline (Currently only support Allen)
midLineCol = sz(2)/2;
% Note 03082022, if only use one seems it bias to one site
midLineCol = [floor(midLineCol), ceil(midLineCol)];
if midLineCol(1) == midLineCol(2)
    midLineCol(2) = midLineCol(2)+1;
end

% Average filter 1: for cortical depth
% Too big will cause distortion for ctxDepth data
avgFilter1 = strel('sphere',3);
avgFilter1 = avgFilter1.Neighborhood;

% Avervage filter 2: for ML & AP
avgFilter2 = strel('sphere',7);
avgFilter2 = avgFilter2.Neighborhood;

%% Outter and innter edges
% Note, a good definition of outter and inner edge can save a lot of
% trouble later on

% Including Claustrum otherwise will have errors
roi = {'Claustrum','Isocortex'};
id = cellfun(@(X) TBS.findAnnoID(annoStruct,X),roi,'Uniformoutput',false);
id = vertcat(id{:});
ctx = ismember(annoMap,id);

% Outer edge --------------------------------------------------------------

% Two voxel wide
SE = strel('diamond',2);

% Layer 1
id = TBS.findAnnoStrID(annoStruct,'Layer 1');
roi = ismember(annoMap,id);

% Make midline empty to compute the outer edge in fold in area
annoMap2 = annoMap > 0;
annoMap2(:,midLineCol,:) = false;

outEdge = imdilate(ctx,SE) & imdilate(roi,SE) & ~annoMap2;

% Inner edge --------------------------------------------------------------
% The outer range of layer 6, touches fiber tract and subcortical nucleoi
% Exclude fold-in area

% One voxel wide
SE = strel('diamond',1);

% Layer 6, includs Claustrum
id = TBS.findAnnoStrID(annoStruct,'Layer 6');
roi = ismember(annoMap,id);
id = TBS.findAnnoID(annoStruct,'Claustrum');
roi = roi | ismember(annoMap,id);

innerEdge = imdilate(ctx,SE) & imdilate(roi,SE) & ~ctx;

% Fiber tract, exclude ee
% (this is needed for excluding edge contacts PIR, PIR contacts with both
% layer 1 and 6)
% 'cranial nerves' touches cortex in olfactoryball
roi = {'Cerebral nuclei','fiber tracts'};
id = cellfun(@(X) TBS.findAnnoID(annoStruct,X),roi,'Uniformoutput',false);
id = vertcat(id{:});
fiberTract = ismember(annoMap,id);
% (otherwise causes errors) exclude two fibers within cortical structure
roi = {'cranial nerves','corpus callosum, extreme capsule'};
id = cellfun(@(X) TBS.findAnnoID(annoStruct,X),roi,'Uniformoutput',false);
id = vertcat(id{:});
fiberTract = fiberTract & ~ismember(annoMap,id);

innerEdge = innerEdge & fiberTract;

% (Check point)
% test = outEdge + innerEdge.*2;
% MIJ.createImage(test);

%% Get reference columns 
% 1. Compute for each hemispher seperately, in case there is cross talk
% 2. Only do min --> max line, the fold in area will cause errors for
% max --> min line

refCol = {};
for i = 1:2
    if i == 1
        col = 1:midLineCol(1);
    elseif i == 2
        col = midLineCol(1)+1:sz(2);
    end
    
    [y,x,z] = TBS.find3(outEdge);
    XYZmin = [x,y,z];
    TF = ismember(x,col);
    XYZmin = XYZmin(TF,:);
    
    [y,x,z] = TBS.find3(innerEdge);
    XYZmax = [x,y,z];
    TF = ismember(x,col);
    XYZmax = XYZmax(TF,:);
    
    % Find the XYZmax for every XYZmin
    [D,I] = pdist2(XYZmax,XYZmin,'euclidean','Smallest',1);
    XYZmax = XYZmax(I',:);
    
    % Minimum 2 voxels away
    TF = D >= 2;
    XYZmin = XYZmin(TF,:);
    XYZmax = XYZmax(TF,:);
   
    % Get the line using two end points
    ixyz = TBS.interp1TwoDot(XYZmin,XYZmax,[0 100],0:100);

    % ROI: Half of the cortex
    roi = false(sz);
    roi(:,col,:) = ctx(:,col,:);
    
    % Voxels within the cortex
    TF = cellfun(@(X) TBS.xyz2v(X,roi),ixyz,'UniformOutput',false); 
    
    % Delete the dots out of cortex
    ixyz = cellfun(@(X,Y) X(Y,:),ixyz,TF,'UniformOutput',false); 
    
    TF = cellfun(@(X) find(X),TF,'UniformOutput',false);
    
    % Delete lines travel out of cortex (With discontinous index)
    TF2 = cellfun(@(X) all(diff(X)==1),TF);
    
    % Minimum 3 pixel for the reference line
    TF = cellfun(@numel,TF);
    TF = TF >= 3;

    ixyz = ixyz(TF & TF2);
    
    refCol{i,1} = ixyz;
end
 
refCol = vertcat(refCol{:});

%% Cortical depth percentage 

% Get depth percentage for each cortical column (0-100%)
% Distance to the first dot
c = cellfun(@(X) pdist2(X,X(1,:)),refCol,'Uniformoutput',false);
% Percentage of the total distance (1-100%)
c = cellfun(@(X) X./X(end).*99,c,'UniformOutput',false);
% Set the first point-1, last -100
c = cellfun(@(X) X + 1,c,'UniformOutput',false);

% Get reference value in cortical column
xyz = vertcat(refCol{:});
c =  vertcat(c{:});
xyz = round(xyz);
[xyz,~,ic] = unique(xyz,'rows');
% Mean if the voxel is founded more than one time
c = TBS.accumarrayMean(ic,c);

ctxDepthPrctile = TBS.xyzv2im(sz,xyz,c);

% Add inner and outter edge, 1 & 100% depth
SE = strel('diamond',1);
roi = imdilate(innerEdge,SE) & ctx;
roi = roi.*100;
ctxDepthPrctile = max(roi,ctxDepthPrctile);
roi = imdilate(outEdge,SE) & ctx;
roi = double(roi);
ctxDepthPrctile = max(roi,ctxDepthPrctile);

% Fill the empty space
ctxDepthPrctile = TBS.fillRefCtx(ctxDepthPrctile,ctx);

% Average filter 
ctxDepthPrctile = TBS.nonzeroAvgFiltCtx(ctxDepthPrctile,avgFilter1,midLineCol);
ctxDepthPrctile(~ctx) = 0;

% (Check point)
MIJ.createImage(ctxDepthPrctile);

%% ML

% ML reference plate location (fuse cortex)--------------------------------
% use >= refDepth% cortical depth 
refPlate = ctx;
refPlate(ctxDepthPrctile < refDepth & ctxDepthPrctile > 0) = false;

% Fuse section for subtraction
TF = annoMap > 0;
TF(ctxDepthPrctile < refDepth & ctxDepthPrctile > 0) = false;

% SE size, 1/2 image width
SE = round(sz(2)/2);
SE = strel('line',SE,0);

% Fuse the gap between cortex
refPlate = imclose(refPlate,SE);
TF = imclose(TF,SE);

SE = strel('sphere',1);
refPlate = ~imerode(refPlate,SE) & refPlate;
TF = ~imerode(TF,SE) & TF;

refPlate = refPlate & TF;

% ML reference plate value ------------------------------------------------

% Inital value: midline 
midline = false(sz);
midline(:,midLineCol,:) = true;

[y,x,z] = TBS.find3(refPlate & midline);
xyzInital = [x y z];

[y,x,z] = TBS.find3(refPlate & ~midline);
xyzRest = [x y z];

refML = TBS.getRefPlateContour(xyzInital,xyzRest,refScale);
% Add 1 to avoid 0
refML(:,end) = refML(:,end) + 1;
refML = TBS.xyzv2im(sz,refML(:,1:3),refML(:,end));

% Assign reference value to column
ctxML = TBS.refPlate2Column(refML,refCol);

% Add reference plate
TF = ctxML == 0;
ctxML(TF) = refML(TF);

% Fill the empty space
ctxML = TBS.fillRefCtx(ctxML,ctx);

% Average filter 
ctxML = TBS.nonzeroAvgFiltCtx(ctxML,avgFilter2,midLineCol);
ctxML(~ctx) = 0;

% (Check point)
MIJ.createImage(ctxML);

%% AP
% AP-direction (AP0): perpendicular to ML; median of all direction, as
% reference line 
% AP-value: distance to the reference line, along ML-contour

% AP-direction, reference line (perpendicular to ML) ======================
% Calculate the axis perpendicular to the ML
% Get AP-reference plate: perpendicular to ML contour
refAP0 = zeros(sz);
for j = 1:2
    % To avoid error, function allow one initial point per z, so do both
    % hemisphere seperately
    if j == 1
        col = 1:midLineCol(1);
    elseif j == 2
        col = (midLineCol(1)+1):sz(2);
    end
    
    im = refML(:,col,:);
    im = TBS.refAPfromML(im,refScale);
    
    refAP0(:,col,:) = im;
end

% Compute the whole volume to get the reference line, smoother ------------
% Assign reference value to column
ctxAP0 = TBS.refPlate2Column(refAP0,refCol);

% Add reference plate
TF = ctxAP0 == 0;
ctxAP0(TF) = refAP0(TF);

% Fill the empty space
ctxAP0 = TBS.fillRefCtx(ctxAP0,ctx);

% Average filter 
ctxAP0 = TBS.nonzeroAvgFiltCtx(ctxAP0,avgFilter2,midLineCol);
ctxAP0(~ctx) = 0;

% (Check point)
MIJ.createImage(ctxAP0);

% %% AP-value, distance based =============================================== 

% AP reference plate location (no fusion) 
refPlate = ctx;
refPlate(ctxDepthPrctile < refDepth & ctxDepthPrctile > 0) = false;

TF = annoMap > 0;
TF(ctxDepthPrctile < refDepth & ctxDepthPrctile > 0) = false;

SE = strel('sphere',1);
refPlate = ~imerode(refPlate,SE) & refPlate;
TF = ~imerode(TF,SE) & TF;

refPlate = refPlate & TF;

% AP reference plate value 
refAP = ctxAP0;
refAP(~refPlate) = 0;

% Reference value, median without min/max value
v = nonzeros(refAP);
v(v == 1 | v == max(v)) = [];

refV = median(v);

% Region posterior to reference lines 
TFpost = ctxAP0 > refV;

% Reference line, 1-voxel anterior
SE = strel('sphere',1);
roi = imdilate(TFpost,SE) & ~TFpost & refAP;

% Voxel in reference line
xyzRef = ctxML;
xyzRef(~roi) = 0;
[y,x,z,v] = TBS.find3(xyzRef);
xyzRef = [x y z v];

% Delete reference line
xyzRest = ctxML;
xyzRest(roi | ~refPlate) = 0;
[y,x,z,v] = TBS.find3(xyzRest);
xyzRest = [x y z v];

% Split into left and right hemisphere
TF = xyzRef(:,1) <= midLineCol(1);
xyzRef = [{xyzRef(TF,:)}; {xyzRef(~TF,:)}];
TF = xyzRest(:,1) <= midLineCol(1);
xyzRest = [{xyzRest(TF,:)}; {xyzRest(~TF,:)}];

% Get ref value from left and right hemisphere seperately
xyz = cellfun(@(X,Y) TBS.getRefPlateContour(X,Y,refScale),...
    xyzRef,xyzRest,'UniformOutput',false);
xyz = vertcat(xyz{:});
% Add 1 to avoid 0
xyz(:,end) = xyz(:,end)+1;

% AP reference plate
refAP = TBS.xyzv2im(sz,xyz(:,1:3),xyz(:,end));

% Exclude any unassigned dots, in case error
refPlate = refPlate & refAP;

% Minus if before the reference line
TF = ctxAP0 < refV;
refAP(TF) = refAP(TF).*(-1);

% Convert all value beyond 0
refAP = refAP - min(refAP,[],'all') + 1;
refAP(~refPlate) = 0;

% Assign reference value to column
ctxAP = TBS.refPlate2Column(refAP,refCol);

% Add reference plate
TF = ctxAP == 0;
ctxAP(TF) = refAP(TF);

% Fill the empty space
ctxAP = TBS.fillRefCtx(ctxAP,ctx);

% Average filter 
ctxAP = TBS.nonzeroAvgFiltCtx(ctxAP,avgFilter2,midLineCol);
ctxAP(~ctx) = 0;

% (Check point)
MIJ.createImage(ctxAP);
