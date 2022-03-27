% Get flatmap
% Cortical flatmap, use ML/AP/depth percentile coordinates system
% Current annotation map: Allen CCF3, 25 um/voxel
% Note: left and right hemisphere is not entirely the same, this different
% pass down to flatmap
% 03202022: change to use pca to determine axes

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

save(fullfile(directory.main,'ctxDepthPrctile.mat'),'ctxDepthPrctile');

% (Check point)
MIJ.createImage(ctxDepthPrctile);

%% Reference plate & axes

% Reference plate about refDepth, fused cortex
% use >= refDepth% cortical depth
refPlate = ctx;
refPlate(ctxDepthPrctile < refDepth & ctxDepthPrctile > 0) = false;

% Fuse section for subtraction
TF = annoMap > 0;
TF(ctxDepthPrctile < refDepth & ctxDepthPrctile > 0) = false;

% SE size for fuse cortex, 1/2 image width
SE = round(sz(2)/2);
SE = strel('line',SE,0);

% Fuse the gap between cortex
refPlate = imclose(refPlate,SE);
TF = imclose(TF,SE);

% Reference plate
SE = strel('sphere',2);
refPlate = ~imerode(refPlate,SE) & refPlate;
TF = ~imerode(TF,SE) & TF;

refPlate = refPlate & TF;

% Get axes using pca ------------------------------------------------------
% Use the right hemisphere to get axis
% axes: 1. AP; 2. ML

[y,x,z] = TBS.find3(refPlate);
xyz = [x y z];

xyz2 = zeros(size(xyz,1),2);
for j = 1:2
    % Right and left hemisphere
    if j == 1
        col = 1:midLineCol(1);
    elseif j == 2
        col = (midLineCol(1)+1):sz(2);
    end
    
    TF = ismember(xyz(:,1),col);
    ixyz = xyz(TF,:);
    
    % Use the left hemisphere to get axes
    % Fliplr for the right hemisphere
    if j == 1
        coeff = pca(ixyz);
    elseif j == 2
        ixyz(:,1) = ixyz(:,1).*-1;
    end
    
    % Get the value for the first two axes: AP, ML
    ixyz = ixyz*coeff(:,1:2);
    
    % Swithch the ML-axes from M->L
    if coeff(1,2) > 0
        ixyz(:,2) = ixyz(:,2).*(-1);
    end
    
    % Set minimum value to 1
    ixyz = ixyz - min(ixyz,[],1)+1;
    
    xyz2(TF,:) = ixyz;
end

% Median value as reference value
refV = median(xyz2,1);

%% AP & ML axis

refML = []; refAP = []; ctxML = []; ctxAP = [];

% Do AP-axis first, the definition of AP is more clear than ML
for a = 1:2
    refPlateV = TBS.xyzv2im(sz,xyz,xyz2(:,a));
    
    % Region posterior to reference lines
    TFpost = refPlateV > refV(a);
    
    % Reference line, 1-voxel anterior
    SE = strel('sphere',1);
    roi = imdilate(TFpost,SE) & ~TFpost & refPlate;
    
    % Voxel in reference line/rest of the regions
    if a == 1        
        [y,x,z] = TBS.find3(roi);
        xyzRef = [x y z];
        
        [y,x,z] = TBS.find3(refPlate & ~roi);
        xyzRest = [x y z];
        
    elseif a == 2
        % AP-axis use ML value contour. Optional, just want to be more accurate
        TF = refAP;
        TF(~roi) = 0;
        [y,x,z,v] = TBS.find3(TF);
        xyzRef = [x y z v];
        
        TF = refAP;
        TF(roi) = 0;
        [y,x,z,v] = TBS.find3(TF);
        xyzRest = [x y z v];
    end
    
    % xyzv of the reference plate
    refAx = TBS.getRefPlateContour(xyzRef,xyzRest,refScale);
    % Add 1 to avoid 0
    refAx(:,end) = refAx(:,end) + 1;
    refAx = TBS.xyzv2im(sz,refAx(:,1:3),refAx(:,end));
    
    % Exclude any unassigned dots, in case error
    refPlate2 = refPlate & refAx;
    
    % Minus if before the reference line
    TF = refPlateV < refV(a);
    refAx(TF) = refAx(TF).*(-1);
    
    % Convert all value beyond 0
    refAx = refAx - min(refAx,[],'all') + 1;
    refAx(~refPlate2) = 0;
    
    % Assign reference value to column
    ctxAx = TBS.refPlate2Column(refAx,refCol);
    
    % Add reference plate
    TF = ctxAx == 0;
    ctxAx(TF) = refAx(TF);
    
    % Fill the empty space
    ctxAx = TBS.fillRefCtx(ctxAx,ctx);
    
    % Average filter
    ctxAx = TBS.nonzeroAvgFiltCtx(ctxAx,avgFilter2,midLineCol);
    ctxAx(~ctx) = 0;
    
    MIJ.createImage(ctxAx);
    
    if a == 2
        refML = refAx; ctxML = ctxAx;
        save(fullfile(directory.main,'ctxML.mat'),'ctxML');
    elseif a == 1
        refAP = refAx; ctxAP = ctxAx;
        save(fullfile(directory.main,'ctxAP.mat'),'ctxAP');
    end
    
end
