% Two main parts:
%   1. register 3D volume to Allen CCF
%       Output: ref2AllenTform.mat
%   2. Construct flatmap using Allen CCF
%       Output: ctxAP/MP/DepthPrctile


% Files needed (under main directory):
%   1. nrrdread
%   2. Miji installed
%   3. nissle & annotation for the reference map (only support same scale)
%   4. annotation text from Allen API matches the annotation map
%
% Note: 
%   1. Miji is needed for alignment
%   2. Use ImageJ for nrrd file, especially the annotation map, may lead to
%   errors, so need to use the nrrdread (by Jeff Mather)
%
% After register all brain seciton to 3D volume

clc

% Settings
% TBS path
addpath('C:\Users\Li Yuan\MATLAB Drive\Script65A');

% Scale factor for the current brain
scaleFactor = 0.055;

% Current mouse
directory.main = 'D:\Li Yuan';
directory.abStitch = 'D:\Li Yuan\65A_AbImage';

sysSetting = TBS.getSysSetting;
imageSetting = TBS.getImageSetting(sysSetting,[]);

mouseID = imageSetting.mouseID;
msBrainName = [mouseID,'_',num2str(scaleFactor),'.tif'];

% Thickness in pixel (scaled)
slideThickness = imageSetting.slideThickness;
resolution = imageSetting.resolution;
slideThickness = slideThickness*(scaleFactor/resolution);

% For interp1 dots
prctileQ = 0:0.05:1;    % Precentage of query
% Exclude the preentage of query if its too different form real data
maxPrctileDiff = 1; % Precentage of the side

% Annotation Info =========================================================

% Setting for reference map
refSetting = TBS.getRefSetting(directory.main);

refScale = refSetting.refScale;

% Reference map: nissle
refMap = TBS.getRefMap('nissl',refSetting);

% Annotation map
annoMap = TBS.getRefMap('anno',refSetting);

% Annotation structure
annoStruct = refSetting.annoStruct;

% Cortex (logical)
id = TBS.findAnnoID(annoStruct,'Isocortex');
ctx = ismember(annoMap,id);

sz = size(refMap);
R = imref3d(sz);

% Current brain, DAPI stack ===============================================
moving = TBS.getStack(msBrainName,[]);

sz = size(moving);
moving = reshape(moving,sz(1),sz(2),4,[]);
% Get DAPI: channel 1
moving = moving(:,:,1,:);
moving = squeeze(moving);

%% Automatic yz-angle finder ==============================================
% Use hip dots to find yz-rotation by finding the symmetrical angles
% Caution! it maynot be totally accurate; therefore, mannual adjustment
% were provided

% Get hipDots for current brain 
% % Hip dots need to be mannually selective
% % Only need to be done once then save the coordinates
% pth = fullfile(directory.main,msBrainName);
% pth = ['path=[',pth,']'];
% MIJ.run('Open...',pth);
% MIJ.run("Split Channels");
% % Then mannually select dots on hip and save the xyz as a matrix
% % in [mouseID,'_hipDot.mat']

% Get y&z rotation from hipDot --------------------------------------------
% Guilded by automatic y&z rotation finder
% 10012021: add mannually adjust

% hipDot, in micron
cd(directory.main);
load([mouseID,'_hipDot.mat']);
tform = TBS.vol2micronTform(imageSetting,scaleFactor);
hipDot = transformPointsForward(tform,hipDot);

% Get left and right groups, on x-axis
hipDot = sortrows(hipDot,1);
c = kmeans(hipDot,2,'Start',hipDot([1 end],:));

% Sort dots along the curve
sortedDot = arrayfun(@(X) hipDot(c == X,:),1:2,'UniformOutput',false);
sortedDot = cellfun(@(X) TBS.sortHipDot(X),sortedDot,'UniformOutput',false);
sortedDot = horzcat(sortedDot{:});

% Smooth the xy of dots along the line (mean filter)
filterSize = 3;
h = fspecial('average',[filterSize 1]);
fh = @(X) [imfilter(X(:,1:2),h,'replicate'),X(:,3)];
sortedDot = cellfun(@(X) fh(X),sortedDot,'Uniformoutput',false);

% Interpolation using percentage of query
tfDotCell = {vertcat(sortedDot{:,1})};
tfDotCell{2,1} = vertcat(sortedDot{:,2});
tfDotCell = cellfun(@(X) TBS.interp1Dot(X,prctileQ,maxPrctileDiff),...
    tfDotCell,'Uniformoutput',false);

% Automatically find y&z angle (can be disable) ---------------------------
% findYZangleHighResolution(dotCell,resolution)
% resolution, resolution for degree
[yAngleAuto,zAngleAuto] = TBS.findYZangleHighResolution(tfDotCell,0.01);

%% Mannual adjustment for y&z angle ========================================

% Whether to include the auto part
if ~exist('yAngleAuto') || isempty(yAngleAuto)
    yAngleAuto = 0;
end

if ~exist('zAngleAuto') || isempty(zAngleAuto)
    zAngleAuto = 0;
end

yAngleMan = -0.8;
zAngleMan = 0;

% % (Check point) check y & z angle rotation --------------------------------
% % Check whether both sides are symmetric
% % Scale tform from micron to refMap scale, otherwise its too big
% S = eye(3).*refScale;
% S(4,4) = 1;
% 
% % tform from 3d volume to LR-symetric (yz transformation)
% tform2 = TBS.roty(yAngleAuto+yAngleMan)*TBS.rotz(zAngleAuto+zAngleMan);
% tform2(4,4) = 1;
% 
% tform2 = tform.T*tform2*S;
% tform2 = affine3d(tform2);
% 
% test = imwarp(im,tform2,'OutputView',R);
% MIJ.createImage(test);

% If it passws the mannual verification -----------------------------------

% Combine y&z-rotation
tform2 = TBS.roty(yAngleAuto+yAngleMan)*TBS.rotz(zAngleAuto+zAngleMan);
tform2(4,4) = 1;
tform = tform.T*tform2;

%% Rough alignment (affine) ===============================================
% Discription: rough align current brain to reference map
% Need to enter the translation, scaling and x-rotation mannually
% Adjustment can be checked in ImageJ

% Mannual adjust 1: scaling
% Including scale to reference size
s = [1.2 1.54 0.98].*refScale;

% Mannual adjust 2: tranlation 
t =  [-38 -50 285];

% Mannual adjust 3: rotation on x-axis 
% During the adjustment, need to pay attention to hip and the point cortex
% seperates
rx = TBS.rotx(-1);
rx(4,4) = 1;

% Get transformaiton matrix
tform2 = eye(3).*s;
tform2(4,1:4) = [t,1];

tform2 = tform2*rx;

% add the previous tform
tform2 = tform*tform2;
tform2 = affine3d(tform2);

% (Check point) affine transformation -------------------------------------
test = imwarp(moving,tform2,'OutputView',R);

test = cat(3,refMap.*2,test);
MIJ.createImage(test);
MIJ.run("Stack to Hyperstack...", "order=xyzct channels=2 slices=528 frames=1 display=Composite");

%% Alignment using coronal plates =========================================
% Alignment here only take into account the coronal plate
% Output as displacement field

% Rough aligned brain
tfMoving = imwarp(moving,tform2,'OutputView',R);

% Whether redo the point-pair selection 
redoTF = false;

transformationSetting = [];
transformationSetting.type = {'polynomial2','pwl'};
transformationSetting.threshold = [{[]};{2}];

% Slide for alignment
zq = [282:5:352,317];

% Select point pair for alignment
% % getAllenD(zq,im,refMap,redoTF,transformationSetting,directory)
% allenD = TBS.getAllenD(zq,im2,refMap,redoTF,transformationSetting,directory);

% Displacment field -------------------------------------------------------

% Table for the displacement field, one slide per row
cd(directory.main)
load('allenD.mat');

sz = size(refMap);

% Get displacement field stack from the table
% interpDcell(Dcell,ax,sz)
Dx = TBS.interpDcell(allenD.D,1,sz);
Dy = TBS.interpDcell(allenD.D,2,sz);

Dz = zeros(sz);
D = cat(4,Dx,Dy,Dz);

% (Check point) Visualize registered volume -------------------------------
test = imwarp(tfMoving,D);

test = cat(3,refMap.*2,test);
MIJ.createImage(test);
MIJ.run("Stack to Hyperstack...", "order=xyzct channels=2 slices=528 frames=1 display=Composite");

return

reg2AllenTform = [];
% Scale for reference map: pixel per micron
reg2AllenTform.refScale = refScale;
% Tform for rough alignment (mannual adjusted)
reg2AllenTform.tform = tform2;
% Displacement field after adjustment
reg2AllenTform.D = D;

save(fullfile(directory.main,'reg2AllenTform.mat'),'reg2AllenTform');

%% Flatmap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(directory.main);
load('ctxAP.mat'); load('ctxDepthPrctile.mat'); load('ctxML.mat');

% Registered current brain ------------------------------------------------
% Transformation matrix and displacement field of registration to CCF
load(fullfile(directory.main,'reg2AllenTform.mat'));
tform = reg2AllenTform.tform;
D = reg2AllenTform.D;

sz = size(annoMap);
R = imref3d(sz);

% Rough aligned brain
tfMoving = imwarp(moving,tform,'OutputView',R);
tfMoving = imwarp(tfMoving,D);

%% SupFig. register flatmap ===============================================
% Delete the region overlap with brain sections

% SupFig, ML --------------------------------------------------------------
% Only include z has current data
tfMoving2 = tfMoving;
tfMoving2(ctxML > 0) = 0;
test = cat(3,ctxML,single(tfMoving2));

MIJ.createImage(test);
MIJ.run("Stack to Hyperstack...", "order=xyzct channels=2 slices=528 frames=1 display=Composite");
% setMinAndMax(-400, 8182);

% SupFig, AP --------------------------------------------------------------
tfMoving2 = tfMoving;
tfMoving2(ctxAP > 0) = 0;
test = cat(3,ctxAP,single(tfMoving2));

MIJ.createImage(test);
MIJ.run("Stack to Hyperstack...", "order=xyzct channels=2 slices=528 frames=1 display=Composite");
% setMinAndMax(6500, 9500);

% SupFig, depth prctile ---------------------------------------------------
% Only include z has current data
tfMoving2 = tfMoving;
tfMoving2(ctxDepthPrctile > 0) = 0;
test = cat(3,ctxDepthPrctile,single(tfMoving2));

MIJ.createImage(test);
MIJ.run("Stack to Hyperstack...", "order=xyzct channels=2 slices=528 frames=1 display=Composite");

%% SupFig. Allen average template/annotation flatmap =========================================

% Settings
% Including Claustrum otherwise will have errors
id = TBS.findAnnoID(annoStruct,'Isocortex');
ctx2 = ismember(annoMap,id);

% Region outline ---------------------------------------------------------
% Annotated region in 3D
refMap = TBS.getAnnoRegion(refSetting);
refMap(~ctx2) = 0;

% Change coordinates before find outline, cleaner
% Flatten the annomap to 3D ML/AP/Depth stack
[y,x,z,v] = TBS.find3(refMap);
ind = sub2ind(size(refMap),y,x,z);
regionOutline = TBS.stack2flatmapIm(ind,v,ctxML,ctxAP,ctxDepthPrctile,...
    @mode,refScale);

% Filled the empty pixel per stack along depth
SE = strel('disk',3);
regionOutline = imclose(regionOutline,SE);

% Mode project the depth
regionOutline = mode(regionOutline,3);

% Mode filter on max projection
regionOutline = modefilt(regionOutline,[7 7]);

% Get outline
SE = strel('disk',1);
regionOutline = imdilate(regionOutline,SE)~= regionOutline |...
    imerode(regionOutline,SE)~= regionOutline;

MIJ.createImage(regionOutline.*255);

% ML-AP view --------------------------------------------------------------

% Average template
avgTemp = TBS.getRefMap('avg',refSetting);
avgTemp(~ctx2) = 0;

% Average template
[y,x,z,v] = TBS.find3(avgTemp);
ind = sub2ind(size(avgTemp),y,x,z);

% stack2flatmapIm(ind,V,ctxML,ctxAP,ctxDepthPrctile,method,refScale);
test = TBS.stack2flatmapIm(ind,v,ctxML,ctxAP,ctxDepthPrctile,@max,refScale);
test = max(test,[],3);
MIJ.createImage(test);

% ML-Depth view -----------------------------------------------------------
% Mannually select the range
lim = 278:282;

[y,x,z,v] = TBS.find3(avgTemp);
ind = sub2ind(size(avgTemp),y,x,z);

test = TBS.stack2flatmapIm(ind,v,ctxML,ctxAP,ctxDepthPrctile,@max,refScale);
test = test(lim,:,:);
test = max(test,[],1);
test = squeeze(test);
MIJ.createImage(test');

%% Additional BC filter using reference map ===============================
% Register soma and axonBC coordiantes to the reference map
cd(directory.main);
load('somaBC.mat'); load('axonBC.mat'); load('codeBook.mat');

load('reg2AllenTform.mat');
tform = reg2AllenTform.tform;
D = reg2AllenTform.D;

% % Register to reference map ===============================================
% % (This part is relatively slow)
% somaBC = TBS.vol2reg(somaBC,reg2AllenTform,scaleFactor,imageSetting,sysSetting);
% axonBC = TBS.vol2reg(axonBC,reg2AllenTform,scaleFactor,imageSetting,sysSetting);
% 
% save(fullfile(directory.main,'somaBC.mat'),'somaBC');
% save(fullfile(directory.main,'axonBC.mat'),'axonBC');

% Exclude barcode outside the brain =======================================
roi = annoMap == 0;
axonBC = TBS.exclBCinROI(roi,axonBC,'xyzRef');
somaBC = TBS.exclBCinROI(roi,somaBC,'xyzRef');

% Exclude floating rolony using ROI =======================================
id = {'Hippocampal region'; 'ventricular systems'; 'fimbria'};
id = cellfun(@(X) TBS.findAnnoID(annoStruct,X),id,'Uniformoutput',false);
id = vertcat(id{:});
roi = ismember(annoMap,id);

% Connect the empty space
% 10062021: use both hip as ROI
SE = strel('sphere',3);
roi = imclose(roi,SE);
roi = TBS.imfillHoles3(roi);

% Delete floating rolonies apear on the same section as rolony in hipI
roiBC = TBS.BCinROI(roi,axonBC,'xyzRef');
axonBC = TBS.exclFloatingUseROI(axonBC,roiBC);
                 
roiBC = TBS.BCinROI(roi,somaBC,'xyzRef');
% 03212021 bug fix, output change into axonBC
axonBC = TBS.exclFloatingUseROI(axonBC,roiBC);

% Exclude soma BC within hipI
somaBC = TBS.exclBCinROI(roi,somaBC,'xyzRef');

% Exclude sporadic BC =====================================================
% CtxI, CtxC, Thal, Str, Mb
% Set minimum threshold for major projection area

% Threshold for sporadic BC per area
regionMinCount = [5 5 5 3 3];

regVoxel = TBS.getRegVoxel(annoMap,annoStruct);

% exclSporadicBC(axonBC,regVoxel,regionMinCount)
axonBC = TBS.exclSporadicBC(axonBC,regVoxel,regionMinCount);

% BCregionCountFilter =====================================================
% (same sas combineBC)
bcSetting = TBS.getBcSetting;
regionMinCount = bcSetting.regionMinCount;

% Region count filter for excluding the BC
% (at least one target region need to reach region min count)
TF = TBS.BCregionCountFilter(axonBC,regionMinCount,sysSetting);

codeBook = codeBook(TF,:);
axonBC = TBS.updateCodeID(axonBC,TF);
somaBC = TBS.updateCodeID(somaBC,TF);

% Eliminated glia cells ===================================================

% Glia radius range
gliaR = bcSetting.gliaR;

somaR = bcSetting.hasSoma.somaR;
minCount = bcSetting.hasSoma.minSomaPixelCount;
load('somaIm.mat');

% Get soma location
% hasSoma(minCount,somaR,somaBC,somaImReg,imageSetting,sysSetting)
TF80 = TBS.hasSoma(minCount,somaR,somaBC,somaIm,imageSetting,sysSetting);

% Soma location for BC with soma
somaLocation = TBS.getSomaLoc(somaBC,somaIm,imageSetting,sysSetting);

% Soma location for BC without soma
axonBCcenter = TBS.getAxonBCcenterInj(axonBC,sysSetting);
somaLocation(~TF80,:) = axonBCcenter(~TF80,:);

% Whether its a glia cell
TF = TBS.isGlia(somaLocation,axonBC,gliaR,bcSetting.minAxonCount);
TF = ~TF;

% Trim out glia BC
codeBook = codeBook(TF,:);
axonBC = TBS.updateCodeID(axonBC,TF);
somaBC = TBS.updateCodeID(somaBC,TF);

save(fullfile(directory.main,'codeBookref.mat'),'codeBook');
save(fullfile(directory.main,'somaBCref.mat'),'somaBC');
save(fullfile(directory.main,'axonBCref.mat'),'axonBC');

% Registered soma location ================================================
% Save the register soma location as variable for analysis

% Soma location
TF80 = TBS.hasSoma(minCount,somaR,somaBC,somaIm,imageSetting,sysSetting);

% Delete the location doesnt have enough pixels
somaLocation = TBS.getSomaLoc(somaBC,somaIm,imageSetting,sysSetting);
somaLocation(~TF80,:) = 0;

% Find the xyzRef info of the soma locaiton pixel
id = vertcat(somaBC.codeID{:});
xyz = vertcat(somaBC.xyz{:});
xyzRef = vertcat(somaBC.xyzRef{:});

parfor i = 1:size(somaLocation,1) % parfor
    ixyz = somaLocation(i,:);
    
    if ~any(ixyz)
        continue
    end
    
    TF = id == i & all(xyz == ixyz,2);
    
    somaLocation(i,:) = xyzRef(TF,:);
end

save(fullfile(directory.main,'somaLocationRef.mat'),'somaLocation');

%% Fig & supFig. Barcode location visualization ===========================
% Note: including both rolony and soma

cd(directory.main);
load('somaBCref.mat'); load('axonBCref.mat');

% Coordinates in 3D, um
xyz = [somaBC.xyzRef; axonBC.xyzRef];
xyz = vertcat(xyz{:});

% Count rolony number in each voxel
xyz = round(xyz);
[xyz,~,ic] = unique(xyz,'rows');
n = accumarray(ic,1);

sz = size(annoMap);
R = imref3d(sz);

%% Fig 1. Barcode location in registered volume ----------------------------
% Register moving volume
tfMoving = imwarp(moving,tform,'OutputView',R);
tfMoving = imwarp(tfMoving,D);

outline = TBS.getBrainOutline(annoMap,5);

test = TBS.xyzv2im(sz,xyz,n);

SE = strel('sphere',1);
SE = SE.Neighborhood;
SE = SE./sum(SE,'all');
tfMoving = imfilter(tfMoving,SE,'same');
test = imfilter(test,SE,'same');

test = cat(3,outline,tfMoving,test);
MIJ.createImage(test);
MIJ.run("Stack to Hyperstack...", "order=xyzct channels=3 slices=528 frames=1 display=Composite");

%% SupFig. Barcode on flatmap ===============================================

cd(directory.main);
if ~exist('ctxAP') || ~exist('ctxML') || ~exist('ctxDepthPrctile')
    load('ctxML.mat'); load('ctxAP.mat'); load('ctxDepthPrctile.mat');
end

ind = sub2ind(sz,xyz(:,2),xyz(:,1),xyz(:,3));

rolonyIm = TBS.stack2flatmapIm(ind,n,ctxML,ctxAP,ctxDepthPrctile,...
    @sum,refScale);

% If need region of interest, draw use the number along AP (didn't use the
% 2nd number)
roiAP = any(rolonyIm,2:3);
roiAP = find(roiAP);
roiAP = [roiAP(1),roiAP(end)]

% Region outline and rolony -----------------------------------------------
test = sum(rolonyIm,3);

MIJ.createImage(log(test));

% slice of z --------------------------------------------------------------
test = permute(rolonyIm,[3 2 1]);

test = test(:,:,330:339);
test = sum(test,3);
MIJ.createImage(log(test));

% Try to look at gfp, cannot see fiber tract, but see scard tissue
