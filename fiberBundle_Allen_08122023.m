% Folder for the Allen nrrd data
directory2 = 'I:\AllenTest\AllenConnectivity';

% Setting for reference map
refSetting = TBS.getRefSetting(directory.main);
refScale = refSetting.refScale;

zFlatmap = 0:2.5:95;

% Minimum intensity to be considered as signal per ML-AP pixel
minInten = 0.5;

% Get cortex pixels on flatmap
ctxTF = ctxML > 0;
ctxTF = TBS.im2flatmapAllen(ctxTF,ctxML,ctxAP,ctxDepthPrctile,refScale,zFlatmap);
ctxTF = any(ctxTF,3);

% layerFlat = TBS.im2flatmapAllen(layer,ctxML,ctxAP,ctxDepthPrctile,refScale);

% Exclude edges
SE = strel('disk',9);
ctxTF = imerode(ctxTF,SE);

%% Plot CT/PT/WT bundle location ==========================================
MIJ.run('Close All');

cd(directory2);
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

%% Visualization of PT/CT bundle location different ---------------------
% stat, cell, with 2-3 column of PT/CT intensity profiling

% Input channel number
z = size(stat{1},2);

n = cellfun(@(X) size(X,1),stat);
n = max(n);
stat2 = [];
for i = 1:numel(stat)
    istat = stat{i};
    % Normalize to the same height
    istat = imresize(istat,[n,z]);
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

%% Bulk tracing projeciton profiling ======================================
% Using flatmap. Currently for 25 um/voxel

z = 0:5:rolonyDepthEdges(end);

% (for AUD)
% Quantification: Lat-left, Med-left, Med-right, Lat-right
% ML-Boundaries
midLine = size(ctxTF,2)/2;
mlBoundary2 = [-mlBoundary 0 mlBoundary].*refScale + midLine;
mlBoundary2(end+1) = size(ctxTF,2);
mlBoundary2 = [0 mlBoundary2];

% -------------------------------------------------------------------------

f = f+1;
iFile = fileName(f);
cd(directory2);
[im,iFile] = TBS.getNrrdAllen(iFile);

im = TBS.correctInjSide(im);

% Get flapmap
flatIm = TBS.im2flatmapAllen(im,ctxML,ctxAP,ctxDepthPrctile,refScale);

% Correct for edge effect
flatIm = flatIm.*ctxTF;

% Depth: 2.5% per voxel, total 40
flatIm = flatIm(:,:,1:z(end)/2.5,:);

% Exclude pixels near saturated pixels, 8.*25 um-200 um range
SE = strel('disk',8);
injSite = any(flatIm == 1,3);
injSite = imdilate(injSite,SE);
flatIm = flatIm.*(~injSite);

% Try to eliminate blood vessles/noise using imreconstruction
TF = sum(flatIm,3);
% Thresholding on sum
TF = TF >= minInten;
TF = flatIm.*TF;
flatIm2 = imreconstruct(TF,flatIm);
% % (Check point)
% MIJ.createImage(cat(3,sum(flatIm,3),sum(flatIm2,3)));

stat = {};
for i = 1:4
    col = [mlBoundary2(i)+1,mlBoundary2(i+1)];
    col = col(1):col(2);
    iIm = flatIm2(:,col,:);
    iIm = sum(iIm,1:2);
    stat{1,i} = squeeze(iIm);
end

tableOut(iFile,:) = stat;
tableOut.Properties.VariableNames = {'LatI','MedI','MedC','LatC'}
% save('tableOut_AUD_bulk_08082023.mat','tableOut')
% save('tableOut_AUD_bulk_L45IT_08082023.mat','tableOut')

% %% ----------------------------------------------------------------------
i = 3;
stat = tableOut{:,i};

% Transparancy in Prism
% (lines beyond 1000 will be no transparancy)
n = cellfun(@sum,stat);
n = n./1000;
n = 1-min(n,1);

stat = cellfun(@(X) reshape(X,2,[]),stat,'Uniformoutput',false);
stat = cellfun(@(X) sum(X,1)',stat,'UniformOutput',false);
stat = horzcat(stat{:});
stat = stat./sum(stat,1).*100;

%% Mannually select ROI for Lat projections ===============================
% fileName = num2cell(fileName)

MIJ.run('Close All');

cd(directory2);
iFile = fileName{i+1};
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

% Projection intensity across layers --------------------------------------
MIJ.run('Close All');
MIJ.createImage(flatIm);
% Pause here: crop using roi
v = MIJ.getCurrentImage();
v = sum(v,1:2);
tableOut.roiC{iFile} = squeeze(v);

% Projection layers -------------------------------------------------------
MIJ.run('Close All');
MIJ.createImage(layerFlat);
% Pause here: crop using roi
v = MIJ.getCurrentImage(); 
% Use mode to get layer number
v(v == 0) = nan;
v = mode(v,1:2);
v = squeeze(v);
% In case there is an error in a pixel
v = medfilt2(v,[3 1]);
tableOut.layer{iFile} = v(1:zFlatmap(end)/2.5);

% save('tableOut_WT_bulk_08092023.mat','tableOut')

%% Mannual quantification of fiber bundles

MIJ.run('Close All');

cd(directory2);
iFile = fileName{i+1};
[iIm,iFile] = TBS.getNrrdAllen(iFile);
iIm = TBS.correctInjSide(iIm);

i = i+1; disp(['Current image: ',num2str(iFile)]);

% ROI selection for fiber intensity measurement
iIm = iIm(85:285,80:236,:);
MIJ.createImage(iIm); MIJ.run('Fire');

% Mannually measure upper lower fiber ratio -------------------------------
% Upper fiber: enter thalamus; lower fiber: enter internal capsul
% Mean * length
userIn = userIn(:,3).*userIn(:,7);
userIn = reshape(userIn,3,2);

userIn = mean(userIn,1);
tableOut.ulInten{iFile} = userIn

%% Stat -------------------------------------------------------------------
% File name of the region
if ~iscell(fileName)
    fileName = cellstr(num2str(fileName));
end
TF = tableOut.Properties.RowNames;
TF = cellfun(@(X) contains(TF,X),fileName,'UniformOutput',false);
TF = horzcat(TF{:});
TF = any(TF,2);

% Have correlation rho >= 0.7
TF2 = tableOut.corr;
TF2 = cellfun(@(X) X(1),TF2);
TF2 = TF2 >= 0.7;
TF = TF & TF2;

% Measured fiber intensity
TF2 = tableOut.ulInten;
TF2 = ~cellfun(@isempty,TF2);
TF = TF & TF2;

tableOut2 = tableOut(TF,:);

% Correlation -------------------------------------------------------------
% Upper/all fiber intensity
x = cellfun(@(X) X(1)/sum(X),tableOut2.ulInten);
% > 60% in contralateral side
y = cellfun(@(X) sum(X(25:end))/sum(X),tableOut2.roiC);   
[x,y].*100
