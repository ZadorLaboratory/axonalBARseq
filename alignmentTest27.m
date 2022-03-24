% Basecall using the finished tileReg
% Get bscall stack of each tile for stitching
% Note: some variable is too big, need to be saved as '-v7.3'
% Note: currently only support <=255 nt long barcode

clear
clc

% 0.1 System settings
% TBS path
addpath('C:\Users\Li Yuan\MATLAB Drive\Script65A');

directory= [];
directory.main = 'D:\Li Yuan\11282019_65A_7th\Batch2\';
directory.stitched = fullfile(directory.main,'stitched');
directory.temporary = fullfile(directory.main,'temporary');
% Folder for stitched chCorrected file for basecalling
directory.soma = fullfile(directory.main,'soma');

sysSetting = TBS.getSysSetting;
% Append for max projection input
sysSetting.tileAppend = sysSetting.localCorrectAppend;

spotMatName = 'spotImage'; % This need to be mannually changed

% 0.2 Imaging Setting (default settings)
imageSetting = TBS.getImageSetting(sysSetting,directory);

% 0.3 bleedThrough Estimation Setting =====================================
% Fixed channel: channel bleedthrough to other channel
% Depended channel: channel get signal from the fixed channel

% Noise for spot picking in the fixed channel
bleedThroughEstSetting.noise = 50;
% Get tile from these locaiton for calculating bleedthrough
bleedThroughEstSetting.location = {'Contra','Thalamus'};
% For tile to be included for calculation, minimum dot number
bleedThroughEstSetting.numPerTilePrctile = 90;
% The bleed through correction will be used in these cycles
bleedThroughEstSetting.correctionForCycle = imageSetting.seqCycles;
bleedThroughEstSetting.fitType = 'poly1';

% Table for bleedthrough estimation
% with info of fixed/depended channel
% seq, sequencing cycle use for identifying spots
bleedThroughEstimation = table();
bleedThroughEstimation.fixedCh{'Ch1Ch2'} = 2;
bleedThroughEstimation.dependedCh{'Ch1Ch2'} = 1;
bleedThroughEstimation.seq{'Ch1Ch2'} = 10;
bleedThroughEstimation.cfun{'Ch1Ch2'} = {};

bleedThroughEstimation.fixedCh{'Ch3Ch4'} = 3;
bleedThroughEstimation.dependedCh{'Ch3Ch4'} = 4;
bleedThroughEstimation.seq{'Ch3Ch4'} = 9;
bleedThroughEstimation.cfun{'Ch3Ch4'} = {};

% 0.4 Channel profile setting ================================================
% For profiling each channel intensity from all cycles
% then normalized them to chScore for each channel per cycle
% So the intensity can be compared acrossed channel and sequencing cycles

chProfileSetting = [];
chProfileSetting.refSeq = 1:3;
% Noise of rolony picking (e.g. lower threshold of ch4) in refSeq
% (If input is 1 noise, it will apply to all 4 channels)
chProfileSetting.noise = 25;
% Minimum of rolony per tile for all channel
chProfileSetting.numOfDot = 20;
% Choose whether going to cap max intensity for chProfiling
% enable expect max intensity
chProfileSetting.enableExpectMaxInten = true;
% For all locations, 'All'
chProfileSetting.location = {'Contra','Thalamus'};
% Fit types for two variables
chProfileSetting.var1fitType = 'poly1';
chProfileSetting.var2fitType = 'poly1';
% Sequencing channel use for profiling, a subset of correctionForCycle
% Exclude seq09 and 10 for VAMP2 due to fixed base
% Exclude seq01 due to high intensity
chProfileSetting.profilingSeq = imageSetting.seqCycles(~ismember(imageSetting.seqCycles,[bleedThroughEstimation.seq{:}]));
chProfileSetting.profilingSeq(1) = [];
% Special case for 06142019_65A due to the change of emission fiter
if contains(directory.main,'06142019_65A') == 1
    chProfileSetting.noise = [25 25 15 15];
end
% Apply the corrections to these cycles
chProfileSetting.correctionForCycle = imageSetting.seqCycles;

% 0.5 FindMaxima Setting ==================================================

% Basecalling noise setting for findMaxima 
bscallNoise = [];
bscallNoise.initial.mean = 1;
bscallNoise.initial.std = -0.9;
bscallNoise.change = zeros(1,4);

% Threshold for soma basecalling
somaBscallThreshold.initial.mean = 1;
somaBscallThreshold.initial.std = 1;
somaBscallThreshold.change = zeros(1,4);

chRatioFilter = [];
chRatioFilter.dot = 0.95;
chRatioFilter.soma = 0.8; % didnt use, 26e

% Scale factor for finding subpixel location
scaleFactor = 3;

% 0.6 Bscall Correction & Analysis Setting ================================
dotMatchingSetting = [];
% Max dist (pixel) for dot matching between pair of sequcing cycles
dotMatchingSetting.maxDist = 5;
% Max interval for empty sequencing cycle
% Missing sequence per image is not included
% (for fast processing during dot match; if it miss multiple cycles it
% likely wont be matched again, so discard it)
dotMatchingSetting.maxSeqInterval = 3;
% minimum bscall cycle for excluding tiles
dotMatchingSetting.minBClen = 9;
dotMatchingSetting.seqNotInMinBscallCycles = [9,10];

%% Get local maximum intensity
% Find the potential dots from the local maximum
% w/ highest intensity in the channel

dotTable = TBS.getDotTable(spotMatName,imageSetting,sysSetting,directory);
save(fullfile(directory.main,'dotTable.mat'),'dotTable','-v7.3');

% %% Bleedthrough & channel corrections

cd(directory.main);
load('dotTable.mat');

% Bscall correciton
% bleed through estimation ================================================
bleedThroughEstimation = TBS.getBleedThroughcfun(bleedThroughEstimation,bleedThroughEstSetting,dotTable);
save(strcat(directory.main,'bleedThroughEstimation.mat'),'bleedThroughEstimation');

% Bleedthrough correciton =================================================
dotTable = TBS.bleedThroughCorrection(dotTable,bleedThroughEstimation,imageSetting);

% Channel profiling and chScore cfun =======================================
% getChScoreCfun(bscallTable,figOutput,chProfileSetting,sysSetting,imageSetting)
% Special case for 06142019_65A due to the change of emission fiter
[chEstimateCfun,chProfileSetting] = TBS.getChScoreCfun(dotTable,true,chProfileSetting,imageSetting);
savefig(fullfile(directory.main,'chProfile.fig'));

if contains(directory.main,'06142019_65A') == 1
    %  specialCaseRound1(dotTable,chProfileSetting,imageSetting,directory)
    chEstimateCfun = specialCaseRound1(dotTable,chProfileSetting,imageSetting,directory);
end
save(strcat(directory.main,'chProfileSetting.mat'),'chProfileSetting');
save(strcat(directory.main,'chEstimateCfun.mat'),'chEstimateCfun');

% Convert channel intensity to chScore =====================================
dotTable = TBS.chScoreConvertion(dotTable,chEstimateCfun,chProfileSetting);

% Bscall intensity ratio correction =======================================
% fh for getting rows pass the ratio
keepRow = @(X) TBS.bscallRatioCorrection(X.chIntensity,chRatioFilter.dot);
dotTable.bscall = cellfun(@(X) X(keepRow(X),:),dotTable.bscall,'Uniformoutput',false);

correctedDotTable = dotTable;
save(strcat(directory.main,'correctedDotTable.mat'),'correctedDotTable','-v7.3');
clearvars dotTable

% %% Image channel correction

load(fullfile(directory.main,'bleedThroughEstimation.mat'));
load(fullfile(directory.main,'chEstimateCfun.mat'));
load(fullfile(directory.main,'chProfileSetting.mat'));

% Assign universal mean and std for intensity calculation
% for the same mouse (mannually selected)
chProfileSetting = TBS.getChScore2Inten(chProfileSetting);

% Modified sysSetting.tileAppend to chCorrected
sysSetting = TBS.imageChCorrectionMain(bleedThroughEstimation,chEstimateCfun,chProfileSetting,imageSetting,sysSetting,directory);

% %% bscallTable
% Rolony for data

% double check the tile append to chCorrectAppend
sysSetting.tileAppend = sysSetting.chCorrectAppend;

load(fullfile(directory.main,'correctedDotTable.mat'))
% strrepRowName(tbl,oldStr,newStr)
correctedDotTable = TBS.strrepRowName(correctedDotTable,sysSetting.localCorrectAppend,sysSetting.tileAppend);

correctedDotTable = TBS.getDotBeyondNoise(correctedDotTable,bscallNoise,imageSetting,directory);
correctedDotTable = TBS.ind2subInTable(correctedDotTable,imageSetting.tileSize);

% %% Subpixel location of the dot center by scaling the image ================
% Range for finding local maxima
SE = strel('diamond',2);
bscallTable = TBS.getScaleBscallTable(correctedDotTable,scaleFactor,SE,imageSetting,sysSetting,directory);

% Delete extra dots, one dot per pixel ====================================
% load(fullfile(directory.main,'bscallTable.mat'));
bscallTable.bscall = cellfun(@(X) TBS.uniqueDotPerPixel(X),...
    bscallTable.bscall,'UniformOutput',false);

save(strcat(directory.main,'bscallTable.mat'),'bscallTable','-v7.3');

clearvars correctedDotTable scaledBscallTable

%% Get fuse image, rough alignment
% Settings
% Scale down for fast processing
scaleRatio = 0.5;

redoTF = true;

% Tile name to check alignment
checkAlignment = false;
if checkAlignment
    imName1 = 'EF65ASlide15L_Seq12_Injection.tif';
    imName2 = 'EF65ASlide15L_Seq16_Injection.tif';
end

% get fuse image & tform for rough alightment -----------------------------
% Note: no alignemnt to get the stitch (just use the grid imaging setting),
% to avoid posibility of error of missing or misalinment; also much faster
fuseStitchTable = fuseImage(scaleRatio,imageSetting,sysSetting,directory);

% Align stitch image ------------------------------------------------------
maxInten = 100;
nTime = 3; % Align # times
fuseAlignTform = TBS.alignFuseImage(redoTF,nTime,maxInten,...
    imageSetting,sysSetting,directory);

% (Optional) check fuse image alignment quality ---------------------------
if checkAlignment    
    load('fuseAlignTform.mat');
    fuseAlignTform = TBS.getAccumulateTform(fuseAlignTform,'movingName','fixName','tform');
    
    % getAlignedImPair(imName1,imName2,tformTable,directory)
    [tfIm1,tfIm2] = TBS.getAlignedImPair(imName1,imName2,fuseAlignTform,...
        directory.temporary);
    
    figure;
    if isempty(maxInten)
        imshowpair(tfIm1,tfIm2);
    else
        imshowpair(min(tfIm1,maxInten),min(tfIm2,maxInten));
    end
end

%% ======================================================
scaleRatio = 0.5;

sysSetting.tileAppend = sysSetting.chCorrectAppend;

profile on

% Get tform across sequencing cycles --------------------------------------
cd(directory.temporary);
load('fuseStitchTable.mat');
load('fuseAlignTform.mat');
fuseAlignTform = TBS.getAccumulateTform(fuseAlignTform,'movingName','fixName','tform');

% Scale back --------------------------------------------------------------
if scaleRatio ~= 1
    % S*T*inv(S)
    scaleTform = TBS.getScaleTform(scaleRatio,2);
    fuseAlignTform.tform = cellfun(@(X) scaleTform*X*inv(scaleTform),...
        fuseAlignTform.tform,'UniformOutput',false);
end

% Add alignment tform to the stitch tform ---------------------------------
% addImTform2TileTform(imTform,tileTform);
fuseStitchTable = TBS.addImTform2TileTform(fuseAlignTform,fuseStitchTable,sysSetting);

% %% Get tformTable =========================================================
% from aligning individual tiles 

redoTF = true;

load(fullfile(directory.main,'bscallTable.mat'));

% Delete wrong images & tiles (identified mannually)
% (This session is here because wrong images were identified using fuse
% image, so just delete the file and rerun this part would be fine)
bscallTable = TBS.delWrongTile(bscallTable);

% getTformTable(spot4Align,fuseStitchTable,redoTF,nTime,imageSetting,sysSetting)
tformTable = TBS.getTformTable(bscallTable,fuseStitchTable,true,imageSetting,sysSetting,directory);

save(fullfile(directory.main,'tformTable.mat'),'tformTable');

% %% accumulateTformTable1: tile alignment -----------------------------------
load(fullfile(directory.main,'tformTable.mat'));

accumulateTformTable1 = TBS.getAccumulateTform(tformTable,'movingName','fixName','tform');

save(strcat(directory.main,'accumulateTformTable1.mat'),'accumulateTformTable1');

% %% Stitching tiles ======================================================
% basing on dot positions

load(fullfile(directory.main,'bscallTable.mat'))
load(fullfile(directory.main,'accumulateTformTable1.mat'))

% Transform dot coordinates for tile stitching 
tfBscallTable = TBS.tranformXYinTable(bscallTable,accumulateTformTable1);

% (Check point) Test accumulate tform ---------------------------------------------------
% im1 = 'EF65ASlide19L_Seq15_Contra_05_chCorrected.tif';
% moving = tfBscallTable.bscall{im1}{:,1:2};
% im2 =  'EF65ASlide19L_Seq09_Contra_05_chCorrected.tif';
% fix = tfBscallTable.bscall{im2}{:,1:2};
% TBS.checkDotAlign(eye(3),moving,fix);

% %% Stitching tiles using transformed dot coordinates
stitchTformTable = TBS.stitchDot(tfBscallTable,accumulateTformTable1,...
    imageSetting,sysSetting);

save(strcat(directory.main,'stitchTformTable.mat'),'stitchTformTable');

% accumulateTformTable2: tile alingment + stitching -----------------------

% combTformInTable(tblA,tblB,fh)
fh = @(A,B) A*B;
accumulateTformTable2 = TBS.combTformInTable(accumulateTformTable1,stitchTformTable,fh);

save(strcat(directory.main,'accumulateTformTable2.mat'),'accumulateTformTable2');

%% Image stitching ======================================================

if ~exist(directory.stitched)
    mkdir(directory.stitched);
end

load(fullfile(directory.main,'accumulateTformTable2.mat'));

paddingSize = imageSetting.tileSize;
overlapTF = true;
% fuseSeqImage(tformTable,redoTF,paddingSize,overlapTF,imageSetting,sysSetting,directory)
% (To redo subset of stitching, set redoTF = false, and deleted the bad
% stitched image from the stitched folder)
paddingTformTable = TBS.fuseSeqImage(accumulateTformTable2,redoTF,...
    paddingSize,overlapTF,imageSetting,sysSetting,directory);

save(strcat(directory.main,'paddingTformTable.mat'),'paddingTformTable');

% accumulateTformTable3: tile alingment + stitching + padding -------------
% Coordinates should be the same as image

fh = @(A,B) A*B;
accumulateTformTable3 = TBS.combTformInTable(accumulateTformTable2,paddingTformTable,fh);

save(strcat(directory.main,'accumulateTformTable3.mat'),'accumulateTformTable3');

%% SomaBscall ===========================================================
% 08292021, ver26e, non-overlap stitching for soma bscalling

% Stitch image for somas
directory.stitched = directory.soma;

load(fullfile(directory.main,'accumulateTformTable3.mat'));

% Only include tiles for injection
imName = accumulateTformTable3.Properties.RowNames;
row = cellfun(@(X) TBS.doSomaBscall(X,sysSetting),imName);
accumulateTformTable3 = accumulateTformTable3(row,:);

% No extra padding
paddingSize = [0 0];
% Delete overlap region
overlapTF = false;
TBS.fuseSeqImage(accumulateTformTable3,true,paddingSize,overlapTF,...
    imageSetting,sysSetting,directory);

% Soma bscalling-----------------------------------------------------------
% (using image from directory.stitched)
somaBscallTable = TBS.somaBscall(somaBscallThreshold,...
    dotMatchingSetting,imageSetting,sysSetting,directory);

save(strcat(directory.main,'somaBscallTable.mat'),'somaBscallTable','-v7.3');

% Change the directory back
directory.stitched = fullfile(directory.main,'stitched');

% %% (Option) check soma correction result ----------------------------------
% tbl = somaBscallTable(5,:);
% im = TBS.getBscallBW(tbl);
% MIJ.createImage(uint8(im));

%% Bscall table ===========================================================
% 08262021,ver 26d,change delete rolony on soma pixel to delNearSomaRolony

bcSetting = TBS.getBcSetting;
somaSectionRnge = bcSetting.somaSectionRnge;    

load(fullfile(directory.main,'somaBscallTable.mat'))
load(fullfile(directory.main,'bscallTable.mat'))
load(fullfile(directory.main,'accumulateTformTable3.mat'))

% Transform coordinates to location in stitch image -----------------------
tfBscallTable = TBS.tranformXYinTable(bscallTable,accumulateTformTable3);

% Combine tiles & seq info
tfBscallTable = TBS.squeezeTile(tfBscallTable,sysSetting);

% Exclude dots near soma --------------------------------------------------
% (otherwise its too slow)

% Get soma sections within the range 
sectionName = somaBscallTable.Properties.RowNames;
sectionNum = TBS.getSectionNumber(sectionName,imageSetting,sysSetting);
TF = sectionNum >= somaSectionRnge(1) & sectionNum <= somaSectionRnge(2);
somaBscallTable = somaBscallTable(TF,:);

if ~isempty(somaBscallTable)
    tfBscallTable = TBS.delNearSomaRolony(tfBscallTable,somaBscallTable,...
        bcSetting,imageSetting);
end

% Match dots --------------------------------------------------------------
bscallOutput = TBS.getBscallTable(tfBscallTable,dotMatchingSetting,imageSetting);

save(strcat(directory.main,'bscallOutput.mat'),'bscallOutput');
beep
return

%% (Optional) check bscalling result

% Image number in the table
iIm = 1;
tbl = somaBscallTable(iIm,:);

% Get image output of bscall
im = TBS.getBscallBW(tbl);
% im = reshape(im,size(im,1),size(im,2),[]);

% Get indivdual bscalling result in the table
I = 100; % row number
TBS.checkBscallResult(tbl,I);

%% Function:    specialCaseRound1
function chEstimateCfun = specialCaseRound1(dotTable,chProfileSetting,imageSetting,directory)
    chProfileSetting.profilingSeq = chProfileSetting.profilingSeq(1:end-2);
    chProfileSetting.correctionForCycle = chProfileSetting.correctionForCycle(1:end-2);
    
    % getChScoreCfun(bscallTable,figOutput,chProfileSetting,imageSetting)
    [chEstimateCfun,~] = TBS.getChScoreCfun(dotTable,true,chProfileSetting,imageSetting);
    savefig(fullfile(directory.main,'chProfile1.fig'));
    
    if contains(directory.main,'Batch1') == 1
        defaultMean = [81.365, 49.315, 55.755,30.51];
    elseif contains(directory.main,'Batch2') == 1
        defaultMean = [45.74 43.265 52.65 43.655];        
    end
    
    addChEstimateCfun = chEstimateCfun(end,:);
    addChEstimateCfun{1}.mean = defaultMean(1); 
    addChEstimateCfun{2}.mean = defaultMean(2); 
    addChEstimateCfun{3}.mean = defaultMean(3); 
    addChEstimateCfun{4}.mean = defaultMean(4); 
        
    chEstimateCfun(16:17,:) = repmat(addChEstimateCfun,2,1);
end

%% Function:    fuseImage
% Discription:  get fuse image for rough alignemnt across seq
function fuseStitchTable = fuseImage(scaleRatio,imageSetting,sysSetting,directory)

% Append of image for stitching
imAppend = sysSetting.chCorrectAppend;

% Output image bit
imageBits = 'uint8';
imageBits = @(X) cast(X,imageBits);

if ~exist(directory.temporary)
    mkdir(directory.temporary);
end

scaleTform = TBS.getScaleTform(scaleRatio,2);

% Output name elements
% 1-sectionElement; 2-seqElement; 3-regionElement
sysSetting.nameElement = 1:3;

tilePos = imageSetting.tilePos;

% Estimate translation for all tile position
tilePosTranslation = cellfun(@(X) TBS.estimateTranslation(X,imageSetting),...
    tilePos,'UniformOutput',false);

% Get image & tile name ---------------------------------------------------
cd(directory.main);
fileName = ls(['**\*',imAppend,'*']);
fileName = cellstr(fileName);

% Get image name
imName = cellfun(@(X) TBS.nameFun(X,sysSetting.nameElement,sysSetting),...
    fileName,'Uniformoutput',false);
imName = cellfun(@(X) [X,sysSetting.delimiter],imName,...
    'UniformOutput',false);
imName = unique(imName);

fuseStitchTable = {};
parfor iIm = 1:size(imName,1) % parfor
    iImName = imName{iIm};
    
    % Get tile name for the current image
    row = contains(fileName,iImName);
    tileName = fileName(row);
    tileNum = TBS.getTileNum(tileName,sysSetting);
    
    % Get estimate translation tform for each tile
    nTile = TBS.estimateTotalNumOfTile(tileNum,tilePos);
    tform = tilePosTranslation{nTile};
    tform = tform(tileNum);
    
    % Add to the stitch table
    fuseStitchTable{iIm,1} = [tileName,tform];
    
    % Stitching -----------------------------------------------------------
    % Get output size (No padding is fine)
    iTilePos = tilePos{nTile};
    sz = TBS.estimateStitchImSize(iTilePos,imageSetting,[0 0]);
    sz = round(sz.*scaleRatio);
    R = imref2d(sz);
    
    % Get directory
    cd(directory.main);
    iDirectory = dir(['**\',iImName,'*',imAppend,'*']);
    iDirectory = iDirectory(1).folder;
    cd(iDirectory);
    
    stack = imageBits(zeros(sz));
    for iTile = 1:size(tileName,1)
        iTileName = tileName{iTile};
        iStack = TBS.getStack(iTileName,[]);
        % Change image bit
        iStack = imageBits(iStack);
        % Max projection
        iStack = max(iStack,[],3);
                
        if scaleRatio ~= 1
            iStack = imresize(iStack,scaleRatio);
        end
        
        % Transform
        iTform = tform{iTile};
        iTform = inv(scaleTform)*iTform*scaleTform;
        iTform = affine2d(iTform);
        
        iStack = imwarp(iStack,iTform,'OutputView',R);
        
        stack = max(stack,iStack);
        
        disp(['Fuse image: ',iTileName]);
    end
    
    % Save image
    cd(directory.temporary);
    iFileName = TBS.imFormat(iImName,sysSetting);
    TBS.saveStack(stack,iFileName);
end

fuseStitchTable = vertcat(fuseStitchTable{:});
fuseStitchTable = table(fuseStitchTable(:,2),'VariableNames',{'tform'},...
    'RowNames',fuseStitchTable(:,1));

cd(directory.temporary);
save('fuseStitchTable.mat','fuseStitchTable');
end
