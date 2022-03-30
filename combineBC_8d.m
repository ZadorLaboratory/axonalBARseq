clc;
clear

% Settings
directory = [];
directory.main = 'D:\Li Yuan';

sysSetting = TBS.getSysSetting;
imageSetting = TBS.getImageSetting(sysSetting,[]);
bcSetting = TBS.getBcSetting;

maxHamming = bcSetting.maxHamming;
degenerateN = bcSetting.degenerateN;

scaleFactor = 0.055; % (change way to set this)
% Soma image for soma location --------------------------------------------
% (median intensity of each seq cycle, for soma)
cd(directory.main);
if ~exist('somaIm.mat')
    somaBCVar = TBS.loadBCvar('somaBscallTable',bcSetting);
    somaIm = TBS.getSomaBCim(somaBCVar.Properties.RowNames,imageSetting.chNum,directory);
    save('somaIm.mat','somaIm','-v7.3');
else
    load('somaIm.mat');
end

% Get soma & axon BC ------------------------------------------------------
somaBCVar = TBS.loadBCvar('somaBscallTable',bcSetting);
axonBCVar = TBS.loadBCvar('bscallOutput',bcSetting);

% Delete somaBC image w/o soma --------------------------------------------
somaSectionRnge = bcSetting.somaSectionRnge;

% Soma section number
sectionName = somaBCVar.Properties.RowNames;
sectionNum = TBS.getSectionNumber(sectionName,imageSetting,sysSetting);
TF = sectionNum >= somaSectionRnge(1) & sectionNum <= somaSectionRnge(2);
somaBCVar = somaBCVar(TF,:);

% Delete rolony close to the soma -----------------------------------------
axonBCVar = TBS.delNearSomaRolony(axonBCVar,somaBCVar,bcSetting,...
    imageSetting);

%% Construct codebook =====================================================

% Use unique axon BC for codebook
bscallCh = vertcat(axonBCVar.bscallCh{:});

[codeBook,~,ic] = unique(bscallCh,'rows');

% Count filter
n = accumarray(ic,1);
row = n >= bcSetting.minAxonCount;
codeBook = codeBook(row,:);

% Exclude repeated BC with extra 0 ----------------------------------------
TF = TBS.hasOverlap(codeBook);

codeBook = codeBook(~TF,:);
disp(['After exclude overlap: ', num2str(sum(TF)),' of total ',...
    num2str(numel(TF))]);

% Merge within hamming distance -------------------------------------------
% 09152021: merge hamming distance as 1
% Keep the BC with the highest count, there can be more than 1
mergeHamming = 1;

% Merge codeBook using axonBC count
for i = 0:mergeHamming
    
    sz = 0;
    while sz ~= size(codeBook,1)
        sz = size(codeBook,1);
        
        % with 0-tolerance
        % codeLookup(BCin,codeBook,maxHamming,minBCcount,tolerate0)
        axonBCLookupTbl = TBS.codeLookup(axonBCVar.bscallCh,codeBook,i,1,true);
        axonBC = TBS.findCode(axonBCVar,axonBCLookupTbl,codeBook);
        axonBC = TBS.correctAxonBC(axonBC,bcSetting);
        
        % Exclude BC with not enough count
        axonN = vertcat(axonBC.codeID{:});
        axonN = accumarray(axonN,1,[size(codeBook,1),1]);
        row = axonN >= bcSetting.minAxonCount;
        codeBook = codeBook(row,:);
        
        if ~all(row)
            disp(['BC does not have enough axonBC count: ',...
                num2str(sum(~row)),', in total: ',num2str(numel(row))]);
            continue
        end
        
        % Merge (exclude the fix sequencing cycle)
        codeBook2 = codeBook;
        codeBook2(:,bcSetting.seqNotInMinBscallCycles) = [];
        
        TF = TBS.mergeAxonBC(codeBook2,axonN,i);
        codeBook = codeBook(TF,:);
        
        disp(['AxonBC merge into ',num2str(sum(TF)),...
            '; in total: ',num2str(numel(TF))]);
    end        
end

% %% Lookup table --------------------------------------------------------
codeBook0 = codeBook;

% Codebook with all barcodes 
% Note 03292022, 0 is treated as mismatch for non-identical matching
% codeLookup(BCin,codeBook,maxHamming,minBCcount,tolerate0)
somaBCLookupTbl = TBS.codeLookup(somaBCVar.bscallCh,codeBook0,maxHamming,1,false);
axonBCLookupTbl = TBS.codeLookup(axonBCVar.bscallCh,codeBook0,maxHamming,1,false);

axonBC = TBS.findCode(axonBCVar,axonBCLookupTbl,codeBook0);
axonBC = TBS.correctAxonBC(axonBC,bcSetting);

axonBC0 = axonBC;

cd(directory.main);
save('codeBook0.mat','codeBook0');
save('axonBC0.mat','axonBC0');
save('somaBCLookupTbl.mat','somaBCLookupTbl');
save('axonBCLookupTbl.mat','axonBCLookupTbl');

%% BC for analysis ========================================================

cd(directory.main);
load('codeBook0.mat'); load('axonBC0.mat');
load('somaBCLookupTbl.mat');

codeBook = codeBook0;
axonBC = axonBC0;

% Exclude degenerate BC ===================================================
% isDegenerateBC(BC,n,tolerate0)
TF = TBS.isDegenerateBC(codeBook,degenerateN,true);
TF = ~TF;
codeBook = codeBook(TF,:);
axonBC = TBS.updateCodeID(axonBC,TF);

disp(['Degenerate BC: ', num2str(sum(~TF)),', in total BC: ',num2str(numel(TF))]);

% Exclude BC have only certain channels ===================================
% pick up Ch3 & 4 only dots in nucleus close to on the surface
% 08312021, hamming > 2 instead 0, cannot match to ch12 or ch34 BC with
% Due to the current illumina kit has bleed through between these two pair
% of channels

minDiffCh = 3;

% Ch 1 & 2
TF = ismember(codeBook,[0 1 2]);
TF = sum(~TF,2) >= minDiffCh;
codeBook = codeBook(TF,:);
axonBC = TBS.updateCodeID(axonBC,TF);

disp(['Exclude BC with only Ch1&2: ', num2str(sum(~TF)),...
    ', in total BC: ', num2str(numel(TF))]);

% Ch3 & 4
TF = ismember(codeBook,[0 3 4]);
TF = sum(~TF,2) >= minDiffCh;
codeBook = codeBook(TF,:);
axonBC = TBS.updateCodeID(axonBC,TF);

disp(['Exclude BC with only Ch3&4: ', num2str(sum(~TF)),...
    ', in total BC: ', num2str(numel(TF))]);

% Exclude BC with low complete matches ====================================
% 0 treated as a match
% countMismatch(tbl,codeBook)
stat = TBS.countMismatch(axonBC,codeBook);

% Complete match has the max axonBC count
[~,I] = max(stat,[],1);
TF = I == 1;

codeBook = codeBook(TF,:);
axonBC = TBS.updateCodeID(axonBC,TF);

disp(['Exclude BC with low complete match: ', num2str(sum(~TF)),...
    ', in total BC: ', num2str(numel(TF))]);

% Count filter ============================================================
% Total axon/soma BC count

somaBC = TBS.findCode(somaBCVar,somaBCLookupTbl,codeBook);

% Filter out barcode with too many axon/somaBC counts
TF = TBS.BCcountFilter(codeBook,somaBC,axonBC,bcSetting);

codeBook = codeBook(TF,:);
axonBC = TBS.updateCodeID(axonBC,TF);

% BCregionCountFilter =====================================================

regionMinCount = bcSetting.regionMinCount;

% Region count filter for excluding the BC
% (at least one target region need to reach region min count)
TF = TBS.BCregionCountFilter(axonBC,regionMinCount,sysSetting);

codeBook = codeBook(TF,:);
axonBC = TBS.updateCodeID(axonBC,TF);

% Exclude BC from 2nd infection ===========================================

% BC from 2nd infection identified mannually
BC2nd = TBS.get2ndInfectionBC;

D = TBS.hammingDist2(BC2nd,codeBook,true);

% Min hamming distance to 2nd infection
% Note: has to do maxHamming * 2, otherwise still pick up 2nd infection
% 09152021, checked again just maxHamming, same issue
row = D(1:(maxHamming*2+1),:);
TF = ~any(row,1);

codeBook = codeBook(TF,:);
axonBC = TBS.updateCodeID(axonBC,TF);

disp(['Exclude 2nd infection BC: ', num2str(sum(~TF)),...
    ' ,in total BC: ',num2str(numel(TF))]);

% Delete same rolony image multiple times =================================
% ver4e, Exclude same rolony from more than one image

% regX/Y, pixel location registered to antibody section
axonBC = TBS.registerXY2Vol(axonBC,scaleFactor,sysSetting,directory);
% xyz, registered coordinates in micron
axonBC = TBS.vol2micron(axonBC,imageSetting,sysSetting);

% Test threshold for overlapping ------------------------------------------
% Threshold in micron
threshold = 0:5:100;
stat = [];
parfor i = 1:numel(threshold) % parfor
    TF = TBS.isSameRolonyDiffIm(axonBC,threshold(i));
    TF = vertcat(TF{:});
    stat(i) = sum(TF)./numel(TF);
end

stat = [threshold; stat.*100]';

overlapStat = stat;

% Get rolony to exclude ---------------------------------------------------

% Use 95% of max as threshold
xq = max(stat(:,2)).*0.95;

% Mannually enter the number accroding ng the stat
threshold = 25; 

disp(['Threshold (um) for excluding same rolony in diff image: ',...
    num2str(threshold)]);

axonBC = TBS.delSameRolonyDiffIm(axonBC,threshold);

% Floating rolony ========================================================

% Soma radius
somaR = bcSetting.hasSoma.somaR;
minCount = bcSetting.hasSoma.minSomaPixelCount;

% Soma pixel registration to find soma slide
somaBC = TBS.findCode(somaBCVar,somaBCLookupTbl,codeBook);
somaBC = TBS.registerXY2Vol(somaBC,scaleFactor,sysSetting,directory);
somaBC = TBS.vol2micron(somaBC,imageSetting,sysSetting);

% BC with soma
TF50 = TBS.hasSoma(minCount,somaR,somaBC,somaIm,imageSetting,sysSetting);

% Soma slide 
somaSection = TBS.getSomaSection(somaBC,somaIm,imageSetting,sysSetting);

% floating rolony slide ---------------------------------------------------
% Region to be included for floating rolony computation
regionName = bcSetting.floatingSlide.regionName;

% Minimum rolony count for floating rolony
floatThreshold = bcSetting.floatingSlid.threhold;

% Note: 0.5 indicates between two sections
floatSection = TBS.findFloatSlide(axonBC,floatThreshold,regionName,...
                bcSetting,imageSetting,sysSetting);

% (Report) Positive & accuracy test for BC with soma ----------------------
I = TF50 & ~isnan(floatSection) & floatSection ~= 0;

err = somaSection(I)-floatSection(I);
err = abs(err) <= 1;

disp('Positive & accuracy rate: ')
disp([sum(I)/sum(TF50), sum(err)/sum(I)]);

% (Report) Positive rate in BC without soma -------------------------------

TF = ~TF50 & ~isnan(floatSection) & floatSection ~=0;

disp('BC without soma: positive and total: ');
disp([sum(TF), sum(~TF50)])

% (Report) Test accuracy in soma-negative cells ---------------------------
% Coordinate: xyz in volume space (self-aligned)

testTF = false; % whether start test

% Cells for testing floating sections
TF = ~TF50 & ~isnan(floatSection) & floatSection ~=0;
I = find(TF);

if testTF
    xyz = vertcat(axonBC.xyz{:});
    id = vertcat(axonBC.codeID{:});
    
    % Subsample to decrease to speed up
    xyz0 = TBS.shuffleRows(xyz);
    xyz0 = xyz0(1:100000,:);
    
    % blue, BC rolony; red, rolony from floating section
    figure;
    for i = I'
        TF = id == i;
        % All rolony distribution in black
        hold off; scatter3(xyz0(:,1),xyz0(:,2),xyz0(:,3),2,'k','filled',...
            'MarkerFaceAlpha',0.1);
        % Rolony from current cell
        hold on; scatter3(xyz(TF,1),xyz(TF,2),xyz(TF,3),10,'b','filled');
        
        title(['BC ',num2str(i)]);
        xlabel('x'); ylabel('y'); zlabel('z'); grid off;
        g = gca; g.YDir = 'reverse'; g.ZDir = 'reverse';
        daspect([1 1 1]); view(0,0); set(gcf,'Position',[100 100 1000 300]);
        
        % Rolony on the floating section 
        iz = floatSection(i);
        % z-location in micron
        iz = (iz-1)*imageSetting.slideThickness+1
        
        TF2 = xyz(TF,3) == iz;
        TF(TF) = TF2;
        scatter3(xyz(TF,1),xyz(TF,2),xyz(TF,3),10,'r','filled');
    end
end

% Delete floating rolonies ------------------------------------------------
% Sections to be deleted
delSection = somaSection;
delSection(~TF50) = floatSection(~TF50);

axonBC = TBS.delFloatingRolony(delSection,axonBC,regionName,...
    imageSetting,sysSetting);

% Region count filter for excluding the BC --------------------------------
TF = TBS.BCregionCountFilter(axonBC,regionMinCount,sysSetting);

codeBook = codeBook(TF,:);
axonBC = TBS.updateCodeID(axonBC,TF);
somaBC = TBS.updateCodeID(somaBC,TF);

% %% Eliminated glia cells ==================================================

% Glia radius range
gliaR = bcSetting.gliaR;

% Get soma location% hasSoma(minCount,somaR,somaBC,somaImReg,imageSetting,sysSetting)
TF50 = TBS.hasSoma(50,somaR,somaBC,somaIm,imageSetting,sysSetting);

% Soma location for BC with soma, xyz in micron
somaLocation = TBS.getSomaLoc(somaBC,somaIm,imageSetting,sysSetting);

% Soma location for BC without soma, median rolony location, xyz in micron
axonBCcenter = TBS.getAxonBCcenterInj(axonBC,sysSetting);
somaLocation(~TF50,:) = axonBCcenter(~TF50,:);

% Whether its a glia cell
TF = TBS.isGlia(somaLocation,axonBC,gliaR,bcSetting.minAxonCount);
TF = ~TF;

% Trim out glia BC
codeBook = codeBook(TF,:);
axonBC = TBS.updateCodeID(axonBC,TF);
somaBC = TBS.updateCodeID(somaBC,TF);

%  Codebook for analysis ==================================================

cd(directory.main);
save('codeBook.mat','codeBook');
save('somaBC.mat','somaBC'); save('axonBC.mat','axonBC');

return

%% (Report) Difference between soma and neigboring slides =================

% Find barcode with 50-500 pixel count soma 
minCount = [50:50:550];

TF0 = {};
for i = minCount
    % hasSoma(minCount,somaR,somaBC,somaImReg,imageSetting,sysSetting)
    TF0{end+1} = TBS.hasSoma(i,somaR,somaBC,somaIm,imageSetting,sysSetting);
end

%% (Report) Estimate degenerate BC ========================================

nIteration = 1000;  % iteration number
nBC = 10000;        % BC number per iteration
nSeq = 17;          % BC length
n = 3:15; % continous N same base sequence

p = TBS.estimateDegenerateBC(nIteration,nBC,nSeq,n);

% Stat output for prism: median and MAD -----------------------------------
p = p.*100;
stat = [];
stat(1,:) = mean(p,1);
stat(2,:) = std(p,1,1);

%% SupFig. Model % BC within hamming distance =============================
% Percentage of barcode have different hamming distance

nIteration = 100;
nBC = 10000;
nSeq = 15;

D = [];
for i = 1:nIteration
    %  randBC(nBC,nSeq,n,minDiffCh)
    BC = TBS.randBC(nBC,nSeq,[],[]);
    
    iD = TBS.hammingDist2(BC,[],false);
    
    % Precentage of barcode with a partner at each hamming distance
    iD = sum(iD > 0,2)/nBC;
    
    D(:,i) = iD;
    disp(['Iteration: ',num2str(i)]);
end

D = D.*100;
stat = [mean(D,2),std(D,0,2)];

%% Fig 1. NT proportion in BC

nCh = imageSetting.chNum;

n = [];
for i = 1:nCh
    iCh = codeBook == i;
    n(i,:) = sum(iCh,1);
end

n = n./sum(n,1);

% Precentage of BC
stat = n'.*100;

% (Report) BC with fixed nt for vamp2 -------------------------------------

% 9th: ch3; 10th: ch2
TF = codeBook(:,9) == 3 & codeBook(:,10) == 2;

stat = sum(TF);
disp(['BC with correct fix NT: ',num2str(stat),newline,...
    'Precentage: ', num2str(stat/size(TF,1))]);

%% Fig 1. Hamming distance

nIteration = 100;
nBC = 10000;
nSeq = 15;

D = [];
for i = 1:nIteration
    %  randBC(nBC,nSeq,n,maxHamming)
    BC = TBS.randBC(nBC,nSeq,[],[]);
    iD = TBS.hammingDist2(BC,[],false);  
        
    D(:,i) = sum(iD,2);
    disp(['Iteration: ',num2str(i)]);
end
stat = D./sum(D,1).*100;
stat = [mean(stat,2),std(stat,0,2)];

codeBook2 = codeBook;
codeBook2(:,bcSetting.seqNotInMinBscallCycles)=[];
Dmat = TBS.hammingDist2(codeBook2,[],true);
Dmat = sum(Dmat,2);
Dmat = Dmat./sum(Dmat).*100;

%% (Check) Codebook with 0
% Get row number for codes with 0, for mannual validation

row = sum(codeBook == 0,2);
disp('Rows in codebook with 0: ')
row = find(row)

%% (Report) Soma and axon BC count

% cd(directory.main);
% load('axonBC.mat'); load('somaBC.mat');

axonN = vertcat(axonBC.codeID{:});
axonN = accumarray(axonN,1);

somaN = vertcat(somaBC.codeID{:});
somaN = accumarray(somaN,1);

% Soma & axon barcode count
% For correlation & histogram
stat = [somaN, axonN];

disp(['Total axonBC count:',num2str(sum(axonN))]);
disp(['Median axonBC count:',num2str(median(axonN))]);

%% SupFig. hamming distance of result to the codeBook
% Calculate the cumulateive precentage of axon/somaBC with different
% hamming distance

tbl = axonBC;
tbl = somaBC;

stat = TBS.countMismatch(tbl,codeBook);
% Precentage
stat = stat./sum(stat).*100;
stat = stat(1:3,:);
stat = cumsum(stat,1);

if isequal(tbl,axonBC)
    stat = [mean(stat,2),std(stat,0,2)];
else
    % Only include BC with soma
    stat = [mean(stat(:,TF50),2),std(stat(:,TF50),0,2)];
end

%% (Check) soma bscalling percentage ======================================
% Discirption: check how many cells were included in the data

cd(directory.main);
% load('codeBook.mat'); load('somaBC.mat'); 
threshold = 300;

% Mannually identify image to check
i = 31;

% imName = somaBC.Properties.RowNames{i};

imName = 'EF65ASlide23L_Injection.tif';

% Median intensity of bscalling channel across cycles
im = somaIm.im{imName};

% Bscalled cells
im2 = TBS.imBCid(somaBC(imName,:));

% Same size of both image
sz = size(im2);
im = im(1:sz(1),1:sz(2));

im = im >= threshold;

im = cat(3,im,im2);

MIJ.createImage(uint8(im).*255);

return

% Open directly in ImageJ to preserve the hyperstack
cd(directory.main);
idir = dir(['**/soma/',imName]);
idir = ['path=[',fullfile(idir.folder,idir.name),']'];
MIJ.run('Open...',idir);

%% (Check) axon bscalling percentage ======================================
% Discirption: check how many cells were included in the data

cd(directory.main);
% load('codeBook.mat'); load('axonBC.mat');

% Mannually identify image to check
i = 47;

% Sequence for stitched image
iSeq = 14;

% Image name
imName = 'EF65ASlide24L_Injection.tif';
imName = 'EF65ASlide36R_Contra.tif';
% imName = axonBC.Properties.RowNames{i}

% Bscalled dots at selective seq
im = TBS.imBCid(axonBC(imName,:),iSeq);
im = im > 0;

SE = strel('disk',1); im = imdilate(im,SE);
im = uint8(im).*255;

% Open directly in ImageJ to preserve the hyperstack
slide = (iSeq-1).*4+(1:4);
cd(directory.main);
idir = dir(['**/stitched/',imName]);
im2 = TBS.getStack(fullfile(idir.folder,idir.name),slide);
im2 = max(im2,[],3);
im2 = uint8(im2);

% Same size and cat
sz = size(im);
im2 = im2(1:sz(1),1:sz(2)); 
im = cat(3,im,im2);

MIJ.createImage(im);

%% (Check) soma location ==================================================

cd(directory.main);
% load('codeBook.mat'); load('somaBC.mat'); 

% Settings

% Soma location, 08112021
% somaLoc = TBS.getSomaLoc(somaBCreg,somaIm,imageSetting,sysSetting);
somaLocation = TBS.getSomaLoc(somaBC,somaIm,imageSetting,sysSetting);

% SomaBC & axon count -----------------------------------------------------
n = vertcat(somaBC.codeID{:});
n = accumarray(n,1);

n2 = vertcat(axonBC.codeID{:});
n2 = accumarray(n2,1);

% Test, sort using soma
[~,I2] = sort(n,'descend'); 

%% Check soma location & bscalling result for indivisual BC ---------------
clc

% Code ID to check (defined mannually)
i = 5352 %I2(4)

% Find the register xyz belongs to the BC
codeID = vertcat(somaBC.codeID{:});
xyz = vertcat(somaBC.xyz{:});
iCenter = somaLocation(i,:);

TF = ismember(xyz,iCenter,'rows') & codeID == i;

% xy coordiantes, image name & bscalling result
sz = cellfun(@(X) size(X,1),somaBC.codeID);
TF = mat2cell(TF,sz,1);
xy = cellfun(@(X1,X2,Y) [X1(Y,:),X2(Y,:)],somaBC.x,somaBC.y,...
    TF,'UniformOutput',false);
xy = vertcat(xy{:});

TF2 = cellfun(@any,TF);
imName = somaBC.Properties.RowNames{TF2};

bc = codeBook(i,:);

disp(['Image name: ',imName,newline,newline,'BC: ', num2str(bc)]);
% disp(['Soma & axon BC count: ',num2str([n(i),n2(i)])]);

% Image name and counts
fh = @(X,Y) disp([X.Properties.RowNames(Y > 0), num2cell(Y(Y>0))]);
TF = cellfun(@(X) sum(X == i),somaBC.codeID); fh(somaBC,TF);
TF = cellfun(@(X) sum(X == i),axonBC.codeID); fh(axonBC,TF);

% Local corrected image ---------------------------------------------------

% Open directly in ImageJ to preserve the hyperstack
cd(directory.main);
idir = dir(['**/soma/',imName]);
idir = ['path=[',fullfile(idir.folder,idir.name),']'];
MIJ.run('Open...',idir);
MIJ.run('Specify...', ['width=15 height=15 x=',num2str(xy(1)),' y=',num2str(xy(2)),' slice=4 centered']);

% Logical soma image ------------------------------------------------------

im = TBS.imBCid(somaBC(imName,:));
im = im == i;
MIJ.createImage(uint8(im).*255);
MIJ.run('Specify...', ['width=15 height=15 x=',num2str(xy(1)),' y=',num2str(xy(2)),' slice=4 centered']);

% %% Background subtracted soma Image ----------------------------------------
% imName = somaBCreg.Properties.RowNames{TF2};
% idir = dir(['**/stitched/',imName]);
% stack = TBS.getStack(fullfile(idir.folder,idir.name),[]);
% sz = size(stack);
% stack = reshape(stack,sz(1),sz(2),[],imageSetting.chNum);
% 
% stack = TBS.somaIntenNormalization(stack);
% stack = TBS.somaBkgrdSubtraction(stack);
% stack = TBS.somaIntenNormalization(stack);
% 
% stack = reshape(stack,sz);
% MIJ.createImage(stack);

%% (Check axon) -----------------------------------------------------------

imName2 = 'EF65ASlide22R_Contra.tif';
im = TBS.imBCid(axonBC(imName2,:),2);
im = im == i;
im = imdilate(im,strel('disk',4));
MIJ.createImage(uint8(im).*255);
MIJ.run('Make Binary');
MIJ.run('Create Selection');
MIJ.run('Make Inverse');

cd(directory.main);
idir = dir(['**/stitched/',imName2]);
idir = ['path=[',fullfile(idir.folder,idir.name),']'];
MIJ.run('Open...',idir);

%% (Check) axon bscalling Accuracy ========================================
% Discirption: check how many cells were included in the data

cd(directory.main);
load('codeBook.mat'); load('axonBC.mat');

% Image to check
i = 10;

% Image name
imName = axonBC.Properties.RowNames{i}

tbl = axonBC(imName,:)

% Open directly in ImageJ to preserve the hyperstack
cd(directory.main);
idir = dir(['**/stitched/',imName]);
idir = ['path=[',fullfile(idir.folder,idir.name),']'];
MIJ.run('Open...',idir);

%% Select specific xy or row to check

col = 2;

j = 12;

% Region center for checking the dots (i.e. injection center)
% Need mannually identify the xy
inputXY = [1650 1500]; 

% Find the dot coordinates around the inputXY
xy = [tbl.x{:}(:,1),tbl.y{:}(:,1)];
[~,I] = pdist2(xy,inputXY,'euclidean','Smallest',20);

% Dots to check
j = I(j); disp(['Dot: ',num2str(j)]);

xy = [tbl.x{:}(:,col),tbl.y{:}(:,col)];
xy = xy(j,:);
bc = tbl.codeID{:}(j); bc = codeBook(bc,:); disp(bc);

MIJ.run('Specify...', ['width=15 height=15 x=',num2str(xy(1)),' y=',num2str(xy(2)),...
    ' slice=',num2str(col),' centered']);

%% Check rolony without soma

% Codebook test1 for statistics
axonBCLookupTbl = TBS.codeLookup(axonBCVar.bscallCh,codeBookTest1,0,1,false);
axonBC = TBS.findCode(axonBCVar,axonBCLookupTbl,codeBookTest1);

axonBC2 = {};
parfor i = 1:size(axonBC,1) % parfor
    axonBC2{i,1} = TBS.trimRepeatID(axonBC(i,:),bcSetting);
end
axonBC = vertcat(axonBC2{:});
