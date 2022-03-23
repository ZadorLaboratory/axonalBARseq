% TileProcessing with pixel/local correction
% Support MATLAB 2021a
%
% Input
%   Directory (need to select folder)
% Output
%   Corrected max projection tiff image 16-bit of all subdirectory folders
%
% Last update: 03132021 by Li Yuan

% clear; 
clc;

% Regular setting =========================================================
% TBS path
addpath('C:\Users\Li Yuan\MATLAB Drive\Script65A');

saveDirectory = 'H:\abImage\06052020_65A_Ab01_Batch1';

directory.main = 'G:\Tape_3\06052020_65A_Ab01_Batch1';

sysSetting = TBS.getSysSetting;

% Get all folder names for multiple cycles processing
% Define prepend and channel name for sequencing/Ab image
seqCycles = TBS.getFolderNumber(directory.main,sysSetting.seqPrepend);

if isempty(seqCycles) == 0  % Seq image
    typePrepend = sysSetting.seqPrepend;
    chName = {'Il-G';'Il-T';'Il-A';'Il-C'};
    % Need do local and pixel correction for sequencing data
    localCorrection = true;    
    % Value for pixel background, min rank # of intensity of the pixel
    minPixelRank = 3;
    
else    % Ab image
    typePrepend = sysSetting.abPrepend;
    seqCycles = TBS.getFolderNumber(directory.main,typePrepend);
    chName = {'P-DAPI-high';'P-GFP-H';'P-RFP';'P-Cy5'};
    % whether to do local and pixel correction
    % Need to skip pixel correction if experiment is saturated
    localCorrection = false;
    % min pixel rank = 0 or empty for no pixel correction
    minPixelRank = 0;
end

nCh  = numel(chName);   % Number of channel

% Image processing settings ===============================================

% Strel Object for finding local maxima
SE = strel('diamond',2);

% Output image class
imClass = 'uint16';

imFormat = sysSetting.imFormat;

%% Image processing

for iSeq = seqCycles
    
    % For spotImage output
    spotImage = {}; rowNames = {};
    
    % Current seq folder for raw data -------------------------------------
    iSeqStr = TBS.seqstr(iSeq);
    iSeqStr = strrep(iSeqStr,'Seq',typePrepend);
    iSeqFolder = fullfile(directory.main,iSeqStr);
    
    % Get all the tile folder under current seq cycle
    cd(iSeqFolder)
    tileFolders = getFolder(['*',iSeqStr,'*']);
    
    % Get save directory --------------------------------------------------
    iSaveDirectory = fullfile(saveDirectory,iSeqStr);
    if exist(iSaveDirectory,'dir') ~=0
        % Delete all the remaining .tif file
        delete(fullfile(iSaveDirectory,['*',imFormat]))
    else
        mkdir(iSaveDirectory)
    end
    
    % Process every tile --------------------------------------------------
    parfor iFolder = 1:numel(tileFolders)  % parfor
        % current tile folder
        iTileName = tileFolders{iFolder};
        iTileFolder = fullfile(iSeqFolder,iTileName);
        
        cd(iTileFolder)
        
        % correct tile name for single tile (e.g. xxx_1)
        iTileName = correctTileName(iTileName,sysSetting);
        
        % **Special case for 65A 
        if contains(iTileName,'EF65ASlide25L_Seq13_Thalamus')
            continue
        elseif contains(iTileName,'EF65ASlide7R_Seq03_ContraContra')
            iTileName = strrep(iTileName,'EF65ASlide7R_Seq03_ContraContra',...
                'EF65ASlide7R_Seq03_Visual');
        elseif contains(iTileName,'EF65ASlide7R_Seq03_Visual')
            iTileName = strrep(iTileName,'EF65ASlide7R_Seq03_Visual',...
                'EF65ASlide7R_Seq03_VisualContra');
        end
        
        % Assemble output name
        pixelCorrectedName = [iTileName,sysSetting.pixelCorrectAppend,imFormat];
        pixelCorrectedName = fullfile(iSaveDirectory,pixelCorrectedName);
        localCorrectedName = [iTileName,sysSetting.localCorrectAppend,imFormat];
        localCorrectedName = fullfile(iSaveDirectory,localCorrectedName);
        
        iSpotImage = {};
        for iCh = 1:nCh
            % Get the stack of the current channel
            iChName = chName{iCh};
            chStack = getChStack(iChName);
            chStack = cast(chStack,imClass);
            
            % **Special case for a dataset for 65A 2nd round ------
            % Due to change to new filter for channel 3 & 4
            if contains(directory.main,'65A_2nd') && ismember(iCh,[3 4]) 
                if  contains(directory.main,'Batch1') && ismember(iSeq,[1 2])
                    chStack = chStack.*2;
                elseif contains(directory.main,'Batch2') && iSeq == 1
                    chStack = chStack.*2;
                end
            end
            
            % No pixel correction =========================================
            if isempty(minPixelRank) || minPixelRank == 0
                
                maxProjection = max(chStack,[],3);
                
                % saveChStack(im,imName,iCh)
                saveChStack(maxProjection,pixelCorrectedName,iCh);
                continue
            end
            
            % Pixel correction ============================================
            % Sort intensity
            [chStack,I] = sort(chStack,3);
            
            % Max projection
            maxProjection = chStack(:,:,end);
            pixelBackground = chStack(:,:,minPixelRank);
            chStack = [];   % Release RAM
            
            % Image output: maxProjPixelCorrected
            maxProjection = maxProjection - pixelBackground;
                        
            % Save
            saveChStack(maxProjection,pixelCorrectedName,iCh);
            
            % Skip local correction
            if localCorrection == false
                continue
            end
            
            % Local correction ============================================
            
            % The max 3 slice are not continue slice
            I = int8(I(:,:,end-2:end));
            
            % The maxIense slide is the neighboring slides or 2/3th slide
            % return to true if either slide is a neighbor of the max
            continueSignal = I(:,:,3) - I(:,:,1:2);
            continueSignal = abs(continueSignal) == 1;
            continueSignal = any(continueSignal,3);
            
            I = [];     % Release RAM
            
            % maxProjLocalCorrected ---------------------------------------
            % Use maxProjPixelCorrected instead of maxProjection, better
            % seperation near the bright soma
            
            localBackground = maxProjection;
            localBackground(continueSignal) = 0;
            
            % imreconstruct(marker,mask)
            localBackground = imreconstruct(localBackground,maxProjection);
            
            % Correct for local background
            maxProjLocalCorrected = maxProjection - localBackground;
            
            % Save
            saveChStack(maxProjLocalCorrected,localCorrectedName,iCh)
            
            % Spot output (cell) ==========================================
            
            % Local maxima needs to be local maxima and continous signal
            localMaxima = imdilate(maxProjection,SE) == maxProjection;
            localMaxima = localMaxima & continueSignal;
            
            % Use sparse logical as output
            localMaxima = sparse(localMaxima);
            iSpotImage{iCh} = localMaxima;
        end
        
        if strcmp(typePrepend,sysSetting.abPrepend)==0
            spotImage{iFolder,1} = iSpotImage;
            rowNames{iFolder,1} = iTileName;
        end
    end
    
    % Save mat output for local correctionn--------------------------------
    if localCorrection == true
        % 03242021: fix error due to special case.
        row = ~cellfun(@isempty,rowNames);
        spotImage = spotImage(row);
        rowNames = rowNames(row);
        
        spotImage = cell2table(spotImage,'RowNames',rowNames);
        save(fullfile(iSaveDirectory,'spotImage.mat'),'spotImage','-v7.3');
    end
    
    disp(['Processed images: ',typePrepend,num2str(iSeq)]);
end

%% Function:    getFolder
% Discription: Get list of folder from the input directory
function folders = getFolder(directory)
% Input: directory (string)
% Output: folders (struct)

% Get the name of all individual folders
folders = dir(directory);
folders = folders([folders.isdir]==1);
folders = {folders.name};
end

%% Function:    correctTileName
% Discription:  correct tile name for single tile
function tileName = correctTileName(tileName,sysSetting)

tileName = strsplit(tileName,sysSetting.delimiter);

% nNameElement - 2: 1 is tile number, the other is tileAppend (not in
% folder name)
if length(tileName)== sysSetting.nNameElement-2
    tileName{sysSetting.tileElement} = '01';
end

tileName = strjoin(tileName,sysSetting.delimiter);
end

%% Function:    getChStack
% Discription:  Get the stack of current channel
function currChstack = getChStack(iChName)
% Input:    iChName, str, (current channel name)
% Output:   currChstack (stack of images belongs to the channel)

currChstack = [];

% Get image names for the current channel
currChImName = dir(['*',iChName,'*']);

for iIm = 1:length(currChImName)
    currChstack(:,:,iIm) = imread(currChImName(iIm).name);
end
end

%% Function:    saveChStack
% Discription:  save stacks for the current channel
function saveChStack(im,imName,iCh)
% Input:    im, mat, image
%           imName, str, image name
%           iCh, num,channel number
% Output:   file save on disk

if iCh == 1
    imwrite(im,imName);
else
    imwrite(im,imName,'WriteMode','append');
end
end
