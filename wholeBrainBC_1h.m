clc

%% Important correction for 65A
% % The batch8 slide22 is misslabel, need to mannually corrrect it and rerun
% % stitch
% cd('F:\Li Yuan\07122020_65A_Ab01_Batch8\Ab01');
% fileName = ls('EF65ASlide24*');
% fileName = cellstr(fileName);
%
% for i = 1:size(fileName,1)
%     oldName = fileName{i};
%     newName = strrep(oldName,'Slide24','Slide22');
%     movefile(oldName,newName);
% end

%% Setting

directory.main = 'D:\Li Yuan';
directory.abStitch = 'D:\Li Yuan\65A_AbImage';


sysSetting = TBS.getSysSetting;
imageSetting = TBS.getImageSetting(sysSetting,[]);


%% Align Seq to Ab image
% Aligned sequencing image from different region to the whole brain
% stitched image
% (Time-consuming step)

redoTF = true;

% channel for alignment in Ab image
abCh = 1;

% Scale factor to speed up
% Cannot be too small, >= 0.2
scaleFactor = 0.25;

% alignSeq2ab(abCh,scaleFactor,redoTF,sysSetting,directory)
seq2abTform = TBS.alignSeq2ab(abCh,scaleFactor,redoTF,sysSetting,directory);

% corrSeq2abTform ==============================================================
load(fullfile(directory.main,'seq2abTform.mat'));

% Combine tform from diff tiles
% combineSeq2abTform(tbl,tileTformName,nameElement,sysSetting,directory)
seq2abTform = TBS.combineSeq2abTform(seq2abTform,...
    'accumulateTformTable3',[1 3],sysSetting,directory);

% Mannual correction (not able to align)
rowName = 'EF65ASlide41R_VisualContra';
seq2abTform(rowName,:) = [];

corrSeq2abTform = seq2abTform;
save(fullfile(directory.main,'corrSeq2abTform.mat'),'corrSeq2abTform');

%% Append aligned seq image
% Stitch the sequencing image of the same brain section and append to the
% whole brain section

% Seq for stitched image
seq = 1:3;

% Output channel in the Ab image
chOut = 5;

load(fullfile(directory.main,'corrSeq2abTform.mat'));

TBS.stitchSeq2Ab(corrSeq2abTform,seq,chOut,imageSetting,directory);

%% Two-step semi-mannual alignment
% transform output were stored in accumulated displacement field, can be
% directly transform
% Note: the reason using accumulated displacement field is the inaccuracy
% of mannual alignment will accumulated by just aligning neigboring
% sections and may results unexpected scaling errors

% Settings 
clc;
% Whether to redo
redoTF = false;

% Scale factor during mannual alignment
scaleFactor = 0.055;

ch4Align= [2 5]; % 65A: For Ab01, use DAPI and vGlut2 for alignment

transformationSetting = [];
transformationSetting.type = {'polynomial2','pwl'};
transformationSetting.threshold = [{[]};{50*scaleFactor}];

% Reference section name (mannually pick the reference image)
refSectionName = '26L';
% Reference section number 
refSectionNum = TBS.getSectionNumber(refSectionName,imageSetting,sysSetting);

varName = {'ch4Align','transformationSetting','scaleFactor','redoTF',...
    'refSectionNum','imageSetting','sysSetting','directory'};
selfAlignmentSetting = TBS.var2struct('base',varName);

% Reference section rotation angle
selfAlignmentSetting.refSectionRotation = -4.5;

startSeq = 106;

selfAlignTform = TBS.selfAlignment(startSeq,selfAlignmentSetting);
save(fullfile(directory.main,'selfAlignTform.mat'),'selfAlignTform');

% correct Self-align selfAlignTform ======================================
% defined mannually

% 1. translation correction
% 2. delete distored sections
corrSelfAlignTform = TBS.correctSelfAlignment(refSectionNum,imageSetting,sysSetting,directory);

% Self-aligned brain volume ===============================================
% Get registrated brain volume from ab images
clc

% Settings
% Channel for output
ch = [1 2 4 5];

% Fix background intensity (mannually defined)
background = [200 110 100 0];

save(fullfile(directory.main,'corrSelfAlignTform.mat'),'corrSelfAlignTform');

TBS.stack2vol(corrSelfAlignTform,ch,background,scaleFactor,imageSetting,sysSetting,directory);

