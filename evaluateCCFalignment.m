% evaluate CCF aliment
% By mannually selecting the lines/edges in the dataset and compared to CCF
% 06122023

clc

% Settings

% Directory with all files
directory.main = 'D:\Li Yuan';

cd(directory.main);

% Output variable for ROI selection
if exist('roiCCFacc.mat')
    load('roiCCFacc.mat');
else
    roiCCFacc = cell(size(tfMoving,3),2);
    % xyzc: xyz coordinates of the roi; c: the roiID on the z (of multiple
    % roi)
    % anno: annotation label
    roiCCFacc = cell2table(roiCCFacc,'VariableNames',{'xyzc','anno'});
end   

% Brain region for ROI selection ------------------------------------------
% Mannual select the regions

% Annotation map & structure
refSetting = TBS.getRefSetting(directory.main);
annoMap = TBS.getRefMap('anno',refSetting);
annoStruct = refSetting.annoStruct;
% voxel size of reference map
refVoxelSz = 1/refSetting.refScale;

% Reference region labels (appear on roiCCFacc.anno)
idLabel = {'Hippocampus'; 'Subcortical';'Striatum';...
    'fasciculus retroflexus'; 'mammillary related'; 'Isocortex'};

id = {};
% Hippocampus (w/o Entorhinal area)
id2 = TBS.findAnnoID(annoStruct,'Hippocampal formation');
TF = TBS.findAnnoID(annoStruct,'Entorhinal area');
TF = ismember(id2,TF);
id2(TF) = [];
id{1} = id2;

% Subcortical: Midbrain + Thalamus + Hypothalamus
id2 = {'Midbrain','Thalamus','Hypothalamus'};
id2 = cellfun(@(X) TBS.findAnnoID(annoStruct,X),id2,'Uniformoutput',false);
id{end+1} = vertcat(id2{:});

% Striatum
id2 = {'Striatum','Lateral amygdalar nucleus','Basolateral amygdalar nucleus'};
id2 = cellfun(@(X) TBS.findAnnoID(annoStruct,X),id2,'Uniformoutput',false);
id{end+1} = vertcat(id2{:});

% fasciculus retroflexus
id{end+1} = TBS.findAnnoID(annoStruct,'fasciculus retroflexus');

% mammillothalamic tract
id{end+1} = TBS.findAnnoID(annoStruct,'mammillary related');

% Cortex
id{end+1} = TBS.findAnnoID(annoStruct,'Isocortex');

refRegion = zeros(size(annoMap),'uint8');
for i = 1:numel(id)
    TF = ismember(annoMap,id{i});
    TF = TF.*i;
    TF = cast(TF,'like',refRegion);
    refRegion = max(refRegion,TF);    
end

% (Check point)
% MIJ.createImage(refRegion);

% Current brain stack------------------------------------------------------
% Just do it simple here, w/o calling all variables
% Transform to registrated voxels

% Channel used for evaluate alignment (i.e. vGlut2)
currentBrainCh = 3;

% Define the first and last z for roi seleciton (mannual)
zLim = 277:359;

% Brain tif file name
msBrainName = '65A_0.055.tif';

load('reg2AllenTform.mat');
% registrated brain
tfMoving = getCurrentBrain(msBrainName,currentBrainCh,reg2AllenTform);

return

%% ROI selection: mannually select the edges of reference region
% (Random selection to minimize bias of the evaluator)
% 1. Randomly select z for ROI selection
% 2. Only half of hemisphere per time, and all flip to M-L direction to
% elimiate bias
% 3. ROI was drawn in MIJ

MIJ.run("Line Width...", "line=1");

% med pixel on x-axis
medX = size(tfMoving,2);
medX = medX/2;
% In case of size difference in left and right half of hemisphere
if mod(medX,1) == 0
    medX = [medX,medX+1];
else
    medX = round(medX);
    medX = [medX,medX];
end

% Randomly select z
I = TBS.shuffleRows(zLim');

for i = I'
       
    im = tfMoving(:,:,i);
    
    % Draw ROI if there is any
    xyzc = roiCCFacc.xyzc{i};
    if isempty(xyzc)
        im2 = zeros(size(im));
    else
        im2 = TBS.xyzv2im(size(im),xyzc(:,1:2),[]);
        % push intensity for visualization
        im2 = im2.*inf;
    end
    im(:,:,2) = im2;
    
    % Randomly choose left and right hemisphere
    isLeftHemi = randi([0 1]);
    
    if isLeftHemi
        im = im(:,1:medX(1),:);
        im = fliplr(im);
    else
        im = im(:,medX(2):end,:);
    end
    
    if ~any(im,'all')
        continue
    end
       
    % ROI selection -------------------------------------------------------
    MIJ.createImage(im);
    MIJ.run("Make Composite","Display Mode = Composite");
    MIJ.run("Set Slice...","1");
    MIJ.run("Grays");
    ij.IJ.setMinAndMax(0, 150);
    MIJ.run("Set Slice...","2");
    
    answer = questdlg('Done with this image?','Log','Yes','Exit','Yes');    
    if strcmp(answer,'Exit')
        return
    end
    
    % Get the xy coordinations of the ROIs
    im = MIJ.getCurrentImage;
    MIJ.run("Close All");
    im = im(:,:,2);
    im = im > 0;
    
    if ~any(im,'all')
        continue
    end
    
    if isLeftHemi
        im = fliplr(im);
    end
    
    [y,x] = find(im);
    
    if ~isLeftHemi
        x = x + medX(2)-1;
    end
    
    xy = [x y];
    
    % Intergrate with the data --------------------------------------------
    anno = roiCCFacc.anno{i};       
    
    if ~isempty(xyzc)
        
        % Delete the xyzc & annotation not selected, in the hemisphere
        TF = ~ismember(xyzc(:,1:2),xy,'rows');
        
        if isLeftHemi
            TF = xyzc(:,1) <= medX(1) & TF;
        else
            TF = xyzc(:,1) >= medX(2) & TF;
        end
        
        xyzc = xyzc(~TF,:);
        anno = anno(~TF,:);
        
        % Delete the overlap ones
        TF = ismember(xy,xyzc(:,1:2),'rows');
        xy = xy(~TF,:);
        xy(:,3:4) = 0;
    end
    
     % Appen the extra rows
    if ~isempty(xy)       
        xyzc = [{xyzc};{xy}]; xyzc = vertcat(xyzc{:});
        
        iAnno2 = repmat({''},size(xy,1),1);
        anno = [{anno};{iAnno2}]; anno = vertcat(anno{:});
    end
    
    roiCCFacc.xyzc{i} = xyzc;
    roiCCFacc.anno{i} = anno;
end

cd(directory.main); save('roiCCFacc.mat','roiCCFacc');

%% Cluster roi lines from the same z 
% Seperate the lines

sz = size(tfMoving,1:2);

for i = 1:size(roiCCFacc,1)
    xyzc = roiCCFacc.xyzc{i};
    anno = roiCCFacc.anno{i};
    
    if isempty(xyzc)
        continue
    end
    
    xy = xyzc(:,1:2);
    
    % Segment the ROI basing on image
    im = TBS.xyzv2im(sz,xy,[]);
    im = bwlabel(im);
    
    % Find xy-coordinates and cluster number
    [y,x,c] = find(im);
    
    % Maintain the row position
    xy2 = [x y];
    [~,ia,ib] = intersect(xy2,xy,'rows');
    xy2 = xy2(ia,:);
    c = c(ia);
    xyzc = xyzc(ib,:);
    anno = anno(ib,:);
    
    % Intergrate with the data --------------------------------------------
    
    if size(xyzc,2) == 4
        c0 = xyzc(:,4);
        
        C = unique([c0 c],'rows');
        
        % Delete the current label links to more than 1 previous label
        [C,~,ic2] = unique(C(:,2));
        n = accumarray(ic2,1);
        TF = n > 1;
        C = C(TF);       
        TF = ismember(c,C);
        anno(TF) = {};        
    end
    
    z = repmat(i,size(x));
    
    roiCCFacc.xyzc{i} = [xy2 z c];
    roiCCFacc.anno{i} = anno;
end

%% Mannual assign roi to the reference region

redoTF = false;

for i = 1:size(roiCCFacc,1)
    xyzc = roiCCFacc.xyzc{i};
    anno = roiCCFacc.anno{i};
    
    if isempty(xyzc)
        continue
    end
    
    [~,~,ic] = unique(xyzc(:,4));
    iIm = refRegion(:,:,i);
    
    for j = 1:max(ic)
        TF = ic == j;
        
        jAnno = anno(TF,:);
        jAnno = unique(jAnno);
        
        if redoTF
            jAnno = {''};
        end
        
        if numel(jAnno) > 1 || isempty(jAnno{:})
            xy = xyzc(TF,1:2);
            iIm2 = TBS.xyzv2im(size(iIm),xy,[]);
            iIm2 = iIm2.*inf;
            iIm2 = cast(iIm2,'like',iIm);
            
            % Select area
            imshowpair(iIm,iIm2);
            set(gcf,'Position',[100 100 1000 800]);
            
            % Press enter once the selection is done
            pause;
            
            % Get the reference region label
            dcm = datacursormode;
            jAnno = dcm.getCursorInfo;
            jAnno = jAnno.Position;            
            jAnno = iIm(jAnno(2),jAnno(1));
            jAnno = idLabel{jAnno};
            disp(['Selected: ',jAnno]);
            anno(TF,:) = {jAnno};
           
            close;
        end        
    end
    
    roiCCFacc.anno{i} = anno;
    cd(directory.main); save('roiCCFacc.mat','roiCCFacc');
end

cd(directory.main); save('roiCCFacc.mat','roiCCFacc');

%% Calculate distance between roi line and reference regions
% If the line inside the ref region, find the shortest distance to out of
% the reference region (outter edge); if the line is outside, find the
% shortest distance to inside of the region (inner edge)

SE = strel('square',3);

xyzc = roiCCFacc.xyzc; xyzc = vertcat(xyzc{:});
anno = roiCCFacc.anno; anno = vertcat(anno{:});
[~,anno] = cellfun(@(X) ismember(X,idLabel),anno);
% Uniqle lines
[~,ia,ic] = unique(xyzc(:,3:4),'rows');

stat = [];
for i = 1:max(ic)
    TF = ic == i;
    ixyzc = xyzc(TF,:);
    ianno = anno(TF);
    
    % Selected reference region
    iIm = refRegion(:,:,ixyzc(1,3));
    iIm = iIm == ianno(1);
    % Use edge pixels to speed up computation
    % Outter edge (outside-false)
    iImF = imdilate(iIm,SE) & ~iIm;
    [y,x] = find(iImF); iImF = [x y];
    % Inner edge (inside-true)
    iImT = ~imerode(iIm,SE) & iIm;
    [y,x] = find(iImT); iImT = [x y];
                                                                                                                                                                  
    ind = sub2ind(size(iIm),ixyzc(:,2),ixyzc(:,1));
    
%     % (Check point)
%     test = false(size(iIm));
%     test(ind) = true;
%     imshowpair(iIm,test);
    
    TF = iIm(ind);
    D = zeros(size(TF));
    % Inside the ref region, to outter edge
    if any(TF)
        D(TF) = pdist2(iImF,ixyzc(TF,1:2),'euclidean','Smallest',1);
    end
    % Outside the ref region, to inner edge
    if any(~TF)
        D(~TF) = pdist2(iImT,ixyzc(~TF,1:2),'euclidean','Smallest',1);
    end
    
    x = median(ixyzc(:,1)); y = median(ixyzc(:,2));
    z = ixyzc(1,3);
    D = median(D);
    
    stat(i,:) = [x y z D];   
end

% Change to micron
stat = stat.*refVoxelSz;

% Split roi stat into regions
stat2 = {};
ia = anno(ia);
for i = 1:max(ia)
    TF = ia == i;
    stat2{i,1} = stat(TF,end);
end

% %% Scatter plot
% figure; scatter3(stat(:,1),stat(:,2),stat(:,3),20,stat(:,4),'filled');
% alpha 0.7; colorbar; 
% xlabel('L-M-L'); ylabel('D-V'); zlabel('P-A');
% g = gca; g.ZDir = 'reverse'; g.YDir = 'reverse'; g.ZTick = [7000 8000]; 
% grid off; view(15,45); daspect([1 1 1]); 
% set(gcf,'Position',[100 100 600 300]);

%% Function: getCurrentBrain
function tfMoving = getCurrentBrain(msBrainName,currentBrainCh,reg2AllenTform)
% Input:    msBrainName, str
%           currentBrainCh, num, the channel for finding lines
%           reg2AllenTform, struct, with tform and D for register to CCF
% Output:   tfMoving, mat, image stack

moving = TBS.getStack(msBrainName,[]);

sz = size(moving);
moving = reshape(moving,sz(1),sz(2),4,[]);

moving = moving(:,:,currentBrainCh,:);
moving = squeeze(moving);


% Transformation matrix and displacement field
tform = reg2AllenTform.tform;
D = reg2AllenTform.D;

% Output size as annotation map
sz = size(D);
R = imref3d(sz(1:3));

tfMoving = imwarp(moving,tform,'OutputView',R);
tfMoving = imwarp(tfMoving,D);
end
