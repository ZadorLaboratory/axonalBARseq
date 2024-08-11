% Rotating brain, highlight data region
% Fig2A & SupVideo1
% 12082023 LY, checked 12182023

% Open the image stack need to be changed in MIJ
% im = MIJ.getCurrentImage();

% ImageJ 16 color LUT
load('LUT16Color.mat');
cmap16 = LUT16Color;

% Scale factor of the original stack
scaleFactor = 2;
nCh = 2;

% Get the limits of z-axis of sequenced area
zLim = vertcat(xyzDot{:});
zLim = zLim(:,3);
zLim = zLim.*scaleFactor;
zLim = floor(min(zLim)): ceil(max(zLim));

% Outline in the sequence region ------------------------------------------
% Change to uint8 for faster imwarp
brainTF = annoMap > 0;
brainTF = imresize3(brainTF,scaleFactor);
brainTF = uint8(brainTF);

outline = im(:,:,end,:);
outline = squeeze(outline);
if max(outline,[],'all') <= 255
  outline = uint8(outline);
end

% Exclud the data region (data-false)
outlineDataF = outline;
outlineDataF(:,:,zLim) = 0;

% Outline of the data region (data-true)
zLim2 = 1:size(outline,3);
TF2 = ismember(zLim2,zLim);
zLim2(TF2) = [];
outlineDataT = outline;
outlineDataT(:,:,zLim2) = 0;

% Reconstruction ----------------------------------------------------------
% do RGB in MATLAB
sz = size(im);
im2 = im(:,:,1,:);
im2 = squeeze(im2);

% The original single cell model has dilation (for FIJI plugin)
SE = strel('sphere',1);
im2 = imerode(im2,SE);
if max(im2,[],'all') <= 255
  im2 = uint8(im2);
end

%% Transformation matrix for 3D-rotation

% % For rotating brain
% tform = {};
% for i = 0:2:358    
%     iTform = TBS.roty(i);
%     iTform(4,4) = 1;
%     tform{end+1} = affine3d(iTform);
% end

% Fig 2A
tform = TBS.rotx(35)*TBS.roty(20)*TBS.rotz(10);
tform(4,4) = 1;
tform = affine3d(tform);
tform = {tform};

imOut = {};
for i = 1:numel(tform)

    iTform = tform{i};
    iIm = imwarp(im2,iTform,'nearest');
    iOutlineDataF = imwarp(outlineDataF,iTform,'nearest');
    iOutlineDataT = imwarp(outlineDataT,iTform,'nearest');
    iBrainTF = imwarp(brainTF,iTform,'nearest');

    nD = size(iIm,3);

    % Reconstruction ------------------------------------------------------
    % Use the closet pixel
    [D,iIm] = findZ1(iIm);

    % Change to RGB
    TF = iIm == 0;
    % Need to +1 for the reach the max-end
    n = size(cmap16,1)+1;
    iIm = ind2rgb(round(iIm./100.*n),cmap16);
    % Set 0 back to black
    TF = repmat(TF,1,1,size(iIm,3));
    iIm(TF) = 0;

    % Reduce intensity along distance
    D = (nD-D)./nD;
    iIm = double(iIm).*D;
    % MIJ.createImage(rgb);

    % Highlight unSeq region ----------------------------------------------

    D = fliplr(1:nD)./nD;
    D = reshape(D,1,1,[]);
    D = single(D);

    iOutlineDataF = single(iOutlineDataF).*D;
    iOutlineDataF = max(iOutlineDataF,[],3);

    iOutlineDataT = single(iOutlineDataT).*D;
    iOutlineDataT = max(iOutlineDataT,[],3);

    % Brain outline-2d ----------------------------------------------------
    % Get the z location of the most anterior voxel
    D = findZ1(iBrainTF);

    % Make anterior end with bigger number
    TF = D == 0;
    iBrainTF = nD-D;
    iBrainTF(TF) = 0;
    SE = strel('disk',1);
    iBrainTF = imdilate(iBrainTF,SE)-imerode(iBrainTF,SE);
    iBrainTF = min(iBrainTF,50);
    % MIJ.createImage(iBrainTF);

    % Save image ----------------------------------------------------------
    iIm = uint8(iIm);
    iOutlineDataF = uint8(iOutlineDataF);
    iOutlineDataT = uint8(iOutlineDataT);
    iBrainTF = uint8(iBrainTF);

    imOut{end+1} = cat(3,iIm,iOutlineDataF,iOutlineDataT,iBrainTF);

    disp(['Done: ',num2str(i)]);
end

return

% Settings for Fig 2A
MIJ.createImage(imOut{:});
MIJ.run("Make Composite","Display mode = Composite");
MIJ.setSlice(1); MIJ.run("Red");    ij.IJ.setMinAndMax(0,200);
MIJ.setSlice(2); MIJ.run("Green");  ij.IJ.setMinAndMax(0,200);
MIJ.setSlice(3); MIJ.run("Blue");   ij.IJ.setMinAndMax(0,200);
MIJ.setSlice(4); MIJ.run("Grays");  ij.IJ.setMinAndMax(0,2000);
MIJ.setSlice(5); MIJ.run("Grays");  ij.IJ.setMinAndMax(0,250);
MIJ.setSlice(6); MIJ.run("Grays");  ij.IJ.setMinAndMax(1,100);

%% Translocation for roating brain

if numel(imOut) == 1
    return
end

R2 = imref2d([640 1300]);
tform = eye(3);
tform(3,1) = 1300/2; tform(3,2) = 640/2;

imOut2 = {};
for j = 1:numel(imOut)
    iIm = imOut{j};

    iTform = eye(3);
    iTform(3,1) = -size(iIm,2)./2;
    iTform(3,2) = -size(iIm,1)./2;
    iTform = iTform*tform;
    iTform = affine2d(iTform);
    
    iIm = imwarp(iIm,iTform,'OutputView',R2);

    imOut2{j} = iIm;
end

imOut2 = cat(3,imOut2{:});

MIJ.createImage(imOut2);
run("Stack to Hyperstack...", "order=xyczt(default) channels=6 slices=180 frames=1 display=Composite");

%% Function:    findZ1
% Description:  fine the col and value of the first non-zero pixel along
% z-axis
function [im2col,im2v] = findZ1(im3)
% Input:    im3, 3d-mat,
% Output:   im2col, 2d-mat, the z-location of the first non-zero pixel
%           im2v, 2d-mat, the value of the first non-zero pixel

% Flat the 3D matrix into two D, one xy per row
sz = size(im3);
im3 = reshape(im3,[],sz(3));

TF = any(im3,2);

im2col = zeros(size(im3,1),1);
im2v = zeros(size(im3,1),1);

parfor i = 1:size(im3,1) % parfor
    if ~TF(i)
        continue
    end

    [~,col,v] = find(im3(i,:),1);

    im2col(i) = col;
    im2v(i) = v;
end

im2col = reshape(im2col,sz(1:2));
im2v = reshape(im2v,sz(1:2));
end

