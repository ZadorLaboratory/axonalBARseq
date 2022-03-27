classdef TBS % Terminal BARseq
    % Add Ab and registration functions
    
    methods (Static)    % General stuff ===================================
        %% Function:    extractStructVar
        % Discription:  extract all variables from structure to workspace
        function extractStructVar(S,ws)
            % Input:        S, structure
            %               ws, 'base' or 'caller'
            % Output:       void, variable in work space
            
            fields = fieldnames(S);
            
            for i = 1:length(fields)
                assignin(ws,fields{i},S.(fields{i}));
            end
        end
        
        %% Function:    var2struct
        % Discription:  pack variable into structure, with the same name
        function S = var2struct(ws,varName)
            % Input:        varName, cell
            %               ws, workspace, caller/base
            % Output:       S, struct
            
            if ~iscell(varName)
                error('Function var2struct: input varName is not cell array')
            end
            
            for i = 1:length(varName)
                S.(varName{i}) = evalin(ws, varName{i});
            end
        end
        
        %% Function:    getStack
        % Discription:  Get stacks of an image
        function stack = getStack(fileName,stackNumber)
            % Input:        fileName, string
            %               chNum, num/mat, specific stack number
            % Output:       stack, M*N*O matrix
            
            % This will have error for large tif file
            if nargin == 1 || isempty(stackNumber)
                info = imfinfo(fileName);
                stackNumber = 1:size(info,1);
            end
            
            for i = 1:length(stackNumber)
                stack(:,:,i) = imread(fileName,stackNumber(i));
            end
        end
        
        %% Funciton:    saveStack
        % Discription:  save image stack
        function saveStack(stack,fileName,appendTF)
            % Input:        stack, 3d mat
            %               fileName, str, w/ or w/o path
            %               appendTF, logical, whether only save as append
            % Output:       void. image will be saved on the disk
            % 07202021: combine with saveStackAppend
            
            if numel(size(stack)) > 3
                stack = reshape(stack,size(stack,1),size(stack,2),[]);
                warning('Stack will be saved as 3-D stack.')
            end
            
            if nargin == 2
                appendTF = false;
            end
            
            for i = 1:size(stack,3)
                % Not append only
                if i == 1 && ~appendTF
                    imwrite(stack(:,:,i),fileName);
                else
                    imwrite(stack(:,:,i),fileName,'WriteMode','append')
                end
            end
        end
        
        %% Function:    cell2stack
        % Discription:  Change cell array (w/ mat) into a stacked mat
        function stack = cell2stack(cellArray,dimension4stackCell)
            % Updated 11082019, change cell2mat to switch and loop, much faster
            % Input:        cellArray, cell array w/ a mat per cell
            %               dimension4stackCell, dimsion for stacking the cells
            %               (the matrix in the cells)
            % Output:       stack, mat
            
            stack = cellArray{1};
            
            for i = 2:length(cellArray)
                stack = cat(dimension4stackCell,stack,cellArray{i});
            end
        end
                
        %% Function:    stack2cell
        % Discription:  change stack (3d-mat) to cell array
        function cellArray = stack2cell(stack)
            % Input:        stack,3d mat
            % Output:       cellArray, cell, 1*1*size(stack,3)
            
            nStack = size(stack,3);
            cellArray = mat2cell(stack,size(stack,1),size(stack,2),ones(nStack,1));
        end
        
        %% Function:    stack2spCell
        % Discription:  change stack 3d mat into 3d cell with sparse matrix
        % Similar as stack to cell, but each cell contains sparse mat
        function stackCell = stack2spCell(stack)
            % Input:        stack, 3d-mat
            % Output:       stackCell, cell with 2d sparse mat
            
            % Right now sparse only support double and logical
            if ~islogical(stack)
                stack = double(stack);
            end
            
            stackCell = TBS.stack2cell(stack);
            % Change into sparse matrix
            stackCell = cellfun(@(X) sparse(X),stackCell,'Uniformoutput',false);
        end
        
        %% Function:    spCell2stack
        % Discription: change cell with sparse matrix per cell to mat
        function stack = spCell2stack(sparseCell)
            % Input:    sparseCell,cell array with sparse mat
            % Output:   stack, mat
            
            sparseCell = cellfun(@full,sparseCell,'Uniformoutput',false);
            stack = cell2mat(sparseCell);
        end
        
        %% Function:    hasThisRow
        % Discription:  whether the table has this row
        function TF = hasThisRow(tbl,rowName)
            % Input:        tbl, table
            %               rowName, str or cell
            % Output:       TF, logical output same size as rowName
            
            if istable(tbl) == 0
                error('Input tbl is not a table');
            end
            
            TF = ismember(rowName,tbl.Properties.RowNames);
        end
        
        %% Function:    strrepRowName
        % Discription: strrep old str to new string in row names
        function tbl = strrepRowName(tbl,oldStr,newStr)
            % Input:    tbl, table
            %           oldStr, newStr, old and new string for the row name
            
            if ~istable(tbl)
                error('Input tbl is not a table');
            end
            
            rowNames = tbl.Properties.RowNames;
            rowNames = cellfun(@(X) strrep(X,oldStr,newStr),rowNames,...
                'Uniformoutput',false);
            tbl.Properties.RowNames = rowNames;
        end
         
        %% Function:    getScaleTform
        % Discription:  get transformation matrix (3x3, 4x4) for scaling
        function scaleTform = getScaleTform(scaleFactor,dim)
            % Input:    scaleFactor, num
            %           dim, dimension for the transformation
            % scaleTform: dim x dim mat
            
            scaleTform = eye(dim+1);
            scaleTform(1,1) = scaleFactor;
            scaleTform(2,2) = scaleFactor;
        end
                
        %% Function:    find3
        % Discription:  find nonzeros element in 3d space
        function [row, col, z, v] = find3(matIn)
            % Input:    matIn, mat, 3d
            
            sz = size(matIn);
            ind = find(matIn);
            
            v = matIn(ind);
            v = double(v);
            [row,col,z] = ind2sub(sz,ind);
        end
        
        %% Function:    im2xyzv
        % Discription:  convert image pixel to xyzv (double)
        function xyzv = im2xyzv(im)
            % Input:    im, mat, image stack
            
            [y,x,z,v] = TBS.find3(im);
            xyzv = [x,y,z,v];
        end
        
        %% Function:   isOutOfRngeDot
        % Discription:  check whether dots are out of pixel area
        function TF = isOutOfRngeDot(sz,dotMat)
            % Input:    sz, mat, image size
            %           dotMat, mat, coordinates, same size as sz
            % Output:   TF, logical vector
            
            % Shift x-y to row-col
            sz(1:2) = fliplr(sz(1:2));
            
            % 02/14/2022: add round
            dotMat = round(dotMat);
            TF = dotMat < 1 | dotMat > sz;
            TF = any(TF,2);
        end
        
        %% Function:    xyzv2im
        % Discription:  convert the sub to image
        function im = xyzv2im(sz,xyz,v)
            % Note, if v is empty, will be logical
            % Update: include out of range and round
            % Input:    sz, size of the image
            %           xyz, mat, xy/xyz of the location
            %           v, vector, value
            
            xyz = round(xyz);
            
            row = TBS.isOutOfRngeDot(sz,xyz);
            if any(row)
                xyz = xyz(~row,:);
                
                if ~isempty(v)
                    v = v(~row,:);
                end
                
                warning(['Function xyzv2im: eliminated out of range dots ',...
                    num2str(sum(row))]);
            end
            
            % Allow both 2 or 3 dimensiton
            if size(xyz,2) == 2
                ind = sub2ind(sz,xyz(:,2),xyz(:,1));
            elseif size(xyz,2) >= 3
                ind = sub2ind(sz,xyz(:,2),xyz(:,1),xyz(:,3));
            end
            
            if isempty(v) || islogical(v)
                im = false(sz);
                im(ind) = true;
            else
                % Add value using the same class type
                classType = class(v);
                classFh = str2func(classType);
                
                im = zeros(sz);
                im = classFh(im);
                im(ind) = v;
            end
        end
        
        %% Function:    xyz2v
        % Discription:  get the value in im
        function v = xyz2v(xyz,im)
            % Input:    xyz, mat, coordinates
            %           im, 3-D or 2-D mat
            % Output:   v, vector,
            
            xyz = round(xyz);
            
            if size(xyz,2) == 2
                ind = sub2ind(size(im),xyz(:,2),xyz(:,1));
            elseif size(xyz,2) == 3
                ind = sub2ind(size(im),xyz(:,2),xyz(:,1),xyz(:,3));
            end
            
            v = im(ind);
        end             
        
        %% Funciton:    getImValidLim
        % Discription:  get the valid pixel limits along an axis
        function axLim = getImValidLim(im,axNum)
            % Input:    im, image stack
            %           axNum, num, of axis
            % Output:   axLim, mat with min and max pixel location
            
            % Find the other axes
            plate = 1:ndims(im);
            plate = plate(plate ~= axNum);
            
            axLim = any(im,plate);
            axLim = [find(axLim,1,'first'),find(axLim,1,'last')];
        end
        
        %% Function:    imtrim
        % Discription:  trim out the zero element from image in specific ax
        function im = imtrim(im,ax)
            % Input:    im, mat, image stack (max 3-d)
            %           ax, mat, axes to be trimmed
            % Output:   im, trimed image stack
            
            for a = ax
                % Get axis limit
                axLim = TBS.getImValidLim(im,ax(a));
                
                % Trim
                switch a
                    case 1
                        im = im(axLim(1):axLim(2),:,:);
                    case 2
                        im = im(:,axLim(1):axLim(2),:);
                    case 3
                        im = im(:,:,axLim(1):axLim(2));
                    otherwise
                        error('Currently only support 3 axis')
                end                
            end
        end
        
        %% Function:    rotz/roty/rotx
        % Discription:  Rotation matrix for rotations around z/y/x-axis
        % Input:    ang, num, angle
        function R = rotz(ang)
            R = eye(3);
            R(1,1) = cosd(ang); R(1,2) = sind(ang);
            R(2,1) = -sind(ang); R(2,2) = cosd(ang);
        end
        
        function R = roty(ang)
            R = eye(3);
            R(1,1) = cosd(ang); R(1,3) = -sind(ang);
            R(3,1) = sind(ang); R(3,3) = cosd(ang);
        end
        
        function R = rotx(ang)
            R = eye(3);
            R(2,2) = cosd(ang); R(2,3) = sind(ang);
            R(3,2) = -sind(ang); R(3,3) = cosd(ang);
        end
                
        %% Function:    imfillHoles3
        % Discription:  fill holes in stacks of images
        function im = imfillHoles3(im)
            % Input/output: im, image stacks
            for z = 1:size(im,3)
                im(:,:,z) = imfill(im(:,:,z),'holes');
            end
        end
        
        %% Function:    outlineSpCell2Full
        % Discription:  get full stack from compressed data
        function stack = outlineSpCell2Full(outlineSpCell)
            % Input:    cell or mat
            % Output:   mat
            
            if iscell(outlineSpCell)
                outlineSpCell = TBS.spCell2stack(outlineSpCell);
            end
            
            stack = TBS.imfillHoles3(outlineSpCell);
        end
        
        %% Function:    RB2paddingSize
        % Discription:  get the padding size from imref2d/3d
        function paddingSize = RB2paddingSize(RB)
            % Input:    RB, imref2d/3d object
            % Output:   paddingSize, 1 x 2/3 vector
            
            if isa(RB,'imref2d')
                paddingSize = -[RB.XWorldLimits(1),RB.YWorldLimits(1)]+0.5;
            elseif isa(RB,'imref3d')
                paddingSize = -[RB.XWorldLimits(1),RB.YWorldLimits(1),...
                    RB.ZWorldLimits(1)]+0.5;
            else
                error('Currently onlly suport imref2d/3d.')
            end
        end
        
        %% Function:    RB2imSize
        % Discription:  get the image size from imref2d/3d
        function imSize = RB2imSize(RB)
            % Input:    RB, imref2d/3d object
            % Output:   image, 1 x 2/3 vector
            
            % row & col
            imSize = [RB.YIntrinsicLimits; RB.XIntrinsicLimits];
            
            if isa(RB,'imref3d')
                imSize = [imSize; RB.ZIntrinsicLimits];
            end
            
            imSize = imSize(:,2)-imSize(:,1);
            imSize = reshape(imSize,1,[]);
        end
        
        %% Function:    repmat2size
        % Discription:  repmat input1 to the same size as input2
        function X = repmat2size(X,Y,dim)
            % Input:    X, mat/cell
            %           Y, cell
            %           dim, dimension for repmat
            % Output:   X, cell
            
            n = cellfun(@(X) size(X,dim),Y);
            
            if isempty(dim) || dim == 1
                X = arrayfun(@(X2,Y2) repmat(X2,Y2,1),X,n,'UniformOutput',false);
            elseif dim == 2
                X = arrayfun(@(X2,Y2) repmat(X2,1,Y2),X,n,'UniformOutput',false);
            end
        end
        
        %% Function:    getParaInObj
        % Discription:  Get the median of parameter in objects from cell
        % array
        function computedPara = getParaInObj(objCell,parameterStr,funStr)
            % Input:        objCell, cell of object
            %               parameterStr, str, string of parameter name
            %               funStr, str, string of function name
            % Output:       computedPara, num, computed parameters
            
            % Get the parapeters
            computedPara = cellfun(@(X) X.(parameterStr),objCell,'UniformOutput',false);
            computedPara = cell2mat(computedPara);
            
            % Apply the function to the parameters
            fh = str2func(funStr);
            if contains(funStr,{'min','max'}) == 1
                computedPara = fh(computedPara,[],'all');
            else
                computedPara = fh(computedPara,'all');
            end
        end
        
        %% Function:    fitWithoutOutlier
        % Discription:  get fit without outlier using the default outlier setting
        function [fitOutlierExclude,outlier] = fitWithoutOutlier(x,y,fitType)
            % Input:    x, mat
            %           y, mat/cell
            %           fitType, str, type for fitting
            % Output:   fitOutlierExclude,fitObject
            %           outlier, mat/cell(if y is cell, same format as y)
            
            sizeY = [];
            % if y is cell array, get the mat size inside each cell
            if iscell(y) % Change required: This need to put outside the funciton
                sizeY = cellfun(@(X) size(X,1),y);
                y = vertcat(y{:});
            end
            
            if iscell(x)
                x = vertcat(x{:});
            end
            
            fitOriginal = fit(x,y,fitType);
            
            % Find outlier, using the diff between modeled and real value
            fdata = feval(fitOriginal,x);
            outlier = fdata - y;
            outlier = isoutlier(outlier);
            
            % Fit without outlier
            fitOutlierExclude = fit(x,y,fitType,'Exclude',outlier);
            
            % Format the outlier into same format as y
            if isempty(sizeY) == 0
                % Format the outliers into cell array as y
                outlier = mat2cell(outlier,sizeY,1);
            end
        end
        
        %% Function:    shuffleRows
        % Discription:  random shuffle rows
        function [A,I] = shuffleRows(A)
            % Input & output: A, cell/mat
            %       I, row number for the shuffling
            
            n = size(A,1);
            I = randperm(n);
            A = A(I,:);
        end
        
        %% Function:    repmat2cell
        % Discription:  repmat content to the same row number as cell array
        function cellOut = repmat2cell(matIn,cellIn)
            % Input:    matIn, cell/mat, content for repmat
            %           cellIn, cell, prvides row format for each cell
            % Output:   cellOut, cell, with repeated rows of matIn per cell
            
            if size(matIn,1)~= size(cellIn,1)
                error('The length of both input are not the same.');
            elseif size(cellIn,2) ~= 1
                error('There are more than one column for 2nd input.')
            end
            
            if ~iscell(matIn)
                matIn = mat2cell(matIn,ones(size(matIn,1),1));
            end
            
            % Allow empty dotCell
            I = ~cellfun(@isempty,cellIn);
            cellOut = cell(size(cellIn));
            
            cellOut(I) = cellfun(@(X,Y) repmat(X,size(Y,1),1),...
                matIn(I),cellIn(I),'UniformOutput',false);
        end
        
        %% Function:    nonzeroAvgFilter
        % Discription:  average filter of non-zero elements
        function im = nonzeroAvgFilter(im,h)
            % 03162022, delete exclude0
            % Input & output:  im, mat, image
            %           h, mat, filter
                        
            % Sum
            im2 = convn(im,h,'same');
            
            % number of voxel/pixel
            fh = @(X) cast(X,class(im));
            n = fh(im ~= 0);
            n = convn(n,h,'same');
            
            im = im2./n;
            
            % Delete pixel with nan
            TF = isnan(im);
            im(TF) = 0;
        end
                
        %% Function:    nonzeroMedFilt
        function im = nonzeroMedFilt(im,h)
            % Note: parfor is used
            % Input:    h, logical matrix, filter voxel
                        
            sz = size(im);
            
            % Loop through all voxel within the filter
            I = find(h);
            
            disp(['Function nonzeroMedFilt in progress, filter size ',...
                num2str(numel(I))]);
                        
            im2 = [];
            parfor i = 1:numel(I) % parfor
                h2 = false(size(h));
                h2(i) = true;
                
                iIm = convn(im,h2,'same');
                im2(:,i) = reshape(iIm,[],1);
                disp(['Voxel: ',num2str(i)]);
            end
            
            % Set 0 to nan, get median without nan
            im2(im2 == 0) = nan;
            im2 = median(im2,2,'omitnan');
            
            im = reshape(im2,sz);            
        end
                
        %% Function:    getMedEdges
        % Discription:  get middle of edges
        function edges = getMedEdges(edges)
            % Input & output: edges, vector
            
            edges = edges(1:end-1)+edges(2:end);
            edges = edges./2;
        end
        
        %% Function:    axLabelSettings
        % Discription:  Set front and size for all axes
        function axLabelSettings(frontName,frontSize)
            % Input:    frontName, str, front name
            %           fontSize, num, front size
            
            ax = {'XAxis'; 'YAxis'; 'ZAxis'};
            
            g = gca;
            for i = 1:numel(ax)
                g.(ax{i}).Label.FontName = frontName;
                g.(ax{i}).Label.FontSize = frontSize;
            end
        end
        
        %% Function:    accumarrayMean
        % Discription:  accumarray-mean function, faster
        function B = accumarrayMean(ind,data)
            % Input:    ind, vector, Output indices
            %           data, vector
            % Output:   B, vector, mean of data for each indices
            
            data = accumarray(ind,data);
            n = accumarray(ind,1);
            
            B = data./n;
        end
        
    end
    
    methods (Static)    % Default settings ================================
        %% Function:    getSystemSetting
        function sysSetting = getSysSetting
            % Output:  sysSetting, struct
            
            sysSetting = [];
            
            % Image setting general
            sysSetting.slidePrepend = 'Slide';
            sysSetting.delimiter = '_';
            sysSetting.imFormat = '.tif';
            
            sysSetting.nNameElement = 5;
            sysSetting.sectionElement = 1;
            sysSetting.seqElement = 2;
            sysSetting.regionElement = 3;
            sysSetting.tileElement = 4;
            sysSetting.nameElement = [sysSetting.sectionElement,...
                sysSetting.regionElement];
            
            % Image prepend: sequencing /immuno data
            sysSetting.seqPrepend = 'Seq';
            sysSetting.abPrepend = 'Ab';
            
            % Image append: correction type
            sysSetting.pixelCorrectAppend = '_pixelCorrected';
            sysSetting.localCorrectAppend = '_localCorrected';
            sysSetting.chCorrectAppend = '_chCorrected';
            sysSetting.temporaryAppend = '_temporary';
            
            % Special setting for experiments
            sysSetting.injectionAppend = 'Injection';
            sysSetting.somaAppend = 'Soma';
        end
        
        %% Function:    getInitialFixSeq
        % Discription:  define fix seq using directory, allow special case
        function initialFixSeq = getInitialFixSeq(directory)
            % Input:    directory, struct/str
            % Output:   initialFixSeq, num
            
            if isobject(directory)
                directory = directory.main;
            end
            
            if contains(directory,{'09062019_65A_3rd','08282019_65A_2nd'})
                initialFixSeq = 8;
            else
                initialFixSeq = 9;
            end
        end
        
        %% Function:    getImageSetting
        function imageSetting = getImageSetting(sysSetting,directory)
            
            imageSetting = [];
            excludeSeq = [];    % Exclude cycles
            
            imageSetting.mouseID = '65A';
            
            imageSetting.sectionLoc = {'L','R'};
            
            imageSetting.chNum = 4;
            imageSetting.tileSize = [1200 1200];
            imageSetting.resolution = 0.55;
            imageSetting.slideThickness = 20;
            
            imageSetting.cameraOffset = 100;
            
            imageSetting.imageBits = 'uint16';
            imageSetting.overlapPrecentage = 15;
            imageSetting = TBS.getOverlapPixel(imageSetting);
            
            % Cell array for tile position
            % ie: 20x whole brain: [11, 20]; 10x whole brain: [6, 10];
            % per row: [nRow, nCol]
            tileSetting = [1 1; 2 2; 3 2; 3 3; 11 20];
            imageSetting.tilePos = TBS.getTilePosition(tileSetting);
            
            if nargin == 1 || isempty(directory)
                return
            end
            
            % Seq setting -------------------------------------------------
            % Get max seq number
            % getFolderNumber(directory,folderAppend)
            maxSeq = TBS.getFolderNumber(directory.main,sysSetting.seqPrepend);
            maxSeq = max(maxSeq);
            
            imageSetting.seqCycles = TBS.getSeqCycles(maxSeq,excludeSeq);
            
            imageSetting.initialFixSeq = TBS.getInitialFixSeq(directory.main);
        end
        
        %% Function:    getBcSetting
        % Discription:  barcode & bscalling setting
        function bcSetting = getBcSetting
            
            bcSetting = [];
            % Maximum intervel for sequencing result
            bcSetting.maxSeqInterval = 3;
            % Sequencing cycles are not included for counting BC length
            bcSetting.seqNotInMinBscallCycles = [9,10];
            % Minimum BC length
            % 08312021: nSeeq - exclude cycle number - maxHamming
            bcSetting.minBClen = 13;
            
            % Minimum digits of nucleotide
            bcSetting.minDiffNt = 3;
            % Max hamming distance for sequantical matching
            bcSetting.maxHamming = 2;
            
            % Number of nt for degenerate BC
            bcSetting.degenerateN = 6;
            
            % Setting for soma areas --------------------------------------
            % Section number of the soma zone (mannually identified)
            bcSetting.somaSectionRnge = [9 70];
            bcSetting.somaSlideRnge = ceil(bcSetting.somaSectionRnge./2);
            
            % min counts for soma pixel (arbitrary)
            bcSetting.nearSomaRolony.minPixelCount = 35;
            % min distance to soma pixel (arbitrary)
            % Note 10 pixel distance doesnt work
            bcSetting.nearSomaRolony.minDistance = 20;
            
            % Setting for count filters -----------------------------------
            % Count filter for soma and axon
            % Unspecific barcode is 9403
            bcSetting.minSomaCount = 0;
            bcSetting.maxSomaCount = 7000;
            bcSetting.minAxonCount = 3;
            bcSetting.maxAxonCount = 1000;
            
            % Min axonBC count for the strongest target region
            bcSetting.regionMinCount = 10;
            
            % Setting for isSoma ------------------------------------------
            % Soma radius range
            bcSetting.hasSoma.somaR = 100;
            % Min pixel count within radius
            bcSetting.hasSoma.minSomaPixelCount = 50;
            
            % Setting for rolony correction -------------------------------
            % max same dot dotID for different dot
            bcSetting.maxDiffDotId = 4;
            
            % Setting for floating rolony slide ---------------------------
            % Region to be included for floating rolony computation
            bcSetting.floatingSlide.regionName = {'Injection','Thalamus','Visual'};
            % Minimum rolony count for floating rolony
            % (not close to any rolony in the neighboring secitons)
            bcSetting.floatingSlid.threhold = 3;
            
            % Settinf for glia --------------------------------------------
            % Radius for glia cells
            bcSetting.gliaR = 200;         
        end
        
        %% Function:    getRefSetting
        % Discription:  get reference settings
        function refSetting = getRefSetting(directory)
            % Input:    directory, str, Directory of reference map & txt
            % Output:   refSetting, struct, reference settings
            
            refSetting = [];
            
            refSetting.directory = directory;
            
            % Reference scale of these maps: micron per pixel
            refSetting.refScale = 1/25;
            
            % reference map: allen ccf3, 25 um, nissle
            refSetting.nisslName = 'ara_nissl_25.nrrd';
            
            % Annotation map, same resolution
            refSetting.annoName = 'annotation_25.nrrd';
            
            % Average template, same resolution
            refSetting.avgName = 'average_template_25.nrrd';
            
            % Annotation structure ----------------------------------------
            % Name for txt file
            refSetting.annoTxtName = 'Allen_api_MouseBrainAtlas.txt';
            
            % Annotation structure
            cd(directory);
            annoStruct = fileread(refSetting.annoTxtName);
            annoStruct = jsondecode(annoStruct);
            annoStruct = annoStruct.msg;
            refSetting.annoStruct = annoStruct;
        end
        
        %% Function:    mouseSlideNum
        % Discription:  special case for slide number for each mouse
        function slideNum = mouseSlideNum(slideNum,mouseID)
            % Input & output: slideNum, vector/num, original slide number
            %            mouseID, str
            
            if strcmp(mouseID,'65A') == 1
                % ie. 12 change to 13
                slideNum(slideNum > 11.5) = slideNum(slideNum > 11.5) + 1;
                % Change to 12
                slideNum(slideNum == 11.5) = slideNum(slideNum == 11.5) + 0.5;
            end
        end
        
         %% Function:    QCtform
        % Discription:  quanlity check for tform using scale
        function TF = QCtform(tform,tolerance)
            % Input:    tform, mat, transformation matrix
            %           tolerance,num, tolerance of error, <1
            % Output:   TF, logical, whether it pass transformation matrix
            
            % Defualt tolerance
            if nargin == 1
                tolerance = 0.8;
            end
            
            % 07032021: exclude reflection
            if tform(1,1) < 0 && tform(2,2) < 0
                TF = false; return
            end
            
            % Check individual axis
            scale = sum(abs(tform(1:2,1:2)),2);
            % 1/toerance: > 1 (upper limit)
            TF1 = all(scale >= tolerance & scale <= 1/tolerance);
            
            % Check area size (determinant)
            % No reflection is allowed (no negative determinant)
            area = det(tform(1:2,1:2));
            TF2 = area >= tolerance & area <= 1/tolerance;
            
            TF = TF1 & TF2;
        end
        
        %% Function:    getTilePosition
        % Discription:  get tile position using number of rows and columns
        function tilePos = getTilePosition(tileSetting)
            % Note, tile position only support snake left-up pattern
            % updated 03132021
            % Input:    tileSetting, mat, one tiling patter per row
            %           1st column, row #; 2nd column, col #
            % Output:   tilePos, cell, one tiling pattern per cell
            
            tilePos = {};
            for i = 1:size(tileSetting)
                nRow = tileSetting(i,1);
                nCol = tileSetting(i,2);
                [nGrid,iTilePos] = indivTilePosition(nRow,nCol);
                tilePos{nGrid} = iTilePos;
            end
            
            % Function:    indivTilePosition ------------------------------
            % Discription:  get individual tile position
            function [nGrid,iTilePos] = indivTilePosition(nRow,nCol)
                % Input:    nRow/nCol, num, row/column number for tiling
                % Output:   nGrid, num, number of tile
                %           iTilePos, mat
                
                nGrid = nRow*nCol; % Number of grid
                
                % Construct current tile positions (snake left-up pattern)
                iTilePos = 1:nGrid;
                iTilePos = reshape(iTilePos,nCol,nRow);
                iTilePos = iTilePos';
                iTilePos(1:2:end,:) = fliplr(iTilePos(1:2:end,:));
                iTilePos = flipud(iTilePos);
            end
            
        end
                
        %% Function:    getChScore2Inten
        % Discription:  default chScore2Inten for 65A
        function chProfileSetting = getChScore2Inten(chProfileSetting)
            % This setting is mannually selected, use for all experiments
            % for 65A, uneven equal intensity will affect the alignment
            % using chCorrected image
            
            chProfileSetting.chScore2Inten.mean = 68.5;
            chProfileSetting.chScore2Inten.std = 39.5;
        end
        
        %% Function:    getChColor
        % Discription:  Get the color of the channel
        function iColor = getChColor(iCh)
            % Input:        iCh, number
            % Output:       iColor, str
            
            switch iCh
                case 1
                    iColor = 'c';
                case 2
                    iColor = 'y';
                case 3
                    iColor = 'm';
                case 4
                    iColor = 'k';
                otherwise
                    error('Function getChColor: Channel number is not supported');
            end
        end
        
        %% Function:    get2ndInfectionBC
        % Discription:  get barcode form secondary infection
        function BC = get2ndInfectionBC
            % Output:   mannual identified BC from secondary infection
            %           one BC per row; one col per seq cycle
            
            BC = [1	2	4	2	4	1	1	1	3	2	1	3	1	4	4	1	1;
                2	3	2	2	1	3	1	4	3	2	4	4	3	3	3	3	1;
                2	3	4	1	2	2	2	2	3	2	1	3	1	4	2	3	4;
                1	2	2	0	2	2	2	0	1	0	2	1	2	0	1	0	0;
                3	4	2	1	3	4	4	1	3	2	2	3	3	2	1	1	2;
                4	4	4	1	2	4	2	1	3	2	1	1	2	3	1	4	3;
                1	2	1	1	2	3	1	3	3	2	2	1	4	1	1	2	2;
                4	1	4	1	3	4	1	3	3	2	1	3	1	2	4	1	3;
                2	1	1	1	4	3	4	1	3	2	1	3	3	4	3	2	1;
                1	2	1	4	4	3	1	3	3	2	3	3	1	3	4	3	4;
                2	3	2	4	2	1	4	3	3	2	2	4	2	1	2	2	3;
                3	2	4	2	4	1	3	2	3	2	3	2	1	1	2	3	3;
                4	1	2	1	3	1	1	2	3	2	3	1	1	2	1	1	1;
                1	2	1	4	2	3	4	4	3	2	1	3	1	4	4	1	4;
                2	4	3	1	2	4	1	3	3	2	1	4	4	1	1	2	3;
                2	2	2	3	3	4	2	2	3	2	2	4	3	1	2	3	2;
                2	1	2	3	4	2	3	3	2	2	0	2	2	0	2	3	1];
        end
        
        %% Function:    delWrongTile
        % Discription:  delete wrong tiles (identified mannually)
        function tbl = delWrongTile(tbl)
            % Input & output: table
            
            % Wrong tile list (defined mannually)
            wrongTile = {'EF65ASlide49R_Seq15_Contra',...
                'EF65ASlide49R_Seq15_Visual',...
                'EF65ASlide49R_Seq15_SupCol',...
                };
            
            % Add delimiter for the tile name
            wrongTile = cellfun(@(X) [X,'_'],wrongTile,'UniformOutput',false);
            
            % Find and delete the row
            rowName = tbl.Properties.RowNames;
            row = contains(rowName,wrongTile);
            tbl = tbl(~row,:);
        end
        
        %% Function:    correctSelfAlignment
        % Discription:  correct self-alignment transformation matrix
        function corrSelfAlignTform = correctSelfAlignment(refSectionNum,imageSetting,sysSetting,directory)
            % Mannual defined correction
            % Input:    refSectionNum, num, reference section number
            %           imageSetting,sysSetting,directory, struct
            % Output:   corrSelfAlignTform, table, corrected self-aligned tform
            
            load(fullfile(directory.main,'selfAlignTform.mat'));
            
            % Section number
            sectionNum = TBS.getSectionNumber(selfAlignTform,imageSetting,sysSetting);
            
            % Translation correction --------------------------------------
            
            % Correction for translation
            corrTranslation = sectionNum - refSectionNum;
            % Mannual define translation for 0.055 scale factor: 0.4
            corrTranslation = abs(corrTranslation).*(-0.4);
            
            D = {};
            parfor i = 1:size(selfAlignTform,1) % parfor
                iTranslation = [corrTranslation(i),0];
                
                % Transformation matrix for translation
                iTform = eye(3);
                iTform(3,1:2) = iTranslation;
                iTform = affine2d(iTform);
                
                sz = size(selfAlignTform.D{i});
                
                % Combine displacement field
                iD = TBS.tform2D(iTform,sz);
                iD = TBS.Doperation(selfAlignTform.D{i},iD);
                
                D{i,1} = iD;
                disp(['Corrected selfAlignTform: ',num2str(i)]);
            end
            
            selfAlignTform.D = D;
            
            % Delete distored sections ------------------------------------
            % Mannually identified
            delSection = [108];
            
            row = ismember(sectionNum,delSection);
            selfAlignTform = selfAlignTform(~row,:);
            
            % Save --------------------------------------------------------
            corrSelfAlignTform = selfAlignTform;
            
            save(fullfile(directory.main,'corrSelfAlignTform.mat'),'corrSelfAlignTform');
        end
        
        %% Function:    getStitchingSeq
        % Discription:  set sequence for stitching tiles due to optical
        % issue of the microscope (uneven illumination)
        function I = getStitchingSeq(tileName,sysSetting,imageSetting)
            % Input:    tileName, cell, tile name
            %           sysSetting/imageSetting, struct
            % Output:   I, vector,
            
            tileNum = TBS.getTileNum(tileName,sysSetting);
            
            % Get tilng grid
            nTile = TBS.estimateTotalNumOfTile(tileNum,imageSetting.tilePos);
            tilePos = imageSetting.tilePos{nTile};
            
            % Exclude non-exist tile
            TF = ismember(tilePos,tileNum);
            tilePos(~TF) = 0;            
            
            [row,col,v] = find(tilePos);
            
            % Set sequence, specific to the microscope
            % Woodbury scope: from right-left, bottom-up
            [~,I] = sort(col,'descend');
            row = row(I); v = v(I);
            [~,I] = sort(row,'descend');
            v = v(I);
                        
            [~,I] = ismember(tileNum,v);
        end
        
        %% Function:    getRegVoxel
        % Discription:  get voxel from brain regions for analysis
        % CtxI, CtxC, Thal, Str, Mb
        function regVoxel = getRegVoxel(annoMap,annoStruct)
            % Input:    annoMap, mat, coronal annotation map
            %           annoStruct, struct, structure of annotation info
            % Output:   regVoxel, cell, logical stack of different brain
            % regions
            
            % Region fro 65A analysis
            regName = {'Isocortex'; 'Thalamus'; 'Striatum'; 'Midbrain'};
            
            % Voxel for each region
            id = cellfun(@(X) TBS.findAnnoID(annoStruct,X),regName,...
                'Uniformoutput',false);
            regVoxel = cellfun(@(X) ismember(annoMap,X),id,...
                'Uniformoutput',false);
                       
            % Get CtxI/C; ipsi is on the right of CCF
            regVoxel = [regVoxel(1);regVoxel];
            TF = TBS.isLeftIm(annoMap);
            regVoxel{1} = regVoxel{1} & ~TF;
            regVoxel{2} = regVoxel{2} & TF;
            
            % % Visualize region for grouping
            % test = regVoxel{1};
            % for i = 2:numel(regVoxel)
            %     iTest = regVoxel{i}.*i;
            %     test = max(test,iTest);
            % end
            % test = uint8(test);
            % MIJ.createImage(test)
        end
        
        %% Function:    isLocalProj
        % Discription:  whether a rolony is consider local projection
        function TF = isLocalProj(mlapdDot,mlapdSoma,localRng)
            % Update 01122022, use soma location and range
            % Input:    mlapdDot, mat, ML/AP/Depth coordinates of rolony
            %           mlapdSoma,
            % Output:   TF, logicalm, whether its local projection
            
            % Only use ML/AP-coordinates
            mlapdDot = mlapdDot(:,1:2);
            mlapdSoma = mlapdSoma(:,1:2);
            
            % Find injection center
            TF = any(mlapdSoma,2);
            mlapdSoma = mlapdSoma(TF,:);
            injCenter = median(mlapdSoma,1);
            
            D = pdist2(mlapdDot,injCenter);
            TF = D < localRng;
        end
        
        %% Function:    loadMLAPDim
        % Discription:  load ctxML/AP/DepthPrctile image without fold in
        % areas
        function [ctxML, ctxAP, ctxDepthPrctile] = loadMLAPDim(directory,refSetting)
            % Input:    directory, str, directory for ctxML/AP/DepthPrctile
            % image
            %           refSetting, struct
            % Output:   ctxML/AP/DepthPrctile, image
            
            cd(directory);
            load('ctxML.mat'); 
            load('ctxAP.mat'); 
            load('ctxDepthPrctile.mat');
            
            % Annotation map
            annoMap = TBS.getRefMap('anno',refSetting);
            
            % Annotation structure
            annoStruct = refSetting.annoStruct;          
            
            % Delete the fold in portion for better visualization
            roi = {'Retrosplenial area, ventral part', ...
                'Anterior cingulate area',...
                'Retrosplenial area, dorsal part'};
            id = cellfun(@(X) TBS.findAnnoID(annoStruct,X),roi,'Uniformoutput',false);
            id = vertcat(id{:});
            roi = ismember(annoMap,id);
            
            % Delete from ML/AP/Depth image
            ctxML(roi) = 0; ctxAP(roi) = 0; ctxDepthPrctile(roi) = 0;
        end
        
    end
    
    methods (Static)    % File, folder and general BARseq function ========
        %% Function:    seqstr
        % Discription:  get string with sequencing prepend
        function seqStr = seqstr(num)
            
            if mod(num,1)~=0
                num = num2str(num);
            else
                num = num2str(num,'%02d');
            end
            
            sysSetting = TBS.getSysSetting;                        
            seqStr = strcat(sysSetting.seqPrepend,num);
        end
        
        %% Function:    getFolderNumber
        % Discription:  Find the number of the folder with the
        % append under the directory
        function folderNum = getFolderNumber(directory,folderPrepend)
            % Input:    directory, str
            %           folderPrepend, str, append of the selected folders
            % Output:   folderNum, vector
            
            % Find all the folders with the append
            folder = fullfile(directory,[folderPrepend,'*']);
            folder = dir(folder);
            folder = folder([folder.isdir]);
            
            % Get the numbers attached to the append
            folderNum = erase({folder.name},folderPrepend);
            folderNum = cellfun(@str2num,folderNum);
        end
        
        %% Funciton:    getSeqCycles
        % Discription:  get sequencing cycle numbers
        function seqCycles = getSeqCycles(maxSeq,excludeSeq)
            % NOte, assume the sequencing cycle start with 1
            seqCycles = 1:maxSeq;
            seqCycles(ismember(seqCycles,excludeSeq)) = [];
        end
        
        %% Function:    nameFun
        % Discription:  extra elements from the name
        function name = nameFun(name,element,sysSetting)
            % Input & output:    name, str,
            %           element, num/mat, element of the name to extract
            %           sysSetting, struct
            
            delimiter = sysSetting.delimiter;
            
            name = strsplit(name,delimiter);
            name = name(element);
            
            % Output as char if there is only one element
            if numel(element) == 1
                name = name{:};
            else
                name = strjoin(name,delimiter);
            end
        end
        
        %% Function:    rowName2ImName
        % Discription:  get the image name the tile belogns to
        function imName = rowName2ImName(tbl,addDelimiter,sysSetting)
            % Input:   tbl, table
            %          addDelimiter, logical, whether add delimiter to the
            %          end
            %          sysSetting,struct
            % Output:  imName, cell, image name the tile belongs to
            
            imName = tbl.Properties.RowNames;
            
            imName = erase(imName,sysSetting.imFormat);
            
            % Extract image name
            nameElement = sysSetting.nameElement;
            imName = cellfun(@(X) TBS.nameFun(X,nameElement,sysSetting),...
                imName,'UniformOutput',false);
            
            % Add delimiter to differenciate Vis vs VisContra
            if addDelimiter
                imName = cellfun(@(X) [X,sysSetting.delimiter],imName,...
                    'UniformOutput',false);
            end
        end
        
        %% Function:    getTileNum
        % Discription:  Get tile number in the image name
        function tileNum = getTileNum(tileName,sysSetting)
            % 03062022, only accept string input
            % Input:    tileName, str
            %           sysSetting, struct
            % Output:   tileNum, num
                        
            tileNum = TBS.nameFun(tileName,sysSetting.tileElement,sysSetting);
            tileNum = str2double(tileNum);
        end
        
        %% Function:    getSeqNum
        % Discription:  extract seq number from name
        function seq = getSeqNum(name,sysSetting)
            % Input:    name, str
            %           sysSetting, sturct
            % Output:   seq, num, sequence number
            
            seq = TBS.nameFun(name,sysSetting.seqElement,sysSetting);
            seq = erase(seq,sysSetting.seqPrepend);
            seq = str2double(seq);
        end
               
        %% Function:    imFormat
        % Discription:  convert to image file name
        function imFormatName = imFormat(name,sysSetting)
            % Input:    name, str, name to convert
            %           sysSetting, struct
            % Otput:    imFormatName, str, with imFormat
            delimiter = sysSetting.delimiter;
            
            % Get rid of extra delimiter
            name = strsplit(name,delimiter);
            
            I = cellfun(@isempty,name);
            name = name(~I);
            
            name = strjoin(name,delimiter);
            
            imFormatName = [name,sysSetting.imFormat];
        end
        
        %% Function:    ensureRowNameAppend
        % Discription:  check and correct whether the row names has the
        % specific append
        function tbl = ensureRowNameAppend(tbl,str)
            % Input:        tbl, table
            %               str, str needs to be added to the row names
            % Output:       tbl, table, with corrected row name
            
            tableRowNames = tbl.Properties.RowNames;
            
            % Row doesnt have the name
            TF = ~contains(tableRowNames,str);
            
            if ~any(TF)
                return
            end
            
            tableRowNames(TF) = cellfun(@(X) [X,str],tableRowNames(TF),'Uniformoutput',false);
            
            tbl.Properties.RowNames = tableRowNames;
        end
              
        %% Funciton:    estimateTotalNumOfTile
        % Discription:  get the estimated total tile number
        function totalTileNum = estimateTotalNumOfTile(tileNum,tilePos)
            % Input:    tileNum, mat, list of tile number
            %           tilePos, cell, tile positions
            % Output:   totalTileNum, num, estimated number
            
            maxTileNum = max(tileNum);
            
            % Find the total tile number available
            possibleTileNum = find(~cellfun(@isempty,tilePos));
            % Find the one closest to it on the bigger side
            totalTileNum = possibleTileNum(maxTileNum <= possibleTileNum);
            totalTileNum = min(totalTileNum);
        end
        
        %% Function:    getAlignmentSeq
        % Discription:  calcualte the sequence for alignment
        function alignmentSeq = getAlignmentSeq(numIn,refNum)
            % Input:    numIn, mat, available section/seq number
            %           refNum, num, reference section number
            % Output:   alignmentSeq, mat, sequence of section number for
            %           the alignment
            
            % Compute the sequence
            if refNum ==  min(numIn)
                alignmentSeq = min(numIn):max(numIn);
            elseif refNum ==  max(numIn)
                alignmentSeq = fliplr(min(numIn):max(numIn));
            else
                % Note: refSecitonNum is repeated once
                alignmentSeq = [refNum:max(numIn),...
                    fliplr(min(numIn):refNum)];
            end
            
            % Delete number when there is no corresponding seciton number
            keepNum = ismember(alignmentSeq,numIn);
            alignmentSeq = alignmentSeq(keepNum);
        end
        
        %% Function:    getAccumulateTform
        % Discription:  accumulate transformation matrix
        function accumulateTform = getAccumulateTform(tformTable,...
                movingVar,fixVar,tformVar)
            % Input:    tformTable, table, with fixVar & tformVar
            %           movingVar, str, variable name of moving column/empty
            %           fixVar, str, variable names of the fix column
            %           tformVar, str, variable names of the tform column
            % Output:   accumulateTform, table
            
            % Assigned row names as moving column
            if isempty(movingVar)
                movingVar = 'movingName';
                tformTable.(movingVar) = tformTable.Properties.RowNames;
            end
            
            % Find the intial fix tform by finding image without fix name
            row = cellfun(@isempty,tformTable{:,fixVar});
            
            % ! Done\'t assume the initial fix tform is eye(3)
            accumulateTform = table();
            accumulateTform.(tformVar) = tformTable.(tformVar)(row);
            accumulateTform.Properties.RowNames = tformTable.(movingVar)(row);
            
            tformTable(row,:) = [];
            
            % Add the remaining image to the table if its fix can be found
            % on the table
            i = 1;
            while i <= size(tformTable,1)
                
                % Current file name
                movingName = tformTable.(movingVar){i};
                
                irow = tformTable.(movingVar);
                irow = contains(irow,movingName);
                iTbl = tformTable(irow,:);
                
                % All fix tiles
                fixName = iTbl.(fixVar);
                [~,lia,lib] = intersect(accumulateTform.Properties.RowNames,fixName);
                
                % If its fix tile is not found in the accumulate tform
                if numel(lib) < numel(fixName)
                    i = i+1;
                    continue
                end
                
                % Accumulate fix tform
                fixTform = accumulateTform.(tformVar)(lia);
                % Current tform
                tform = iTbl.(tformVar)(lib);
                
                % Add pervious tform
                % for example: B = AT; BR = A(TR);
                tform = cellfun(@(X,Y) X*Y, tform,fixTform,'Uniformoutput',false);
                
                % Get median tform
                tform = reshape(tform,1,1,[]);
                tform = cell2mat(tform);
                tform = median(tform,3);
                
                % Add to output table
                accumulateTform(movingName,tformVar) = {tform};
                
                % Delete the aligned image
                tformTable(irow,:) = [];
                
                % Restart the loop
                if i == size(tformTable,1) && size(tformTable,1) > 0
                    i = 1;
                end
            end
        end
        
        %% Function:    doSomaBscall
        % Discription:  whether do soma bascall for current image
        function TF = doSomaBscall(imName,sysSetting)
            % Input:        iImName, str
            %               sysSetting,structure
            % Do somaBscall if this is an injection site
            if contains(imName,sysSetting.injectionAppend)==1
                TF = true;
            else
                TF = false;
            end
        end
        
        %% Function:    getSlideNumLoc
        % Discription:  get slide# and section location from image name
        function [slideNum,sectionLocation] = getSlideNumLoc(imName,...
                sysSetting,possibleSectionLoc)
            % Input:    imName, str
            % Output:   slideNum & sectionLocation, num
            
            % Get section location and section name
            imName = strsplit(imName,{sysSetting.delimiter,sysSetting.slidePrepend});
            imName = imName{2};
            
            I = cellfun(@(X) contains(imName,X),possibleSectionLoc);
            sectionLocation = find(I);
            
            slideNum = erase(imName,possibleSectionLoc{I});
            slideNum = str2double(slideNum);
        end
                
        %% Function:    getSectionNumber
        % Discription:  Number all brain secitons for analysis
        function sectionNum = getSectionNumber(input,imageSetting,sysSetting)
            % Input:    input, image names
            %           use row names for table
            %           imageSetting/sysSetting, struct
            % Output:   sectionNum, mat
            
            mouseID = imageSetting.mouseID;
            
            % Prepare input into cell array -------------------------------
            if istable(input) == 1
                input = input.Properties.RowNames;
            elseif ischar(input) == 1
                input = cellstr(input);
            end
            
            % Add slide prepend if there is none --------------------------
            rows = ~contains(input,sysSetting.slidePrepend);
            input(rows) = cellfun(@(X) [sysSetting.slidePrepend,X],...
                input(rows),'Uniformoutput',false);
            
            % Get possible section locaiton for the animal ----------------
            sectionLoc = imageSetting.sectionLoc;
            nSectionLoc = numel(sectionLoc);
            
            [slideNum,sectionLocation] = cellfun(@(X) ...
                TBS.getSlideNumLoc(X,sysSetting,sectionLoc),input);
            
            % Fix slide number according for this mouse
            slideNum = TBS.mouseSlideNum(slideNum,mouseID);
            
            % Calculate section number
            sectionNum = (slideNum-1).*nSectionLoc + sectionLocation;
        end
        
    end
    
    methods (Static)    % Fuse image ======================================
        %% Function:    getOverlapPixel
        % Discription:  calculate overlap pixel number
        function imageSetting = getOverlapPixel(imageSetting)
            % Input & output: imageSetting,struct
            
            tileSize = imageSetting.tileSize;
            overlapPrecentage = imageSetting.overlapPrecentage;
            
            % In case sometimes the input is fraction instead of precentage
            if overlapPrecentage >= 1
                overlapPrecentage = overlapPrecentage.*0.01;
            end
            
            overlapPixel = tileSize.*overlapPrecentage;
            overlapPixel = round(overlapPixel);
            imageSetting.overlapPixel = overlapPixel;
        end
        
        %% Funciton:    getImageTileNames
        % Discription:  get the tile name for the image
        function tileName = getImageTileNames(tileName,imageName,delimiter)
            % Input:    tileName, cell, all the tile names
            %           imageName, str, the input image name
            % Output:   tileName, cell, tile names for the current image
            
            imageName = strsplit(imageName,delimiter);
            
            % Add delimiter if for the element comes later, for perfect
            % match. ie, seperate Contra/Visual from VisulContra
            imageName(2:end) = cellfun(@(X) [delimiter,X,delimiter],...
                imageName(2:end),'UniformOutput',false);
            
            % Find tile names have all the elements
            imageName2 = repmat(imageName,size(tileName,1),1);
            tileName2 = repmat(tileName,1,size(imageName,2));
            
            row = cellfun(@(X,Y) contains(X,Y),tileName2,imageName2);
            row = all(row,2);
            
            tileName = tileName(row);
        end
        
        %% Function:    getMovingTileNum
        % Discription:  get the matching/moving tile number for the fix tile
        function movingTileNum = getMovingTileNum(fixTileNum,tilePos,tileNum)
            % Input:    fixTileNum, num,
            %           tilePos, mat, tile positions for imaging grid
            %           tileNum, mat, tile number with images available
            % Output:   matchingTileNum, moving tile number for the fixed tile
            
            % Get the neighboring tile number
            movingTileNum = tilePos == fixTileNum;
            SE = strel('diamond',1);
            movingTileNum = imdilate(movingTileNum,SE);
            % Delete the fix tile itself
            movingTileNum(tilePos == fixTileNum) = false;
            
            % Get tile number within the available tile number
            movingTileNum = tilePos(movingTileNum);
            movingTileNum = intersect(movingTileNum,tileNum);
        end
        
        %% Funciton:    getInitialFixTileNum
        % Discription:  get the middle tile as initial tile number
        function fixTileNum = getInitialFixTileNum(tilePos,tileNum)
            % Input:    tilePos, mat, tile positions in the imaging grid
            %           tileNum, mat, tile number with image available
            % Output:   fixTileNum, num, intial fix tile
            
            % Start with the middle tile
            fixTileNum = round(size(tilePos)./2);
            fixTileNum = tilePos(fixTileNum(1),fixTileNum(2));
            
            % If the tile number doesn't have an image available
            % Find directly connected neigboring tiles
            if ~ismember(fixTileNum,tileNum)
                
                fixTileNum2 = TBS.getMovingTileNum(fixTileNum,tilePos,...
                    tileNum);
                
                if isempty(fixTileNum2)
                    error('Function getInitialFixTileNum: currently only support direct connected tiles.')
                end
                
                fixTileNum = fixTileNum2(1);
            end
        end
        
        %% Function:    prepareFixMovingIm4Alignment
        % Discription: Get only the overlap area of moving and fix
        function [fixIm,movingIm] = prepareFixMovingIm4Alignment(fixIm,movingIm,...
                fixTileNum,movingTileNum,tilePos,overlapPixel)
            % Input:    fixIm/movingIm, fixed/moving image
            %           fixTileNum/movingTileNum, num, fix/moving tile number in tile
            %           position
            %           tilePos, mat, tile position of the fuse image
            %           overlapPixel, num, overlap pixel number
            % Output:   fixIm/movingIm, image only with the overlap area
            
            [fixRow,fixCol] = find(tilePos == fixTileNum);
            [movingRow,movingCol] = find(tilePos == movingTileNum);
            
            % Set the non-overlap area to 0
            % This is more accurate than only include the overlap part (?)
            if fixRow < movingRow
                fixIm(1:end-overlapPixel(1),:,:) = 0;
                movingIm(overlapPixel(1)+1:end,:,:) = 0;
            elseif fixRow > movingRow
                fixIm(overlapPixel(1)+1:end,:,:) = 0;
                movingIm(1:end-overlapPixel(1),:,:) = 0;
            elseif fixCol < movingCol
                fixIm(:,1:end-overlapPixel(2),:) = 0;
                movingIm(:,overlapPixel(2)+1:end,:) = 0;
            elseif fixCol > movingCol
                fixIm(:,overlapPixel(2)+1:end,:) = 0;
                movingIm(:,1:end-overlapPixel(2),:) = 0;
            end
            
        end
        
        %% Function:    stitchTileChTform
        % Discription:  get transformation matrix by aligning in specific
        % channels
        function tform = stitchTileChTform(iCh,alignmentSetting,tform,...
                fixTileNum,tileNum,tileName,iTilePos,overlapPixel)
            % Input:        alignCh, num, channel(s) for alignment
            %               translaitonOnly, logical, whether only do translaiton for
            %               the alignment
            %               tform, cell, transformation matrixes
            %               fixTileNum, num, fix tile number
            %               tileNum, mat, tile number available
            %               tileName, cell, corresponding name of the tiles
            %               iTilePos, mat, tile positions in imaging grid
            %               overlapPixel, mat, pixel number in overlap areas
            % Output:       tform, cell, updated transformation matrixes
            
            % Special case for only one tile
            if numel(tform) == 1
                tform = {eye(3)};
                return
            end
            
            if isempty(tform)
                tform = cell(numel(iTilePos),1);
            end
            
            alignmentCh = alignmentSetting.ch{iCh};
            
            if ~contains('method',fieldnames(alignmentSetting))
                alignmentMethod = [];
            elseif numel(alignmentSetting.method) < iCh
                alignmentMethod = [];
            else
                alignmentMethod = alignmentSetting.method{iCh};
            end
            
            % Get intial fixed tile =======================================
            if isempty(fixTileNum)
                fixTileNum = TBS.getInitialFixTileNum(iTilePos,tileNum);
            end
            tform{fixTileNum} = eye(3);
            
            fixTileName = tileName{tileNum == fixTileNum};
            disp(['Intial fixed tile: ',fixTileName]);
            
            % Align the rest of the tiles =================================
            % Tile hasn't used a fixed tile
            availableFix = tileNum;
            
            while ~isempty(availableFix)
                
                fixTileName = tileName{tileNum == fixTileNum};
                fixTform = tform{fixTileNum};
                
                % Get moving tile number for the current fix
                movingTileNum = TBS.getMovingTileNum(fixTileNum,iTilePos,tileNum);
                movingTileNum = reshape(movingTileNum,1,[]);
                
                % Loop through each matching tile =========================
                for iMovingTileNum = movingTileNum
                    
                    % If the tile already aligned
                    if ~isempty(tform{iMovingTileNum})
                        continue
                    end
                    
                    % Get fix & moving images -----------------------------
                    movingTileName = tileName{tileNum == iMovingTileNum};
                    
                    movingIm = TBS.getStack(movingTileName,alignmentCh);
                    fixIm = TBS.getStack(fixTileName,alignmentCh);
                    
                    % Use method if there is one
                    if ~isempty(alignmentMethod)
                        fixIm = alignmentMethod(fixIm);
                        movingIm = alignmentMethod(movingIm);
                    end
                    
                    % Set non-overlap region to 0 for alginment
                    [fixIm,movingIm] = TBS.prepareFixMovingIm4Alignment...
                        (fixIm,movingIm,fixTileNum,iMovingTileNum,iTilePos,overlapPixel);
                    
                    % Alignment: Allow use multiple channel for alignment -
                    % Get tform for each stack
                    iTform = zeros(3,3,size(fixIm,3));
                    for iStack = 1:size(fixIm,3)
                        if alignmentSetting.translaitonOnly == true % Translation only
                            iStackTform = imregcorr(movingIm(:,:,iStack),...
                                fixIm(:,:,iStack),'translation','Window',false);
                        else % default similarity
                            iStackTform = imregcorr(movingIm(:,:,iStack),...
                                fixIm(:,:,iStack),'Window',false);
                        end
                        
                        iTform(:,:,iStack) = iStackTform.T;
                    end
                    
                    % Delete tfrorm outlier -------------------------------
                    if size(iTform,3)~=1
                        outlier = isoutlier(iTform,3);
                        outlier = any(outlier,[1 2]);
                        iTform = iTform(:,:,~outlier);
                    end
                    
                    if isempty(iTform)
                        warning(['Did not matched: ',movingTileName])
                        continue
                    end
                    
                    % Get median as tform ---------------------------------
                    iTform = median(iTform,3);
                    
                    % QC: a very small scale factor
                    if ~TBS.QCtform(iTform)
                        warning(['Did not matched: ',movingTileName])
                        continue
                    end
                    
                    % Add the tform to the fix tform ----------------------
                    tform{iMovingTileNum} = iTform*fixTform;
                    
                    disp(['Aligned: ',movingTileName]);
                end
                
                % Prepare for the next cycle ------------------------------
                % Delete the current ref (fix) form the ref list
                availableFix(availableFix == fixTileNum) = [];
                
                % Find the next fix tile number
                fixTileNum = find(~cellfun(@isempty,tform));
                fixTileNum = intersect(fixTileNum,availableFix);
                if isempty(fixTileNum) 
                    break
                else % Pick the middle of it
                    I = ceil(numel(fixTileNum)/2);
                    fixTileNum = fixTileNum(I);
                end
            end
            
        end
        
        %% Function:    stitchTileInfo
        % Discription:  get tile stiching infomation
        function tblOut = stitchTileInfo(tileName,tileNum,imageSetting,...
                alignmentSetting,fixTileNum,directory)
            % Output:   tblOut, table, stitch info for the current image
            
            tilePos = imageSetting.tilePos;
            overlapPixel = imageSetting.overlapPixel;
            ch = alignmentSetting.ch;
            
            % Get tilePos using tile number & tilePos ---------------------
            % for tile number < grid it created
            nTile = TBS.estimateTotalNumOfTile(tileNum,tilePos);
            iTilePos = imageSetting.tilePos{nTile};
            
            % Confirm the directory ---------------------------------------
            cd(directory.imFolder)
            iFolder = dir(['**\',tileName{1},'*']);
            iFolder = iFolder.folder;
            cd(iFolder);
            
            % Get tform ---------------------------------------------------
            tform = cell(numel(iTilePos),1);
            
            % Go through all alignment channel combinations
            % Trie new channels if it didnt aligned previously
            for iCh = 1:length(ch)
                
                tform = TBS.stitchTileChTform(iCh,alignmentSetting,...
                    tform,fixTileNum,tileNum,tileName,iTilePos,overlapPixel);
                
                % Stop if all the tiles were aligned
                if all(~cellfun(@isempty,tform(tileNum)))
                    break
                end
            end
            
            % Update table output -----------------------------------------
            % 03142021: delete append
            tblVar = {'directory','imageTile','tilePos','tform'};
            tblOut = table();
            tblOut(1,tblVar) = [{iFolder},{tileName},{iTilePos},{tform}];
        end
        
        %% Function:    stitchImInfo (main)
        % Discription:  get information for stitching images
        function tformStitchTable = stitchImInfo(tformStitchTable,...
                alignmentSetting,imageSetting,sysSetting,directory)
            % Input:    tformStitchTable, table
            %           alignmentSetting, struct, setting for alignment
            %           imageSetting/sysSetting/dirctory, struct,
            % Output:   tformStitchTable, table, with stitching info. also saved on
            % disk
            
            % Variable name for saving on disk
            tformVarName = alignmentSetting.tformVarName;
            
            % Append of image for alignment
            imAppend = alignmentSetting.imAppend;
            
            % Whether skip the aligned image
            if contains('skipAlignedIm',fieldnames(alignmentSetting))
                skipAlignedIm = alignmentSetting.skipAlignedIm;
            else
                skipAlignedIm = false;
            end
            
            % All tif with the append under the current folder ------------
            cd(directory.imFolder);
            fileName = ls(['**\*',imAppend,'*']);
            fileName = cellstr(fileName);
            
            % Get unique image names
            imageName = cellfun(@(X) TBS.nameFun(X,sysSetting.nameElement,...
                sysSetting),fileName,'Uniformoutput',false);
            imageName = unique(imageName);
            
            % Get stiching info for these images --------------------------
            % 03282021: switch to new output for parfor
            tformStitchTable2 = table(); rowNames = {};
            parfor iIm = 1:length(imageName) %parfor
                % Current image name
                iImageName = imageName{iIm};
                
                % Skip aligned images
                if skipAlignedIm == true && TBS.hasThisRow(tformStitchTable,iImageName)
                    continue
                end
                
                % Pick out tile name for the current image
                tileName = TBS.getImageTileNames(fileName,iImageName,sysSetting.delimiter);
                % Get tile number for this section
                tileNum = TBS.getTileNum(tileName,sysSetting);
                
                fixTileNum = [];
                tblOut = TBS.stitchTileInfo(tileName,tileNum,imageSetting,...
                    alignmentSetting,fixTileNum,directory);
                
                rowNames{iIm,1} = iImageName;
                tformStitchTable2(iIm,:) = tblOut;
            end
            
            if numel(tformStitchTable2) ~= 0
                row = ~cellfun(@isempty,rowNames);
                rowNames = rowNames(row);
                tformStitchTable2 = tformStitchTable2(row,:);
                
                varNames = {'directory','imageTile','tilePos','tform'};
                tformStitchTable(rowNames,varNames) = tformStitchTable2;
                
                save(fullfile(directory.stitched,[tformVarName,'.mat']),'tformStitchTable');
            end
        end
        
        %% Function:    estimateStitchImSize
        % Discription:  estimate size for fused image
        function outputSize = estimateStitchImSize(tilePos,...
                imageSetting,paddingSize)
            % Input:    tilePos, mat
            %           imageSetting, struct
            %           paddingSize, 1x2 mat
            % Output:   outputSize, 1x2 mat
            
            tileSize = imageSetting.tileSize;
            overlapPixel = imageSetting.overlapPixel;
            
            outputSize = tileSize.*size(tilePos);
            outputSize = outputSize-overlapPixel.*(size(tilePos)-1);
            outputSize = outputSize + paddingSize;
        end
        
        %% Function:    getPaddingTform
        % Discription:  get padding transformation matrix for fuse image
        % (can take tform into consideraton)
        function paddingTform = getPaddingTform(outputSize,tform)
            % Input:    outputSize, mat
            %           tform, cell,
            % Output:   paddingTform, 3 x 3 mat
            
            % Get center of the image
            paddingTform = eye(3);
            paddingTform(3,1:2) = fliplr(outputSize)./2;
            
            if isempty(tform)
                return
            end
            
            % Get the median translation of all tform to make sure the center of
            % tform is in the center of image (necessary to reduce errors)
            tformTranslation = tform(~cellfun(@isempty,tform));
            tformTranslation = cellfun(@(X) X(3,1:2),tformTranslation,...
                'UniformOutput',false);
            tformTranslation = cell2mat(tformTranslation);
            tformTranslation = median(tformTranslation,1);
            paddingTform(3,1:2) = paddingTform(3,1:2) - tformTranslation;
        end
        
        %% Function:    getTrimSetting
        % Discription:  get setting for trimming the fuse image
        % (do pad-trim because it's not straight forward to estimate the
        % size of output)
        function [imSize,trimTform] = getTrimSetting(fuseIm)
            % Input:    fuseIm, image stack
            % Output:   imSize, row vector, trimed image size
            %           trimTform, affine2d object, tform for translation
            
            % Find valide pixels
            fuseIm = any(fuseIm,3);
            rowLim = TBS.getImValidLim(fuseIm,1);
            colLim = TBS.getImValidLim(fuseIm,2);
            
            % Trimed image size
            imSize = [rowLim(end)-rowLim(1), colLim(end)-colLim(1)]+1;
            
            % Tform for triming
            trimTform = eye(3);
            trimTform(3,1:2) = [colLim(1), rowLim(1)].*-1+1;
            trimTform = affine2d(trimTform);
        end
        
        %% Function:    fuseImage (main)
        % Discription:  fuse image according to the transformation matrix
        % Max projeciton in the overlap region
        function tformStitchTable = fuseImage(tformStitchTable,...
                alignmentSetting,imageSetting,sysSetting,directory)
            % Input:    tformStitchTable, table
            %           alignmentSetting, struct, setting for alignment
            %           imageSetting/sysSetting/dirctory, struct,
            % Output:   tformStitchTable, table, with stitching info. also saved on
            % disk
            
            paddingSize = alignmentSetting.paddingSize;
            tformVarName = alignmentSetting.tformVarName;
            
            % Image process (i.e. max project-@max)
            if ~contains('method',fieldnames(alignmentSetting))
                alignmentMethod = [];
            elseif isempty(alignmentSetting.method)
                alignmentMethod = [];
            else
                alignmentMethod = alignmentSetting.method{1};
            end
            
            % Fuse image --------------------------------------------------
            for iIm = 1:size(tformStitchTable,1)
                
                % Fuse name
                fuseImName = tformStitchTable.Properties.RowNames{iIm};
                fuseImName = [fuseImName,sysSetting.imFormat];
                
                % Only stitch absent images
                if exist(fullfile(directory.stitched,fuseImName))
                    continue
                end
                
                % Get tilePos, tform, tile names
                iTilePos = tformStitchTable.tilePos{iIm};
                tileName = tformStitchTable.imageTile{iIm};
                tform = tformStitchTable.tform{iIm};
                % (Need tile number because length of tileName and tform
                % are not the same)
                tileNum = TBS.getTileNum(tileName,sysSetting);
                
                % Estimate stitched image size ----------------------------
                % (Use this instead of imref2d to have a unique imSize
                % during the loop)
                outputSize = TBS.estimateStitchImSize(iTilePos,...
                    imageSetting,paddingSize);
                R = imref2d(outputSize);
                
                % Get padding value to center the image
                paddingTform = TBS.getPaddingTform(outputSize,tform);
                
                % Get directory of tiles ----------------------------------
                % Modified the disk name in case there is a change
                iDirectory = tformStitchTable.directory{iIm};
                iDirectory(1) = directory.imFolder(1);
                cd(iDirectory);
                
                % Transform image -----------------------------------------
                nCh = [];
                for iTile = 1:length(tform)
                    iTform = tform{iTile};
                    
                    if isempty(iTform)
                        continue
                    end
                    
                    iTileName = tileName{tileNum == iTile};
                    iTileImage = TBS.getStack(iTileName,[]);
                    
                    % Process image according to method
                    if ~isempty(alignmentMethod)
                        iTileImage = alignmentMethod(iTileImage);
                    end
                    
                    % Initiate new empty image
                    % With current channel size and bit
                    if isempty(nCh)
                        nCh = size(iTileImage,3);
                        fuseIm = zeros([outputSize,nCh]);
                        
                        imageBits = class(iTileImage);
                        fuseIm = cast(fuseIm,imageBits);
                    end
                    
                    % Add padding to the tform to center the image
                    iTform = iTform * paddingTform;
                    iTform = affine2d(iTform);
                    
                    % Transform image
                    iTileImage = imwarp(iTileImage,iTform,'Outputview',R);
                    
                    % Max projection to fuse image
                    fuseIm = max(fuseIm,iTileImage);
                    
                    disp(['Fused tile: ',iTileName]);
                end
                
                % Trim out the extra edges of empty pixels ----------------
                [imSize,trimTform] = TBS.getTrimSetting(fuseIm);
                
                % Trim and save the image
                fuseIm = imwarp(fuseIm,trimTform,'Outputview',imref2d(imSize));
                TBS.saveStack(fuseIm,fullfile(directory.stitched,fuseImName));
                
                % Save trimming info --------------------------------------
                % Combine padding and triming tform
                paddingTform = paddingTform*trimTform.T;
                
                % Add the padding and imSize to table output
                tformStitchTable.paddingTform{iIm} = paddingTform;
                tformStitchTable.imSize{iIm} = imSize;
                
                save(fullfile(directory.stitched,[tformVarName,'.mat']),'tformStitchTable');
                
                disp(['Save image:',fuseImName]);
            end
            
        end
        
    end
    
    methods (Static)    % Dot table =======================================
        %% Function:    delLowBscallInten
        % Discription:  delete row with bscall intensity < max inten or = 0
        function dotTable = delLowBscallInten(dotTable)
            % Input & output:   dotTable, table, with chIntensity &
            % bscallCh
            
            chIntensity = dotTable.chIntensity;
            
            maxInten = max(chIntensity,[],2);
            
            sz = size(chIntensity);
            ind = sub2ind(sz,(1:sz(1))',dotTable.bscallCh);
            
            % Get rows with maxInten in the bascall channel
            row = chIntensity(ind) == maxInten & maxInten > 0;
            dotTable = dotTable(row,:);
        end
        
        %% Function:    getTileDotTable
        % Discription:  get spot as local maxima & max intensity
        function dotTable = getTileDotTable(iSpotMat,iTileImage)
            % Input:        iSpotMat, cell, spot mat for current tile
            %               one cell per channel (position of local maxima)            %
            %               iTileImage, mat, image stack
            % Output:       dotTable, table, one row for one dot
            
            % Row: index; column: channel
            iTileImage = reshape(iTileImage,[],size(iTileImage,3));
            
            index = {};
            bscallCh = {};
            chIntensity = {};
            for iCh = 1:numel(iSpotMat)
                % Get all spots from image processing
                iChSpot = iSpotMat{iCh};
                iChSpot = reshape(iChSpot,[],1);
                
                % Dot index
                iIndex = find(iChSpot);
                
                % Dot intensity in all channels
                iChIntensity = iTileImage(iIndex,:);
                
                % column for bscall channel number
                iBscallCh = repmat(iCh,size(iIndex));
                
                index{iCh,1} = iIndex;
                bscallCh{iCh,1} = iBscallCh;
                chIntensity{iCh,1} = iChIntensity;
            end
            
            index = vertcat(index{:});
            bscallCh = vertcat(bscallCh{:});
            chIntensity = vertcat(chIntensity{:});
            
            bscallCh = uint8(bscallCh);
            
            dotTable = table(index,bscallCh,chIntensity);
        end
        
        %% Function:    getDotTable (main)
        % Discription:  get dot bscall according to the spotMat.mat
        function dotTable = getDotTable(spotMatName,imageSetting,sysSetting,directory)
            % Input:        spotMatName,str
            %               imageSetting & sysSetting & directory, struct
            % Output:       bscallTable, table, one row for one dot
            %               VariableNames: index,bscallIntensity,chIntensity
            
            spotMatName = strcat(spotMatName,'.mat');
            chNum = imageSetting.chNum;
            
            dotTable = table();
            for iSeq = imageSetting.seqCycles
                % Go to seq folder
                iSeqFolder = TBS.seqstr(iSeq);
                iSeqFolder = fullfile(directory.main,iSeqFolder);
                cd(iSeqFolder)
                
                % Local maxima for each channel (in .mat)
                spotMat = load(spotMatName);
                spotMat = spotMat.spotImage;
                
                % Tile name; if the row name doesnt have the tile
                % append & format, added it
                spotMat = TBS.ensureRowNameAppend(spotMat,sysSetting.tileAppend);
                spotMat = TBS.ensureRowNameAppend(spotMat,sysSetting.imFormat);
                
                iDotTable = cell(size(spotMat));
                for iTile = 1:size(spotMat,1)  % parfor
                    iTileName = spotMat.Properties.RowNames{iTile};
                    
                    iTileImage = TBS.getStack(iTileName,1:chNum);
                    
                    % Get spot matrix from spotMat
                    iSpotMat = spotMat{iTile,1};
                    if numel(iSpotMat) == 1
                        error('???????');  %???????
                        iSpotMat = iSpotMat{:};
                    end
                    
                    % VariableNames: index,bscallCh,chIntensity
                    iTileDotTable = TBS.getTileDotTable(iSpotMat,iTileImage);
                    
                    % Delete nonspecific bascall results
                    iTileDotTable = TBS.delLowBscallInten(iTileDotTable);
                    
                    % Get the new table output of current tile
                    iDotTable{iTile} = iTileDotTable;
                end
                
                iDotTable = cell2table(iDotTable,...
                    'RowNames',spotMat.Properties.RowNames,...
                    'VariableNames',{'bscall'});
                dotTable = [dotTable; iDotTable];
                
                disp(strcat('Found Spot: ', iSeqFolder))
            end
        end
        
    end
    
    methods (Static)    % Bleedthrough Correciton =========================
        %% Function:    getBleedThroughcfun
        % Discription: get bleed through estimation, coeffients for bleedthrouh
        % channels
        function bleedThroughEstimation = getBleedThroughcfun(...
                bleedThroughEstimation,bleedThroughEstSetting,bscallTable)
            % Input & output: bleedThroughEstimation, table
            %               bleedThroughEstSetting, struct
            %               bscallTable, table, bscall result, one tile per row
            
            % Brain regions use for bleed through correction
            location = bleedThroughEstSetting.location;
            
            % fittype for calcuating bleedthrough
            fitType = bleedThroughEstSetting.fitType;
            
            % Cycles will be corrected for bleedthough
            correctionForCycle = bleedThroughEstSetting.correctionForCycle;
            
            if ~strcmp(fitType,'poly1')
                error('Currently only support poly1 as fittype.')
            end
            
            for iCh = 1:size(bleedThroughEstimation,1)
                
                % Get of fix and dependend channel
                fixedCh = bleedThroughEstimation.fixedCh{iCh};
                dependedCh = bleedThroughEstimation.dependedCh{iCh};
                
                % Get the seq cycle for bleedthrough estimation
                iSeq = bleedThroughEstimation.seq{iCh};
                iSeq = TBS.seqstr(iSeq);
                
                % Get row names belong to the sequence & location
                rowNames = bscallTable.Properties.RowNames;
                rowNames = rowNames(contains(rowNames,iSeq));
                rowNames = rowNames(contains(rowNames,location));
                
                % Get chIntensity of the bscall channel
                % from teh rows
                chIntensity = bscallTable{rowNames,1};
                row = cellfun(@(X) X.bscallCh == fixedCh,...
                    chIntensity,'Uniformoutput',false);
                chIntensity = cellfun(@(X,Y) X{Y,'chIntensity'},...
                    chIntensity,row,'Uniformoutput',false);
                
                % Get dots w/ intensity beyond noise in the fixed channel
                fh = @(X) X(:,fixedCh) >= bleedThroughEstSetting.noise;
                chIntensity = cellfun(@(X) X(fh(X),:),chIntensity,'Uniformoutput',false);
                
                % Only include tile with enough dots
                nDot = cellfun(@(X) size(X,1),chIntensity);
                row = nDot >= prctile(nDot,bleedThroughEstSetting.numPerTilePrctile);
                chIntensity = chIntensity(row);
                
                coeffval = [];
                % Not that straight forward to use cellfun, so use loop instead
                parfor iTile = 1:size(chIntensity,1) % parfor
                    
                    iChIntensity = chIntensity{iTile};
                    % Need double for fval
                    iChIntensity = double(iChIntensity);
                    
                    fixedChInten = iChIntensity(:,fixedCh);
                    dependedChInten = iChIntensity(:,dependedCh);
                    
                    %[fitObject,outlier] = fitWithoutOutlier(x,y,fitType)
                    [fitOutlierExclude,~] = TBS.fitWithoutOutlier(fixedChInten,dependedChInten,fitType);
                    
                    % Coeffient value of the fit
                    coeffval(iTile,:) = coeffvalues(fitOutlierExclude);
                end
                
                % Take median as coefficients
                coeffval = median(coeffval);
                
                % Assemble cfun
                ffun = fittype(fitType);
                cfun = cfit(ffun,coeffval(1),coeffval(2));
                
                % Get the median of coeff from different images
                bleedThroughEstimation.cfun{iCh}(correctionForCycle,:)=...
                    repmat({cfun},length(correctionForCycle),1);
            end
            
            disp('Done: getBleedThroughcfun.')
        end
        
        %% Function:    chBleedthroughCorrection
        % Description:  Depended channel intnesity correction
        %               Correct the bleed through from channelFix
        function result = chBleedthroughCorrection(result,fixedCh,dependedCh,cfun)
            % Input & output: result,mat
            %               cfun, equantion for bleed through estimation
            
            chFix = result(:,fixedCh);
            chDepend = result(:,dependedCh);
            
            % for the result from feval is > 0 for ch == 0
            zeroRow = chDepend == 0;
            
            chDepend = chDepend - feval(cfun,chFix);
            
            chDepend(zeroRow) = 0;
            chDepend = max(chDepend, 0);
            
            result(:,dependedCh) = chDepend;
        end
        
        %% Function:    bscallBleedthroughCorrection
        % Discription:  correct bscall result bleedthrough in the table
        function dotTable = bscallBleedthroughCorrection(dotTable,iSeq,bleedThroughEstimation)
            % Input & output:    dotTable, table
            %           iSeq, num
            %           bleedThroughEstimation, struct
            
            chIntensity = dotTable.chIntensity;
            chIntensityClass = class(chIntensity);
            % Convert to double for cfun
            chIntensity = double(chIntensity);
            
            for iCh = 1:size(bleedThroughEstimation,1)
                fixedCh = bleedThroughEstimation.fixedCh{iCh};
                dependedCh = bleedThroughEstimation.dependedCh{iCh};
                cfun = bleedThroughEstimation.cfun{iCh}{iSeq};
                
                % chBleedthroughCorrection(result,fixedCh,dependedCh,cfun)
                chIntensity = TBS.chBleedthroughCorrection(chIntensity,fixedCh,dependedCh,cfun);
            end
            
            % Convert back to orginal class
            chIntensity = cast(chIntensity,chIntensityClass);
            dotTable.chIntensity =  chIntensity;
            
            % Delete nonspecific bscalls
            dotTable = TBS.delLowBscallInten(dotTable);
        end
        
        %% Function:    bleedThroughCorrection
        % Discription:  bleedthrough correcting the bscall table
        %               dots with bscallCh not the max intensity will be deleted
        function dotTable = bleedThroughCorrection(dotTable,bleedThroughEstimation,imageSetting)
            % Input:        dotTable, table, one tile per row
            %               bleedThroughEstimation, table, with bleedthrough estimation
            %               cfun
            %               imageSetting, struct
            %               sysSetting, struct
            % Output:       corrected bscallTable
            
            rowNames = dotTable.Properties.RowNames;
            
            for iSeq = imageSetting.seqCycles
                seqName = TBS.seqstr(iSeq);
                
                iRowName = contains(rowNames,seqName);
                iBscallTable = dotTable{iRowName,1};
                
                % bscallBleedthroughCorrection(bscallTable,iSeq,bleedThroughEstimation)
                fh = @(X) TBS.bscallBleedthroughCorrection(X,iSeq,bleedThroughEstimation);
                iBscallTable = cellfun(@(X) fh(X),iBscallTable,'Uniformoutput',false);
                
                dotTable{iRowName,1} = iBscallTable;
            end
            
            disp('Done: bleedThroughCorrection')
        end
        
    end
    
    methods (Static)    % Channel profiling ===============================
        %% Function:    tileChProfilingRef
        % Discription:  select tile for channel profiling
        function [tileName,expectDotNum,expectMaxInten] = ...
                tileChProfilingRef(refSeq,chProfileSetting,bscallTable,imageSetting)
            % Input:        refSeq, num/mat, sequencing cycle number
            %               chProfileSetting, struct
            %               bscallTable, table,
            %               imageSetting, struct
            % Output:       selectiveTile, tile name of selected tiles
            %               expectDotNum, number of dots expected for the
            %               tile
            %               expectMaxInten, expected max intensity
            
            nCh = imageSetting.chNum;
            
            % Min number of dots per tile per channel
            numOfDot = chProfileSetting.numOfDot;
            
            % Location & sequencing cycle of the tile
            location = chProfileSetting.location;
            % Allow multiple reference seq
            refSeq = num2cell(refSeq);
            refSeq = cellfun(@(X) TBS.seqstr(X),refSeq,'Uniformoutput',false);
            
            % Noise for thresholding dots
            noise = chProfileSetting.noise;
            if numel(noise) == 1
                noise = repmat(noise,nCh,1);
            end
            
            % Get tile meet the reqirement
            rowNames = bscallTable.Properties.RowNames;
            rowNames = rowNames(contains(rowNames,refSeq));
            rowNames = rowNames(contains(rowNames,location));
            bscallTable = bscallTable{rowNames,1};
            
            sizeMat = []; expectMaxInten = {};
            for iCh = 1:nCh
                % Bascall channel intensity
                row = cellfun(@(X) X.bscallCh == iCh,bscallTable,'Uniformoutput',false);
                iChStats = cellfun(@(X,Y) X.chIntensity(Y,iCh),bscallTable,row,'Uniformoutput',false);
                
                % Find dots beyond noise, as rolonies
                iChStats = cellfun(@(X) X(X >= noise(iCh)),iChStats,'Uniformoutput',false);
                
                % Pick out outliers (float)
                iChStats = cellfun(@(X) X(~isoutlier(single(X))),iChStats,'Uniformoutput',false);
                
                % Find how many such dots per tile
                sizeMat(:,iCh) = cellfun(@(X) size(X,1),iChStats);
                
                % Output as cell array because there are empty cells
                expectMaxInten(:,iCh) = cellfun(@max,iChStats,'UniformOutput',false);
            end
            
            % Trim reference channel profiling ----------------------------
            % Discard tiles with rolony number < numOfDot from either Ch
            delRow = sizeMat < numOfDot;
            delRow = any(delRow,2);
            
            % Discard tiles with uneven rolony number
            unevenCh = mean(sizeMat,2) - std(sizeMat,0,2).*3;
            unevenCh = unevenCh < 0;
            
            delRow = delRow | unevenCh;
            
            rowNames = rowNames(~delRow);
            % Delete seq info
            for i = 1:numel(refSeq)
                rowNames = strrep(rowNames,refSeq{i},'*');
            end
            tileName = rowNames;
            
            expectDotNum = sizeMat(~delRow,:);
            
            expectMaxInten = expectMaxInten(~delRow,:);
            expectMaxInten = cell2mat(expectMaxInten);
        end
        
        %% Function:    chProfilingRef
        % Discription:  consolidate channel reference value for tiles
        function [selectiveTile,expectDotNum,expectMaxInten] = ...
                chProfilingRef(selectiveTile,expectDotNum,expectMaxInten)
            % Input & output:   seletiveTile, cell, tile names
            %       expectDotNum, mat, expect dot number for each tile
            %       expectMaxInten, mat, expect max intensity
            
            [selectiveTile,~,ic] = unique(selectiveTile);
            
            expectDotNum2 = []; expectMaxInten2 = [];
            for i = 1:max(ic)
                row = ic == i;
                iexpectDotNum = expectDotNum(row,:);
                iexpectMaxInten = expectMaxInten(row,:);
                
                % DotNum: median; maxInten: max
                expectDotNum2(i,:) = round(median(iexpectDotNum,1));
                expectMaxInten2(i,:) = max(iexpectMaxInten,[],1);
            end
            
            expectDotNum = expectDotNum2;
            expectMaxInten = expectMaxInten2;
        end
        
        %% Function:    chProfiling
        % Discription:  get mean and std of channel intensity of location
        %               requirements: discard inten beyond expectedMaxInten
        %               Calculate the mean&std of the expected dot number
        %               w/ top intensity
        function chProfile = chProfiling(iSeq,chProfileSetting,bscallTable,imageSetting)
            % Input:        iSeq, num/mat, sequencing cycle number
            %               chProfileSetting, struct
            %               bscallTable, table,
            %               sysSetting & imageSetting, struct
            % Output:       chProfile, table, variableNames, maxInten &
            %               minIntens one row is one channel
            
            nCh = imageSetting.chNum;
            
            % Check whether iSeq beyond the min refSeq
            iSeqAfterRefSeq = iSeq >  min(chProfileSetting.refSeq);
            
            % iSeqStr
            iSeq = TBS.seqstr(iSeq);
            
            % Get the selective tiles of current sequencing cycle
            selectiveTile = chProfileSetting.selectiveTile;
            rowNames = strrep(selectiveTile,'*',iSeq);
            
            % Get the N number of dots from chProfileSetting
            expectDotNum = chProfileSetting.expectDotNum;
            % Expected max intensity
            expectMaxInten = chProfileSetting.expectMaxInten;
            
            % Check & get row available
            % hasThisRow(tableInput,rowName)
            hasThisRow = TBS.hasThisRow(bscallTable,rowNames);
            
            if ~any(hasThisRow)
                chProfile = table([],[],'VariableNames',{'Var1','Var2'});
                return
            end
            
            rowNames = rowNames(hasThisRow);
            expectDotNum = expectDotNum(hasThisRow,:);
            expectMaxInten = expectMaxInten(hasThisRow,:);
            
            bscallTable = bscallTable{rowNames,1};
            
            chVar1 = []; chVar2 = [];
            for iCh = 1:nCh
                % Bascall channel intensity
                row = cellfun(@(X) X.bscallCh == iCh,bscallTable,'Uniformoutput',false);
                iChStats = cellfun(@(X,Y) X.chIntensity(Y,iCh),bscallTable,row,'Uniformoutput',false);
                
                % Cannot pick out outlier here since there is no noise
                % threholding
                
                % Cap the max inten if iSeq beyond min refSeq -------------
                if iSeqAfterRefSeq && chProfileSetting.enableExpectMaxInten
                    iExpectMaxInten = num2cell(expectMaxInten(:,iCh));
                    % Discard the dots above expected max intensity
                    iChStats = cellfun(@(X,Y) X(X<=Y),iChStats,iExpectMaxInten,'Uniformoutput',false);
                end
                
                % Get the top N number of dots for represent the intensity-
                iChStats = cellfun(@(X) sort(X,'Descend'),iChStats,'Uniformoutput',false);
                
                iExpectDotNum = num2cell(expectDotNum(:,iCh));
                iChStats = cellfun(@(X,Y) X(1:Y),iChStats,iExpectDotNum,'Uniformoutput',false);
                
                iChStats = cellfun(@double,iChStats,'Uniformoutput',false);
                chVar1(:,iCh)= cellfun(@mean,iChStats);
                chVar2(:,iCh) = cellfun(@std,iChStats);
            end
            
            chProfile = table({chVar1},{chVar2},'VariableNames',{'Var1','Var2'});
        end
        
        %% Function:    getChProfileVars
        % Discription:  get the channel profile in mat format
        function chProfileVar = getChProfileVars(chProfileIn)
            % Input:        chProfile, table, with var1/var2
            %               one row is one channel
            % Output:       chProfileVar, cell
            %               one seq cycle per row (excluded cycle is not
            %               included), column channel
            
            chProfileVar = {};
            
            nCh = size(chProfileIn{1},2);
            for iCh = 1:nCh
                chProfileVar(:,iCh) = cellfun(@(X) X(:,iCh),...
                    chProfileIn,'Uniformoutput',false);
            end
        end
        
        %% Function:    getChEstimateCfun
        % Discription:  get the cfun for chScore for all channels and cycles
        function chEstimateCfun = getChEstimateCfun(var1Fit,var2Fit,x)
            % Input:        var1Fit, cell, estimation of the max
            %               var2Fit, cell, estimation of the min
            %               (cell of cfit object, one channel per cell)
            %               x, mat
            % Output:       chEstimateCfun, cell of cfit object for each
            % channel (col) in each sequencing cycles (row)
            
            % chScore equation
            cfun = fittype('(x - mean)/std');
            
            chEstimateCfun = {};
            for iCh = 1:length(var1Fit)
                % Get the estimated mean/std for seq cycles of current Ch
                var1 = feval(var1Fit{iCh},x);
                var2 = feval(var2Fit{iCh},x);
                % Create chScore cfun for each cycle of the current ch
                chEstimateCfun(:,iCh) = arrayfun(@(X,Y) cfit(cfun,X,Y),var1,var2,'Uniformoutput',false);
            end
        end
        
        %% Function:    chProfileFig
        function chProfileFig(chProfileVar1,chProfileVar2,var1Fit,var2Fit,var1Outlier,var2Outlier,chProfileSetting)
            % Input:        chProfileMean, chProfileStd, mat
            %               one column per channel, one row per seq cycle
            %               var1Fit & var2Fit, cell of cfit object, one cfun for
            %               every channel
            %               var1Outlier & var2Outlier, cell of logical
            %               chProfileSetting, struct
            % Output:       figure
            
            x = chProfileSetting.profilingSeq';
            
            figure
            for iCh =  1:size(var1Fit,2)
                
                iColor = TBS.getChColor(iCh);
                
                subplot(1,2,1)  % Plot of var1 ---------------------------
                % chProfileSubplot(scatterX,scatterY,outlier,fitLine,iColor)
                chProfileSubplot(x,chProfileVar1(:,iCh),...
                    var1Outlier{iCh},var1Fit{iCh},iColor);
                title('Mean');
                setFig(chProfileSetting);
                
                subplot(1,2,2)  % Plot of var2 ----------------------------
                chProfileSubplot(x,chProfileVar2(:,iCh),...
                    var2Outlier{iCh},var2Fit{iCh},iColor)
                title('Std');
                setFig(chProfileSetting);
            end
            
            function chProfileSubplot(scatterX,scatterY,outlier,fitLine,iColor)
                
                scatterY = cellfun(@(X,Y) X(~Y),scatterY,outlier,'UniformOutput',false);
                
                % std
                errY = cellfun(@mean,scatterY);
                err = cellfun(@std,scatterY);
                
                % Plot ---------------------------------------------------
                hold on; errorbar(scatterX,errY,err,'.','Color',iColor,'MarkerSize',10);
                
                % Fit curve
                hold on; plot(fitLine,iColor);
            end
            
            function setFig(chProfileSetting)
                xlabel('Seq'); ylabel('Value');legend('off');
                
                g = gca;
                g.XLim = [0,max(chProfileSetting.correctionForCycle)];
                g.YLim = g.YLim + [-2 2];
            end
        end
        
        %% Function:    getChScoreCfun (main)
        % Discription:  get estimate cfun for chScore for every channel
        function [chEstimateCfun,chProfileSetting] = getChScoreCfun(bscallTable,figOutput,chProfileSetting,imageSetting)
            % Input:        bscallTable, table, row name is tile name, with chIntensity
            % of all channels
            %               figOutput, logical, whether have figure output or now
            %               chProfileSetting, struct
            %               sysSetting, struct
            %               imageSetting, struct
            
            profilingSeq = chProfileSetting.profilingSeq;
            
            % Get reference value -----------------------------------------
            % Allow multiple reference seq
            refSeq = chProfileSetting.refSeq;
            if numel(refSeq) > numel(profilingSeq)
                refSeq = refSeq(1:numel(profilingSeq));
            end
            
            % tileChProfilingRef(iSeq,chProfileSetting,bscallTable,imageSetting)
            [selectiveTile,expectDotNum,expectMaxInten] = ...
                TBS.tileChProfilingRef(refSeq,chProfileSetting,bscallTable,imageSetting);
            
            % Consolidate ref ch profilling for each tile
            [selectiveTile,expectDotNum,expectMaxInten] = ...
                TBS.chProfilingRef(selectiveTile,expectDotNum,expectMaxInten);
            
            chProfileSetting.selectiveTile = selectiveTile;
            chProfileSetting.expectDotNum = expectDotNum;
            chProfileSetting.expectMaxInten = expectMaxInten;
            
            % Get channel profile: mean and std ---------------------------
            chProfile = table();
            for iSeq = profilingSeq
                % chProfiling(iSeq,chProfileSetting,bscallTable,sysSetting,imageSetting)
                chProfile(end+1,:) = TBS.chProfiling(iSeq,...
                    chProfileSetting,bscallTable,imageSetting);
            end
            
            % split the table; Var1: mean; Var2: std
            chProfileVar1 = TBS.getChProfileVars(chProfile{:,1});
            chProfileVar2 = TBS.getChProfileVars(chProfile{:,2});
            
            % Get cfun of the mean and std --------------------------------
            chProfileVar1Fit = {}; chProfileVar2Fit = {};
            chProfileVar1Outlier = {}; chProfileVar2Outlier = {};
            for iCh = 1:size(chProfileVar1,2)
                % Expand x for y in cell
                x = TBS.repmat2size(profilingSeq',chProfileVar1(:,iCh),1);
                
                % [fitObject,outlier] = fitWithoutOutlier(x,y,fitType)
                [chProfileVar1Fit{iCh},chProfileVar1Outlier{iCh}] = ...
                    TBS.fitWithoutOutlier(x,chProfileVar1(:,iCh),chProfileSetting.var1fitType);
                [chProfileVar2Fit{iCh},chProfileVar2Outlier{iCh}] = ...
                    TBS.fitWithoutOutlier(x,chProfileVar2(:,iCh),chProfileSetting.var2fitType);
            end
            
            % Output figure if figOutput is true
            if figOutput == true
                TBS.chProfileFig(chProfileVar1,chProfileVar2,chProfileVar1Fit,chProfileVar2Fit,chProfileVar1Outlier,chProfileVar2Outlier,chProfileSetting);
            end
            
            % Estimate the mean and std using cfun, for each channel per
            % cycle cfun of chScore
            % getChEstimateCfun(var1Fit,var2Fit,x)
            chEstimateCfun = TBS.getChEstimateCfun(chProfileVar1Fit,chProfileVar2Fit,chProfileSetting.correctionForCycle);
            
            % Calculate the chScore2Inten to convert chScore to intensity
            % Use median of vars from the 1st cycle
            % getParaInObj(objCell,parameterStr,funStr)
            chProfileSetting.chScore2Inten.mean = TBS.getParaInObj(chEstimateCfun(1,:),'mean','median');
            chProfileSetting.chScore2Inten.std = TBS.getParaInObj(chEstimateCfun(1,:),'std','median');
            
            disp('Done: getchScoreCfun')
        end
        
    end
    
    methods (Static)    % dotTable correction =============================
        %% Function:    chScore2intnesity
        % Discription:  change chScore to intensity
        function intensityMat = chScore2intnesity(chScoreMat,chProfileSetting,outputType)
            % Input:        chScoreMat, mat of chScore
            %               chProfileSetting, struct
            %               outputType, str, class type of the output
            % Output:       intensityMat, mat
            
            meanVal = chProfileSetting.chScore2Inten.mean;
            stdVal = chProfileSetting.chScore2Inten.std;
            
            % chScore to intensity convertion
            intensityMat = chScoreMat.*stdVal + meanVal;
            intensityMat = cast(intensityMat,outputType);
        end
        
        %% Function:    chScoreConvertion
        % Discription:  convert the channel intensity of bscall into chScore
        function bscallTable = chScoreConvertion(bscallTable,chScoreCfun,chProfileSetting)
            % Input:        bscallTable, table
            %               chScoreCfun, cell of cfit object
            %               row- seq, column- channel
            % Output:       bscallTable, intensity changed into chScore, value below
            % mean (chScore < 0) = 0, rows were deleted with the max intensity not in
            % the bscall channel
            
            for iSeq = 1:size(chScoreCfun,1)
                
                iSeqStr = TBS.seqstr(iSeq);
                
                % Get row names belong to the sequence
                rowNames = bscallTable.Properties.RowNames;
                rowNames = rowNames(contains(rowNames,iSeqStr));
                
                % Get max and bscall channel from the rows
                ibscallTable = bscallTable{rowNames,1};
                chIntensity = cellfun(@(X) X.chIntensity,ibscallTable,'Uniformoutput',false);
                chIntensity = cellfun(@single,chIntensity,'Uniformoutput',false);
                
                % Convert the channel intensity into chScore
                chScore = {};
                for iCh = 1:size(chScoreCfun,2)
                    chScore(:,iCh) = cellfun(@(X) feval(chScoreCfun{iSeq,iCh},X(:,iCh)),chIntensity,'Uniformoutput',false);
                end
                
                % Convert chScore to intensity
                for iTile = 1:size(chIntensity,1)
                    iChIntensity = [chScore{iTile,:}];
                                        
                    % chScore2intnesity(chScoreMat,chProfileSetting,outputType)
                    iChIntensity = TBS.chScore2intnesity(iChIntensity,chProfileSetting,'uint16');
                    
                    % zero intensity remains 0
                    zeroInten = chIntensity{iTile,:} ==0;
                    % Bug03242022, fixed and tested
                    iChIntensity(zeroInten) = 0;
                    
                    ibscallTable{iTile}.chIntensity = iChIntensity;
                end
                
                bscallTable{rowNames,1} = ibscallTable;
            end
            
            % Delete nonspecific bscall
            bscallTable{:,1} = cellfun(@(X) TBS.delLowBscallInten(X),...
                bscallTable{:,1},'UniformOutput',false);
            
            disp('Done: chScoreConvertion');
        end
        
        %% Function:    bscallRatioCorrection
        % Discription:  discard the bscalling rows cannot reach the
        % requiremen: 2ndMaxChIntensity/maxChIntensity > chRatioFilter
        function keepRow = bscallRatioCorrection(chIntensity,chRatioFilter)
            % Input:        chIntensity, mat; chRatioFilter, num
            % Output:       keepRow, logical
            
            % Get the max and 2nd max intensity
            chIntensity = sort(chIntensity,2,'descend');
            % Decrease variable size (07182020)
            chIntensity = chIntensity(:,1:2);
            % chRatioFilter usually is double, so convert to float if the
            % input is not
            chIntensity = single(chIntensity);
            
            % Check whether pass the intensity fiter
            keepRow = chIntensity(:,1).*chRatioFilter > chIntensity(:,2);
        end
        
    end
    
    methods (Static)    % Image channel correction ========================
        %% Function:    imageChCorrection
        % Discription:  correct image for bleedthrough and intensity
        % difference
        function imageChCorrection(imName,iSeq,bleedThroughEstimation,...
                chScoreCfun,chProfileSetting,sysSetting)
            % Input:        imageName, str
            %               iSeq, number
            %               bleedThroughEstimation,table, w/ bleedthrough cfun
            %               chScoreCfun, cell of cfit obj
            %               sysSetting, struct
            % Output:       void, image on the disk, w/ chCorrection append
            
            nCh = size(chScoreCfun,2);
            
            stack = TBS.getStack(imName,1:nCh);
            
            % Get class type & size of the image
            classType = class(stack);
            
            % Change to ch column for easy processing
            imSize = size(stack);
            stack = reshape(stack,[],nCh);
            stack = single(stack);
            
            % Bleed through correction ------------------------------------
            cfunCh1Ch2 = bleedThroughEstimation.cfun{'Ch1Ch2'}{iSeq};
            cfunCh3Ch4 = bleedThroughEstimation.cfun{'Ch3Ch4'}{iSeq};
            % chBleedthroughCorrection(result,fixedCh,dependedCh,cfun)
            stack = TBS.chBleedthroughCorrection(stack,...
                bleedThroughEstimation.fixedCh{'Ch1Ch2'},...
                bleedThroughEstimation.dependedCh{'Ch1Ch2'},...
                cfunCh1Ch2);
            stack = TBS.chBleedthroughCorrection(stack,...
                bleedThroughEstimation.fixedCh{'Ch3Ch4'},...
                bleedThroughEstimation.dependedCh{'Ch3Ch4'},...
                cfunCh3Ch4);
            
            % Intensity correction ----------------------------------------
            zeroInten = stack == 0;
            
            % Convert to chScore
            for iCh = 1:nCh
                stack(:,iCh) = feval(chScoreCfun{iSeq,iCh},stack(:,iCh));
            end
            
            % Convert chScore to intensity, also change class to make it smaller
            stack = TBS.chScore2intnesity(stack,chProfileSetting,classType);
            
            % Recover the 0-intensity pixel
            stack(zeroInten) = 0;
            
            % Save the image
            outputName = strrep(imName,sysSetting.localCorrectAppend,sysSetting.chCorrectAppend);
            stack = reshape(stack,imSize);
            TBS.saveStack(stack,outputName);
        end
        
        %% Function:    imageChCorrectionMain (main)
        % Discription:  main method for image channel correction
        % All the correction were done on localCorrected image
        function sysSetting = imageChCorrectionMain(bleedThroughEstimation,...
                chScoreCfun,chProfileSetting,imageSetting,sysSetting,directory)
            % Input:        bleedThroughEstimation, table
            %               chScoreCfun, cell of cfun
            %               chProfileSetting, struct
            %               imageSetting & sysSetting & directory, struct
            % Output:       void, on disk, with chCorrection append
            
            for iSeq = imageSetting.seqCycles
                iSeqFolder = TBS.seqstr(iSeq);
                iSeqFolder = fullfile(directory.main,iSeqFolder);
                cd(iSeqFolder)
                
                imName = ls(strcat('*',sysSetting.localCorrectAppend,sysSetting.imFormat));
                imName = cellstr(imName);
                
                parfor iIm = 1:size(imName,1) % parfor
                    iImName =  imName{iIm};
                    TBS.imageChCorrection(iImName,iSeq,bleedThroughEstimation,...
                        chScoreCfun,chProfileSetting,sysSetting);
                end
                
                disp(strcat('Image channel correction, done: ',iSeqFolder))
            end
            
            sysSetting.tileAppend = sysSetting.chCorrectAppend;
            warning('Function imageChCorrectionMain: sysSetting.tileAppend is changed into chCorrectAppend.')
        end
        
    end
    
    methods (Static)    % Get dot bscall by noise =========================
        %% Function:    getNoise
        % Discription:  Calculate noise for findMaxima
        function noise = getNoise(iSeq,initialNoise,changeNoise)
            % Input:        iSeq, num, sequence number
            %               initialNoise, num/mat, initial noise from setting
            %               changeNoise, num/mat, noise change per cycle, from setting
            % Output:       noise, num/mat, calculated noise for the sequencing cycle
            
            noise = initialNoise+round(changeNoise.*(iSeq-1));
        end
        
        %% Function:    getInitialNoise
        % Discription:  get initial noise for the noise structure basing on
        % the chScore2Inten of the current experiment
        function noiseStructure = getInitialNoise(noiseStructure,imageSetting,directory)
            load(fullfile(directory.main,'chProfileSetting.mat'));
            noiseStructure.initial = round(chProfileSetting.chScore2Inten.mean*noiseStructure.initial.mean ...
                + chProfileSetting.chScore2Inten.std*noiseStructure.initial.std);
            noiseStructure.initial = repmat(noiseStructure.initial,1,imageSetting.chNum);
        end
        
        %% Function:    getDotBeyondNoise
        % Discription:  get dot index and bscallCh by thresholding the
        % intensity of local maxima
        function bscallTable = getDotBeyondNoise(bscallTable,noiseStructure,imageSetting,directory)
            % Input & output:  bscallTable, table,with index, chIntensity, bscallCh
            %           noiseStructure,imageSetting,sysSetting,directory:
            %           struct
            % Note: chIntensity will get eliminated in the output
            
            noiseStructure = TBS.getInitialNoise(noiseStructure,imageSetting,directory);
            
            rowNames = bscallTable.Properties.RowNames;
            
            for iSeq = imageSetting.seqCycles
                seqName = TBS.seqstr(iSeq);
                
                % getNoise(iSeq,initialNoise,changeNoise)
                noise = TBS.getNoise(iSeq,noiseStructure.initial,noiseStructure.change);
                
                % Get rows for the current seq
                iRowName = contains(rowNames,seqName);
                
                if ~any(iRowName)
                    continue
                end
                
                iBscallTable = bscallTable.bscall(iRowName);
                
                % Get beyond noise signal for all channels
                iBscallTable = cellfun(@(X) chBeyondNoise(X,noise),...
                    iBscallTable,'UniformOutput',false);
                
                % delLowBscallInten(dotTable)
                iBscallTable = cellfun(@(X) TBS.delLowBscallInten(X),...
                    iBscallTable,'Uniformoutput',false);
                
                % Get table output
                fh = @(X) X(:,{'index','bscallCh'});
                iBscallTable = cellfun(@(X) fh(X),iBscallTable,...
                    'Uniformoutput',false);
                
                bscallTable.bscall(iRowName) = iBscallTable;
            end
            
            % function: chBeyondNoise -------------------------------------
            % Discription: check whether channel intensity are beyond noise
            function bscallTable = chBeyondNoise(bscallTable,noise)
                bscallTable.chIntensity = bscallTable.chIntensity >= noise;
            end
        end
        
        %% Function:    ind2subInTable
        % Discription:  change index to xy in table table
        % Similar to in2sub
        function bscallTable = ind2subInTable(bscallTable,imSize)
            % Input:        bscallTable, table, one bscall table/index per row
            %               imSize, mat, image size
            % Output:       bscallTable, with index changed to xy coordinates
            
            % Check whether did the ind2sub convertion before
            if ~any(contains(bscallTable{1,1}{:}.Properties.VariableNames,...
                    'index','IgnoreCase',true))
                warning('Function ind2subInTable: variable x & y already existed.')
                return
            end
            
            for iTile = 1:size(bscallTable,1)
                iBscallTable = bscallTable.bscall{iTile};
                index = iBscallTable.index;
                
                % Delete index
                iBscallTable.index = [];
                
                [y,x] = ind2sub(imSize,index);
                iBscallTable.x = x;
                iBscallTable.y = y;
                
                % Update the table
                bscallTable.bscall{iTile} = iBscallTable;
            end
        end
        
        %% Function:    getScaleBscallTable
        % Discription:  scale the bscall table (not the somaBscall table)
        %               and pick the maxima more precisely
        function bscallTableOut = getScaleBscallTable(bscallTable,...
                scaleFactor,SE,imageSetting,sysSetting,directory)
            % Input:        bscallTable, table,
            %               scaleFactor, num, scale up to find subpixel
            %               location
            %               SE, strel object, for finding local maxima
            %               imageSetting & sysSetting & directory, stuct
            % Output:       bscallTableOut, table, bascall result w/
            % subpixel
            
            nCh = imageSetting.chNum;
            
            % tform for the scale factor
            scaleTform = eye(2).*scaleFactor;
            
            % ensureRowNameAppend(tbl,str)
            bscallTable = TBS.ensureRowNameAppend(bscallTable,sysSetting.imFormat);
            rowNames = bscallTable.Properties.RowNames;
            
            bscallTableOut = {};
            parfor iTile = 1:size(rowNames,1) % parfor
                iTileName = rowNames{iTile};
                
                cd(directory.main);
                iDirectory = dir(fullfile('*',iTileName));
                
                cd(iDirectory.folder);
                stack = TBS.getStack(iTileName,1:nCh);
                
                % Current bscall result
                bscall = bscallTable.bscall{iTileName};
                
                sz = size(stack);
                
                % Get local max area
                % Convert uint8 to float
                bscall.bscallCh = single(bscall.bscallCh);
                pickedDot = [bscall.x, bscall.y, bscall.bscallCh];
                pickedDot = TBS.xyzv2im(sz,pickedDot,[]);
                
                % Resize ROI and image to subpixel
                pickedDot = imresize(pickedDot,scaleFactor);
                stack = imresize(stack,scaleFactor,'bilinear');
                
                % Only pick the local maximum in selected region
                stack(~pickedDot) = 0;
                stack = imdilate(stack,SE) == stack & pickedDot;
                
                % Find xy coordinates, and channel (z)
                [y,x,z,~] = TBS.find3(stack);
                xyz = [x y z];
                
                % Scale back in xy axis
                xyz(:,1:2) = xyz(:,1:2)*inv(scaleTform);
                xyz = single(xyz);
                
                iBscallTableOut = array2table(xyz,'VariableNames',...
                    {'x','y','bscallCh'});
                iBscallTableOut.bscallCh = uint8(iBscallTableOut.bscallCh);
                bscallTableOut{iTile,1} = iBscallTableOut;
                
                disp(['Scale bscallTable, Done: ',iTileName]);
            end
            
            bscallTableOut = table(bscallTableOut,'RowNames',...
                rowNames,'VariableNames',{'bscall'});
        end
        
        %% Function:    uniqueDotPerPixel
        % Discription:  trim extra dot and get unique dot per pixel
        function tbl = uniqueDotPerPixel(tbl)
            % (will slightly bias to smaller x)
            % Input & output:    tbl, table, with xy coordinates and bscall result
            
            [~,I] = sort(tbl.x,'ascend');
            tbl = tbl(I,:);
            
            xy = [tbl.x, tbl.y];
            xy = round(xy);
            [~,ia,~] = unique(xy,'rows');
            
            tbl = tbl(ia,:);
        end
        
    end
        
    methods (Static)    % Fuse image alignment ============================
        %% Function:    estimateTranslation
        % Discription:  estimate translation matrix using tile position
        function tformCell = estimateTranslation(tilePos,imageSetting)
            % Input:    tilePos, mat, tile position
            %           imageSettng,struct
            % Output:   tformCell, cell, one translation matrix per cell
            %           sorted using tile number
            
            if isempty(tilePos)
                tformCell = [];
                return
            end
            
            tileSize = imageSetting.tileSize;
            overlapPrecentage = imageSetting.overlapPrecentage;
            
            overlapPixel = tileSize.*overlapPrecentage.*0.01;
            
            [row,col,v] = find(tilePos);
            
            y = row*tileSize(1) - (row-1)*overlapPixel(1) - tileSize(1);
            x = col*tileSize(2) - (col-1)*overlapPixel(2) - tileSize(2);
            
            tformCell = arrayfun(@(X,Y) [1 0 0; 0 1 0; X Y 1],x,y,'UniformOutput',false);
            
            [~,I] = sort(v,'ascend');
            tformCell = tformCell(I);
        end
        
        %% Function:    alignFuseImage
        % Discription:  aligned fuse image for rough alignment
        function fuseAlignTform = alignFuseImage(redoTF,nTime,maxInten,...
                imageSetting,sysSetting,directory)
            % Input:    imageSetting/sysSetting/directory, struct
            %           redoTF, logical, whether redo
            %           maxInten, max intensity for fused image for
            %           alignment
            %           nTime, num, time for repeat alignment
            % Output:   fuseAlignTform, table, with alignment
            % transformation matrix
            
            % Output file name
            fileName = 'fuseAlignTform.mat';
            
            directory = directory.temporary;
            cd(directory)
            % Load the exists variable
            if exist(fileName) && ~redoTF
                load(fileName);
            else
                fuseAlignTform = table();
            end
            
            % Image name under the folder
            cd(directory)
            imName = ls(['*',sysSetting.imFormat]);
            imName = cellstr(imName);
            
            % Alignment sequence
            % (here use initial)
            seqCycles = imageSetting.seqCycles;
            alignmentSeq = TBS.getAlignmentSeq(seqCycles,1);
            
            % (Use moving seq to find fix seq)
            for i = 2:numel(alignmentSeq)
                movingSeq = alignmentSeq(i);
                movingSeq = TBS.seqstr(movingSeq);
                
                % Images belongs to current moving seq
                iImName = imName(contains(imName,movingSeq));
                
                tblOut = {};
                parfor iIm = 1:size(iImName,1) % parfor
                    movingName = iImName{iIm};
                    
                    if ~redoTF && ~isempty(fuseAlignTform) && ...
                            any(contains(fuseAlignTform.movingName,movingName))
                        continue
                    end
                    
                    % Find fix forward
                    jTblOut = {};
                    for j = fliplr(1:i-1)
                        
                        % Allow multiple alignment for the same image
                        if ~isempty(jTblOut) && j < i - nTime
                            break
                        end
                        
                        fixSeq = alignmentSeq(j);
                        fixSeq = TBS.seqstr(fixSeq);
                        
                        fixName = strrep(movingName,movingSeq,fixSeq);
                        
                        % Check whether the image is available
                        if ~contains(fixName,imName)
                            continue
                        end
                        
                        cd(directory);
                        fixIm = TBS.getStack(fixName,[]);
                        movingIm = TBS.getStack(movingName,[]);
                        
                        % 06072021, set max intensity
                        % (to decrease impact of strong nonspecific signal)
                        if ~isempty(maxInten)
                            fixIm = min(fixIm,maxInten);
                            movingIm = min(movingIm,maxInten);
                        end
                        
                        % Alignment
                        tform = imregcorr(movingIm,fixIm,'Window',false);
                        tform = tform.T;
                        
                        % QC: a very small scale factor
                        if ~TBS.QCtform(tform)
                            warning(['Did not matched: ',fixName])
                            continue
                        end
                        
                        jTblOut = [jTblOut; {movingName},{fixName},{tform}];
                    end
                    
                    tblOut{iIm,1} = jTblOut;
                    
                    disp(['Aligned fuse image: ',movingName]);
                end
                
                if isempty(tblOut)
                    continue
                end
                
                tblOut = vertcat(tblOut{:});
                
                varNames = {'movingName','fixName','tform'};
                
                if isempty(fuseAlignTform) % when i == 2
                    n = size(tblOut,1);
                    fuseAlignTform = cell2table([tblOut(:,2),...
                        cell(n,1),repmat({eye(3)},n,1)],...
                        'VariableNames',varNames);
                end
                
                fuseAlignTform = [fuseAlignTform; cell2table(...
                    tblOut,'VariableNames',varNames)];
                
                save(fullfile(directory,fileName),'fuseAlignTform');
            end
        end
        
        %% Function:    addImTform2TileTform
        % Discription:  add image tform to tile tform
        function tileTform = addImTform2TileTform(imTform,tileTform,sysSetting)
            % Input & output:   imTform, table, tform for image
            %                   tileTform, table, tform for tile
            
            % Get the image name from row name
            % Bug fix 04092021_erase tif
            imName = imTform.Properties.RowNames;
            imName = erase(imName,sysSetting.imFormat);
            imName = cellfun(@(X) [X sysSetting.delimiter],imName,...
                'UniformOutput',false);
            
            tileName = tileTform.Properties.RowNames;
            for i = 1:size(imName,1)
                iImName = imName{i};
                iImTform = imTform.tform{i};
                
                row = contains(tileName,iImName);
                
                iTileTform = tileTform.tform(row);
                
                iTileTform = cellfun(@(X) X*iImTform,iTileTform,...
                    'UniformOutput',false);
                
                tileTform.tform(row) = iTileTform;
            end
        end
        
        %% Function:    getAlignedImPair
        % Discription:  transformed images 1 & 2 for visualize the
        % aligmneent
        function [tfIm1,tfIm2] = getAlignedImPair(imName1,imName2,tformTable,directory)
            % Input:    imName1/imName2, str, image names
            %           tformTable, table, row: image names, tform: transformation
            %           matrix
            %           directory, str, directory of the images
            % Output:   tfIm1/tfIm2, transformed image, size as image2
            
            cd(directory);
            im1 = TBS.getStack(imName1,[]);
            im2 = TBS.getStack(imName2,[]);
            
            % Do max projection if it is stack
            im1 = max(im1,[],3);
            im2 = max(im2,[],3);
            
            % Use image2 size as output
            R = size(im2);
            R = imref2d(R);
            
            % Transformation
            % Note: Y is mat
            fh = @(X,Y) imwarp(X,affine2d(Y),'OutputView',R);
            tform = tformTable.tform{imName1};
            tfIm1 = fh(im1,tform);
            tform = tformTable.tform{imName2};
            tfIm2 = fh(im2,tform);
        end
        
    end
        
    methods (Static)    % Dot alignment (DA) ==============================
        %% Function:    DA_getDotPair
        % Discription:  get index of moving-fix pair, using translaiton
        function [movingI,fixI] = DA_getDotPair(moving,fix,dim,binSz,dotLim,maxDist)
            % Note: error range is presetted
            %       3rd dimension for hamming distance only
            % Input:    moving, mat, coordinates of moving dots
            %           fix, mat, coordinates of fixed dots
            %           dim, num, dimension of moving & fix
            %           binSz, num, binning size, for finding transloation
            %           dotLim, num, max dots for alignment
            %           maxDist, num, max distance for translocation
            % Output:   movingI/fixI, index of moving and fix dot pairs
            
            % Trim out extra dots to speed up
            [movingI,fixI] = TBS.DA_withinRng(moving,fix,maxDist);
            
            if size(fixI,1) < 3 || size(movingI,1) < 3
                movingI = []; fixI = []; return
            end
            
            % Trim out crowded area if there are too many dots ============
            % to decrease the effect of croweded area, also speed up
            if ~isempty(dotLim)
                fix2 = fix(fixI,1:dim);
                moving2 = moving(movingI,1:dim);
                
                if dim == 2
                    fixI2 = TBS.DA_trimXY(fix2,dotLim);
                    movingI2 = TBS.DA_trimXY(moving2,dotLim);
                elseif dim == 3
                    fixI2 = TBS.DA_trimXY3(fix2,dotLim);
                    movingI2 = TBS.DA_trimXY3(moving2,dotLim);
                end
                
                fixI = fixI(fixI2);
                movingI = movingI(movingI2);
            end
            
            fix2 = fix(fixI,1:dim);
            moving2 = moving(movingI,1:dim);
            
            % Difference along axes =======================================
            diffX = moving2(:,1)-fix2(:,1)';
            diffY = moving2(:,2)-fix2(:,2)';
            
            % (dim3) only compare similarity
            if dim == 3
                diffZ = moving2(:,3) == fix2(:,3)';
            end
            
            % Find most frequence distance in x & y-axes ------------------
            
            % Exclude distance beyond maxDist
            I = abs(diffX) <= maxDist &  abs(diffY) <= maxDist;
            
            % (dim3) exclud unidentical 3rd dimension
            if dim == 3
                I = I & diffZ;
            end
            
            if ~any(I,'all')
                movingI = []; fixI = []; return;
            end
            
            diffXY2 = [diffX(I), diffY(I)];
            
            % Binning (currently only support int16)
            diffXY2 = int16(diffXY2./binSz).*(binSz);
            
            % Most frequent diff
            [C,~,ic] = unique(diffXY2,'rows');
            
            [M,F] = mode(ic,1);
            if F <= 1
                movingI = []; fixI = []; return
            end
            
            diffXY2 = C(M,:);
            
            % Find dot pairs ----------------------------------------------
            % Expand the range tolaration for rotaiton adjustment
            errRng = binSz*2;
            cmpFh = @(X,Y) X >= (Y-errRng) & X <= (Y+errRng);
            
            % Find the pair of moving and fix within the range
            I = cmpFh(diffX,diffXY2(1)) & cmpFh(diffY,diffXY2(2));
            
            % (dim3) exclud unidentical 3rd dimension
            if dim == 3
                I = I & diffZ;
            end
            
            % Find the index for each pair in moving and fix
            [movingI2, fixI2] = find(I);
            
            fixI = fixI(fixI2);
            movingI = movingI(movingI2);
            
            % Check whether the dots can be used for affine3d -------------
            fix2 = fix(fixI,1:2);
            moving2 = moving(movingI,1:2);
            
            try
                tform = fitgeotrans(moving2,fix2,'affine');
            catch
                fixI = []; movingI = [];
            end
            
            % (CheckPoint) Look at selective dot pairs
            % TBS.checkDotAlign(eye(3),moving2,fix2);
        end
        
        %% Funciton:    DA_pair2Tform
        % Discription:  get transformation matrix using fitgeotrans
        function tform = DA_pair2Tform(moving2,fix2)
            % Input:    moving, mat, coordinates of moving dots
            %           fix, mat, coordinates of fixed dots
            % Output:   tform, 3 x 3 transformation matrix
            
            if isempty(moving2) || isempty(fix2)
                tform = []; return
            end
            
            moving2 = moving2(:,1:2);
            fix2 = fix2(:,1:2);
            
            % affine & projective -----------------------------------------
            % 06242021 alignmentTest 26
            tform = [];
            
            % min 10 dor for projective (can change)
            fh = @(X,Y) size(unique(round(X),'rows'),1) > Y;
            
            if fh(moving2,10) && fh(fix2,10)
                try
                    tform = fitgeotrans(moving2,fix2,'projective');
                end
            end
            
            % QC: if it doesnt pass QC, use affine
            if ~isempty(tform) && ~TBS.QCtform(tform.T)
                tform = [];
            else
                return
            end
            
            try
                tform = fitgeotrans(moving2,fix2,'affine');
            end
            
            % QC: if it doesnt pass QC, use affine
            if ~isempty(tform) && ~TBS.QCtform(tform.T)
                tform = [];
            end
        end
        
        %% Function:    DA_withinRng
        % Discription:  Trim out far away dots to speed up
        function [movingI,fixI] = DA_withinRng(moving,fix,maxDist)
            % Input:    moving/fix, mat, coordinates
            %           maxDist, num, tolerance
            % Output:   movingI/fixI, vector, index of dots within the
            % range
            
            % Find the lim of coordinates comparison
            minFh = @(X) min(X(:,1:2),[],1);
            maxFh = @(X) max(X(:,1:2),[],1);
            
            minLim = max(minFh(moving),minFh(fix));
            maxLim = min(maxFh(moving),maxFh(fix));
            
            % Add tolerance
            minLim = minLim - maxDist;
            maxLim = maxLim + maxDist;
            
            % Select rows within the limits
            rowFh = @(X) find(all(X(:,1:2)>= minLim & X(:,1:2)<= maxLim,2));
            fixI = rowFh(fix); movingI = rowFh(moving);
        end
        
        %% Function:    DA_trimXY
        % Discription: pick the N dots with longest nearest distance if
        % there are too many dots
        function I = DA_trimXY(xy,dotLim)
            % Input:    xy, mat, dot coordinates
            %           dotLim, num, max amount of dots
            % Output:   I, vector, index
            
            n = size(xy,1);
            
            if n <= dotLim
                I = 1:n; return
            end
            
            % Distance to the nearest neighbor
            dist = pdist2(xy(:,1:2),xy(:,1:2),'euclidean','Smallest',2);
            dist = dist(2,:);
            dist = dist';
            [~,I] = sort(dist,'descend');
            I = I(1:dotLim);
        end
        
        %% Function:    DA_trimXY3
        % Discription: trimXY in 3 dimension input (channel/seq)
        function I = DA_trimXY3(xy,dotLim)
            % Input:    xy, mat, dot coordinates
            %           dotLim, num, max amount of dots
            % Output:   I, vector, dot index
            
            [C,~,ic] = unique(xy(:,3));
            
            I = {};
            for i = 1:size(C,1)
                row = ic == i;
                row = find(row);
                
                ixy = xy(row,:);
                iI = TBS.DA_trimXY(ixy,dotLim);
                I{i,1} = row(iI);
            end
            
            I = vertcat(I{:});
        end
        
        %% Function:    alignDot
        % Discription:  alignDot basing on XY (mainly translocation)
        function tform = alignDot(moving,fix,dim,binSz,dotLim,maxDist)
            % (3rd dimension for hamming distance only)
            % Input:    moving, mat, coordinates of moving dots
            %           fix, mat, coordinates of fixed dots
            %           dim, num, dimension of moving & fix
            %           binSz, num, binning size, for finding transloation
            %           dotLim, num, max dots for alignment
            %           maxDist, num, max distance for translocation
            % Output:   tform, 3 x 3 transformation matrix
            
            % Index of moving-fix pair
            [movingI,fixI] = TBS.DA_getDotPair(moving,fix,dim,binSz,dotLim,maxDist);
            
            % Get tform (presepective/affine) -----------------------------
            moving2 = moving(movingI,:);
            fix2 = fix(fixI,:);
            
            tform = TBS.DA_pair2Tform(moving2,fix2);
        end
        
        %% Function:    DA_splitTile
        % Discription:  split dots into n tile/segment
        function cellOut = DA_splitTile(matIn,lim,nTile)
            % Note, no overlap between tiles
            % Input:    matInt, mat, dot coordinates
            %           lim, mat, 1st row, min; 2nd row, max; 1st col, x; 2nd col, y
            %           n, mat, tile/seciton to split along x & y
            % Output:   cellOut, cell, dot index in each segment
            
            if numel(nTile) == 1
                nTile = repmat(nTile,1,2);
            end
            
            % Intervel for x & y per tile
            interval = diff(lim)./nTile;
            
            % x & y lim for each tile
            xLim = arrayfun(@(X) lim(1,1)+[(X-1) X].*interval(1),1:nTile(1),...
                'UniformOutput',false);
            yLim = arrayfun(@(X) lim(1,2)+[(X-1) X].*interval(2),(1:nTile(2))',...
                'UniformOutput',false);
            xLim = repmat(xLim,nTile(2),1);
            yLim = repmat(yLim,1,nTile(1));
            
            % Get dot index within the segment
            row = cellfun(@(X,Y) matIn(:,1) >= X(1) & matIn(:,1)<X(2) &...
                matIn(:,2) >= Y(1) & matIn(:,2)<Y(2),xLim,yLim,'UniformOutput',false);
            cellOut = cellfun(@find,row,'UniformOutput',false);
        end
        
        %% Function:    DA_alignDotTiling
        % Discription:  alignDot with split them into tile
        function tform = DA_alignDotTiling(moving,fix,tform,nTile,...
                dim,binSize,dotLim,maxDist)
            % Note: Prealignment is necessary, due to the tiling
            % Input:    moving/fix, mat, coordinates of moving/fix dots
            %           tform, affine2d/projective2d, transformation mat
            %           dim, num, dimension of moving & fix
            %           nTile, mat/num, number of tile for splitting dots
            %           binSz, num, binning size, for finding transloation
            %           dotLim, num, max dots for alignment
            %           maxDist, num, max distance for translocation
            % Output:   tform, 3 x 3 transformation matrix
            
            tfMoving = transformPointsForward(tform,moving(:,1:2));
            if size(moving,2) == 3
                tfMoving(:,3) = moving(:,3);
            end
            
            % Split dots into tiles ---------------------------------------------------
            % Limition (min & max of the field of view)
            lim = [min(tfMoving,[],1); max(tfMoving,[],1)];
            lim = lim(:,1:2);
            lim(1,:) = lim(1,:)- maxDist;
            lim(1,:) = lim(1,:)+ maxDist;
            % Index and coordinates for each tile
            moving2I = TBS.DA_splitTile(tfMoving,lim,nTile);
            fix2I = TBS.DA_splitTile(fix,lim,nTile);
            
            % Delete tiles has too few dots
            fh = @(X) size(X,1)<10;
            row = cellfun(@(X,Y) fh(X) | fh(Y), moving2I,fix2I);
            row = ~row;
            if sum(row,'all') <= numel(row)/4
                tform = false; return
            end
            moving2I = moving2I(row);
            fix2I = fix2I(row);
            
            moving2 = cellfun(@(X) tfMoving(X,:),moving2I,'UniformOutput',false);
            fix2 = cellfun(@(X) fix(X,:),fix2I,'UniformOutput',false);
            
            % Alignment of tiles (for rotation) ---------------------------------------
            % TBS.DA_getDotPair(moving,fix,dim,binSz,dotLim,maxDist)
            [movingI,fixI] = cellfun(@(X,Y) TBS.DA_getDotPair(X,Y,dim,binSize,dotLim,maxDist),...
                moving2,fix2,'UniformOutput',false);
            
            % Get original index
            movingI = cellfun(@(X,Y) X(Y),moving2I,movingI,'Uniformoutput',false);
            fixI = cellfun(@(X,Y) X(Y),fix2I,fixI,'Uniformoutput',false);
            
            movingI = vertcat(movingI{:});
            fixI = vertcat(fixI{:});
            
            % Get tform using original coordinates
            moving2 = moving(movingI,:);
            fix2 = fix(fixI,:);
            
            % 25e: delete outlier pairs, from outlier tiles
            D = sqrt(sum((moving2-fix2).^2,2));
            D = zscore(D);
            row = D <= 5;
            moving2 = moving2(row,:);
            fix2 = fix2(row,:);
            
            tform = TBS.DA_pair2Tform(moving2,fix2);
        end
        
        %% Function:    alignDotRot
        % Discription:  align dots for rotaiton
        function tform = alignDotRot(moving,fix,tform,nTile,maxNtile,...
                maxIteration,dim,binSize,dotLim,maxDist)
            % Input:    moving/fix, mat, coordinates of moving/fix dots
            %           tform, affine2d/projective2d, transformation mat
            %           dim, num, dimension of moving & fix
            %           nTile, mat/num, number of tile for splitting dots
            %           maxIteration, num, max iteration for looping
            %           binSz, num, binning size, for finding transloation
            %           dotLim, num, max dots for alignment
            %           maxDist, num, max distance for translocation
            % Output:   tform, 3 x 3 transformation matrix
            
            if ~isobject(tform)
                tform = projective2d(tform);
            end
            
            nTile0 = nTile;
            
            for i = 1:maxIteration
                previousTform = tform.T;
                
                % Alignment with tiling -----------------------------------
                tform2 = TBS.DA_alignDotTiling(moving,fix,tform,nTile,...
                    dim,binSize,dotLim,maxDist);
                
                % Temporary change the tiling setting while no tform has been found
                % (!! it's still possible it will run into error)
                if isempty(tform2) && nTile < maxNtile
                    nTile = nTile+1;
                    disp(['Function alignDotRot nTile+1: ',num2str(nTile)]);
                    continue
                    % Beyond maxNtile
                elseif isempty(tform2) || (islogical(tform2) && ~tform2)
                    tform = [];
                    disp('Function alignDotRot not aligned.');
                    return
                end
                
                nTile = nTile0; % go back to original setting
                tform = tform2;
                
                % Compare tform difference --------------------------------
                % Tolerance, 'rotation' 0.001; translation: 1
                diffTform = abs(tform.T-previousTform);
                diffTform = diffTform(:,1:2);
                diffTform = any(diffTform(1:2,:)>= 10.^-3,'all') |...
                    any(diffTform(3,:)>= 1);
                
                if ~diffTform
                    break
                end
            end
        end
        
    end
    
    methods (Static)    % Rolony alignment ================================
        %% Function:    findClosestFix
        % Discription:  find the closest fix tile basing on the
        % pre-alignment transformation matrix
        function [fixName,fixTform] = findClosestFix(tbl,movingTform)
            % Input:    tbl, table, candidate fix tiles with tform and name
            %           movingTform, mat, moving tile prealigned tform
            % Output:   fixName, str
            %           fixTform, mat, fix tile prealigned tform
            
            fixTform = tbl.tform;
            
            % Change this to double to avoid warning
            classFh = class(fixTform{1});
            classFh = @(X) cast(X,classFh);
            movingTform = classFh(movingTform);
            
            % Find transformation matrix with the shortest distance
            dist = cellfun(@(X) pdist2(X(3,1:2),movingTform(3,1:2)),...
                fixTform);
            [~,I] = min(dist);
            
            fixName = tbl.Properties.RowNames{I};
            
            % Fix tform in fuse image
            fixTform = tbl.tform{fixName};
        end
        
        %% Function:    checkDotAlign
        % Discription:  figure output to check the alignment result
        function checkDotAlign(tform,moving,fix)
            % Input:    tform, transformation matrix
            %           moving/fix, mat, coordinates for fix and moving dot
            if numel(tform) == 9
                tform = projective2d(tform);
            end
            
            moving = moving(:,1:2); fix = fix(:,1:2);
            
            tfMoving = transformPointsForward(tform,moving);
            
            hold off; scatter(fix(:,1),fix(:,2),5,'filled','MarkerFaceAlpha',0.6);
            hold on; scatter(tfMoving(:,1),tfMoving(:,2),5,'filled','MarkerFaceAlpha',0.5);
            
            xlabel('X (pixel)'); ylabel('Y (pixel)');
            g = gca; g.YDir = 'reverse';
            daspect([1 1 1]); pbaspect([1 1 1])
        end
        
        %% Function:    getTformTable (main)
        % Discription:  align tiles basing on strong dots
        % Note: prealignment is needed!!
        function tformTable = getTformTable(bscallTable,fuseStitchTable,...
                redoTF,imageSetting,sysSetting,directory)
            % Input:    dotTable,table, coordinates of dots
            %           fuseStitchTable,table, prealignment tform
            %           imageSetting,sysSettign, sturct
            %           redoTF, logical, whether redo alignment
            % Output:   tformTable, tform for each tile
            
            % Alignment sequence
            % (Allow Align from the middle of sequencing cycle)
            seqCycles = imageSetting.seqCycles;
            initialFixSeq = imageSetting.initialFixSeq;
            alignmentSeq = TBS.getAlignmentSeq(seqCycles,initialFixSeq);
            
            tileName = bscallTable.Properties.RowNames;
            nameElement = 1:3;
            
            cd(directory.main)
            if ~redoTF && exist('tformTable.mat')
                tformTable = load('tformTable.mat').tformTable;
            else
                tformTable = table();
            end
            
            for i = 1:numel(alignmentSeq)
                movingSeq = alignmentSeq(i);
                movingSeq = TBS.seqstr(movingSeq);
                
                % Images belongs to current moving seq
                iImName = tileName(contains(tileName,movingSeq));
                
                tblOut = {};
                % Its faster Not use parfor, delay from some tiles
                for iIm = 1:size(iImName,1) % parfor
                    movingName = iImName{iIm};
                    
                    % Skip the tile didn't stitched
                    if ~TBS.hasThisRow(fuseStitchTable,movingName)
                        continue
                        % Skip the tile already aligned
                    elseif ~isempty(tformTable) && ...
                            any(contains(tformTable.movingName,movingName))
                        continue
                    end
                    
                    movingTform = fuseStitchTable.tform{movingName};
                    
                    if i == 1
                        tblOut(iIm,:) = [{movingName},{[]},{movingTform}];
                        continue
                    end
                    
                    % Moving dot coordinates
                    moving = bscallTable.bscall{movingName};
                    
                    if isempty(moving)
                        continue
                    end
                    
                    moving = [moving.x,moving.y];
                    moving = transformPointsForward(affine2d(movingTform),moving);
                    
                    % alignment sequence for finding fix seq, 06132021
                    % Note03252022, this is not optimal after getAlignmentSeq
                    % but it works
                    if alignmentSeq(i) >  alignmentSeq(i-1)
                        jAlignmentSeq = seqCycles(1:alignmentSeq(i-1));
                        jAlignmentSeq = fliplr(jAlignmentSeq);
                    else
                        jAlignmentSeq = seqCycles(alignmentSeq(i-1):end);
                    end
                    
                    % Test 25e, first align to the inital fix
                    if jAlignmentSeq(1)~= alignmentSeq(1)
                        jAlignmentSeq = [alignmentSeq(1),jAlignmentSeq];
                    end
                    
                    % Find fix and align
                    for j = jAlignmentSeq
                        
                        fixSeq = TBS.seqstr(j);
                        
                        fixName = strrep(movingName,movingSeq,fixSeq);
                        fixName = TBS.nameFun(fixName,nameElement,sysSetting);
                        % Add delimiter to differicinate Vis from VisC etc
                        fixName = [fixName,sysSetting.delimiter];
                        
                        % Potential fix tiles
                        row = contains(fuseStitchTable.Properties.RowNames,fixName);
                        
                        if ~any(row)
                            continue
                        end
                        
                        % Find the fix tile with the closest distance to the moving
                        [fixName,fixTform] = TBS.findClosestFix(...
                            fuseStitchTable(row,:),movingTform);
                        
                        if ~TBS.hasThisRow(bscallTable,fixName)
                            continue
                        end
                        
                        % Fix dot coordinates
                        fix = bscallTable.bscall{fixName};
                        
                        if isempty(fix)
                            continue
                        end
                        
                        fix = [fix.x,fix.y];
                        fix = transformPointsForward(affine2d(fixTform),fix);
                        
                        fix = single(fix); moving = single(moving);
                        
                        % alignDot(moving,fix,dim,binSz,dotLim,maxDist)
                        % Parameter need to be detertermin mannually
                        % tform: projective/affine object
                        tform = TBS.alignDot(moving,fix,2,5,4000,200);
                        
                        if isempty(tform)
                            continue
                        end
                        
                        % alignDotRot(moving,fix,tform,nTile,maxNtile,...
                        % maxIteration,dim,binSize,dotLim,maxDist)
                        tform = TBS.alignDotRot(moving,fix,tform,2,4,100,2,5,1000,200);
                        
                        if isempty(tform)
                            continue
                        end
                        
                        tform = tform.T;
                        
                        % (Check point) For manual testing alignment ---------
                        % TBS.checkDotAlign(tform,moving,fix);
                        
                        tform = movingTform*tform*inv(fixTform);
                        
                        tblOut(iIm,:)= [{movingName},{fixName},{tform}];
                        
                        disp(['Aligned tile: ',movingName]);
                        break
                    end
                end
                
                if isempty(tblOut)
                    continue
                end
                
                % Delete empty rows
                row = cellfun(@isempty,tblOut(:,1));
                tblOut = tblOut(~row,:);
                
                varNames = {'movingName','fixName','tform'};
                
                tformTable = [tformTable; cell2table(tblOut,'VariableNames',varNames)];
            end
        end
        
    end
    
    methods (Static)    % Tile transform ==================================
        %% Function:    intersectTable
        % Discription:  get sorted table with intersect rows
        function [tblA, tblB] = intersectTable(tblA,tblB)
            % Input & output:    tblA/B, table
            
            rowNameA = tblA.Properties.RowNames;
            rowNameB = tblB.Properties.RowNames;
            
            [~,ia,ib] = intersect(rowNameA,rowNameB);
            
            tblA = tblA(ia,:);
            tblB = tblB(ib,:);
        end
        
        %% Function:    combTformInTable
        % Discription:  combine tform from table A & B
        function tblA = combTformInTable(tblA,tblB,fh)
            % Input:    tblA/B, intput A & B
            %           fh, function handle
            % Output:   tblA, with combined tform frm tblA & B
            
            [tblA, tblB] = TBS.intersectTable(tblA,tblB);
            
            tblA.tform = cellfun(@(A,B) fh(A,B),tblA.tform,...
                tblB.tform,'Uniformoutput',false);
        end
        
        %% Function:    tranformXYinTable
        % Discription:  transform coordinates in table
        function xyTable = tranformXYinTable(xyTable,tformTable)
            % Input:    xyTable, table, with xy coordinates
            %           tformTable, table, transformaiton matrix
            % Output:   xyTable, table, with transformed coordinates
            
            [xyTable,tformTable] = TBS.intersectTable(xyTable,tformTable);
            
            for i = 1:size(xyTable,1)
                
                tform = tformTable.tform{i};
                tform = projective2d(tform);
                
                iTable = xyTable{i,1}{:};
                xy = [iTable.x,iTable.y];
                
                xy = transformPointsForward(tform,xy);
                
                % Replace the transformed XY
                iTable.x = xy(:,1);
                iTable.y = xy(:,2);
                
                xyTable{i,1} = {iTable};
            end
            
            disp('Done: tranformXYinTable.');
        end
        
        %% Function:    squeezeTile
        % Discription:  Convert seq cycle into seq for each image
        function tblOut = squeezeTile(tfBscallTable,sysSetting)
            % Input:    tfBscallTable, table
            %           sysSetting, struct
            % Output:   tblOut, table
            
            % Image the tile belongs to
            imName = TBS.rowName2ImName(tfBscallTable,true,sysSetting);
            
            % All available image
            imNameC = unique(imName);
            
            tblOut = {};
            for iIm = 1:size(imNameC,1)
                iImName = imNameC{iIm};
                
                row = contains(imName,iImName);
                tileName = tfBscallTable.Properties.RowNames;
                tileName = tileName(row);
                
                iTblOut = tfBscallTable.bscall(tileName);
                
                % Sequencing cycle of the tile
                seq = cellfun(@(X) TBS.getSeqNum(X,sysSetting),...
                    tileName);
                seq = TBS.repmat2size(seq,iTblOut,1);
                
                iTblOut = vertcat(iTblOut{:});
                seq = vertcat(seq{:});
                
                iTblOut.seq = uint8(seq);
                
                tblOut(iIm,:)=[{iImName},{iTblOut}];
            end
            
            % Convert to image name
            tblOut(:,1) = cellfun(@(X) TBS.imFormat(X,sysSetting),...
                tblOut(:,1),'UniformOutput',false);
            
            tblOut = table(tblOut(:,2),'VariableNames',{'bscall'},...
                'RowNames',tblOut(:,1));
            
            disp('Done: squeezeTile.')
        end
        
        %% Function:    delDotPixelInSoma
        % Discription:  delete dot pixel with soma basecall
        function tfBscallTable = delDotPixelInSoma(tfBscallTable,somaBscallTable)
            % Input & output: tfBscallTable, table, with xy coordiantes of
            % dot bscall
            %           somaBscallTable, table, with xy coordinates of soma
            %           bscall pixel
            
            for i = 1:size(somaBscallTable,1)
                imName = somaBscallTable.Properties.RowNames{i};
                
                % Pixel with soma bscall
                somaPixel = [somaBscallTable.x{i},somaBscallTable.y{i}];
                
                iBscallTable = tfBscallTable.bscall{imName};
                
                % Pixel with dot bscall
                dotPixel = [iBscallTable.x,iBscallTable.y];
                dotPixel = round(dotPixel);
                
                % Delete dot pixel with soma bscall
                row = ismember(dotPixel,somaPixel,'rows');
                iBscallTable = iBscallTable(~row,:);
                
                tfBscallTable.bscall{imName} = iBscallTable;
            end
            
            disp('Done: delDotPixelInSoma.')
        end
        
    end
    
    methods (Static)    % Stitch dots =====================================
        %% Function:    SD_tfSpot4fixTile
        % Discription: consolidate tfSpot for fix tiles
        % (function for stitch dot)
        function [fixNameC,tfSpot] = SD_tfSpot4fixTile(fixName,tfSpot)
            % Input & output: fixName, cell,
            %           tfSpot, cell, one mat per cell
            
            [fixNameC,~,ic] = unique(fixName);
            
            tfSpot2 = {};
            for i = 1:numel(fixNameC)
                row2 = ic == i;
                tfSpot2{i,1} = tfSpot(row2);
            end
            
            tfSpot = cellfun(@(X) vertcat(X{:}),tfSpot2,'UniformOutput',false);
        end
        
        %% Function:    SD_sortSpotNum
        % Discription:  sort spot number from more to less
        % (function for stitch dot)
        function [fixName,tfSpot] = SD_sortSpotNum(fixName,tfSpot)
            % Input & output: fixName, tfSpot, cell
            
            n = cellfun(@(X) size(X,1),tfSpot);
            [~,I] = sort(n,'descend');
            fixName = fixName(I); tfSpot = tfSpot(I);
        end
        
        %% Function:    stitchDot(main)
        % Discription:  find transformation matrix for each tile
        function stitchTformTable = stitchDot(tfDotTable,accumulateTformTable,...
                imageSetting,sysSetting)
            
            % Get image name
            % rowName2ImName(tbl,delimiter,sysSetting)
            rowNames = TBS.rowName2ImName(tfDotTable,true,sysSetting);
            
            imNames = unique(rowNames);
            
            % Initial alignment seq -middle sequencing cycle
            initialSeq = imageSetting.initialFixSeq;
            initialSeq = TBS.seqstr(initialSeq);
            
            stitchTformTable = {};
            parfor iIm = 1:size(imNames,1) % parfor
                iImName = imNames{iIm};
                
                % Tiles belongs to the image
                tileName = contains(rowNames,iImName);
                tileName = tfDotTable.Properties.RowNames(tileName);
                
                % Find initial fix for every tile -------------------------
                % Tile from the inital alignment
                row = contains(tileName,initialSeq);
                row = tileName(row);
                initialTileTable = accumulateTformTable(row,:);
                
                % tform for the tile (with fuse image alignment)
                tileTform = accumulateTformTable.tform(tileName);
                
                % Find the closest initial tile
                [tileFixName,~] = cellfun(@(X) TBS.findClosestFix(initialTileTable,X),...
                    tileTform,'UniformOutput',false);
                
                % Get the spot coordinates for each inital fix ------------
                % Dot alignment dim = 3 (align with channel)
                tfSpot = tfDotTable{tileName,1};
                tfSpot = cellfun(@(X) [X.x, X.y, single(X.bscallCh)],...
                    tfSpot,'UniformOutput',false);
                
                % Unique fixName and cooresponding tfSpot coordinates
                [fixName,tfSpot] = TBS.SD_tfSpot4fixTile(tileFixName,tfSpot);
                
                % Sort the tile from more --> less dots
                [fixName,tfSpot] = TBS.SD_sortSpotNum(fixName,tfSpot);
                
                % Loop from more to less spot tiles -----------------------
                tblOut = table();
                tblOut.tform{fixName{1}}= eye(3);
                fix = tfSpot{1};
                
                iTile = 2;
                while iTile <= size(fixName,1)
                    
                    movingName = fixName{iTile};
                    
                    % Skip if it's aligned
                    if TBS.hasThisRow(tblOut,movingName)
                        iTile = iTile+1;
                        continue
                    end
                    
                    moving = tfSpot{iTile};
                    
                    if isempty(moving)
                        iTile = iTile + 1;
                        continue
                    end
                    
                    % alignDot(moving,fix,dim,binSz,dotLim,maxDist)
                    % tform: projective/affine object
                    tform = TBS.alignDot(moving,fix,3,5,1000,60);
                    
                    if isempty(tform)
                        iTile = iTile + 1;
                        continue
                    end
                    
                    tform = tform.T;
                    
                    % (Optional) For manual testing alignment
                    % TBS.checkDotAlign(tform,moving,fix);
                    
                    % Add transformed dot to fix
                    tfMoving = transformPointsForward(projective2d(tform),...
                        moving(:,1:2));
                    tfMoving = [tfMoving, moving(:,3)];
                    fix = [fix; tfMoving];
                    
                    tblOut.tform{movingName} = tform;
                    
                    disp(['Stitched tiles: ',movingName]);
                    
                    % Go to the beginining of the loop to go down again
                    iTile = 2;
                end
                
                % Delete fix tile non-aligned
                row = cellfun(@(X) TBS.hasThisRow(tblOut,X),tileFixName);
                tileFixName = tileFixName(row);
                tileName = tileName(row);
                
                tform = tblOut.tform(tileFixName);
                stitchTformTable{iIm,1} = [tileName, tform];
            end
            
            stitchTformTable = vertcat(stitchTformTable{:});
            stitchTformTable = table(stitchTformTable(:,2),...
                'VariableNames',{'tform'},'RowNames',stitchTformTable(:,1));
        end
        
    end
    
    methods (Static)    % Fuse sequence images (FSI)=======================
        %% Function:    FSI_estimateStitchImSize
        % Discription:  estimate size for image with multiple cycles
        % (function for fusing dot image)
        function sz = FSI_estimateStitchImSize(nameQ,tableIn,paddingSize,imageSetting,sysSetting)
            % Tile size can be difference across cycles, will use the max
            % to define output image size
            % Input:    nameQ, str, name on query
            %           tableIn, table, table for providing row names
            %           paddingSize, num, padding size for output
            %           imageSetting/sysSetting, struct
            % Output:   sz, mat, image size
            
            % rowName2ImName(tbl,delimiter,sysSetting)
            name = TBS.rowName2ImName(tableIn,true,sysSetting);
            
            tileName = contains(name,nameQ);
            tileName = tableIn.Properties.RowNames(tileName);
            
            % Extract tile numbers
            % bug fix? 08292021
            tileNum = cellfun(@(X) TBS.nameFun(X,sysSetting.tileElement,...
                sysSetting),tileName,'UniformOutput',false);
            tileNum = cellfun(@str2num,tileNum);
            
            % Estimate tile position and size
            tilePos = imageSetting.tilePos;
            iTilePos = TBS.estimateTotalNumOfTile(tileNum,tilePos);
            iTilePos = tilePos{iTilePos};
            
            sz = TBS.estimateStitchImSize(iTilePos,imageSetting,paddingSize);
        end
        
        %% Function:    fuseSeqImage (main)
        function paddingTformTable = fuseSeqImage(tformTable,redoTF,...
                paddingSize,overlapTF,imageSetting,sysSetting,directory)
            % Cannot directly use fuse image due to multiseq
            % Input:    tformTable,table, with transformation matrix
            %           redoTF, logical, whether redo the stitched image
            %           paddingSize, mat, padding area for output
            %           overlapTF, logical whether the fuse image has max
            %           projection in overlap
            %           imageSetting/sysSetting/directory,struct
            % Output:   paddingTformTable, table, padding tform for each
            % tile
            
            % If there is no output, there is no be trimming
            if ~nargout
                disp('Function fuseSeqImage: No trimming for stitching image.');
            end
            
            if ~exist(directory.stitched)
                mkdir(directory.stitched);
            end
            
            outputMatName = 'paddingTformTable.mat';
            
            % Function handle to change image bit
            imageBits = imageSetting.imageBits;
            imageBits = @(X) cast(X,imageBits);
            
            nCh = imageSetting.chNum;
            
            % Sequencing cycle for output
            seqCycle = imageSetting.seqCycles;
            seqCycle = min(seqCycle):max(seqCycle);
            
            % Add padding tform to the current tform ----------------------
            paddingTform = eye(3);
            paddingTform(3,1:2) = paddingSize./2;
            
            tformTable.tform = cellfun(@(X) X*paddingTform,...
                tformTable.tform,'UniformOutput',false);
            
            % Get image name ----------------------------------------------
            % rowName2ImName(tbl,delimiter,sysSetting)
            rowNames = TBS.rowName2ImName(tformTable,true,sysSetting);
            
            imNames = unique(rowNames);
            
            % Stitch image -----------------------------------------------
            paddingTformTable2 = {};
            for iIm = 1:size(imNames,1)
                iImName = imNames{iIm};
                
                % Image file nmae
                fileName = TBS.imFormat(iImName,sysSetting);
                
                % Skip stitched image
                if exist(fullfile(directory.stitched,fileName)) && ~redoTF
                    continue
                end
                
                % Estimate the size of stitched image
                sz = TBS.FSI_estimateStitchImSize(iImName,tformTable,...
                    paddingSize,imageSetting,sysSetting);
                R = imref2d(sz);
                
                im = imageBits(zeros([sz nCh numel(seqCycle)]));
                tblOut = {};
                parfor iSeq = seqCycle % parfor
                    iSeqstr = TBS.seqstr(iSeq);
                    
                    % Tiles belongs to the image
                    tileName = contains(rowNames,iImName);
                    tileName = tformTable.Properties.RowNames(tileName);
                    tileName = tileName(contains(tileName,iSeqstr));
                    
                    % Skip the sequence if no tile has been found
                    if isempty(tileName)
                        continue
                    end
                    
                    % 08292021, stitch sequence
                    % Allow delete later overlap regions due to possible
                    % optical distortion
                    I = TBS.getStitchingSeq(tileName,sysSetting,imageSetting);                    
                    tileName = tileName(I,:);
                    
                    % Add additional stack for tile overlap pixels
                    stack = imageBits(zeros([sz nCh+1]));
                    
                    for iTile = 1:size(tileName,1)
                        iTileName = tileName{iTile};
                        
                        cd(directory.main);
                        iDirectory = dir(fullfile('**',iTileName));
                        cd(iDirectory.folder);
                        
                        iStack = TBS.getStack(iTileName,[]);
                        % add stack for tile position
                        iStack(:,:,end+1) = inf; 
                        
                        tform = tformTable.tform{iTileName};
                        % Need to set to double to avoid warning
                        tform = double(tform);
                        tform = projective2d(tform);
                        
                        iStack = imwarp(iStack,tform,'OutputView',R);
                        
                        % 08292021, allow delete later overlap region
                        if ~overlapTF
                            overlapRegion = any(stack,3);
                            overlapRegion = repmat(overlapRegion,1,1,size(iStack,3));
                            iStack(overlapRegion) = 0;
                        end
                        
                        % Combine with the current image
                        stack = max(iStack,stack);
                    end
                    
                    stack = stack(:,:,1:nCh); 
                    
                    im(:,:,:,iSeq) = stack;
                    tblOut{iSeq,1} = tileName;
                    
                    disp(['Stitched image ',fileName,', ', iSeqstr]);
                end
                
                tblOut = vertcat(tblOut{:});
                
                % Trim image --------------------------------------------
                [sz,trimTform] = TBS.getTrimSetting(im);
                
                % No triming (when no output)
                if nargout > 0
                    R = imref2d(sz);
                    im = imwarp(im,trimTform,'OutputView',R);
                end
                
                % Combind padding tform to trimTform
                trimTform = paddingTform*trimTform.T;
                
                tblOut(:,2) = repmat({trimTform},size(tblOut));
                paddingTformTable2{iIm,1} = tblOut;
                
                cd(directory.stitched);
                TBS.saveStack(im,fileName);
            end
            
            % Convert to table --------------------------------------------
            if isempty(paddingTformTable2)
                error('No tiles has been stitched.')
            end
            paddingTformTable2 = vertcat(paddingTformTable2{:});
            
            cd(directory.main)
            if exist(outputMatName) && ~redoTF
                load(outputMatName);
                paddingTformTable(paddingTformTable2(:,1),:) = ...
                    paddingTformTable2(:,2);
            else
                paddingTformTable = table(paddingTformTable2(:,2),...
                    'VariableNames',{'tform'},'RowNames',paddingTformTable2(:,1));
            end
        end
        
    end
    
    methods (Static)    % Soma bascalling =================================
        %% Function:    BClongEngouth
        % Discription:  test whether BC is long enough
        function TF = BClongEngouth(bcTF,excludeSeq,minBClen)
            % Input:    bcTF, logical, BC pass some specific criteria
            %           excludeSeq, mat/num, seq cycle to be excluded
            %           minBClen, num, minimum BC length
            % Output:   TF, logical, whether BC is long enough
            
            if ~islogical(bcTF)
                bcTF = bcTF~=0;
                warning('Function BClongEngouth: converted the input to logical');
            end
            
            TF = squeeze(bcTF);
            TF(:,excludeSeq) = [];
            TF = sum(TF,2) >= minBClen;
        end
        
        %% Function:    somaIntenNormalization
        % Discription:  normalized image for soma bscalling, using
        % reference intensity
        function im = somaIntenNormalization(im)
            % Input & output: im, 4-D image stack
            
            imBit = class(im);
            imBitFh = @(X) cast(X,imBit);
            
            im = single(im);
            
            refInten = max(im,[],3);
            refRatio = 1./refInten.*refInten(:,:,:,1);
            
            im = im.*refRatio;
            
            im = imBitFh(im);
        end
        
        %% Function:    somaBkgrdSubtraction
        % Discription:  remove soma background by subtracting previous
        % cycle
        function im = somaBkgrdSubtraction(im)
            % Note: this is for removing phasing or other nonspecific
            % signals
            % Input & output:   im, image stack,
            %           nCh, channel number
            %           nSeq, sequence number
            
            imBit = class(im);
            imBitFh = @(X) cast(X,imBit);
            
            % For floatting calculation
            im = single(im);
            
            % Background as pervious seq cycle
            bkgrd = im(:,:,:,1:end-1);
            
            % 50% max intensity as background, other channel remain the same
            maxBkgrd = max(bkgrd,[],3);
            % Find max inten pixel
            maxBkgrd = bkgrd == maxBkgrd;
            maxBkgrd = single(maxBkgrd);
            % Current subtraction setting: 50%
            maxBkgrd = maxBkgrd.*(-0.5) + 1;
            bkgrd = bkgrd.*maxBkgrd;
            
            bkgrd = imBitFh(bkgrd); im = imBitFh(im);
            
            % Background subtraction
            im(:,:,:,2:end) = im(:,:,:,2:end) - bkgrd;
        end
        
        %% Function:    somaBscall (main)
        % Discription:  soma base calling using pixel
        % Pixel with >= min BC len pass the threshold, will bscall all seq
        % bscall channel will be filtered
        function somaBscallTable = somaBscall(somaBscallThreshold,...
                dotMatchingSetting,imageSetting,sysSetting,directory)
            % Input:    somaBscallThreshold, struct
            %           dotMatchingSetting, struct
            %           imageSetting/sysSetting/directory, sturct
            % Output:   somaBscallTable, table, with pixel index and bscall
            % channel
            
            % Filter size for image smoothing (micron)
            filterSize = 5;
            % Convert to pixel distance
            filterSize = filterSize/imageSetting.resolution;
            
            nCh = imageSetting.chNum;
            
            seqCycle = imageSetting.seqCycles;
            seqCycle = min(seqCycle):max(seqCycle);
            nSeq = numel(seqCycle);
            
            % Sequencing cycle not included in the length count
            excludeSeq = dotMatchingSetting.seqNotInMinBscallCycles;
            % Minimum barcode lenght
            minBClen = dotMatchingSetting.minBClen;
            
            % Get noise in 3-D (2-D, ch; 3-D, seq) ------------------------------------
            noiseStructure = TBS.getInitialNoise(somaBscallThreshold,imageSetting,directory);
            noise = [];
            for iSeq = seqCycle
                % Noise for current cycle
                iNoise = TBS.getNoise(iSeq,noiseStructure.initial,noiseStructure.change);
                noise = cat(3,noise,iNoise);
            end
            
            % Soma bscall -------------------------------------------------
            cd(directory.stitched)
            imName = ls(['*',sysSetting.imFormat]);
            imName = cellstr(imName);
            
            % Image for soma bscall
            row = cellfun(@(X) TBS.doSomaBscall(X,sysSetting),imName);
            imName = imName(row);
            
            tblOut = {};
            parfor iIm = 1:size(imName,1) % parfor
                iImName = imName{iIm};
                
                cd(directory.stitched)
                stack = TBS.getStack(iImName,[]);
                sz = size(stack,[1 2]);
                
                % Use original image (no filter) to pick positive pixel
                % 06212021 alignmentTest25f
                stack0 = stack;
                
                % Somoothing using resize
                % Also speed up
                stack = imresize(stack,1/filterSize);
                
                % stack, change into 4-D
                stack = reshape(stack,size(stack,1),size(stack,2),nCh,nSeq);
                
                % Adjust to similar brightness across cycle
                % For background subtraction (previous cycle)
                stack = TBS.somaIntenNormalization(stack);
                
                % Background subtraction
                stack = TBS.somaBkgrdSubtraction(stack);
                
                % Adjust to similar brightness after subtraction
                % Normalize seq1 and other cycles
                stack = TBS.somaIntenNormalization(stack);
                
                % Scale back
                stack = imresize(stack,sz);
                
                % Convert to 3d: column-channel, 3rd dimension-sequence
                stack = reshape(stack,[],nCh,numel(seqCycle));
                
                % Pixel pass the threshold --------------------------------
                % Use the original image (stack0)
                stack0 = reshape(stack0,[],nCh,numel(seqCycle));
                abvNoise = stack0 >= noise;
                
                % Pixel with BC long enough -------------------------------
                row = any(abvNoise,2);
                row = TBS.BClongEngouth(row,excludeSeq,minBClen);
                ind = find(row);
                
                chIntensity = stack(row,:,:);
                
                [v,bscallCh] = max(chIntensity,[],2);
                
                % 08222021, more than one max channel
                % (a lof of z)
                TF = chIntensity == v;
                TF = sum(TF,2) > 1;
                bscallCh(TF) = 0;
                
                bscallCh = squeeze(bscallCh);
                bscallCh = uint8(bscallCh);
                
                [y,x] = ind2sub(sz(1:2),ind);
                x = single(x); y = single(y);
                
                % 08222021, trim again after allow 0 for chIntensity
                TF = bscallCh ~= 0;
                TF = TBS.BClongEngouth(TF,excludeSeq,minBClen);
                bscallCh = bscallCh(TF,:);
                x = x(TF,:); y = y(TF,:);
                
                tblOut{iIm,1} = [{iImName},{bscallCh},{x},{y}];
                
                disp(strcat('Get soma bscall: ', iImName));
            end
            
            tblOut = vertcat(tblOut{:});
            
            varName = {'bscallCh','x','y'};
            somaBscallTable = cell2table(tblOut(:,2:end),'VariableNames',...
                varName,'RowNames',tblOut(:,1));
            
            disp('Done: somaBscall.');
        end
        
    end
    
    methods (Static)    % Dot bscalling ===================================
        %% Function:    splitBscallTableIntoCell
        % Discription:  split bscall table into cell basing on seq cycles
        %               one row with one column/seq
        function bscallCell = splitBscallTableIntoCell(bscallTable)
            % Input:    bscallTable, table, with seq info in z
            % Output:   bscallCell, cell
            
            nSeq = max(bscallTable.seq);
            
            bscallCell = arrayfun(@(X) bscallTable(bscallTable.seq == X,...
                {'x','y','bscallCh'}),1:nSeq,'UniformOutput',false);
        end
        
        %% Function:    getDotID
        % Discription:  get dotID for bscall table with coordinates
        function bscallTable = getDotCode(bscallTable)
            % Input:    bscallTable, table, with xy (x, column)
            % Output:   bscallTable with dotID
            
            bscallTable = sortrows(bscallTable,'x');
            
            n = size(bscallTable,1);
            bscallTable.id = (1:n)';
        end
        
        %% Function:    getInterSeqPair
        % Discription:  get seq tile pair for dot matching
        function interSeqPair = getInterSeqPair(seqCycle,maxSeqInterval)
            % Input:    seqCycle, row vector, sequencing cycle
            %           maxSeqInterval, number, max sequencing interval
            % Output:   interSeqPair, cell, seq cycle pair for inter tile alignment
            
            interSeqPair = {};
            for i = 1:maxSeqInterval
                
                iSeqCycle = seqCycle(1:end-i)';
                
                interSeqPair{i,1} = arrayfun(@(X) [X,X+i],iSeqCycle,'UniformOutput',false);
            end
            
            interSeqPair = vertcat(interSeqPair{:});
        end
        
        %% Function:    getDotPair
        % Discription:  match dots (id) across seq cycles
        function pairID = getDotPair(moving,fix,maxDist)
            % Input:    moving/fix, table, with dot ID and coordinates
            %           maxDist, num, max distance for each pair
            % Output:   pairID, single mat, column 1: fix id; column 2: moving id
            
            pairID = single([]);
            while ~isempty(moving) && ~isempty(fix)
                % (to speed up) exclude fix with min distance larger than
                % maxDist -------------------------------------------------
                D = pdist2([moving.x,moving.y],[fix.x,fix.y],...
                    'euclidean','Smallest',1);
                row = D <= maxDist;
                fix = fix(row,:);
                
                if isempty(fix)
                    return
                end
                
                % Find the nearest fix to every moving --------------------
                [D,I] = pdist2([fix.x,fix.y],[moving.x,moving.y],...
                    'euclidean','Smallest',1);
                D = D'; I = I';
                
                % pair to be determine: Distance, fix id, moving id
                pairTBD = [D,fix.id(I),moving.id];
                pairTBD = single(pairTBD);
                
                % Delete the pair above max distance
                row = D <= maxDist;
                pairTBD = pairTBD(row,:);
                % Delete moving dots with nearest fix above max distance
                moving = moving(row,:);
                
                if isempty(moving)
                    return
                end
                
                % Find pair -----------------------------------------------
                [~,I] = sort(pairTBD(:,1),'ascend');
                % Fix id, moving id
                pairTBD = pairTBD(I,2:3);
                % Fix pair with smallest distance
                [~,ia,~] = unique(pairTBD(:,1));
                iPairID = pairTBD(ia,:);
                
                % Delete the fix and moving ID in the matched pair --------
                row = ~ismember(fix.id,iPairID(:,1));
                fix = fix(row,:);
                row = ~ismember(moving.id,iPairID(:,2));
                moving = moving(row,:);
                
                pairID = [pairID; iPairID];
            end
        end
        
        %% Function:    blwMaxInterval
        % Discription:  barcode with interval below max interval
        function TF = blwMaxInterval(id,maxSeqInterval)
            % Max interval is not included missing seqing cycle
            % Note03252022: function name change from abvMaxInterval to
            % below to fit the function
            % Input & output: id, mat, barcode matrix, one row per bc, one
            % column per seq
            %           maxSeqInterval, num, maximum interval with no
            %           bascalling
            
            % without bscall
            noBscall = id == 0;
            
            % Delete col for missing cycles (no bscall result)
            col = all(noBscall,1);
            noBscall = noBscall(:,~col);
            
            % bc matrix at least have two basecalling seq
            if size(noBscall,2) <= maxSeqInterval+2
                TF = true(size(id,1),1); return
            end
            
            % Size beyond the max interval: maxInterval+1
            if mod(maxSeqInterval,2)==1
                isEven = true;
                SE = ones(1,maxSeqInterval);
            else
                isEven = false;
                SE = ones(1,maxSeqInterval+1);
            end
            
            % Padding false for edge effect, size of padding
            padding = false(size(noBscall,1),ceil(numel(SE)/2));
            noBscall = [padding,noBscall,padding];
            
            % Find rows with interval >= max interval
            row = imerode(noBscall,SE);
            
            % Even SE: two continous odd
            if isEven
                row = row(:,1:end-1) & row(:,2:end);
            end
            
            row = any(row,2);
            
            TF = ~row;
        end
        
        %% Function:    mergeIDlogical
        % Discription:  merge BC ID basing on logical
        function id = mergeIDlogical(id,tolerate0)
            % Input:    id, mat, barcode matrix,  one row per bc, one
            % column per seq
            %           tolerate0, logical, whether tolerate 0 in bc matrix
            % Output:   id, mat, merged barcode matrix
            
            % size check for merge
            nID = size(id,1);
            
            nSeq = size(id,2);
            
            i = 1;
            while i <= nSeq
                id = sortrows(id,i);
                
                % Identical BC row ----------------------------------------
                TF = id(1:end-1,:) == id(2:end,:);
                
                if tolerate0 % Ignore 0
                    TF = TF | id(1:end-1,:)==0 | id(2:end,:)==0;
                end
                
                % Identical to the next row
                TF = all(TF,2);
                
                % Merge BC ------------------------------------------------
                if any(TF)
                    row = find(TF);
                    
                    if tolerate0
                        mergeBC = max(id(row,:),id(row+1,:));
                        id(row,:) = mergeBC;
                    end
                    
                    % Delete repeat rows
                    id(row+1,:)=[];
                    continue
                end
                
                if i < nSeq
                    i = i+1;
                    continue
                elseif size(id,1) == nID
                    break
                elseif size(id,1) < nID
                    i = 1; nID = size(id,1);
                end
            end
        end
        
        %% Function:    id2bscallTbleVar
        % Discription:  get bscall table variable using row number (dot id)
        function var = id2bscallTbleVar(id,bscallTbl,varName)
            % Input:    id, vector, row number/dot ID
            %           bscallTbl, table, with bscall info
            %           varName, str, variable name of table
            % Output:   var, vector, of dot ID
            
            % Same class as the table variable
            classType = class(bscallTbl.(varName));
            classFh = @(X) cast(X,classType);
            
            var = zeros(size(id));
            var = classFh(var);
            
            % Only fill the non-zero rows
            row = id~=0;
            id2 = id(row);
            
            var(row) = bscallTbl.(varName)(id2);
        end
        
        %%  Function:   hasOverlap
        % Discription:  find BC with identical BC with more 0
        function TF = hasOverlap(BCin)
            % Using ismember instead of for-loop to speed up
            % Input:    BCin, mat, one row per BC
            % Output:   TF, logcial, BC with identical BC and extra 0
            
            % Find 0 position pattern
            pos0i = BCin == 0;
            [pos0i,~,ic] = unique(pos0i,'row');
            
            % BC with extra 0 will be identical to other BC without
            % 0-position column, not the other way around
            TF = {};
            parfor i = 1:size(pos0i,1) % parfor
                
                iPos0 = pos0i(i,:);
                
                % Row with the 0 pattern
                rowi = ic == i;
                iBCin = BCin(rowi,:);
                
                % Delete 0 column
                iBCin = iBCin(:,~iPos0);
                BCin2 = BCin(:,~iPos0);
                
                [Lia,Lcob] = ismember(BCin2,iBCin,'rows');
                
                % Count overlap
                Lcob = Lcob(Lia);
                n = accumarray(Lcob,1,[size(iBCin,1),1]);
                % Find more than 1 overlap (ttself is a match)
                n = n > 1;
                
                rowi(rowi) = n;
                
                TF{1,i} = rowi;
            end
            
            % Combine ovelrap from each pattern
            TF = horzcat(TF{:});
            TF = any(TF,2);
        end
        
        %% Function:    getBscallTable (main)
        % Discription:  get bscall table, with barcode,x,y
        function bscallTable = getBscallTable(tfBscallTable,...
                dotMatchingSetting,imageSetting)
            % Input:    tfBscallTable, table, bscall info of dot
            %           dotMatchingSetting/imageSetting,struct
            % Output:   bscallTable, table, each cell, one row per bc
            
            seqCycle = imageSetting.seqCycles;
            maxSeqInterval = dotMatchingSetting.maxSeqInterval;
            maxDist = dotMatchingSetting.maxDist;
            excludeSeq = dotMatchingSetting.seqNotInMinBscallCycles;
            minBClen = dotMatchingSetting.minBClen;
            
            bscallTable = {};
            for iIm = 1:size(tfBscallTable,1)
                imName = tfBscallTable.Properties.RowNames{iIm};
                
                % One seq cycle per cell
                bscallCell = TBS.splitBscallTableIntoCell(tfBscallTable.bscall{iIm});
                
                % Dot ID is the row number of the sequencing cycle
                bscallCell = cellfun(@(X) TBS.getDotCode(X),bscallCell,'UniformOutput',false);
                
                % Seq pair for alingment ----------------------------------
                inteSeqPair = TBS.getInterSeqPair(seqCycle,maxSeqInterval);
                
                % Delete empty seq
                emptySeq = cellfun(@isempty,bscallCell);
                emptySeq = find(emptySeq);
                row = cellfun(@(X) any(ismember(X,emptySeq)),inteSeqPair);
                inteSeqPair = inteSeqPair(~row);
                
                % Sort alingment sequence,08282021
                % position 2 need to be ascend
                % Position 1 need to be descend to decrease multiple
                % possible match in ismember(iPairID(:,1),id(:,iPair(1)))
                % when adding a new digit, start from the closet nt 
                I = cellfun(@(X) X(1),inteSeqPair);
                [~,I] = sort(I,'descend');
                inteSeqPair = inteSeqPair(I);
                
                I = cellfun(@(X) X(2),inteSeqPair);
                [~,I] = sort(I,'ascend');
                inteSeqPair = inteSeqPair(I);
                
                % Match dots ----------------------------------------------
                % pair of dot ID for fix & moving for each seq pair
                % TBS.getDotPair(moving,fix,maxDist);
                pairID = cellfun(@(X) TBS.getDotPair(bscallCell{X(2)},...
                    bscallCell{X(1)},maxDist),inteSeqPair,'UniformOutput',false);
                
                % ID for each barcode -------------------------------------
                % 04262021: allow missing cycle at the begnining
                id = single([]);
                id(:,inteSeqPair{1}) = pairID{1};
                for i = 2:size(inteSeqPair,1)
                    
                    % sequce for matching
                    iPair = inteSeqPair{i};
                    
                    % Dot ID for this sequence pair
                    iPairID = pairID{i};
                    
                    % 08292021,sort using basecalling digit before matching
                    % When there are two identical id, pick the one with
                    % longer BC 
                    % (descend and ascend will give different result)
                    n = sum(id ~= 0,2);
                    [~,I] = sort(n,'descend');
                    id = id(I,:);
                    
                    % Find the fix id in existed BC id
                    % 08272021 Not multiple match for A, inflate a lot in
                    % later seq, Use other seq to find the match
                    [Lia,Lcob] = ismember(iPairID(:,1),id(:,iPair(1)));
                    % 08282021, only do 1:iPair(1), leave the middle nt
                    iId = id(Lcob(Lia),1:iPair(1));
                    iId(:,iPair(2)) = iPairID(Lia,2);
                    
                    % For empty first N-seq, add 0 ------------------------
                    % 08272021, bug fix
                    if iPair(1) <= maxSeqInterval+1
                        iId2 = [];
                        iId2(:,iPair) = iPairID(~Lia,:);                          
                        iId = [iId; iId2];                        
                    end
                    
                    % Add the new matching to the output ------------------
                    % Pad array if the size are not the same
                    if size(id,2) < size(iId,2)
                        id(1,size(iId,2)) = 0;
                    end
                    id = [id; iId];
                    
                    % Clean output: delete interval > max interval --------
                    % 08282021, delete after matching all the 4th position,
                    % there is still 3-nt space in the middle
                    TF = TBS.blwMaxInterval(id,maxSeqInterval+1);
                    id = id(TF,:);
                    
                    % Merge rows ------------------------------------------
                    % 10052021 ver27: add hasOveralp
                    TF = TBS.hasOverlap(id);
                    id = id(~TF,:);
                    
                    % mergeBClogical(bc,tolerate0);
                    id = TBS.mergeIDlogical(id,true);
                end
                
                % Delete the maxSeqInterval+1 interval for the last digit
                TF = TBS.blwMaxInterval(id,maxSeqInterval);
                id = id(TF,:);
                
                % Exclude short barcodes ----------------------------------
                % BClongEngouth(bcTF,excludeSeq,minBClen);
                row = TBS.BClongEngouth(id > 0,excludeSeq,minBClen);
                id = id(row,:);
                
                % Bscall info (Ch, x, y) ----------------------------------
                id2 = mat2cell(id,size(id,1),ones(size(id,2),1));
                
                % id2bscallTbleVar(id,bscallTbl,varName)
                bscallCh = cellfun(@(X,Y) TBS.id2bscallTbleVar(X,Y,'bscallCh'),...
                    id2,bscallCell,'UniformOutput',false);
                x = cellfun(@(X,Y) TBS.id2bscallTbleVar(X,Y,'x'),...
                    id2,bscallCell,'UniformOutput',false);
                y = cellfun(@(X,Y) TBS.id2bscallTbleVar(X,Y,'y'),...
                    id2,bscallCell,'UniformOutput',false);
                
                bscallCh = horzcat(bscallCh{:});
                x = horzcat(x{:});
                y = horzcat(y{:});
                
                bscallTable(iIm,:) = [{imName},{id},{bscallCh},{x},{y}];
                
                disp(['Got BC: ',imName]);
            end
            
            varName = {'id','bscallCh','x','y'};
            bscallTable = cell2table(bscallTable(:,2:end),...
                'VariableNames',varName,'RowNames',bscallTable(:,1));
            
            disp('Done: getBscallTable.')
        end
        
    end
    
    methods (Static)    % Dot bscall correction ===========================
        %% Function:    deletePureChDot
        % Discription:  delete dots with (almost) pure bscalling reuslts
        function tbl = deletePureChDot(tbl,minDiffNt,tolerate0)
            % Input & Output:   tbl, table, with id, bscallCh, x & y
            %       minDiffNt, num, min digits of diff channel nt
            %       tolerate0, logical, count 0 or not
            
            bscallCh = tbl.bscallCh{:};
            nSeq = size(bscallCh,2);
            
            % Find max number of nt in all channels
            fh = @(X) sum(bscallCh == X,2);
            mostCh = [fh(1),fh(2),fh(3),fh(4)];
            mostCh = max(mostCh,[],2);
            
            % Including 0
            if tolerate0
                mostCh = mostCh+ fh(0);
            end
            
            % With too many identical nt
            % (i.e. minDiffNt = 2, for nSeq = 17, BC cannot have >15
            % identical nt plus 0)
            TF = mostCh > (nSeq - minDiffNt);
            
            if ~any(TF)
                return
            end
            
            tbl{1,:} = cellfun(@(X) X(~TF,:),tbl{1,:},'Uniformoutput',false);
        end
        
        %% Function:   delPureChDotInTable
        % Discription: delete pure channel dot for each row in table
        function tbl = delPureChDotInTable(tbl,bcSetting)
            % Input & output:   tbl, table of bscall results
            %           bcSetting, struct
            
            minDiffNt = bcSetting.minDiffNt;
            
            % Delete pure channel dots
            for i = 1:size(tbl,1)
                % deletePureChDot(tbl,minDiffNt,tolerate0)
                tbl(i,:) = TBS.deletePureChDot(tbl(i,:),minDiffNt,true);
            end
            
            % Delete empty row
            row = cellfun(@isempty,tbl.bscallCh);
            tbl = tbl(~row,:);
        end
        
        %% Function:    getBscallBW
        % Discription:  for mannually checking bscalling result (logical)
        % Note, output size may not be the same as the image
        function im = getBscallBW(tbl)
            % Input:    table, with bscall result info
            % Output:   image, logical stack`
            
            v = tbl.bscallCh{:};
            v = single(v);
            
            nSeq = size(v,2);
            nCh = max(v,[],'all');
            
            xy = cat(3,tbl.x{:},tbl.y{:});
            % Soma bscall
            if size(xy,2) < nSeq && size(xy,2)==1
                xy = repmat(xy,1,nSeq,1);
            end
            
            % Image size -----------------------------------------
            sz = max(xy,[],[1 2]);
            sz = reshape(sz,1,2);
            % Get row and colum for sz
            sz(1:2) = fliplr(sz(1:2));
            sz = round(sz);
            sz = [sz,nCh];
            
            for iSeq = 1:nSeq
                ixyz = xy(:,iSeq,:);
                ixyz = squeeze(ixyz);
                
                ixyz = [ixyz,v(:,iSeq)];
                
                % Delete empty pixel for the cycle
                row = any(ixyz == 0,2);
                ixyz = ixyz(~row,:);
                
                iIm = TBS.xyzv2im(sz,ixyz,[]);
                
                im(:,:,:,iSeq) = iIm;
            end
            
            % Convert to 3-D
            im = reshape(im,size(im,1),size(im,2),[]);
        end
        
        %% Function:    checkBscallResult
        % Discription:  check single bscalling result
        function checkBscallResult(tbl,I)
            % Input:    table, bscalling result
            %           I, num, row number of result in the table
            % Output:   display the result
            
            bc = tbl.bscallCh{:}(I,:);
            
            x = tbl.x{:}(I,:); x = round(x);
            y = tbl.y{:}(I,:); y = round(y);
            xy = [x' y'];
            xy = xy(1,:);
            
            disp(['Seq 01 xy: ', mat2str(xy)]);
            disp(['BC: ', mat2str(bc)]);
        end
        
    end
    
    methods (Static)    % Align Seq to Ab =================================
        % Discription: Aligned sequencing images of individual regions into
        % the whole brain sticted image
        
        %% Function:    alignSeq2ab (main)
        % Discription:  align sequencing tile to stitched antibody image
        function seq2abTform = alignSeq2ab(abCh,scaleFactor,...
                redoTF,sysSetting,directory)
            % Note: use pixelCorrected image
            % Input:    abCh, num, channel for alignment in ab image
            %           nAlingment, num, # times for alignment
            %           scaleFactor, num, scale factor to speed up
            %           redoTF, logical, wheteher redo alignment
            %           imageSetting/sysSetting/directory, struct
            % Output:   seq2abTform, table, transformation matrix for
            % stitched sequencing image to antibody image
            
            outputName = 'seq2abTform.mat';
            
            % Output table
            cd(directory.main);
            if exist(outputName)
                load(outputName);
            else
                seq2abTform = table('Size',[0 2],'VariableTypes',...
                    {'cell','cell'},'VariableNames', {'fixName','tform'});
            end
            
            % Scale tform
            scaleTform = TBS.getScaleTform(scaleFactor,2);
            
            % File name of stitched Ab image
            cd(directory.abStitch);
            fileName = ls(['*',sysSetting.imFormat]);
            fileName = cellstr(fileName);
            
            for i = 1:size(fileName) % parfor
                fixName = fileName{i};
                fix = TBS.getStack(fullfile(directory.abStitch,fixName),abCh);
                % Adjust intensity
                fix = imadjust(fix);
                
                % To speed up, scale down
                fix = imresize(fix, scaleFactor);
                fix = im2uint8(fix);
                
                % Find sequencing image from initial fix tile -------------
                imName = TBS.nameFun(fixName,1,sysSetting);
                
                cd(directory.main);
                iFolder = dir(['**\',imName,'*']);
                iFolder = iFolder(1).folder;
                
                % Get initial fix seq
                initialFixSeq = TBS.getInitialFixSeq(iFolder);
                initialFixSeq = TBS.seqstr(initialFixSeq);
                
                cd(directory.main);
                imName = dir(['**\',imName,'*',initialFixSeq,'*',...
                    sysSetting.pixelCorrectAppend,'*']);
                
                % Loop through individual regions -------------------------
                tbl = {};
                for j = 1:size(imName,1) % parfor
                    movingName = imName(j).name;
                    
                    % Skip alignmed image
                    if ~redoTF && TBS.hasThisRow(seq2abTform,movingName)
                        continue
                    end
                    
                    disp(['Aligning to whole brain section: ',movingName]);
                    
                    cd(imName(j).folder);
                    moving = TBS.getStack(movingName,[]);
                    % Use min projection for alignment
                    moving = min(moving,[],3);
                    
                    % To speed up, scale down
                    moving = imresize(moving,scaleFactor);
                    moving = uint8(moving);
                    
                    % Alignment
                    tform = imregcorr(moving,fix,'Window',false);
                    tform = tform.T;
                    
                    % Quality check of tform
                    if ~TBS.QCtform(tform)
                        continue
                    end
                    
                    % tform correction for imregcorr ------------------------------
                    % Correct with sum of size of fix and moving, dont know why
                    % but it works
                    if tform(3,1)<= -size(moving,2)
                        tform(3,1) =  tform(3,1)+size(moving,2)+size(fix,2);
                    end
                    
                    if tform(3,2)<= -size(moving,1)
                        tform(3,2) =  tform(3,2)+size(moving,1)+size(fix,1);
                    end
                    
                    % % (Check point) ---------------------------------------------------
                    %                     R = imref2d(size(fix));
                    %                     tfMoving = imwarp(moving,affine2d(tform),'OutputView',R);
                    %                     imshowpair(fix,tfMoving.*20);
                    
                    % Convert back to original scale
                    tform = scaleTform*tform*inv(scaleTform);
                    
                    tbl(j,:) = [{movingName},{fixName}, {tform}];
                end
                
                % Delete empty rows
                row = cellfun(@isempty,tbl(:,1));
                tbl = tbl(~row,:);
                
                seq2abTform(tbl(:,1),:) = tbl(:,2:3);
                
                % (The process is too long, so save it in case bug)
                save(fullfile(directory.main,'seq2abTform.mat'),'seq2abTform');
            end
        end
        
        %% Function:    combineSeq2abTform
        % Discription:  combine tform from the same image (median)
        function tbl = combineSeq2abTform(tbl,tileTformName,nameElement,...
                sysSetting,directory)
            % Input & output:   tbl, table, with fixName & tform
            %               tileTformName, str, var name for stitching tile
            %               into image
            %               accumulateTformTable3, tbl
            %               sysSetting, struct
            
            % Get tileTform (acumulateTform3) -----------------------------
            tileTformName = erase(tileTformName,'.mat');
            
            % Get all the tile tform
            cd(directory.main);
            fileName = dir(['**',filesep,tileTformName,'.mat']);
            for i = 1:size(fileName,1)
                iFile = fullfile(fileName(i).folder,fileName(i).name);
                iFile = load(iFile);
                iFile = iFile.(tileTformName);
                if i == 1
                    tileTform = iFile;
                else
                    tileTform = [tileTform; iFile];
                end
            end
            
            % Change row name to pixel append
            tileTform = TBS.strrepRowName(tileTform,...
                sysSetting.chCorrectAppend,sysSetting.pixelCorrectAppend);
            
            % Get image name for each region ------------------------------
            fileName = tbl.Properties.RowNames;
            
            fileName = cellfun(@(X) TBS.nameFun(X,nameElement,sysSetting),...
                fileName,'UniformOutput',false);
            
            [fileName,~,ic] = unique(fileName);
            
            tbl2 = {};
            for i = 1:size(fileName,1)
                row = ic == i;
                
                iTbl = tbl(row,:);
                fixName = iTbl.fixName(1);
                movingName = iTbl.Properties.RowNames;
                
                % Exclude rows without tileTform
                row = cellfun(@(X) TBS.hasThisRow(tileTform,X),movingName);
                movingName = movingName(row);
                iTbl = iTbl(row,:);
                
                % tform in stitched imageg
                iTileTform = tileTform.tform(movingName);
                
                % Combine with the tile stitch tform
                % (need to reverse stitch tform before add the new tform)
                tform = iTbl.tform;
                tform = cellfun(@(X,Y) inv(X)*Y,iTileTform,tform,...
                    'UniformOutput',false);
                % Convert to the same data type
                tform = cellfun(@double,tform,'UniformOutput',false);
                
                % Use median for the tform
                tform = reshape(tform,1,1,[]);
                tform = cell2mat(tform);
                
                % Outlier exclusion
                I = isoutlier(tform,3);
                I = sum(I,1:2);
                I = I < 2;
                tform = tform(:,:,I);
                
                tform = median(tform,3);
                
                tbl2(i,:) = [fileName(i), fixName, {tform}];
            end
            
            tbl = cell2table(tbl2(:,2:3),'VariableNames',...
                {'fixName','tform'},'RowNames',tbl2(:,1));
        end
        
        %% Function:    stitchSeq2Ab
        % Discription:  stitch seq image to ab image
        function stitchSeq2Ab(seq2abTform,seq,chOut,imageSetting,directory)
            % Input:    seq2abTform, table, seq & ab image name, tform
            %           seq, mat, seq for stitched seq image
            %           chOut, num, channel for stiched seq in ab image
            %           imageSetting/directory, struct
            % Output:   append stitched sequencing image in chOut
            
            disp('Append Sequencing image to Ab image...')
            
            % Class of output image
            bitFh = imageSetting.imageBits;
            bitFh = @(X) cast(X,bitFh);
            
            % Channel number of seq image
            nCh = imageSetting.chNum;
            
            % Sequencing image stack number for ab output
            slide = (seq-1).*nCh +(1:nCh)';
            slide = reshape(slide,1,[]);
            
            [fixName,~,ic] = unique(seq2abTform.fixName);
            
            % Loop through the stiched image
            for i = 1:numel(fixName)
                iFixName = fixName{i};
                
                % Find the rows with the current fix image
                row = ic == i;
                row = find(row);
                
                cd(directory.abStitch);
                fix = TBS.getStack(iFixName,[]);
                
                sz = size(fix,[1 2]);
                R = imref2d(sz);
                
                moving = bitFh(zeros(sz));
                % Go through individual region
                for j = 1:numel(row)
                    iMovingName = seq2abTform.Properties.RowNames{row(j)};
                    iTform = seq2abTform.tform{row(j)};
                    iTform = projective2d(iTform);
                    
                    cd(directory.main);
                    iDir = dir(['**',filesep,iMovingName,'*']);
                    iDirectory = iDir.folder;
                    
                    iMovingName = iDir.name;
                    iDirectory = fullfile(iDirectory,iMovingName);
                    
                    % Max proj of seq01 image
                    iMoving = TBS.getStack(iDirectory,slide);
                    iMoving = max(iMoving,[],3);
                    
                    % Transformation
                    iMoving = imwarp(iMoving,iTform,'OutputView',R);
                    
                    % Add to the stitched image
                    moving = max(moving,iMoving);
                    
                    disp(['Done: ',iMovingName]);
                end
                
                fix(:,:,chOut) = moving;
                
                TBS.saveStack(fix,fullfile(directory.abStitch,iFixName));
            end
            
        end
        
    end
    
    methods (Static)    % Registration (w/ non-rigid) =====================
        %% Function:    scatteredInterpolantRng
        % Discription:  speed up scatteredInterpolant by excluding extra
        % points
        function Vq = scatteredInterpolantRng(X,Y,V,Xq,Yq,rng)
            % Input:    X,Y, mat/vector, coordinates
            %           V, mat/vector, can be multiple value Xq/Yq,
            %           mat/vector, point of inquary
            %           rng, number, range of reference points from inquary
            %           points to be consider
            % Output:   Vq, vector, value of inqury points
            
            fh = @(X) reshape(X,[],1);
            X = fh(X); Y = fh(Y);
            Xq = fh(Xq); Yq = fh(Yq);
            
            V = reshape(V,size(X,1),[]);
            
            % Exclude out of range dots
            if ~isempty(rng)
                % (To speed up) decrease point number by rounding
                % pdist2 is the limited step for a lot of inqury points
                XYq = [Xq,Yq];
                XYq = round(XYq);
                XYq = unique(XYq,'rows');
                
                D = pdist2(XYq,[X Y],'euclidean','Smallest',1);
                row = D' <= rng;
                
                % Delete XY far from inqury points (outside of range)
                X = X(row);
                Y = Y(row);
                V = V(row,:);
            end
            
            % Get inqury value
            Vq = [];
            for i = 1:size(V,2)
                % (return to object)
                iVq = scatteredInterpolant(X,Y,V(:,i));
                
                Vq(:,i) = iVq(Xq,Yq);
            end
            
        end
        
        %% Function:    interp2nonzeros
        % Discription:  get grid data by interperlating nonzero datapoints
        function im = interp2nonzeros(im)
            % Input & output:   im, m*n*2 mat, displacement field
            
            sz = size(im);
            x = 1:sz(2); y = 1:sz(1);
            [x,y] = meshgrid(x,y);
            xy = cat(3,x,y);
            
            % Within range pixels
            I = im == 0;
            I = ~any(I,3);
            
            im = reshape(im,[],2);
            xy = reshape(xy,[],2);
            
            % (return to object)
            Vx = scatteredInterpolant(xy(I,1),xy(I,2),im(I,1));
            Vy = scatteredInterpolant(xy(I,1),xy(I,2),im(I,2));
            
            % Get value across the grid
            Vx = Vx(x,y);
            Vy = Vy(x,y);
            
            im = cat(3,Vx,Vy);
        end
                
        %% Function:    tform2D
        % Discription:  tform object to displacement field (D)
        function D = tform2D(tform,sz)
            % 03272022, checked with imwarp
            % Input:    tform, MATLAB tform object, can work with
            %           non-rigid transformation object
            %           sz, mat, size of output (i.e. fix image)
            % Output:   D, mat, distplacement field
            
            % Grid of pixel location in moving
            x = 1:sz(2); y = 1:sz(1);
            [X,Y] = meshgrid(x,y);
            U = cat(3,X,Y);
            
            % Transformation
            R = imref2d(sz);
            V = imwarp(U,tform,'OutputView',R);
            
            % Interp using within range pixels
            V = TBS.interp2nonzeros(V);
            
            % Compute difference: displacement field
            D = V - U;                      
        end
        
        %% Function:    Doperation
        % Discription:  Combine displacement field
        function D = Doperation(D1,D2)
            % Input:    D1, (m,n,2) mat, 1st displacement field
            %           D2, (m,n,2) mat, 2nd displacement field
            % Output:   D, (m,n,2) mat, output displacment field
            
            % Grid of pixel location
            sz = size(D2);
            x = 1:sz(2); y = 1:sz(1);
            [x,y] = meshgrid(x,y);
            xy = cat(3,x,y);
            
            % D1 & D2-transformed pixel location
            % For example: original pixel location, xy; V1, pixel location
            % after D1
            V1 = D1 + xy;
            V2 = D2 + xy;
            
            % Using D2-D1 coordinates find D1-orignal coordinates
            Vx = interp2(V1(:,:,1),V2(:,:,1),V2(:,:,2),'cubic',0);
            Vy = interp2(V1(:,:,2),V2(:,:,1),V2(:,:,2),'cubic',0);
            
            V = cat(3,Vx,Vy);
            V = TBS.interp2nonzeros(V);
            
            % Calculate displacement field
            D = V - xy;            
        end
       
        %% Function:    transformPointsForwardD
        % Discription:  transform points forward using displacement field
        function X = transformPointsForwardD(D,U,rng)
            % Input:    D, mat, displacement field
            %           U, mat, input dot coordinates
            %           rng, num, range for interpolation
            % Output:   X, mat, transformed dot coordinates
            
            % Mesh of pixel location
            sz = size(D);
            x = 1:sz(2); y = 1:sz(1);
            [x,y] = meshgrid(x,y);
            xy = cat(3,x,y);
            
            % Pixel location in transformed image
            % D: where current location came from
            % xy: current pixel location, value, original pixel location
            % Note, this & inverseD may be counter intuitive, do a
            % checkboard and tform to check
            V = xy + D;
            
            % Get grid of before-after transform
            % (Use range to speed up)
            Xq = U(:,1); Yq = U(:,2);
            X = TBS.scatteredInterpolantRng(V(:,:,1),V(:,:,2),xy,Xq,Yq,rng);
        end
        
        %% Function:    transformPointsInverseD
        % Discription:  transform points inverse using displacement field
        function U = transformPointsInverseD(D,X)
            % Input:    D, mat, displacement field
            %           U, mat, input dot coordinates
            %           interpRng, num, interperlation range, to speed up
            % Output:   X, mat, transformed dot coordinates
            
            % Mesh of pixel location: xy
            sz = size(D);
            x = 1:sz(2); y = 1:sz(1);
            [x,y] = meshgrid(x,y);
            xy = cat(3,x,y);
            
            % Pixel location in transformed image
            % D: where current location came from
            % xy: current pixel location, value, original pixel location
            % Note, this & forwardD may be counter intuitive, do a
            % checkboard and tform to check
            V = xy + D;
            
            % Get the original pixel location by interpolation
            Xq = X(:,1); Yq = X(:,2);
            Vx = interp2(V(:,:,1),Xq,Yq,'cubic');
            Vy = interp2(V(:,:,2),Xq,Yq,'cubic');
            
            U = [Vx, Vy];
        end
        
        %% Function:    thresholdD
        % Discription:  thresholding selective displacement field
        function D = thresholdD(D, threhold)
            % Input & output:   D, cell, each cell contains m*n*2,
            % displacement field
            %           threshold, cell, threshold for each
            %           transformation method
            
            for i = 1:numel(D)
                iThreshold = abs(threhold{i});
                
                if isempty(iThreshold)
                    continue
                end
                
                iD = D{i};
                iD(iD < -iThreshold) = -iThreshold;
                iD(iD > iThreshold) = iThreshold;
                
                D{i} =  iD;
            end
        end
        
    end
    
    methods (Static)    % Self alignment (SA) =============================
        % Discription:  Align whole brain stitched section into a 3-D volume
        
        %% Function:    getSectionNameNum
        % Discription:  get sorted section name and number of the mouse
        function [sectionName,sectionNum] = getSectionNameNum(...
                directory,imageSetting,sysSetting)
            % Input:    directory, str, directory of stitched image files
            %           imageSetting/sysSetting, struct
            % Output:   sectionName, cell,
            %           sectionNum, mat, coresponding section number
            
            imFormat = sysSetting.imFormat;
            
            cd(directory);
            sectionName = ls(['*',imFormat]);
            sectionName = cellstr(sectionName);
            
            sectionNum = TBS.getSectionNumber(sectionName,imageSetting,sysSetting);
            
            [sectionNum,I] = sort(sectionNum,'ascend');
            sectionName = sectionName(I,:);
        end
        
        %% Function:    getFitgeotransFh
        % Discription:  get functional handles for each transformation
        function fitgeotransFh = getFitgeotransFh(transformationType)
            % Input:    transformationType, cell, transformation types
            % Output:   fitgeotransFh, cell, one function handle per cell
            
            % Function handels for each step
            fitgeotransFh = {};
            for i = 1:numel(transformationType)
                
                if contains(transformationType{i},'polynomial')
                    n = erase(transformationType{i},'polynomial');
                    n = str2double(n);
                    fitgeotransFh{i} = @(X,Y) fitgeotrans(X,Y,'polynomial',n);
                    
                elseif contains(transformationType{i},'lwm')
                    n = erase(transformationType{i},'lwm');
                    n = str2double(n);
                    fitgeotransFh{i} = @(X,Y) fitgeotrans(X,Y,'lwm',n);
                    
                else
                    fitgeotransFh{i} = @(X,Y) fitgeotrans(X,Y,transformationType{i});
                end
            end
        end
        
        %% Function:    cp2D
        % Discription:  control points to displacement field
        function D = cp2D(mp,fp,sz,transformationType)
            % Input:    mp/fp, mat, moving & fix point coordinates
            %           sz, mat, output size
            %           transformationType, cell, transformation types for
            %           fitgeotrans
            % Output:   D, cell, displacement field
            
            if ischar(transformationType)
                transformationType = {transformationType};
            end
            
            % Function handles for transformation types
            fitgeotransFh = TBS.getFitgeotransFh(transformationType);
            
            % Get displacement field from each step
            D = {};
            for i = 1:numel(fitgeotransFh)
                
                % Delete nan rows
                row = isnan([mp,fp]);
                row = ~any(row,2);
                mp = mp(row,:); fp = fp(row,:);
                
                tform = fitgeotransFh{i}(mp,fp);
                iD = TBS.tform2D(tform,sz);
                
                if isa(tform,'affine2d') || isa(tform,'projective2d')
                    mp = transformPointsForward(tform,mp);
                else
                    mp = TBS.transformPointsForwardD(iD,mp,50);
                end
                
                D{i} = iD;
            end
        end
        
        %% Function:    cp2Dn
        % Discription:  get displacement field from sequencial (n)
        % transformation using control points
        function D = cp2Dn(mp,fp,sz,transformationSetting)
            % Input:    mp/fp, mat, moving & fix point coordinates
            %           transformationType, cell, transformation types for fitgeotrans
            %           threshold, cell, threshold for each
            %           sz, mat, output size
            % Output:   D, m*n*2 mat, displacement field
            
            transformationType = transformationSetting.type;
            threshold = transformationSetting.threshold;
            
            % Get the displacement field for transformation
            D = TBS.cp2D(mp,fp,sz,transformationType);
            
            % Thresholding selective displacement field
            D = TBS.thresholdD(D,threshold);
            
            % Combine displacment field
            while numel(D) >= 2
                D{1} = TBS.Doperation(D{1},D{2});
                D(2) = [];
            end
            D = D{:};
        end
        
        %% Function:    cpselectIter
        % Discription:  control point selection with iteration
        function [mp,fp,D] = cpselectIter(moving,fixed,mp,fp,...
                transformationSetting,maxIteration)
            % Note, currently only support two step tform, can be expanded
            % Input:    moving/fixed, moving and fixed image
            %           transformationSetting, struct,
            %           maxIteration, max interaciton of point selection
            % Output:   mp/fp, n x 2 mat of coordinates
            %           D, mat, displacement field
            
            % Setting -----------------------------------------------------
            
            % Function handle to change to image class
            imageBitFh = class(fixed);
            imageBitFh = str2func(imageBitFh);
            
            sz = size(fixed);
            
            % Initial result ----------------------------------------------
            if isempty(mp) && isempty(fp)
                
                % Select moving and fix points
                [mp,fp] = cpselect(moving,fixed,'Wait',true);
            end
            
            D = [];
            try
                % Get displacement field, using selective method
                D = TBS.cp2Dn(mp,fp,sz,transformationSetting);
            end
            
            % Pick more points if there are not enough point number -------
            while isempty(D)
                disp('Point number are not enough for current method.')
                [mp,fp] = cpselect(moving,fixed,mp,fp,'Wait',true);
                
                try
                    % Get displacement field, using selective method
                    D = TBS.cp2Dn(mp,fp,sz,transformationSetting);
                end
            end
            
            % Transform image
            tfMoving = imwarp(moving,D);
            
            % Alignment iterations ----------------------------------------
            i = 1;
            while i <= maxIteration
                % transformed moving points
                tfMP = TBS.transformPointsForwardD(D,mp,50);
                
                % Delete NaN if there is any
                row = isnan(tfMP);
                row = ~any(row,2);
                tfMP = tfMP(row,:); fp = fp(row,:);
                
                % Overlap the annotation for better visualizaiton
                % Truecolor image (tfMoving-M, fixed-G)
                overlayTFmoving = cat(3,tfMoving,fixed,tfMoving);
                
                fixedGreen = zeros([size(fixed),3]);
                fixedGreen = imageBitFh(fixedGreen);
                fixedGreen(:,:,2) = fixed;
                % Select additional points basing on current selection
                [tfMP,fp] = cpselect(overlayTFmoving,fixedGreen,tfMP,fp,'Wait',true);
                
                % Get moving points using reverse transform ---------------
                mp = TBS.transformPointsInverseD(D,tfMP);
                
                % Update tform & transformed image ------------------------
                D = TBS.cp2Dn(mp,fp,sz,transformationSetting);
                tfMoving = imwarp(moving,D);
                
                figure; imshowpair(fixed,tfMoving);
                
                % Check with user whether this is good enough -------------
                answer = questdlg('Add more points','Log','Yes','No','Yes');
                if strcmp(answer,'No')
                    break
                end
                
                i = i+1;
            end
            
            close all;
        end
        
        %% Function:    SA_getImage4Alignment
        % Discription:  get images for alignment, intensity adjusted
        function im = SA_getImage4Alignment(directory,imName,...
                ch4Align,scaleFactor)
            % Input:    directory, str, directory of the image
            %           imName, str, image name
            %           ch4Align, mat/num, channels use for alignment
            %           scaleFactor, num, ratio for scale down
            % Output:   im, mat, scaled and intensity adjusted uint8 image
            
            cd(directory);
            idir = dir(['**',filesep,'*',imName,'*']);
            
            % Get images for alignment
            cd(idir.folder);
            im = TBS.getStack(idir.name,ch4Align);
            im = max(im,[],3);
            
            % Scale down for fast computation
            im = imresize(im,scaleFactor);
            
            % Subtract background
            minInten = prctile(nonzeros(im),5);
            im = im - minInten;
            
            % Convert to uint8 to speed up
            im = imadjust(im);
            im = im2uint8(im);
        end
        
        %% Function:    selfAlignment (main)
        % Discription:  alignment brain sections within a brain
        function selfAlignTform = selfAlignment(startSeq,selfAlignmentSetting)
            % Input:        startSeq, num, starting seq number
            %               outputStruct, struct, for output
            %               selfAlignmentSetting, struct, settings for the alignment
            % Output:       selfAlignTform, (the tform and D are scaled)
            
            % Settings
            % ch4Align, imageSetting, directory, sysSetting,
            % redoTF, refSectionNum, refSectionRotation
            % scaleFactor, alignmentSeq
            % Extract variables from the setting
            TBS.extractStructVar(selfAlignmentSetting,'caller');
            
            % Additional fixed settings
            [sectionName,sectionNum] = TBS.getSectionNameNum(...
                directory.abStitch,imageSetting,sysSetting);
            
            % Alignment seq (modify this part if its necessary)
            alignmentSeq = TBS.getAlignmentSeq(sectionNum,refSectionNum);
            
            % reference image rotation (do empty to skip)
            initialTform = eye(3);
            if ~isempty(refSectionRotation)
                initialTform = TBS.rotz(refSectionRotation);
            end
            
            % Scaling for speeding up
            scaleTform = TBS.getScaleTform(scaleFactor,2);
            scaleTform = affine2d(scaleTform);
            
            % Can start in different point in the alignmentSeq (default: 1)
            if isempty(startSeq)
                startSeq = 1;
            end
            
            % Whether skip the rest image
            skipTheRest = false;
            
            % Padding size ------------------------------------------------
            % Padding size doesnt affect the downstream result, just used
            % for visualizaiton during alignment
            
            % Ref image
            movingName = sectionName{sectionNum == refSectionNum};
            moving = TBS.SA_getImage4Alignment(directory.abStitch,...
                movingName,ch4Align,1);
            
            % 20% for each edge
            paddingSize = round(size(moving,[1 2]).*0.2);
            paddingTform = eye(3);
            paddingTform(3,[2 1]) = paddingSize./2;
            
            % Add to initial tform
            initialTform = inv(scaleTform.T)*initialTform*paddingTform*scaleTform.T;
            initialTform = affine2d(initialTform);
            
            % Output size
            sz = size(moving,[1 2]) + paddingSize;
            
            Rsz = sz*scaleFactor;
            Rsz = round(Rsz);
            
            % Load annotAlignTform ----------------------------------------
            cd(directory.main);
            if exist('selfAlignTform.mat')
                load('selfAlignTform.mat');
            else
                selfAlignTform = table('Size',[0 5],'VariableTypes',...
                    repmat({'cell'},1,5),'VariableNames',...
                    {'mp','D','fixName','fp','size'});
            end
            
            % Alignment ===================================================
            
            cd(directory.abStitch);
            for iIm = alignmentSeq(startSeq:end)
                tic
                % Get current image with scaling
                movingName = sectionName{sectionNum == iIm};
                moving = TBS.SA_getImage4Alignment(directory.abStitch,...
                    movingName,ch4Align,scaleFactor);
                
                % For reference image -------------------------------------
                if iIm == alignmentSeq(1)
                    % Add to output
                    selfAlignTform.D{movingName} = TBS.tform2D(initialTform,Rsz);
                    selfAlignTform.size{movingName} = sz;
                    continue
                end
                
                % Get fixed image -----------------------------------------
                % Find a fixed image for the current image or use the old in table
                if redoTF || ~TBS.hasThisRow(selfAlignTform,movingName)
                    for j = 1:10
                        % Get fix number, j digit before the current number
                        fixNum = find(alignmentSeq==iIm);
                        fixNum = alignmentSeq(fixNum-j);
                        
                        fixName = sectionName{sectionNum == fixNum};
                        if ~TBS.hasThisRow(selfAlignTform,fixName)
                            continue
                        end
                        
                        % Get fixed image
                        fix = TBS.SA_getImage4Alignment(...
                            directory.abStitch,fixName,ch4Align,scaleFactor);
                        
                        % Check with user to see whether this image is good to use
                        figure; imshow(fix);
                        answer = questdlg('This as fixed image?','Log','Yes','No','Yes');
                        if strcmp(answer,'Yes')
                            close all; break
                        end
                    end
                    
                else
                    fixName = selfAlignTform{movingName,'fixName'}{:};
                    fix = TBS.SA_getImage4Alignment(directory.abStitch,...
                        fixName,ch4Align,scaleFactor);
                end
                
                % Transform the fixed image
                fixD = selfAlignTform.D{fixName};
                
                % Transform fixed using tformCoarse/Fine
                fix = imwarp(fix,fixD);
                
                % Step 1: coarose alignment -------------------------------------------
                % mp: moving point; fp: fixed point
                
                % Whether use exist mp & fp (moving/fixed points)
                if ~redoTF && TBS.hasThisRow(selfAlignTform,movingName)
                    mp = selfAlignTform{movingName,'mp'}{:};
                    fp = selfAlignTform{movingName,'fp'}{:};
                    
                    % Scale to current scale
                    mp = transformPointsForward(scaleTform,mp);
                    fp = transformPointsForward(scaleTform,fp);
                    
                    % Transform the fixed using its tform
                    fp = TBS.transformPointsForwardD(fixD,fp,50);
                else
                    mp = []; fp = [];
                end
                
                % Select/modify points, or get displacement field from
                % skipped image
                if ~skipTheRest
                    % cpselectModify(moving,fixed,mp,fp,transformationType,maxIteration)
                    [mp,fp,D] = TBS.cpselectIter(moving,fix,mp,fp,transformationSetting,20);
                else
                    % Get displacement field from point pairs
                    D = TBS.cp2Dn(mp,fp,Rsz,transformationSetting);
                end
                
                % Transform to original scale -----------------------------
                % The fixed point here transformed to the original fixed
                fp = TBS.transformPointsInverseD(fixD,fp);
                
                % Transform fixed & moving to orginial scale
                fp = transformPointsInverse(scaleTform,fp);
                mp = transformPointsInverse(scaleTform,mp);
                
                % Update variables ----------------------------------------
                selfAlignTform(movingName,:)= table({mp},{D},...
                    {fixName},{fp},{sz});
                
                % return var to workspace
                assignin('base','selfAlignTform',selfAlignTform);
                
                disp(['Done: ',movingName]); close all; toc
                
                if skipTheRest
                    continue
                end
                
                answer = questdlg('Skip the rest images?','Log','Yes','No','No');
                if strcmp(answer,'Yes')
                    skipTheRest = true;
                end
            end
        end
        
        %% Function:    stac2vol (main)
        function stack2vol(tbl,ch,background,scaleFactor,imageSetting,sysSetting,directory)
            % Input:    tbl, table, image info with displacement field
            %           ch, mat, output channels
            %           background, mat, background intensity
            %           imageSetting/sysSetting/directory, struct
            % Output:   scalled image on disk, uint8
            
            % Settings
            outputName = [imageSetting.mouseID,'_',num2str(scaleFactor),'.tif'];
            
            background = reshape(background,1,1,[]);
            
            % Output image size
            sz = cellfun(@(X) size(X,1:2),tbl.D,'Uniformoutput',false);
            sz = vertcat(sz{:});
            sz = max(sz);
            
            % Sort transformation matrix table
            sectionNum = TBS.getSectionNumber(tbl,imageSetting,sysSetting);
            [sectionNum,I] = sort(sectionNum,'ascend');
            
            tbl = tbl(I,:);
            
            cd(directory.main);
            for i = min(sectionNum):max(sectionNum)
                row = sectionNum == i;
                
                % Allow empty seciton
                if ~any(row)
                    iIm = zeros(sz,numel(ch));
                    iIm = uint8(iIm);
                    
                    if i == 1
                        TBS.saveStack(iIm,outputName);
                    else
                        TBS.saveStack(iIm,outputName,true);
                    end
                    
                    continue
                end
                
                imName = tbl.Properties.RowNames{row};
                
                iIm = TBS.getStack(fullfile(directory.abStitch,imName),ch);
                
                % Scale image
                iIm = imresize(iIm,scaleFactor);
                
                % Background subtraction
                classFh = class(iIm);
                classFh = @(X) cast(X,classFh);
                
                iIm = iIm - classFh(background);
                % Change to uint8 after background subtraction
                iIm = uint8(iIm);
                
                % Transformation
                D = tbl.D{imName};
                iIm = imwarp(iIm,D);
                
                if i == min(sectionNum)
                    TBS.saveStack(iIm,outputName);
                else
                    TBS.saveStack(iIm,outputName,true);
                end
                
                disp(['Transform image: ',imName]);
            end
        end
        
    end
    
    methods (Static)    % ComebineBC ======================================
        %% Function:    isDegenerateBC
        % Discription:  whether a BC is degenerate BC
        function TF = isDegenerateBC(BC,n,tolerate0)
            % Input:    BC, barcode matrix
            %           n, number, digit for continous nt
            %           tolerate0, logical, whether tolerate 0
            % Output:   TF, logical, whether is degenerate BC
            
            % 09012021: max digit of continous BC
            % exlcuded digit need to +1
            n = n+1;
            
            nCh = max(BC,[],'all');
            
            % For even number, imerode doesnt work
            % when n is even, erode with n-1 strel object, and find the >2
            % continous true along column
            if mod(n,2) == 0
                n2 = n-1;
            else
                n2 = n;
            end
            
            SE = strel('line',n2,0);
            
            % Padding false for edge effect, size of padding
            padding = false(size(BC,1),ceil(n2/2));
            
            % Whether it is degenerate BC for each channel
            TF = false(size(BC,1),1);
            for iCh = 1:nCh
                row = BC == iCh;
                
                if tolerate0
                    row = row | BC == 0;
                end
                
                row = [padding,row,padding];
                
                % Have continous N nucleotides
                row = imerode(row,SE);
                
                % Two continous true for even number nt
                if mod(n,2) == 0
                    row = row(:,1:end-1) & row(:,2:end);
                end
                
                row = any(row,2);
                
                TF = TF | row;
            end
        end
        
        %% Function:    estimateDegenerateBC
        % Discription:  model possibility for degenerate BC
        function p = estimateDegenerateBC(nIteration,nBC,nSeq,n)
            % Input:    nInteration, num, number of interation
            %           nBC, number of BC per iteraiton
            %           nSeq, number of BC length
            %           excludeSeq, seq exclude from calculation
            %           n, number of continous same nt
            % Output:   p, mat, possibility of degenate BC
            %           one iteration per row; one n per column
            
            p = [];
            parfor i = 1:nIteration % parfor
                
                BC = TBS.randBC(nBC,nSeq,[],[]);
                
                % Fraction of degenerate BC for different length
                % isDegenerateBC(BC,n,tolerate0)
                p2 = arrayfun(@(X) TBS.isDegenerateBC(BC,X,true),...
                    n,'UniformOutput',false);
                p2 = cellfun(@(X) sum(X)/size(X,1),p2);
                
                p(i,:) = p2;
                disp(['Function estimateDegenerateBC iteration: ',...
                    num2str(i)]);
            end
        end
        
        %% Function:    hammingDist2
        % Discription:  calculate hamming distance to BC1 of BC2
        function D = hammingDist2(BC1,BC2,tolerate0)
            % Note: 08042021 no row sum of BC2
            % tried split into blocks, not faster
            % Input:    BC1/2, mat, BC matrix
            %           tolerate0, logical
            % Output:   D, vector, count of hamming distance == n, from 0
            % -> nSeq. Itself is excluded if its self-self distance
            
            % For calculating self-self distance
            if isempty(BC2)
                isSelf = true;
                BC2 = BC1;
            elseif all(size(BC1) == size(BC2)) && all(BC1 == BC2,'all')
                isSelf = true;
            else
                isSelf = false;
            end
            
            % (to speed up)
            BC1zero = BC1 == 0;
            
            % Output size for each BC: nSeq+1
            sz = [size(BC1,2)+1,1];
            
            D = [];
            parfor i = 1:size(BC2,1) % parfor
                iBC2 = BC2(i,:);
                
                % Whether BC are the same
                row = BC1 == iBC2;
                
                if tolerate0
                    row = row | BC1zero | iBC2 == 0;
                end
                
                % Distance
                iD = sum(~row,2);
                
                % Set min from 0 to 1 by adding 1
                iD = iD + 1;
                
                % Counts for each hamming dist from 0 to nSeq
                iD = accumarray(iD,1,sz);
                
                D(:,i) = iD;
            end
            
            % Delete self-self comparison
            if isSelf
                D(1,:) = D(1,:)-1;
            end
        end
        
        %% Function:    randBC
        % Discription:  get random unique BC
        function bc = randBC(nBC,nSeq,n,maxHamming)
            % Input:    nBC, num, number of barcode
            %           nSeq, num, number of nt per barcode
            %           n, number, digit of degenerated barcode
            %           maxHamming,num,
            % Output:   bc, barcode matrix
            
            
            nBC2 = nBC * 1.5;
            
            % Random unique BC matrix
            bc = randi(4,nBC2,nSeq);
            bc = uint8(bc);
            % Bug fix 09162021
            bc = unique(bc,'stable','row');
            
            % Exclude degenerate BC
            if ~isempty(n)
                TF = TBS.isDegenerateBC(bc,n,true);
                bc = bc(~TF,:);
            end
            
            % Exclude BC have only certain channels
            if ~isempty(maxHamming)
                % Ch 1 & 2
                TF = ismember(bc,[0 1 2]);
                TF = sum(~TF,2) > maxHamming;
                bc = bc(TF,:);
                
                % Ch3 & 4
                TF = ismember(bc,[0 3 4]);
                TF = sum(~TF,2) > maxHamming;
                bc = bc(TF,:);
            end
                        
            % Get the required number of BC
            bc = bc(1:nBC,:);
        end
        
        %% Function:    estimateHammingDist
        % Discription:  model posibility for hamming distance
        function D = estimateHammingDist(nIteration,method,nBC,nSeq,n,maxHamming)
            % Input:    nInteration, num, number of interation
            %           fh, function handel/method for each stimulation
            %           nBC, number of BC per iteraiton
            %           nSeq, number of BC length
            %           n, number of continous same nt
            % Output:   p, mat, possibility of hamming distance
            %           one iteration per column; row: 0-nSeq
            
            disp('Function estimateHammingDist on progress...');
            
            D = [];
            for i = 1:nIteration
                
                BC = TBS.randBC(nBC,nSeq,n,maxHamming);
                
                iD = TBS.hammingDist2(BC,[],true);
                
                switch method
                    case 'min'                        
                        iD = (0:nSeq)'.*(iD > 0);
                        iD(iD == 0) = nan;
                        iD = min(iD,[],1,'omitnan');        
                        iD = accumarray(iD',1,[nSeq,1]);
                end
                                
                D(:,i) = iD;
            end
        end
        
        %% Function:    getCombineVar
        % Discription:  combine all variable under the directory/sub-folder
        function varOut = getCombineVar(directory,fileName)
            % Input:    directory, str
            %           fileName, str
            % Output:   varOut, combined variable
            
            fileName = erase(fileName,'.mat');
            
            cd(directory);
            fileDir = dir(['**\',fileName,'.mat']);
            
            varOut = [];
            for i = 1:size(fileDir,1)
                iDirectory = fileDir(i).folder;
                iName = fileDir(i).name;
                
                iFile = load(fullfile(iDirectory,iName));
                iFile = iFile.(fileName);
                
                if isempty(iFile)
                    continue;
                end
                
                if isempty(varOut)
                    varOut = iFile;
                else
                    varOut = [varOut; iFile];
                end
            end
            
            disp('Function getCombineVar: done.');
        end
        
        %% Function:    loadBCvar
        % Discription:  load axon/soma table form all experiments
        function varOut = loadBCvar(str,bcSetting)
            % Input:    str, name of the variable
            % Output:   varOut, table
            
            directory = evalin('base','directory');
            directory = directory.main;
            
            varOut = TBS.getCombineVar(directory,str);
            
            varOut = TBS.delPureChDotInTable(varOut,bcSetting);
            
            % Delete barcode doesnt have enough length --------------------
            % Added 08252021
            excludeSeq = bcSetting.seqNotInMinBscallCycles;
            minBClen = bcSetting.minBClen;
            
            for i = 1:size(varOut,1)
                TF = varOut.bscallCh{i};
                TF = TF ~= 0;
                
                TF = TBS.BClongEngouth(TF,excludeSeq,minBClen);
                
                varOut{i,:} = cellfun(@(X) X(TF,:),varOut{i,:},'Uniformoutput',false);
            end
            
        end
                
        %% Function:    delNearSomaRolony
        % Discription:  delete rolony near soma
        function axonBC = delNearSomaRolony(axonBC,somaBCVar,bcSetting,imageSetting)
            % Input & output: axonBC, table, axon bscall result
            %               somaBCVar, table
            %               bcSetting/imageSetting,struct
            
            % Setings
            % min counts for soma pixel
            minPixelCount = bcSetting.nearSomaRolony.minPixelCount;
            
            % min pixel distance to soma pixel
            minDistance = bcSetting.nearSomaRolony.minDistance;
            minDistance = minDistance./imageSetting.resolution;
            minDistance = round(minDistance);
            
            % max hamming distance
            maxHamming = bcSetting.maxHamming;
            
            disp('Function delNearSomaRolony on progress...')
            
            % Functional handle for getting xy from table
            fh = @(X,Y,Z) [varXY(X,Y,'x',Z),varXY(X,Y,'y',Z)];
            
            % Strel for conv2
            SE = strel('disk',minDistance);
            SE = SE.Neighborhood;
            
            for i = 1:size(somaBCVar,1)
                
                imName = somaBCVar.Properties.RowNames{i};
                % Number of col/seq
                % varXY(tbl,rowName,varName,col)
                nIter = size(varXY(axonBC,imName,'x'),2);
                
                if ~nIter
                    continue
                end
                
                % Soma xy coordinates
                somaXY = fh(somaBCVar,imName,1);
                
                % Get image of somaXY
                % Do conv2 to get the sum of nearby pixel
                sz = max(somaXY);
                % Extra area for SE
                sz = fliplr(sz)+ minDistance + 3;
                
                % Count of soma pixel within the range
                somaTF = TBS.xyzv2im(sz,somaXY,[]);
                somaTF = conv2(somaTF,SE,'same');
                % is near the soma, not including minmum
                somaTF = somaTF > minPixelCount;
                
                % Whether the rolony within the range across sequencing cycles
                % Count how many times it within the range of soma
                TF = [];
                for j = 1:nIter % parfor
                    % Rolony xy of the current cycle
                    axonXY = fh(axonBC,imName,j);
                    axonXY = round(axonXY);
                    
                    % Rolony within the image
                    iTF = TBS.isOutOfRngeDot(sz,axonXY);
                    iTF = ~iTF;
                    
                    % Nearby soma pixel number (intensity)
                    ind = sub2ind(sz,axonXY(iTF,2),axonXY(iTF,1));
                    ind = somaTF(ind);
                    
                    iTF(iTF) = ind;
                    
                    TF(:,j) = iTF;
                end
                
                % For multiple column, sum of cycles close to the soma 
                if size(TF,2)> 1                   
                    % sum of cycles near the soma
                    TF = sum(TF,2);
                    TF = TF > maxHamming;
                end
                
                % Delete the rows in axonBC
                axonBC{imName,:} = cellfun(@(X) X(~TF,:),...
                    axonBC{imName,:},'UniformOutput',false);
            end
            
            % Function: varOut --------------------------------------------
            function varOut = varXY(tbl,rowName,varName,col)
                if ismember(varName,tbl.Properties.VariableNames)
                    varOut = tbl.(varName){rowName};
                else
                    varOut = tbl{rowName,1}{:}.(varName);
                end
                
                if nargin == 4 && ~isempty(col)
                    varOut = varOut(:,col);
                end                
            end
            
        end
        
        %% Function:    registerXY2Vol
        % Discription:  register rolony/soma xy to volume space
        function varOut = registerXY2Vol(varIn,scaleFactor,sysSetting,directory)
            
            disp('Function registerXY2Vol on progress...');
            
            % Settings
            cd(directory.main)
            corrSeq2abTform = load('corrSeq2abTform.mat');
            corrSelfAlignTform = load('corrSelfAlignTform.mat');
            corrSeq2abTform = corrSeq2abTform.corrSeq2abTform;
            corrSelfAlignTform = corrSelfAlignTform.corrSelfAlignTform;
            
            % Ensure both table has same format of row names
            imFormat = sysSetting.imFormat;
            varIn = TBS.ensureRowNameAppend(varIn,imFormat);
            corrSeq2abTform = TBS.ensureRowNameAppend(corrSeq2abTform,imFormat);
            
            % Scaling factor for brain seciton to volume
            scaleTform = TBS.getScaleTform(scaleFactor,2);
            scaleTform = affine2d(scaleTform);
            
            varOut = {};
            parfor i = 1:size(varIn,1) % parfor
                imName = varIn.Properties.RowNames{i};
                
                if ~TBS.hasThisRow(corrSeq2abTform,imName) ||...
                        isempty(varIn{imName,1}{:})
                    continue
                end
                
                % tform: seq image to whole brain section
                iSeq2AbTform = corrSeq2abTform.tform{imName};
                if ~isobject(iSeq2AbTform)
                    iSeq2AbTform = projective2d(iSeq2AbTform);
                end
                
                % Ab image name
                iAbName = corrSeq2abTform.fixName{imName};
                
                if ~TBS.hasThisRow(corrSelfAlignTform,iAbName)
                    continue
                end
                
                % Displacement field for self-alignment
                iAbD = corrSelfAlignTform.D{iAbName};
                
                % BC coordinates
                x = varIn.x{imName};
                y = varIn.y{imName};
                
                % Collapse to single column
                x = median(x,2);
                y = median(y,2);
                xy = [x, y];
                
                % Transform 1: region to section --------------------------
                xy = transformPointsForward(iSeq2AbTform,xy);
                
                % Transform 2: section to brain volume --------------------
                % Scale down for displacement field
                xy = transformPointsForward(scaleTform,xy);
                
                % Transform to brain volumn space
                classFh = class(xy);
                classFh = @(X) cast(X,classFh);
                
                xy = double(xy);
                xy = TBS.transformPointsForwardD(iAbD,xy,50);
                
                % Scale back
                xy = transformPointsInverse(scaleTform,xy);
                
                xy = classFh(xy);
                varOut(i,:) = [{imName},{xy(:,1)},{xy(:,2)}];
            end
            
            % Delete empty rows
            row = cellfun(@isempty,varOut(:,1));
            varOut = varOut(~row,:);
            
            % Delete non-registrated rows
            % To use varIn for output to keep other variables
            rowNames = varOut(:,1);
            row = ismember(varIn.Properties.RowNames,rowNames);
            varIn = varIn(row,:);
            
            % Update xy in table
            varIn.xReg(rowNames) = varOut(:,2);
            varIn.yReg(rowNames) = varOut(:,3);
            varOut = varIn;
        end
        
        %% Function:    vol2micron
        % Discription:  change coordinates from volume space to micron
        function tbl = vol2micron(tbl,imageSetting,sysSetting)
            % Input & output: tbl, table, with coordinates
            %           imageSetting/sysSetting,struct
            
            slideThickness = imageSetting.slideThickness;
            resolution = imageSetting.resolution;
            
            % Get section number via image name
            secitonName = tbl.Properties.RowNames;
            sectionNumber = TBS.getSectionNumber(secitonName,...
                imageSetting,sysSetting);
            
            % z in micron
            z = (sectionNumber-1).*slideThickness+1;
            
            % Collapse to single column
            x = tbl.xReg; y = tbl.yReg;
            x = cellfun(@(X) median(X,2),x,'Uniformoutput',false);
            y = cellfun(@(X) median(X,2),y,'Uniformoutput',false);
            
            % xy in micron
            xy = cellfun(@(X,Y) [X,Y],x,y,'Uniformoutput',false);
            xy = cellfun(@(X) X.*resolution,xy,'UniformOutput',false);
            
            n = cellfun(@(X) size(X,1),xy);
            z = arrayfun(@(X,Y) repmat(X,Y,1),z,n,'UniformOutput',false);
            
            xyz = cellfun(@(X,Y) [X,Y],xy,z,'Uniformoutput',false);
            
            tbl.xyz = xyz;
        end
                
        %% Function:    vol2reg
        % Discription:  change coordinates from volume space to reference
        % map
        function tbl = vol2reg(tbl,reg2AllenTform,scaleFactor,imageSetting,sysSetting)
            % Note: resolution is the same as the reference map
            % Input & output: tbl, table, with coordinates
            %           reg2AllenTform, struct, info register to reference
            %           map like Allen 
            %           scaleFactor, num, scaling factor from original
            %           image to volume
            %           imageSetting/sysSetting,struct
            
            disp('Function vol2reg on progress...');
            
            % Tform from aligned 3D volumn to Allen reference 
            % (i.e. 25 um resolution)
            tform = reg2AllenTform.tform;
            % Displacement field after tform
            D = reg2AllenTform.D;
            
            % Scaling factor for brain seciton to volume
            scaleTform = TBS.getScaleTform(scaleFactor,2);
            scaleTform = affine2d(scaleTform);
            
            % Get section number via image name
            secitonName = tbl.Properties.RowNames;
            sectionNumber = TBS.getSectionNumber(secitonName,...
                imageSetting,sysSetting);
            
            % x/yReg is in original image resolution (i.e. 0.55 um/pixel)
            % Collapse to single column
            x = tbl.xReg; y = tbl.yReg;
            x = cellfun(@(X) median(X,2),x,'Uniformoutput',false);
            y = cellfun(@(X) median(X,2),y,'Uniformoutput',false);
            
            % volume use section number as z
            % (sz will be used for mat2cell later)
            sz = cellfun(@(X) size(X,1),x);
            z = arrayfun(@(X,Y) repmat(X,Y,1),sectionNumber,sz,'UniformOutput',false);
            
            xyz = [vertcat(x{:}),vertcat(y{:}),vertcat(z{:})];
            
            % Change xy to volume rosolution
            xyz(:,1:2) = transformPointsForward(scaleTform,xyz(:,1:2));
            
            % Transform coordinates from volume to reference map ----------
            xyz = transformPointsForward(tform,xyz);
            
            % Get displacement field from each coordinates use
            % interpolation
            xyz2 = [];
            for i = 1:size(xyz,2)
                xyz2(:,i) = interp3(D(:,:,:,i),xyz(:,1),xyz(:,2),...
                    xyz(:,3),'spline');
            end
            
            % Transform using displacement field
            xyz = xyz - xyz2;
            
            xyz = mat2cell(xyz,sz,3);
            
            % Reference xyz
            tbl.xyzRef = xyz;
        end
              
        %% Function:    mergeAxonBC 
        % Discription:  merge codeBook BC using axon count
        function TF = mergeAxonBC(codeBook,axonBCn,mergeHamming)
            % Input:    codeBook, mat, template barcode
            %           axonBCn, table, axonBC count
            %           mergeHamming, num, max hamming distance for merging
            % Output:   TF, logical, included BC after merge
            
            disp('Function mergeAxonBC on progress...');
                                    
            % (To speed up)
            codeBookZero = codeBook == 0;
            
            % Identify similar BC
            similarBC = {};
            parfor i = 1:size(codeBook,1) % parfor
                iBC = codeBook(i,:);
                
                % BC within hamming distance
                % Tolerate to 0 during the merge
                hammingDist = codeBook == iBC | codeBookZero | iBC == 0;
                hammingDist = sum(~hammingDist,2);
                
                row = hammingDist <= mergeHamming;
                similarBC{i,1} = find(row);
            end
            
            % Axonal pixel count for potential BC
            n = cellfun(@(X) axonBCn(X),similarBC,'UniformOutput',false);
            
            % Find rows equal to max count
            maxN = cellfun(@max,n);
            
            TF = axonBCn == maxN;
        end        
        
        %% Function:    updateCodeID
        % Discription:  update codeID in axon/somaBC
        function tblBC = updateCodeID(tblBC,TF)
            % Input & output:    tblBC, table, with codeID
            %           TF, logcial, whether the BC presented in the new
            %           codebook, size as the old codebook
            
            newId = find(TF);
            
            id = tblBC.codeID;
            [~,id] = cellfun(@(X) ismember(X,newId),id,'UniformOutput',false);
            
            tblBC.codeID = id;
            
            % Delete the excluded BC rows
            for i = 1:size(tblBC,1)
                TF = id{i}~=0;
                tblBC{i,:} = cellfun(@(X) X(TF,:),tblBC{i,:},'UniformOutput',false);
            end            
        end
        
        %% Function:    codeLookupIdem
        % Discription:  get look up table with perfect match with 0
        % tolerance (faster using ismember)
        function [BCq,templateBC,unmatchBC] = codeLookupIdem(BCin,codeBook)
            % Input:    BCin, mat, one row per BC
            %           codeBook, mat, one row per BC
            % Output:   BCq, mat, BC of inquery, matched BC in BCin
            %           templateBC, mat, mached codebook BC cooresponding
            %           to BCq
            %           unmatchBC, mat, BC without hammingDist=0 match
            
            % (To speed up) Perfect match with or without 0
            [Lia,Lcob] = ismember(BCin,codeBook,'rows');
            BCq0 = BCin(Lia,:);
            templateBC0 = codeBook(Lcob(Lia),:);
            
            % Delete the perfectly matched rows
            BCin = BCin(~Lia,:);
            
            % Find complete match with 0 ===================================
            
            % 0-pattern/position in BC, for BCin and codebook
            pos0i = BCin == 0;
            [pos0i,~,ic] = unique(pos0i,'row');
            
            pos0j = codeBook == 0;
            [pos0j,~,jc] = unique(pos0j,'row');
            
            % Loop through 0-pattern in BCin (i) and codebook (j)
            BCq = {}; codeID = {};
            parfor i = 1:size(pos0i,1) % parfor
                % 0 columns in BCinput
                iPos0 = pos0i(i,:);
                
                % Row with the 0 pattern
                rowi = ic == i;
                rowi = find(rowi);
                iBCin = BCin(rowi,:);
                
                iBCq = {}; iCodeID = {}; 
                for j = 1:size(pos0j,1)
                    
                    % Non-0 position in both matrix
                    jPos0 = pos0j(j,:);
                    jPos0 = ~iPos0 & ~jPos0;
                    
                    % Row with the 0 pattern
                    rowj = jc == j;
                    rowj = find(rowj);
                    jCodeBook = codeBook(rowj,:); 
                    
                    % Input with with overlap columns
                    jBCin = iBCin(:,jPos0);
                    jCodeBook = jCodeBook(:,jPos0);
                    
                    % Prefect match codebook to BCin
                    [Lia,Lcob] = ismember(jCodeBook,jBCin,'rows');
                    
                    if ~any(Lia)
                        continue
                    end
                    
                    % Matched codebook
                    iCodeID{j,1} = rowj(Lia);
                    
                    % Matched BCin
                    iBCq{j,1} = rowi(Lcob(Lia));
                end
 
                iBCq = vertcat(iBCq{:});
                iCodeID = vertcat(iCodeID{:});
                
                BCq{i,1} = iBCq;
                codeID{i,1} = iCodeID;
            end
            
            BCq = vertcat(BCq{:});
            codeID = vertcat(codeID{:});
            
            % Umatched BC, exclude matched rows (single & multiple match)
            unmatchBC = 1:size(BCin,1);
            unmatchBC = ~ismember(unmatchBC,BCq);
            unmatchBC = BCin(unmatchBC,:);
            
            % Find single matched BCin, and cooresponding codeBook
            [~,~,ic] = unique(BCq);
            n = accumarray(ic,1);
            row = find(n == 1);
            row = ismember(ic,row);
            
            BCq = BCq(row);
            codeID = codeID(row);
            
            BCq = BCin(BCq,:);
            templateBC = codeBook(codeID,:);
            
            % Combine with tthe initial perfect match
            BCq = [BCq0; BCq];
            templateBC = [templateBC0; templateBC];
            
            disp('Got identical BC lookup matrix.')
        end
        
        %% Function:    codeLookup
        % Discription:  Lookup table between input BC and code/template BC
        function lookupTbl = codeLookup(BCin,codeBook,maxHamming,...
                minBCcount,tolerate0)
            % 08252021, added seperate method for matching identical BC
            % Input:    BCin, mat/cell, input BC
            %           codeBook, mat, one template BC per row
            %           maxHamming, num, maximum hamming distance
            %           minBCcount, num, min BC count for input BC
            %           tolerate0, logical, whether tolerate 0 in
            %           non-identical hamming distance calculation
            % Output:   lookupTbl, table, look up table between input BC
            % and template BC
            
            disp('Function codeLookup on progress...');
            
            if iscell(BCin)
                BCin = vertcat(BCin{:});
            end
            
            % (To speed up) exlude BC with less than min BC count ---------
            % (Its not practical to include all 1 counts)
            [BCin,~,ic] = unique(BCin,'rows');
            
            % BC count
            n = accumarray(ic,1);
            
            % Count threshold
            row = n >= minBCcount;
            BCin = BCin(row,:);
            
            % Match identical BC, with 0 tolerance ------------------------

            % unmatchBC, exclude single & multiple matched BC
            % [BCq,templateBC,unmatchBC] = codeLookupIdem(BCin,codeBook)
            [BCq,templateBC,BCin] = TBS.codeLookupIdem(BCin,codeBook);
            
            if maxHamming == 0               
                lookupTbl = table();
                lookupTbl.templateBC = templateBC;
                lookupTbl.BCq = BCq;
                return
            end
                        
            % Match non-identical BC --------------------------------------
            
            % (To speed up)
            codeBookZero = codeBook == 0;
            
            codeID2 = {}; tic
            parfor i = 1:size(BCin,1) % parfor
                iBC = BCin(i,:);
                
                % Hamming distance                
                % 08312021, add option of 0-tolerance
                if tolerate0
                    hammingDist = codeBook == iBC | codeBookZero | iBC == 0;
                else
                    hammingDist = codeBook == iBC;
                end
                
                hammingDist = sum(~hammingDist,2);
                
                minHamming = min(hammingDist);
                
                if minHamming > maxHamming
                    continue
                end
                
                % Row smaller than hamming distance
                row = hammingDist == minHamming;
                
                codeID2{i,1} = find(row);
            end            
            
            % Only include rows with one match
            n = cellfun(@numel,codeID2);
            row = n == 1;
            BCq2 = BCin(row,:);
            codeID2 = codeID2(row,:);
            codeID2 = cell2mat(codeID2);
            templateBC2 = codeBook(codeID2,:);
            
            % Output table
            lookupTbl = table();
            lookupTbl.templateBC = [templateBC; templateBC2];
            lookupTbl.BCq = [BCq; BCq2];
            
            disp(['Got BC lookup table: ',num2str(toc)]); 
        end
        
        %% Function:    findCode
        % Discription:  find code id using lookup table
        function BCtable = findCode(BCtable,lookupTbl,codebook)
            % Input & output:  BCtable, table, BC table with bscallCh
            %           lookupTbl, table, lookup table for input BC &
            %           template BC
            %           codeBook, mat,
            
            BCq = vertcat(BCtable.bscallCh{:});
            
            templateBC2 = lookupTbl.templateBC;
            BCq2 = lookupTbl.BCq;
            
            % templateBC id
            [TF,Locb] = ismember(BCq,BCq2,'rows');
            [~,Locb2] = ismember(templateBC2,codebook,'rows');
            Locb(TF) = Locb2(Locb(TF));
            
            % Add to the table
            sz = cellfun(@(X) size(X,1),BCtable.bscallCh);
            BCtable.codeID = mat2cell(Locb,sz,1);
            
            % Delete rows not matched to template BC ----------------------
            TF = cellfun(@(X) X~=0, BCtable.codeID,'UniformOutput',false);
            
            for i = 1:size(BCtable,2)
                BCtable{:,i} = cellfun(@(X,Y) X(Y,:),BCtable{:,i},TF,...
                    'Uniformoutput',false);
            end
            
            disp('Done: findCode');
        end
        
        %% Function:    delMultiBCdot
        % Discription:  delete axon dots link to multiple BC
        function tbl = delMultiBCdot(tbl,bcSetting)
            % Input & output: tbl, table, axonBC result of single image
            %           dotMatchingSetting, struct
            
            excludeSeq = bcSetting.seqNotInMinBscallCycles;
            minBClen = bcSetting.minBClen;
            maxSeqInterval = bcSetting.maxSeqInterval;
            
            % Find dotId link to > 1 templateBC ---------------------------
            dotId = tbl.id{:};
            bcId = tbl.codeID{:};
            
            TF = false(size(dotId));
            for i = 1:size(dotId,2)
                
                iDotId = dotId(:,i);
                
                % Unique pair of dotID-templateBCID
                C = [iDotId,bcId];
                C = unique(C,'rows');
                
                % Dot ID has > 1 templateBC
                C = nonzeros(C(:,1));
                [C,~,ic] = unique(C);
                n = accumarray(ic,1);
                C = C(n > 1);
                
                TF(:,i) = ismember(iDotId,C);
            end
            
            % Delete the dot in results -----------------------------------
            % Set bscallCh,x,y, id to 0
            for i = 1:size(tbl,2)
                if size(tbl{1,i}{:},2) == size(TF,2)
                    tbl{1,i}{:}(TF) = 0;
                end
            end
            
            % Trim rows for identical dotIDs ------------------------------
            [~,ia,~] = unique(tbl.id{:},'rows');
            tbl{1,:} = cellfun(@(X) X(ia,:),tbl{1,:},'UniformOutput',false);
            
            % BC is not long enough ---------------------------------------
            TF = tbl.bscallCh{:} ~= 0;
            TF = TBS.BClongEngouth(TF,excludeSeq,minBClen);
            
            % BC has too big intervel -------------------------------------
            TF2 = TBS.blwMaxInterval(tbl.bscallCh{:},maxSeqInterval);
            
            % Rows fits interval and length requirements
            TF = TF & TF2;
            
            % Delete the rows for all variable
            tbl{1,:} = cellfun(@(X) X(TF,:),tbl{1,:},'UniformOutput',false);
        end
        
        %% Function:    trimRepeatID
        % Discription:  trim repeated dots/ID
        function tbl = trimRepeatID(tbl,bcSetting)
            % Input & output: table, axonBC result
            
            % max same dotID for two different dots
            maxDiffDotId = bcSetting.maxDiffDotId;
            
            dotId = tbl.id{:};
            nSeq = size(dotId,2);
            
            % Sort using bscall number: high to low
            % no basecall: dotId = 0
            n = sum(dotId~=0,2);
            [~,I] = sort(n,'descend');
            BC1 = dotId(I,:);
            
            % 08132021:Change from unique to hamming distance
            dotId2 = [];  hammingMat = false;
            while ~isempty(hammingMat)
                
                % Hamming distance matrix ---------------------------------
                % Not tolerate 0 besides itself
                BC2 = permute(BC1,[3 2 1]);
                hammingMat = BC1 == BC2 & BC1 ~= 0;
                hammingMat = sum(hammingMat,2);
                hammingMat = squeeze(hammingMat);
                
                % Set itself to nSeq
                eyeTF = logical(eye(size(hammingMat)));
                hammingMat(eyeTF) = nSeq;
                
                % Find the first hamming distance -------------------------
                % Faster use loop than unique & sortrows
                hamming1 = zeros(size(hammingMat,1),1);
                for i = 1:size(hammingMat,2)
                    [~,~,v] = find(hammingMat(:,i),1,'first');
                    hamming1(i) = v;
                end
                
                row = hamming1 == nSeq;
                dotId2 = [dotId2; BC1(row,:)];
                
                % Only keep rows with first hamming distance within range
                row = hamming1 <= maxDiffDotId;
                BC1 = BC1(row,:);
            end
            
            % Get the remaining var in the table -------------------------
            TF = ismember(dotId,dotId2,'rows');
            
            tbl{1,:} = cellfun(@(X) X(TF,:),tbl{1,:},'UniformOutput',false);
        end
        
        %% Function:    correctAxonBC
        % Discription:  correct and trim axonBC table
        function axonBC = correctAxonBC(axonBC,bcSetting)
            
            stat = sum(cellfun(@(X) size(X,1),axonBC.bscallCh));
            
            % Delete dots link to multiple BC
            for i = 1:size(axonBC,1)
                axonBC(i,:)= TBS.delMultiBCdot(axonBC(i,:),bcSetting);
            end
            
            % Trim dots link to same BC
            axonBC2 = {};
            parfor i = 1:size(axonBC,1) % parfor
                axonBC2{i,1} = TBS.trimRepeatID(axonBC(i,:),bcSetting);
            end

            axonBC = vertcat(axonBC2{:});
            
            % Delete empty image
            TF = cellfun(@isempty,axonBC.bscallCh);
            axonBC = axonBC(~TF,:);
            
            stat(2) = sum(cellfun(@(X) size(X,1),axonBC.bscallCh));
            stat = stat(2)/stat(1);
            disp(['AxonBC after correction: ',num2str(stat)]);
        end
        
        %% Function:    countMismatch
        % Discription:  count mismatch axon/somaBC for individual BC
        function stat = countMismatch(tbl,codeBook)
            % Note: after axonBC correction, multipleBC dot were set to 0,
            % so those nt won't be counted as mismatch
            % Input:    tbl, table, axon/somaBC
            %           codeBook, mat, one row per BC; col, seq cycles
            % Output:   stat, mat, mismatch counts; one BC per col; row:
            % 0-nSeq mismatch
            
            % code ID
            id = vertcat(tbl.codeID{:});
            % basecall channels
            bscallCh = vertcat(tbl.bscallCh{:});
            
            nSeq = size(bscallCh,2);
            
            stat = [];
            parfor i = 1:size(codeBook,1) % parfor
                iBC = codeBook(i,:);
                
                row = id == i;
                iBscallCh = bscallCh(row,:);
                
                % Number of mismatch
                D = iBscallCh == iBC | iBscallCh == 0 | iBC == 0;
                D = sum(~D,2);
                D = D + 1;
                
                % Count mismatch
                D = accumarray(D,1,[nSeq+1,1]);
                stat(:,i) = D;
            end
        end
        
        %% Function:    BCcountFilter 
        % Discription:  axon/soma coutn filter for template BC
        function TF = BCcountFilter(codebook,somaBC,axonBC,bcSetting)
            % Input:    codeBook, mat, one BC per row
            %           axonBC/somaBC, table, with bscallCh of axon/soma
            %           bcSetting, struct
            % Output:   TF,logical, BC pass the filter
            
            TBS.extractStructVar(bcSetting,'caller');
            
            TF = true(size(codebook,1),1);
            
            % Function handle for BC count
            fh = @(X) accumarray(vertcat(X{:}),1);
            
            % Trim out BC with too few axon count ------------------------
            n = fh(axonBC.codeID);
            if ~isempty(minAxonCount)
                row = n < minAxonCount;
                TF(row,:) = false;
            end
            
            if ~isempty(maxAxonCount)
                row = n > maxAxonCount;
                TF(row,:) = false;
            end
            
            % Trim out BC with too much soma count ------------------------
            n = fh(somaBC.codeID);
            
            if ~isempty(maxSomaCount)
                row = n < minSomaCount;
                TF(row,:) = false;
            end
            
            if ~isempty(maxSomaCount)
                row = n > maxSomaCount;
                TF(row,:) = false;
            end
                        
            disp(['Codebook count filter done: ',num2str(sum(TF)), ...
                ' of total ',num2str(size(TF,1))])
        end       
            
        %% Function:    BCregionCountFilter
        % Discription:  filter of axonBC count per region
        function TF = BCregionCountFilter(axonBC,regionMinCount,sysSetting)
            % Input:    axonBC, table,
            %           regionMinCount, num, minimum barcode count for the strongest
            %           projection region
            %           sysSetting, struct
            % Output:   TF, logical, hether pass the filter
            
            % Total BC count
            n = vertcat(axonBC.codeID{:});
            n = max(n);
            
            % Regions
            regionName = axonBC.Properties.RowNames;
            regionName = cellfun(@(X) TBS.nameFun(X,2,sysSetting),regionName,'UniformOutput',false);
            [regionName,~,ic] = unique(regionName);
            
            % AxonBC count for each image
            axonN = cellfun(@(X) accumarray(X,1,[n 1]),axonBC.codeID,'UniformOutput',false);
            
            % AxonBC count for each region
            axonN2 = [];
            for i = 1:numel(regionName)
                row = ic == i;
                iN = axonN(row);
                iN = horzcat(iN{:});
                axonN2(:,i) = sum(iN,2);
            end
            
            % Thresholding
            TF = axonN2 >= regionMinCount;
            TF = any(TF,2);
            
             disp(['BCregionCountFilter filter done: ',num2str(sum(TF)), ...
                ' of total ',num2str(size(TF,1))])
        end              
                
        %% Function:    isSameRolonyDiffIm
        % Discription:  find the same rolony in different images
        function TF = isSameRolonyDiffIm(axonBCreg,threshold)
            % Input:    axonBCreg, table, with register coordinates
            %           threshold, num, threshold to tell whether its the
            %           same rolony, in micron
            % Output:  TF, logical, row to be delete, same rolony in later
            % image in the table
            
            sz = cellfun(@(X) size(X,1),axonBCreg.codeID);
            % Assign each image an number for differenciate whehter the
            % same rolony apread twice
            im = arrayfun(@(X,Y) repmat(X,Y,1),(1:numel(sz))',...
                sz,'UniformOutput',false);
            
            id = vertcat(axonBCreg.codeID{:});
            xyz = vertcat(axonBCreg.xyz{:});
            im = vertcat(im{:});
            
            TF = false(size(id));
            for i = 1:max(id)
                
                % coordiantes and image number for rolony of the BC
                row = id == i;
                ixyz = xyz(row,:);
                iIm = im(row,:);
                
                if size(ixyz,1) < 2
                    continue
                end
                
                % Get all-to-all distance of rolony of the same barcode
                [D,I] = pdist2(ixyz(:,1:2),ixyz(:,1:2),'euclidean','Smallest',2);
                D = D(2,:); I = I(2,:); % Exclude itself
                D = D'; I = I';
                
                % Find rolony within the range, from the same seciton but different
                % image (get the one from later image to delete)
                % Same section: z; distance within threhold
                iTF = D < threshold & ixyz(:,3) == ixyz(I,3) & iIm > iIm(I);
                
                TF(row) = iTF;
            end
            
            TF = mat2cell(TF,sz);
        end
        
        %% Function:    delSameRolonyDiffIm
        % Discription:  delete extra rolony imaged more than once
        function axonBCreg = delSameRolonyDiffIm(axonBCreg,threshold)
            % Input & output: axonBCreg, table, with register coordinates
            %           threshold, num, threshold to tell whether its the
            %           same rolony, in micron
            
            % Get the rolony to be deleted
            TF = TBS.isSameRolonyDiffIm(axonBCreg,threshold);
            
            for i = 1:size(axonBCreg,1)
                if ~any(TF{i})
                    continue
                end
                axonBCreg{i,:} = cellfun(@(X) X(~TF{i},:),axonBCreg{i,:},...
                    'Uniformoutput',false);
            end
        end        
        
        %% Function:    plotSomaBCcount
        % Discription:  Plot soma BC counts for different radius
        function plotSomaBCcount(D,threshold)
            % Input:    D, cell, one cell per BC, distance of soma pixel to center
            %           threshold, num, min soma BC count
            % Output:   figure
            
            % Soma radius range
            radius = 10:10:100;
            
            % Edges for histogram
            edges = 0:5:1000;
            
            % Hisotgram for soma range and pixel count within range -------
            h = [];
            for i = 1:numel(radius)
                % Pixel count within radius
                imat = cellfun(@(X) sum(X <= radius(i)),D);
                h(:,i) = histcounts(imat,edges,'Normalization','cumcount');
            end
            
            % Plot
            edges = (edges(2:end)+edges(1:end-1))./2;
            c = jet(numel(radius));
            figure; hold on;
            for i = 1:size(h,2)
                plot(edges,h(:,i),'Color',c(i,:));
            end
            
            disp(['Min soma pixel count within radius: ',num2str(threshold)]);
            xline(threshold,'k--','LineWidth',2);
            
            % Colorbar
            colormap(jet); caxis([min(radius) max(radius)]);
            c = colorbar; c.Label.String = 'Radius (\mum)';
            
            xlabel('Soma pixel within radius'); ylabel('BC counts');
            g = gca; g.XLim = [0 200];
            pbaspect([1 1 1]);
            set(gcf,'Position',[100 100 300 300])
        end
        
        %% Function:    getSomaBCim
        % Discription:  get soma image finding soma location, with birghtes
        % pixel intensity
        function somaBCim = getSomaBCim(somaImName,nCh,directory)
            % Input:    somaImName, cell, name of soma image
            %           imageSetting, struct,
            % Output:   somaBCim, table, image of median of max intensity
            % across cycles
            
            disp('Function getSomaBCim on progress...');
            
            somaBCim = {};
            parfor i = 1:numel(somaImName) % parfor
                
                % Get image
                imName = somaImName{i};
                
                cd(directory.main);
                idir = dir(['**/soma/',imName]);
                im = TBS.getStack(fullfile(idir.folder,imName),[]);
                
                % max projection for each sequencng cycle
                sz = size(im);
                im = reshape(im,sz(1),sz(2),nCh,[]);
                im = max(im,[],3);
                im = squeeze(im);
                
                % Median pixel intensity across sequencing cycle
                somaBCim{i,1} = median(im,3);
            end
            
            somaBCim = cell2table(somaBCim,'VariableNames', {'im'},...
                'RowNames',somaImName);
        end
        
        %% Function:    getSomaSection
        % Discription:  get soma section number (sum of intensity)
        function somaSectionNum = getSomaSection(somaBC,somaIm,imageSetting,sysSetting)
            % Note: 08112021, use soma image
            % Input:    somaBC, table, soma pixel coordinates and BCid
            %           somaIm, table, image of median of max intensity
            % across cycles
            %           imageSetting/sysSetting, struct
            % Output:   somaLoc, mat, soma coordinates in micron, one BC per row
            
            disp('Function getSomaSection on progress...');
            
            % Organize rows as somaBCreg
            [TF,Lcob] = ismember(somaBC.Properties.RowNames,somaIm.Properties.RowNames);
            somaIm = somaIm(Lcob(TF),:);
            
            % Section number from the image table
            sectionName = somaIm.Properties.RowNames;
            sectionNumber = TBS.getSectionNumber(sectionName,imageSetting,sysSetting);
            
            % Image size
            sz = cellfun(@size,somaIm.im,'Uniformoutput',false);
            
            % xy-coordiantes and value
            xy = cellfun(@(X,Y) [X,Y],somaBC.x,somaBC.y,'UniformOutput',false);
            ind = cellfun(@(X,Y) sub2ind(X,Y(:,2),Y(:,1)),sz,xy,'Uniformoutput',false);
            v = cellfun(@(X,Y) X(Y),somaIm.im,ind,'UniformOutput',false);
            
            codeID = vertcat(somaBC.codeID{:});
            
            somaSectionNum = [];
            parfor i = 1:max(codeID) % parfor
                % Get pixels belongs to the BC
                row = cellfun(@(X) X == i,somaBC.codeID,'UniformOutput',false);
                
                % Pixel value in each image
                iv = cellfun(@(X,Y) X(Y),v,row,'Uniformoutput',false);
                
                % Sum of intensity in each image
                iv = cellfun(@sum,iv);
                
                % Find the section with max value
                [~,I] = max(iv,[],'omitnan');
                
                somaSectionNum(i,1) = sectionNumber(I);
            end
        end
        
        %% Function:    getSomaLoc
        % Discription:  soma location in xyz coordinates
        function somaLoc = getSomaLoc(somaBCreg,somaIm,imageSetting,sysSetting)
            % Note: 08112021, use soma image
            % Input:    somaBC, table, soma pixel coordinates and BCid
            %           somaIm, table, soma image
            %           imageSetting/sysSetting, struct
            % Output:   somaLoc, mat, soma coordinates in micron, one BC per row
            
            disp('Function getSomaLoc on progress...');
            
            % Section number from the image table
            sectionName = somaIm.Properties.RowNames;
            sectionNumber = TBS.getSectionNumber(sectionName,imageSetting,sysSetting);
            
            % Soma section number for each barcode
            sectionNumberBC = TBS.getSomaSection(somaBCreg,somaIm,imageSetting,sysSetting);
            
            codeID = vertcat(somaBCreg.codeID{:});
            
            somaLoc = [];
            for i = 1:max(codeID) % NO parfor, image table is too big
                
                % Find the image name with soma location
                row = sectionNumber == sectionNumberBC(i);
                
                imName = somaIm.Properties.RowNames{row};
                
                % Soma image
                im = somaIm.im{imName};
                
                % Find the coordinates of the soma pixel
                xy = [somaBCreg.x{imName},somaBCreg.y{imName}];
                xyz = somaBCreg.xyz{imName};
                
                TF = somaBCreg.codeID{imName} == i;
                xy = xy(TF,:);
                xyz = xyz(TF,:);
                
                if ~any(TF)
                    somaLoc(i,:) = nan(1,3); continue;
                end
                
                % Get the center of highest 90% of the soma pixel
                TF = TBS.xyzv2im(size(im),xy,[]);
                v = im(TF);
                I = v >= prctile(v,90);
                jCenter = median(xy(I,:),1);
                
                % Find the pixel closet to xy
                [~,I] = pdist2(xy,jCenter,'euclidean','Smallest',1);
                
                somaLoc(i,:) = xyz(I,:);
            end            
        end
                
        %% Function:    hasSoma
        % Discription:  determine whether the BC has soma pixels
        function TF = hasSoma(minCount,somaR,somaBC,somaImReg,imageSetting,sysSetting)
            % Input:    minCount, num, min BC count within the radius
            %           somaR, num, soma radius (from the center)
            %           somaBC, table
            %           somaIm, table, image of median of max intensity
            % across cycles
            %           imageSetting/sysSetting, sturct
            % Output:   TF, logical, whether the BC has enough pixel count
            % within the selective radius
            
            disp('Function hasSoma on progress...');
            
            % Get coordinates for soma pixel, in micron
            xyz= vertcat(somaBC.xyz{:});
            codeID = vertcat(somaBC.codeID{:});
            
            % Get soma location
            somaLocation = TBS.getSomaLoc(somaBC,somaImReg,imageSetting,sysSetting);
            
            % Get pixel to center distance for each BC --------------------
            D = {};
            parfor i = 1:max(codeID) % parfor
                
                % Get soma pixels of the BC
                row = codeID == i;
                ixyz = xyz(row,:);
                
                iSomaLocation = somaLocation(i,:);
                % Avoid pdist2 warning
                iSomaLocation = double(iSomaLocation);
                ixyz = double(ixyz);
                
                % Distance of each pixel to the center
                iD = pdist2(iSomaLocation,ixyz);
                D{i,1} = iD';
            end
            
            % Plot relationship between soma radius and counts
            TBS.plotSomaBCcount(D,minCount);
            
            % Whether BC passed the requirement
            % Pixel count within radius, beyond threhold
            n = cellfun(@(X) sum(X <= somaR),D);
            TF = n >= minCount;
            
            disp(['Has soma: ', num2str(sum(TF)),', in total ',...
                num2str(numel(TF))]);            
        end
        
        %% Function:    sBCcount
        % Discription:  barcode count for section or slide
        function [BCcount,sectionRegion] = sBCcount(tbl,str,imageSetting,sysSetting)
            % Input:    tbl, table, one image per row, with codeID
            %           str, 'section'/'slide', barcode count for section
            %           or slide
            %           imageSetting/sysSetting, struct
            % Output:   BCcount, mat, BC count per BC per section per
            % region. 2nd-D is the section number
            %           sectionRegion, cell, region name along the 3rd-D
            
            sectionLoc = imageSetting.sectionLoc;
            mouseID = imageSetting.mouseID;
            
            imName = tbl.Properties.RowNames;
            % Section/slide number
            switch str
                case 'section'
                    sNum = TBS.getSectionNumber(imName,imageSetting,sysSetting);
                    
                case 'slide'
                    sNum = cellfun(@(X) TBS.getSlideNumLoc(X,sysSetting,sectionLoc),imName);
                    
                    % Fix slide number according for this mouse
                    sNum = TBS.mouseSlideNum(sNum,mouseID);
            end
            
            % BC count for each image -------------------------------------
            % number of BC
            nBC = vertcat(tbl.codeID{:});
            nBC = max(nBC);
            
            BCid = tbl.codeID;
            BCid = cellfun(@(X) accumarray(X,1,[nBC 1]),BCid,'UniformOutput',false);
            
            % BC count per BC per section per region ----------------------
            % Get region
            sectionRegion = cellfun(@(X) TBS.nameFun(X,2,sysSetting),...
                imName,'UniformOutput',false);
            sectionRegion = erase(sectionRegion,sysSetting.imFormat);
            [sectionRegion,~,ic] = unique(sectionRegion);
            
            % One BC per row; one colume per section; 3rd D per region
            BCcount = zeros(nBC,max(sNum),numel(sectionRegion));
            for i = 1:numel(sectionRegion)
                % Row number of this region
                row = ic == i;
                row = find(row);
                
                % Section and slide number
                iSNum = sNum(row);
                
                [iSNum,~,ic2] = unique(iSNum);
                
                % BCcount of the region in this section/slide
                iBCid = {};
                for j = 1:numel(iSNum)
                    
                    row2 = ic2 == j;
                    row2 = row(row2);
                    iBCid{1,j} = BCid(row2);
                end
                
                iBCid = cellfun(@(X) horzcat(X{:}),iBCid,'UniformOutput',false);
                iBCid = cellfun(@(X) sum(X,2),iBCid,'UniformOutput',false);
                
                % col: slide/section number
                iBCid = horzcat(iBCid{:});
                
                BCcount(:,iSNum,i) = iBCid;
            end
        end
        
        %% Function:    floatingStat
        % Discription:  flating rolony stats, compare counts between soma
        % section and its neigbors
        function stat = floatingStat(axonBC,somaSection,nearbyN,imageSetting,sysSetting)
            % Input:    axonBC, table,
            %           nearbyN, num, nearby N-slide from the soma
            %           somaSlide, vector, section number with soma 
            %           imageSetting/sysSetting, struct
            % Output:   stat, mat, difference of BC count between soma and
            % its neigboring sections
            
            % BC count per slide ('slide'/'section')
            % One BC per row; col, slide/section; 3rd, region
            [BCcount,sectionRegion] = TBS.sBCcount(axonBC,'section',imageSetting,sysSetting);
            
            % Section number difference to soma
            % col: slide number; one row per BC
            x = repmat(1:size(BCcount,2),size(BCcount,1),1);
            x = x - somaSection;
            
            % Note, 0.5 will be rounded to interger
            x = round(x);
            
            stat = [];
            for i = 1:size(BCcount,3)
                
                for j = 1:size(BCcount,1)
                    % Difference to the soma section
                    jx = x(j,:);
                    jBCcount = BCcount(j,:,i);
                    
                    % Get axonBC count on soma and nearby slides
                    jSoma = jx == 0;
                    jNearby = abs(jx)<= nearbyN;
                    jNearby = jNearby & ~jSoma;
                    
                    jSoma = jBCcount(jSoma);
                    % Cannot use median, it is easy to be 0
                    jNearby = mean(jBCcount(jNearby));
                    
                    % Cannot use /, the jNearby can be 0
                    stat(j,i) = jSoma-jNearby;
                end
            end            
        end
        
        %% Function:    delCloseDot
        % Discription:  delete dots close to gether
        function im = delCloseDot(im,rng)
            % Note, the script likely has bias
            % Input & output:   im
            %           rng, num, range in pixel
            
            % Coordinates
            [y, x] = find(im);
            xy = [x y];
            
            while ~isempty(xy)
                % Calculate all to all distance
                D = pdist(xy);
                D = squareform(D);
                
                % Dot within the range
                D = D < rng & D > 0;
                
                % Delete thr first one
                [I,~] =  find(D,1,'first');     
                
                if isempty(I)
                    break
                end
                
                xy(I,:) = [];                
            end
            
            im2 = TBS.xyzv2im(size(im),xy,[]);
            % Can keep the value
            im(~im2) = 0;
        end
                
        %% Function:    findFloatSlide
        % Discription:  find slides with floating rolony basing on rolony
        % number and distribution
        function floatSlide = findFloatSlide(tbl,minDotCount,regionName,...
                bcSetting,imageSetting,sysSetting)
            % Input:    tbl, table, with codeID and xyz
            %           minDotCount, num, minimum axonBC count on the
            %           seciton for computation
            %           regionName, cell, region included into the
            %           computation
            %           bcSetting/imageSetting/sysSetting, struct,
            % Output:   floatSlide, vector, slide number with floating
            % rolonies
            
            disp('Function findFloatSlide on progress...')
            
            % Setting                                              
            % Slide thickeness: n secitons (for binning)
            sectionThickness = imageSetting.slideThickness;
            
            % Minimum pixel range for dots within a slide
            % To calculate the coverage of rolony
            minRng = round(100./sectionThickness);
            
            % Filter diameter (pixel)
            % 09122021, tested
            dia = 14;
            
            % Filter 1: sphere filter
            % Exclude dots with any dots within the filter
            SE1 = strel('sphere',dia);
            SE1 = SE1.Neighborhood;
            % Leave the middle 3 slide empty
            % calculation due to strel fnction may have different size
            dia2 = ceil(size(SE1,1)/2);
            SE1(:,:,dia2+[-1:1]) = 0;
            
            % Trim out the extra to speed up
            TF = TBS.getImValidLim(SE1,1);
            SE1 = SE1(TF(1):TF(2),:,:);
            TF = TBS.getImValidLim(SE1,2);
            SE1 = SE1(:,TF(1):TF(2),:);
            
            % Filter 2: disk filter
            % Filter for rolony cover area
            SE2 = strel('disk',dia);
            
            
            % Slide range for soma
            somaSectionRnge = bcSetting.somaSectionRnge;            
            % Function handle of slide within the soma range
            withinRng = @(X) X >= min(somaSectionRnge) & X <= max(somaSectionRnge);
            
            % Only include image of selective region
            regionName = cellfun(@(X) ['_',X,'.tif'],regionName,'UniformOutput',false);
            
            imName = tbl.Properties.RowNames;
            row = contains(imName,regionName);
            tbl = tbl(row,:);
            imName = imName(row);
                       
            % Section number for each rolony
            sectionNum = TBS.getSectionNumber(imName,imageSetting,sysSetting);
            sectionNum = TBS.repmat2cell(sectionNum,tbl.codeID);
                        
            id = vertcat(tbl.codeID{:});
            xyz = vertcat(tbl.xyz{:});
            sectionNum = vertcat(sectionNum{:});
            
            nBC = max(id);
            
            % NaN: not able to find floating slide
            % 0: no rolony identified in the region/section
            % number: floating slide identified
            floatSlide = nan(nBC,1);
            parfor i = 1:nBC % parfor
                
                row = id == i;
                
                % Registered coordinates
                ixyz = xyz(row,:);
                
                % Slide number
                iSectionNum = sectionNum(row,:);
                
                % No rolony in the regions
                % Or no rolony in the soma range
                if ~any(row) || ~any(withinRng(iSectionNum))
                    floatSlide(i) = 0;
                    continue
                end       
                
                % Get 3D mat of rolony location ---------------------------
                ixyz = round(ixyz./sectionThickness);
                ixyz(:,3) = iSectionNum;
                
                sz = max(ixyz,[],1);
                sz(1:2) = fliplr(sz(1:2));                
                % Edge effect for conv (leave a little bit more, + 3)
                sz(1:2) = sz(1:2) + dia + 3;
                
                im = TBS.xyzv2im(sz,ixyz,[]);
                
                % 1. Exclude dots with rolony on the top/bottom -----------  
                TF = imdilate(im,SE1);
                im1 = im;
                im1(TF) = false;
                
                % Delete slide with dots below minimum
                slide = sum(im1,1:2) >= minDotCount;   
                slide = find(slide);
                
                % Exclude dots outside of the range
                % Note, cannot delete dot before exclusion filter
                row = withinRng(slide);
                slide = slide(row);
                
                im1 = im1(:,:,slide);
                
                if isempty(im1)
                    continue
                end
                
                % 2. Check overlap ----------------------------------------
                % Use disk filter to measure sparese/covered area
                % (to speed up) select slide and do 2D filter instead of 3D
                                
                im2 = zeros(size(im1));
                for j = 1:size(im1,3)
                    
                    % Combine neigboring slides (-1 & +1)
                    % i.e. the soma can be cut in the middle of both slides
                    row = slide(j);
                    row = row + (-1:1);
                    row = ismember(slide,row);
                    
                    jim = im1(:,:,row);
                    jim = max(jim,[],3);
                                        
                    % Delete dots close together, get cluster/dot number
                    % Need to have enough cluster to be countted for
                    % sparseness
                    jim = TBS.delCloseDot(jim,minRng);
                    
                    if sum(jim,'all') < minDotCount
                        continue
                    end
                    
                    im2(:,:,j) = imdilate(jim,SE2);
                end
                
                % No slide has enough rolony to be sparse
                if ~any(im2,'all')
                    continue
                end
                
                % Area counts size/sparse
                im2 = im2 ~= 0;
                D = squeeze(sum(im2,1:2));   
                
                [D,I] = sort(D,'descend');
                slide = slide(I);
                
                % Combine neigborhing slides (get median)
                if numel(D) > 1 &&...
                    abs(diff(slide(1:2))) == 1 && D(1)== D(2)   
                
                    slide(1) = mean(slide(1:2));
                    slide(2) = [];  D(2) = [];
                end
               
                % Get the slide with most sparse area
               floatSlide(i) = slide(1);
            end
        end
        
        %% Function:    delFloatingRolony
        % Discription:  delete floating rolony 
        function tbl = delFloatingRolony(somaSection,tbl,regionName,imageSetting,sysSetting)
            % Input & output: somaSlide, slide number the soma is localized, one row/BC
            %           tbl, table, with codeID
            %           regionName, cell, region to delete floating
            %           rolonies
            %           imageSetting/sysSetting, struct
            
            % Get section number to be deleted ----------------------------
            % For interger, delete 3 slides, -1:1
            % For non-interger, delete two, ceil & floor
            somaSection2 = {};
            for i = 1:size(somaSection,1)
                iNum = somaSection(i);
                
                if mod(iNum,1)==0
                    somaSection2{i,1} = [-1:1]' + iNum;
                else
                    somaSection2{i,1} = [floor(iNum);ceil(iNum)];
                end       
            end
            
            % barcode id
            id = TBS.repmat2cell([1:size(somaSection,1)]',somaSection2);
            id = vertcat(id{:});
            
            somaSection2 = vertcat(somaSection2{:});  
            
            % Delete axon rolony from all images of the slide -------------
            
            imName = tbl.Properties.RowNames;
            sectionNum = TBS.getSectionNumber(imName,imageSetting,sysSetting);
            
            % Only include image of selective region
            regionName = cellfun(@(X) ['_',X,'.tif'],regionName,'UniformOutput',false);            
            
            for i = 1:size(sectionNum,1)
                
                if ~contains(imName{i},regionName)
                    continue
                end
                                
                row = sectionNum(i);
                
                % Find the codeID with soma on this slide
                row = somaSection2 == row;
                row = id(row);
                
                if isempty(row)
                    continue
                end
                
                iId = tbl.codeID{i};
                TF = ~ismember(iId,row);
                tbl{i,:} = cellfun(@(X) X(TF,:),tbl{i,:},'UniformOutput',false);
            end            
        end
        
        %% Function:    getAxonBCcenterInj
        % Discription:  get the center of axonBC in injection site
        function axonBCcenter = getAxonBCcenterInj(axonBC,sysSetting)
            % Input:    axonBC, table, registrated axonBC coordinates
            %           sysSetting, struct
            % Output:   axonBCcenter, mat, xyz one BC per row
            
            imName = axonBC.Properties.RowNames;
            
            % Image for soma bscalling (i.e. injection)
            row = cellfun(@(X) TBS.doSomaBscall(X,sysSetting),imName);
            axonBC = axonBC(row,:);
            
            id = vertcat(axonBC.codeID{:});
            xyz = vertcat(axonBC.xyz{:});
            
            % Use median as the center
            axonBCcenter = [];
            for i = 1:max(id)
                row = id == i;
                
                if ~any(row)
                     axonBCcenter(i,:) = nan(1,3); continue
                end
                
                ixyz = xyz(row,:);
                
                axonBCcenter(i,:) = median(ixyz,1);
                
                % (Check point) ===========================================  
%                 
%                 % Plot dot location
%                 % All dots
%                 hold off; scatter3(xyz(:,1),xyz(:,2),xyz(:,3),0.5,'k','filled','MarkerEdgeAlpha',0.1);
%                 % Dots of the current barcoded cell
%                 hold on; scatter3(ixyz(:,1),ixyz(:,2),ixyz(:,3),7,'r','filled','MarkerEdgeAlpha',0.8);
%                 % Center
%                 scatter3(iCenter(:,1),iCenter(:,2),iCenter(:,3),50,'k');
%                 daspect([1 1 1]); pbaspect([1 1 1]); grid off; view([0 0]);
%                 g = gca; g.YDir = 'reverse'; g.ZDir = 'reverse';
%                 g.XLim = [500 10500]; g.YLim = [0 5000]; title(i);
            end
        end
        
        %% Function:    isGlia
        % Discription:  whether the barcoded cell is glia, only has rolony
        % within the defined radius
        function TF = isGlia(somaLocation,axonBCreg,gliaR,minAxonCount)
            % Input:    somaLocation, mat, in registrated coordinates
            %           axonBCreg, table, with registrated axonBC coordinates
            %           gliaR, num, defined glia radius
            %           minAxonCount, num, minimum count of axon rolony
            % Output:   TF, logical, whether the BC is regonized as a glia
            
            % Get axonal rolony soma distance 
            codeID = vertcat(axonBCreg.codeID{:});
            xyz = vertcat(axonBCreg.xyz{:});
            
            D = {};
            parfor i = 1:max(codeID) % parfor
                iSomaLoccation = somaLocation(i,:);
                
                row = codeID == i;
                
                ixyz = xyz(row,:);
                
                ixyz = double(ixyz);
                D{i,1} = pdist2(iSomaLoccation,ixyz);
            end
            
            % Soma count of axonal rolony beyond glia radius
            TF = cellfun(@(X) X >= gliaR, D,'Uniformoutput',false);
            % Count filter for axonal rolony
            TF = cellfun(@sum,TF);
            TF = TF < minAxonCount;
            
            % Rows without somaLocation
            row = isnan(somaLocation);
            row = any(row,2);
            TF(row) = false;
            
            disp(['Found glia: ',num2str(sum(TF)), ' in ',num2str(size(TF,1))]);
        end
                
        %% Function:    regionCountFilter
        % Discription:  exclude projected region below threhold
        function axonBC = regionCountFilter(axonBC,minCount,selectRegion,...
                imageSetting,sysSetting)
            
            % BC count per section
            [BCcount,sectionRegion] = TBS.sBCcount(axonBC,'section',imageSetting,sysSetting);
            
            % Only include the selective region
            if ~isempty(selectRegion)
                selectRegion = ismember(sectionRegion,selectRegion);
                
                sectionRegion = sectionRegion(selectRegion);
                BCcount = BCcount(:,:,selectRegion);
            end
            
            sectionRegion = cellfun(@(X) ['_',X,'.'],sectionRegion,'UniformOutput',false);
            
            % BC count of the region
            n = sum(BCcount,2);
            n = squeeze(n);
            n = n >= minCount;
            
            for i = 1:numel(sectionRegion)
                imName = axonBC.Properties.RowNames;
                row = contains(imName,sectionRegion{i});
                row = find(row);
                
                % codeID beyond threshold
                id = n(:,i);
                id = find(id);
                
                % Delete id less than threshold
                TF = axonBC.codeID(row);
                TF = cellfun(@(X) ismember(X,id),TF,'Uniformoutput',false);
                
                for j = 1:size(TF,1)
                    irow = row(j);
                    axonBC{irow,:} = cellfun(@(X) X(TF{j},:),...
                        axonBC{irow,:} ,'UniformOutput',false);
                end
            end
        end
        
        %% Function:    exclBCinROI
        % Discription:  excle barcodes within ROI
        function BCtable = exclBCinROI(roi,BCtable,varName)
            % Input & output:  roi, logical stack, region of interest
            %           BCtable, table, with barcode coordinates
            %           varName, str, variableName in BCtable to determine 
            %           whether BC is within ROI
            
            xyz = BCtable.(varName);
            
            sz = cellfun(@(X) size(X,1),xyz);
            
            xyz = vertcat(xyz{:});
            
            xyz = round(xyz);
            
            % xyz within range
            TF = TBS.isOutOfRngeDot(size(roi),xyz);
            TF = ~TF;
            
            % xyz outside ROI, xyz to keep
            ind = sub2ind(size(roi),xyz(TF,2),xyz(TF,1),xyz(TF,3));
            TF(TF) = ~roi(ind);
            
            disp(['After exclude out of ROI: ', num2str(sum(TF)),...
                ' of total ',num2str(numel(TF))]);
            
            % Exclude the rows in the table -------------------------------
            TF = mat2cell(TF,sz);
            for i = 1:size(BCtable,1)
                BCtable{i,:} = cellfun(@(X) X(TF{i},:),BCtable{i,:},...
                    'Uniformoutput',false);
            end            
        end
        
        %% Function:    BCinROI
        % Discription:  find barcode with rolony inside ROI
        function roiBC = BCinROI(roi,BCtable,varName)
            % Input:    roi, logical stack, region of interest
            %           BCtable, table, with barcode coordinates
            %           varName, str, variableName in BCtable to determine
            %           whether BC is within ROI
            % Output:   roiBC, table, with barcode id for each image
            
            xyz = BCtable.(varName);
            
            sz = cellfun(@(X) size(X,1),xyz);
            
            xyz = vertcat(xyz{:});
            
            xyz = round(xyz);
            
            % xyz within range
            TF = TBS.isOutOfRngeDot(size(roi),xyz);
            TF = ~TF;
            
            % xyz inside ROI
            ind = sub2ind(size(roi),xyz(TF,2),xyz(TF,1),xyz(TF,3));
            TF(TF) = roi(ind);
            
            TF = mat2cell(TF,sz);
            
            id = cellfun(@(X,Y) X(Y),BCtable.codeID,TF,...
                'UniformOutput',false);
            
            roiBC = cell2table(id,'VariableName',{'codeID'},...
                'RowNames',BCtable.Properties.RowNames);
        end
        
         %% Function:   exclFloatingUseROI
        % Discription:  exclude floating rolony using BC identified by ROI,
        % delete all rolony on the seciton
        function BCtable = exclFloatingUseROI(BCtable,roiBC)
            % Input & output: BCtable, table, with barcode coordinates
            %           roiBC, table, with barcode id for each image
            
            TF2 = {};
            for i = 1:size(roiBC,1)
                imName = roiBC.Properties.RowNames{i};
                
                % rows doesnt belong to BC in ROI
                TF = ~ismember(BCtable.codeID{imName},roiBC.codeID{imName});
                
                BCtable{imName,:} = cellfun(@(X) X(TF,:),...
                    BCtable{imName,:},'UniformOutput',false);
                
                TF2{i,1} = TF;
            end
            
            TF = vertcat(TF2{:});
            disp(['After exclude flaoting rolony using ROI: ',...
                num2str(sum(TF)),' of total ',num2str(numel(TF))]);
        end
        
        %% Function:    exclSporadicBC
        % Discription:  exclude sporadic BC in different regions
        function axonBC = exclSporadicBC(axonBC,regVoxel,regionMinCount)
            % Input & output:   axonBC, table
            %               regVoxel, cell, one cell for one logical stack for a brain region
            %               regionMinCount, vector, threshold for sporadic barcode for
            %               each region
            
            xyz = axonBC.xyzRef;
            id = axonBC.codeID;
            sz = cellfun(@(X) size(X,1),id);
            
            xyz = vertcat(xyz{:});
            id = vertcat(id{:});
            
            % index for the coordinate
            xyz = round(xyz);
            ind = sub2ind(size(regVoxel{1}),xyz(:,2),xyz(:,1),xyz(:,3));
            
            [~,~,ic] = unique(id);
            
            TF = [];
            for i = 1:numel(regVoxel)
                iReg = regVoxel{i};
                
                % whether coordinates are in the region
                inReg = iReg(ind);
                
                % Count rolony within the region for each BC
                n = accumarray(ic,inReg);
                
                % BC with less than min count in the region
                n = n < regionMinCount(i);
                I = find(n);
                
                % coordinates from BC to be exclude for the current region
                iTF = inReg & ismember(id,I);
                
                TF(:,i) = iTF;
            end
            
            TF = any(TF,2);
            
            disp(['After exclude sporadic BC: ', num2str(sum(~TF)),...
                ' of total ',num2str(numel(TF))]);
            
            TF = mat2cell(TF,sz);
            
            % Delete rolony doesnt match the criteria ---------------------
            for i = 1:size(axonBC,1)
                axonBC{i,:} = cellfun(@(X) X(~TF{i},:),axonBC{i,:},...
                    'UniformOutput',false);
            end
            
        end
        
        %% Function:    imBCid (for validation)
        % Discription:  get image of BC location and ID
        function im = imBCid(tbl,iSeq)
            % Input:    tbl, table, bscall result & coordinates
            %           iSeq, number, sequence or coordinates to check
            % Output:   im, mat, image of BC id in the location
            
            if nargin == 1 || isempty(iSeq)
                iSeq = 1;
            end
            
            x = tbl.x; y = tbl.y;
            if iscell(x) && iscell(y)
                x = x{:}; y = y{:};
            end
            
            x = x(:,iSeq); y = y(:,iSeq);
            
            sz = [max(y),max(x)];
            sz = ceil(sz);
            
            if contains('codeID',tbl.Properties.VariableNames)
                v = tbl.codeID{:};
            else
                v = [];
            end
            
            im = TBS.xyzv2im(sz,[x y],v);
        end
        
    end
    
    methods (Static)    % Reference map ===================================
        %% Function:    getRefMap
        % Discription:  get reference map (nissl, avgTemplate, annotation)
        function refMap = getRefMap(str,refSetting)
            % Input:    str, string, name of the reference map
            %           refSetting, struct, reference settings
            % Output:   refMap, mat, image stack
            
            directory = refSetting.directory;
            
            switch str
                case 'nissl'
                    refName = refSetting.nisslName;
                case 'avg'
                    refName = refSetting.avgName;
                case 'anno'
                    refName = refSetting.annoName;
            end
            
            % Get imgage (need nrrdread under the driectory)
            cd(directory)
            refMap = nrrdread(refName);
            
            if strcmp('nissl',str)
                % (don't know why cannot directly convert to uint8)
                % Annomap cannot be change due to the ID
                refMap = uint16(refMap);
                refMap = im2uint8(refMap);
            end
            
            % Sagital to coronal plate
            refMap = permute(refMap,[1 3 2]);            
        end
                
        %% Function:    findAnnoLevel
        % Discription:  find annotation structure with the selective level
        function list = findAnnoLevel(annoStruct,level)
            % Input:    annoStruct, sturct, of annotation
            %           level, num, level of the region
            %           currently (only support one)
            % Output:   list, struct, of areas with the level
            
            list = [];
            while 1 > 0
                annoStruct = vertcat(annoStruct.children);
                if isempty(annoStruct)
                    break
                end
                
                iLevel = [annoStruct.st_level];
                
                TF = iLevel == level;
                list{end+1,1} = annoStruct(TF,:);
            end
            
            list = vertcat(list{:});
        end
        
        %% Function:    getAnnoRegion
        % Discription:  getAnnotation regions
        function refMap = getAnnoRegion(refSetting)
            % Input:    refSetting, sturcture, reference settings
            % Output:   refMap, mat, reference map stack
            
            annoStruct = refSetting.annoStruct;
            
            % Get level 11 id & their parents id
            list = TBS.findAnnoLevel(annoStruct,11);
            list = [[list.parent_structure_id];[list.id]];
            list = list';
            
            refMap = TBS.getRefMap('anno',refSetting);
            
            % Combined the regions
            [Lia,Locb] = ismember(refMap,list(:,2));
            refMap(Lia) = list(Locb(Lia),1);
        end
        
        %% Function:    annoRegionOutline
        % Discription:  get outline of annotated regions
        function regionOutline = annoRegionOutline(annoRegion)
            % Input:    annoRegion, mat, 2d/3d works, value indicating
            % regions
            % Output:   regionOutline, logical, outline of the regions
            
            % Draw outline using single pixel
            if ndims(annoRegion) ~= 3
                SE = strel('disk',1);
            else
                SE = strel('sphere',1);
            end
            
            regionOutline = imdilate(annoRegion,SE)~= annoRegion |...
                imerode(annoRegion,SE)~= annoRegion;
        end
        
        %% Function:    getBrainOutline
        % Discirption:  get brain outline uses annotation map
        % (for 3d-rotating brain)
        function outline = getBrainOutline(annoMap,filterSz)
            % Input:    annoMap, mat, annotation map
            %           filterSz, num, filter size for avaerge filter
            % Output:   outline, mat, brain outline
            
            % Get outline of the brain
            outline = annoMap > 0;
            outline = imdilate(outline,ones(3)) & outline == 0;
            outline = uint8(outline).*255;
            
            % 3D-average filter
            h = fspecial3('average',filterSz);
            outline = imfilter(outline,h,'replicate');
        end

        %% Function:    findAnnoID
        % Discription:  find all ID for the brain region
        function id = findAnnoID(annoStruct,str)
            % Input:    annoStruct, sturct, of annotation
            %           str, char, region name
            %           currently (only support one)
            % Output:   id, vector, all id belongs to the region
                   
            % Find the struct (including children) of the region
            TF = 0;
            while ~any(TF)
                annoStruct = vertcat(annoStruct.children);
                iName = {annoStruct.name};
                TF = contains(iName,str);
            end
            
            annoStruct = annoStruct(TF);
            
            % Get all the id of it and its children
            id = annoStruct.id;
            TF = 1;
            while TF
                annoStruct = vertcat(annoStruct.children);
                if isempty(annoStruct)
                    break
                end
                id = [id; [annoStruct.id]'];
            end
        end
        
        %% Function:    findAnnoStrID
        % Discription:  find IDs of the region with the keywords
        function id = findAnnoStrID(annoStruct,str)
            % Input:    annoStruct, sturct, of annotation
            %           str, string, key work of annotation name
            % Output:   id, vector, id of the area with the key words
            
            list = [];
            while 1 > 0
                annoStruct = vertcat(annoStruct.children);
                if isempty(annoStruct)
                    break
                end
                
                list{end+1,1} = annoStruct;
            end
            
            list = vertcat(list{:});
            
            % Find region name contains the string
            TF = cellfun(@(X) contains(X,str,'IgnoreCase',true),{list.name});
            % Get the id of these regions
            id = [list(TF).id]';
        end
        
        %% Function:    isLeftIm
        % Discription:  is left of imge (for left/right hemispher)
        function TF = isLeftIm(im)
            % Input:    im, mat, image stack
            % Output:   TF, logical stack
            
            TF = false(size(im));
            TF(:,1:round(size(im,2)/2),:) = true;
        end
                       
    end
    
    methods (Static)    % Align to Allen map ==============================
        %% Function:    vol2micronTform
        % Discription:  convert coordinates from 3d volume scale to micron
        function tform = vol2micronTform(imageSetting,scaleFactor)
            % Input:    imageSetting, struct,
            %           scaleFactor, num, scale factor for 3d volume
            % Output:   tform, affine3d object
            
            resolution = imageSetting.resolution;
            slideThickness = imageSetting.slideThickness;
            
            resolution = scaleFactor/resolution;
            
            S = eye(4);
            S(1,1) = 1/resolution;
            S(2,2) = 1/resolution;
            S(3,3) = S(3,3).*slideThickness;
            
            % 03072022, not sure why there are two translation...?
            % Translation1: - 1
            T1 = eye(4);
            T1(4,1:3) = -1;
            
            % Translation2: + 1
            T2 = eye(4);
            T2(4,1:3) = 1;
                       
            tform = T1*S*T2;
            
            tform = affine3d(tform);
        end
        
        %% Function:    sortHipDot
        % Discription:  sort dots along a sequence, front-anterior end
        function sortedDot = sortHipDot(dotXYZ)
            % Input:    dotXYZ, mat, N*3
            % Output:   sortedDot,cell, with mat of sorted dots
            %               2x2, upper & bottom half of L/R side
            
            maxZ = max(dotXYZ(:,3));
            
            % Get the upper/lower half
            centerYRow = dotXYZ(:,3) == maxZ;
            centerY = mean(dotXYZ(centerYRow,2));
            
            sortedDot = {};
            for jHalf = 1:2
                if jHalf == 1
                    row = dotXYZ(:,2) <= centerY;
                    iDir = 'ascend';
                else
                    row = dotXYZ(:,2) >= centerY;
                    iDir = 'descend';
                end
                
                jDot = dotXYZ(row,:);
                jDot = sortrows(jDot,3,iDir);
                sortedDot{jHalf,1} = jDot;
            end
            
            % Delete the repeated copy of middle
            if sum(centerYRow) == 1
                iDist = [sortedDot{1}(end-1,:);sortedDot{2}(2,:)];
                iDist = pdist2(sortedDot{1}(end,:),iDist);
                if iDist(1) > iDist(2)
                    sortedDot{1}(end,:) = [];
                else
                    sortedDot{2}(1,:) = [];
                end
            end
            
        end
        
        %% Function:    getDistPrctile
        % Discription: get the prctile location of the dots in the
        % whole length
        function distPrctile = getDistPrctile(xyz)
            % Input:    xyz, mat, coordinates, sorted
            % Output:   distPrctile, vector
            
            iDist = diff(xyz,1,1);
            iDist = sqrt(sum(iDist.^2,2));
            iDist = cumsum(iDist);
            distPrctile = iDist./iDist(end);
            % Add the first dot
            distPrctile = [0; distPrctile];
        end
        
        %% Function:    interp1Dot
        % Discription:  get dots along the line using interp1
        function xyzQ = interp1Dot(xyz,prctileQ,maxPrctileDiff)
            % Input:    xyzIn, mat, coordinates of dots
            %           prctileQ, query prectile
            %           maxPrctileDiff, max prctile different from raw data
            % Output:   xyzQ, mat, query coordinates
            
            distPrctile = TBS.getDistPrctile(xyz);
            
            % Pick the 1st one if there is a duplication
            [distPrctile,ia,~] = unique(distPrctile);
            xyz = xyz(ia,:);
                        
            % Exclude reference precentage too far from the real data
            if ~isempty(maxPrctileDiff)
                prctileQ = reshape(prctileQ,[],1);
                [D, ~] = pdist2(distPrctile,prctileQ,'euclidean',...
                    'Smallest',1);
                row = D <= maxPrctileDiff;
                prctileQ = prctileQ(row,:);
            end
            
            % Inpterp though all axes
            for a = 1:size(xyz,2)
                xyzQ(:,a) = interp1(distPrctile,xyz(:,a),prctileQ,'spline');
            end
        end
        
         %% Function:   totalDist
        % Discription:  total distance of the coordinates
        function dist = totalDist(xyz)
            % Input:    xyz, mat
            % Ouput:    dist, num
            dist = sum(sqrt(sum(diff(xyz).^2,2)),'all');
        end
                              
        %% Function:    findYZangle
        % Discription:  find z & y rotation angle using hipDot
        function [yAngle,zAngle] = findYZangle(dotCell,yAngleInput,zAngleInput)
            % Input:    dotCell, cell, 1x2, dot coordinates per side
            %           yAngleInput/zAngleInput, row vector, angles for rotation
            % Output:   yAngle/zAngle, num, degree of rotation along y/z-axis
            
            % (Speed up by not using loop)
            % x-axis, yAngle; y-axis: zAngle
            % Get tform
            [gridY,gridZ] = meshgrid(yAngleInput,zAngleInput);
            tform = arrayfun(@(Y,Z) TBS.roty(Y)*TBS.rotz(Z),gridY,gridZ,'Uniformoutput',false);
            
            % Transform dots from each side
            dotLR = {};
            dotLR(:,:,1) = cellfun(@(X) dotCell{1}*X,tform,'Uniformoutput',false);
            dotLR(:,:,2)= cellfun(@(X) dotCell{2}*X,tform,'Uniformoutput',false);
            
            % Only include Y & Z, convert to single to speed up
            dotLR = cellfun(@(X) X(:,2:3),dotLR,'Uniformoutput',false);
            dotLR = cellfun(@single,dotLR,'Uniformoutput',false);
            
            % Calculate the scale ratio according to distance between both side -------
            % ie, longer side with <1 scale ratio
            scaleRatio = cellfun(@(X) TBS.totalDist(X),dotLR);
            scaleRatio = scaleRatio./max(scaleRatio,[],3);
            % Flip both side
            scaleRatio = cat(3,scaleRatio(:,:,2),scaleRatio(:,:,1));
            
            % Resample same number of dots in both side -------------------------------
            % Use the side with less does for reference percentile
            % scaleRatio == 1: full scale
            [~,I] = sort(scaleRatio,3,'descend');
            refPrctile = dotLR(I);
            refPrctile = refPrctile(:,:,1);
            refPrctile = cellfun(@(X) TBS.getDistPrctile(X),refPrctile,'Uniformoutput',false);
            
            % Scale similar equal distance using scale ratio
            refPrctile = cellfun(@(X,Y) X.*Y,repmat(refPrctile,1,1,2),...
                num2cell(scaleRatio),'Uniformoutput',false);
            
            dotLR = cellfun(@(X,Y) TBS.interp1Dot(X,Y,[]),dotLR,refPrctile,...
                'Uniformoutput',false);
            
            % Calculate cost
            fh = @(X1,X2) sqrt(sum((X1-X2).^2,2));
            cost = cellfun(@(X1,X2) fh(X1,X2),dotLR(:,:,1),dotLR(:,:,2),...
                'Uniformoutput',false);
            cost = cellfun(@std,cost);
            
            % Find localMin
            localMin = imerode(cost,ones(3)) == cost;
            localMin = localMin.*cost;
            [row,col,v] = find(localMin);
            
            if numel(v) > 1
                [~,I] = min(v);
                row = row(I); col = col(I); v = v(I);
                warning('Funciton findYZangle: more than one local minimum was founded.')
            end
            
            yAngle = yAngleInput(col); zAngle = zAngleInput(row);
            
            hold on; imagesc(yAngleInput,zAngleInput,cost);
        end
        
        %% Function:    findYZangleHighResolution
        % Discription:  find z & y rotation angle using hipDot
        function [yAngle,zAngle] = findYZangleHighResolution(dotCell,resolution)
            % This function is pretty fast, so just loop through high resolution
            % Input:    dotCell, cell, 1x2, dot coordinates per side
            %           resolution, num, resolution for degree
            % Output:   yAngle/zAngle, num, degree of rotation along y/z-axis
            
            disp('Start finding Z & Y angles...'); tic;
            
            % Angle interval
            aInterval = 3;
            
            % intial center degree
            yAngle = 0; zAngle = 0;
            
            figure;
            while aInterval >= resolution
                disp(['Finding YZ-angle: Resolution ',num2str(aInterval)]);
                yAngleInput = [-45:1:45].*aInterval+yAngle;
                zAngleInput = [-45:1:45].*aInterval+zAngle;
                
                [yAngle,zAngle] = TBS.findYZangle(dotCell,yAngleInput,zAngleInput);
                
                aInterval = aInterval/3;
            end
            
            title(['Cost of finding rotation angle',newline,'in y&z-Axis']);
            xline(yAngle,'k:'); yline(zAngle,'k:');
            xlabel('y-axis'); ylabel('z-axis');
            daspect([1 1 1]);
            set(gcf, 'Position',  [100, 100, 240, 240]);
        end
        
        %% Function:    getAllenD
        % Discription:  get displacement field to allen reference map
        function allenD = getAllenD(zq,im,refMap,redoTF,transformationSetting,directory)
            % Input:    zq, vector, slide number for registration
            %           im, mat, stack of current brain, same size as ref
            %           refMap, mat, stack of reference map
            %           redoTF, logical, whether to redo
            %           transformationSetting, struct, setting for
            %           alignment
            %           directory, struct
            % Output:   allenD, table, with point pairs, displacement field
            % for transformaiton; one slide per row
            
            % Load allenTform 
            cd(directory.main);
            if exist('allenD.mat')
                load('allenD.mat');
            else
                allenD = table('Size',[0 5],'VariableTypes',...
                    repmat({'cell'},1,5),'VariableNames',...
                    {'mp','D','fixName','fp','size'});
            end
                 
            sz = size(refMap,1:2);
            
            tic
            for i = zq
                movingName = i;
                
                % fix & moving
                fix = refMap(:,:,i);
                fix = imadjust(fix);
                moving = im(:,:,i);                
               
                % mp: moving point; fp: fixed point                
                % Whether use exist mp & fp (moving/fixed points)
                if ~redoTF && size(allenD,1)>i &&~isempty(allenD.mp{movingName})
                    mp = allenD{movingName,'mp'}{:};
                    fp = allenD{movingName,'fp'}{:};
                else
                    mp = []; fp = [];
                end
                
                % cpselectModify(moving,fixed,mp,fp,transformationType,maxIteration)
                [mp,fp,D] = TBS.cpselectIter(moving,fix,mp,fp,transformationSetting,20);
                
                % Update variables
                allenD(movingName,:)= table({mp},{D},{[]},{fp},{sz});
                
                % return var to workspace
                assignin('base','allenD',allenD);
                
                save(fullfile(directory.main,'allenD.mat'),'allenD');
                
                disp(['Done: ',num2str(movingName)]); close all; toc
            end
        end
        
        %% Function:    interpDcell
        % Discription:  interpolation displacement field in the cell, for
        % each axis, wih inter/extrapolation
        function D = interpDcell(Dcell,ax,sz)
            % Input:    Dcell, cell, with displacement field, one row per
            % slide
            %           ax, number, axis number for the displacement field
            %           sz, vector, size of the output displacement field
            % Output:   D, mat, displacement field of the axis
            
            D = zeros(sz);
            
            % Get displacement field from the table
            Z = [];
            for i = 1:size(Dcell,1)
                iD = Dcell{i};
                
                if isempty(iD)
                    continue
                end
                
                D(:,:,i) = iD(:,:,ax);
                
                % Slide with dispalcement field
                Z(end+1) = i;
            end
                        
            % Interpolation of dispalcement field between first and last --
            Z2 = Z(1):Z(end);
            
            % Interpolate 1-D for every pixel
            % (speed is ok)
            for i = 1:sz(1)
                for j = 1:sz(2)
                    
                    % Stack with displacement field
                    iD = D(i,j,Z);
                    iD = reshape(iD,size(Z));
                    D(i,j,Z2)= interp1(Z,iD,Z2,'spline');
                end
            end
            
            % Extrapolate the both end on z-axis --------------------------
            % Same as the first and last displacement field
            iD = D(:,:,Z(1));
            D(:,:,1:Z(1))= repmat(iD,1,1,Z(1));
            iD = D(:,:,Z(end));
            D(:,:,Z(end):end)= repmat(iD,1,1,sz(3)-Z(end)+1);            
        end
        
    end
    
    methods (Static)    % Construct flatmap, 03/2022 ======================
        %% Function:    interp1TwoDot
        % Discription:  interpolate point between two dots
        function xyzCell = interp1TwoDot(xyz1,xyz2,x,xq)
            % Input:    xyz1/2, mat, xyz coordinates of the begining/end
            % dot
            %           x, vector, sample points
            %           xq, vector, query points
            % Outut:    xyzCell, cell, coordinates for the qeury dots
            % between begining & end
            
            if numel(x)~=2
                error('Two sample points for interp1.');
            end
            
            xyzCell = {};
            for i = 1:size(xyz1,1)
                                                               
                % Loop throught each axis
                ixyz = [];
                for j = 1:size(xyz1,2)
                    jAx = interp1(x,[xyz1(i,j),xyz2(i,j)],xq,'spline');
                    ixyz(:,j) = jAx';
                end
                
                xyzCell{i,1} = ixyz;
            end            
        end
                
        %% Function:    fillRefCtx
        % Discription:  fill the empty space between reference values
        function refCtx = fillRefCtx(refCtx,ctx)
            % Note, use 7 pixel ball
            % quickly tried scatteredInterpolant to fill cortex, 
            % no significant different by eye
            % Input & output: reference cortex, mat
            %           ctx, logical stack, cortical area
            
            disp('Function fillRefCtx in progress...')
            
            % Step per filling
            SE = strel('sphere',3);
            SE = SE.Neighborhood;
            
            checkTF = inf;
            while checkTF > 0
                
                im = TBS.nonzeroAvgFilter(refCtx,SE);
                
                % Only include the empty area
                TF = ctx & ~refCtx;
                
                refCtx(TF) = im(TF);
                
                % Remaining voxel
                n = sum(~refCtx & ctx,'all');
                
                if checkTF == n
                    return
                end
                
                checkTF = n;
                
                disp(['Remaining voxel to fill: ',num2str(n)]);
            end
        end
        
        %% Function:    nonzeroAvgFiltCtx
        % Discription:  apply average filter on nonzero value to each side
        % of the cortex seperately
        function im = nonzeroAvgFiltCtx(im,h,midLineCol)
            % Note, because the fold in regions are connected, so average
            % filter need to be done seperately
            % Assumption: the midline is in the middle pixel 
            % Input & output:   im, mat, cortical value to be filtered
            %           h, logical mat, filter
            %           midLineCol, vector, column location for midline
            
            % Midline to split both hemisphere
            sz = size(im);
            
            for j = 1:2
                
                if j == 1
                    col = 1:midLineCol(1);
                elseif j == 2
                    col = midLineCol(1)+1:sz(2);
                end
                
                jIm = im(:,col,:);
                
                jIm = TBS.nonzeroAvgFilter(jIm,h);
                
                im(:,col,:) = jIm;
            end            
        end
       
        %% Function:    getRefPlateContour
        % Discription:  Calculate distance along contour on reference plate
        function xyz = getRefPlateContour(xyzInital,xyzRest,refScale)
            % Input:    xyzInital, mat, initial dots coordinates, ML 
            %           xyzRest, mat, dots to be assigned AP-vlaue
            % Output:   xyz, mat, dots with AP-value
            
            disp('Function getRefPlateContour on progress...');
                       
            % Set reference point value to 0 (5th column)
            xyzInital(:,end+1) = 0;
            xyz = {xyzInital};
            
            ixyz = xyz{end}; 
            while ~isempty(xyzRest) && ~isempty(ixyz)
                ixyz = xyz{end};
                
                % Find the closest alligned dot, within 2-voxel distance
                [D,I] = pdist2(ixyz(:,1:3),xyzRest(:,1:3),'euclidean','Smallest',5);
                TF = D <= 2;
                
                if size(ixyz,2) > 4
                    % Find the reference dot (close range) with closest ML-value
                    v = ixyz(:,4);
                    v = v(I);
                    v(~TF) = nan;
                    v = abs(v - xyzRest(:,4)');
                    TF = v == min(v,[],1,'omitnan');
                end
                
                % If there is multiple hit, find the one with shortest distance
                D(~TF) = nan;
                TF = D == min(D,[],1,'omitnan');
                
                % I, index in ixyz; col, index in xyzRest
                I(~TF) = 0;
                [~, col,I] = find(I);
                
                % Distance in micron
                D = ixyz(I,1:3)- xyzRest(col,1:3);
                D = sqrt(sum(D.^2,2));
                D = D./refScale;
                
                % Add distance to the reference dot
                D = D + ixyz(I,end);
                
                % Take mean distance if there is multiple dot
                [col,~,ic] = unique(col);
                D = TBS.accumarrayMean(ic,D);
                                
                xyz{end+1} = [xyzRest(col,:),D];
                
                % Delete the assigned dots
                xyzRest(col,:) = [];
            end
            
            xyz = vertcat(xyz{:});
        end
                                    
        %% Function:    refPlate2Column
        % Discription:  Assign colum with reference value
        function refPlate = refPlate2Column(refPlate,xyz)
            % Input:    refPlate, mat, volumn with reference value
            %           xyz, cell, coordinates of cortical colume
            % Output:   mat, volume with value assigned to colume
            
            disp('Function refPlate2Column on progress...');
            
            % Get reference value
            % (Flat to mat for speed)
            sz = cellfun(@(X) size(X,1),xyz);
            v = vertcat(xyz{:});
            v = TBS.xyz2v(v,refPlate);
            v = mat2cell(v,sz,1);
            % Mean of all the non-zero value
            v = cellfun(@(X) mean(nonzeros(X)),v);
            
            % Delete column with no reference value
            TF = ~isnan(v);
            xyz = xyz(TF);
            v = v(TF);
            
            % Assign to all the voxel of the column
            v = TBS.repmat2cell(v,xyz);
            
            xyz = vertcat(xyz{:});
            v = vertcat(v{:});
            
            % Get the mean if multiple value assign to a voxel
            [xyz,~,ic] = unique(xyz,'rows');
            % (this is faster than doing accumarry & @mean)     
            
            % Mean value for each xyz
            v = TBS.accumarrayMean(ic,v);
            
            refPlate = TBS.xyzv2im(size(refPlate),xyz,v);
        end
                                     
    end
       
    methods (Static)    % Data processing using flatmap ===================
        %% Function:    stack2flatmapIm
        % Description:  convert the image stack to flatmap-stack
        function flatStack = stack2flatmapIm(ind,V,ctxML,ctxAP,...
                ctxDepthPrctile,method,refScale)
            % Input:    ind, vector, index of inquery points
            %           V, vector, value of inquery points
            %           ctxML, stack, lookup table with ML-value for each voxel
            %           ctxAP, stack, lookup table with AP-value for each voxel
            %           ctxDepthPrctile, stack, lookup table with depth prectile for
            %           each voxel
            %           method, str, method to combine voxels
            %           refScale, num, scale for reference map, pixel per micron
            % Output:   flatStack, mat, faltmap stack using ML/AP/Depth coordinates
            %           same scale as reference map
            
            % Change ML to flatmap x-coordinates on image (>0)
            % Find left image-right hemiphere, change to minus
            TF = TBS.isLeftIm(ctxML);
            ctxML(~TF) = ctxML(~TF).*(-1);
            ctxML = ctxML - min(ctxML,[],'all')+1;
            
             % Change to same scale as reference map
            ctxML = ctxML.*refScale;
            ctxAP = ctxAP.*refScale;
            
            % Transfer depth precetile to length, assume depth is 1000 um
            ctxDepthPrctile = ctxDepthPrctile./100.*1000.*refScale;         
                        
            sz = [max(ctxAP,[],'all'),max(ctxML,[],'all'),max(ctxDepthPrctile,[],'all')];
            sz = round(sz);
            
            % Get ML,AP,Depth prectile value
            ML = ctxML(ind);
            AP = ctxAP(ind);
            depth = ctxDepthPrctile(ind);
            
            % Get voxel for the flatmap
            xyz = [ML,AP,depth];
            xyz = round(xyz);
            [xyz,~,ic] = unique(xyz,'rows');
            % voxel intensity compute using the method
            v = accumarray(ic,V,[],method);            
            
            flatStack = TBS.xyzv2im(sz,xyz,v);            
        end
        
        %% Function:    getMLAPD
        % Discription:  get ML/AP/depth coordinates
        function mlapd = getMLAPD(xyz,ctxML,ctxAP,ctxDepthPrctile)
            % Input:    xyz, mat, corrdinates in reference framework
            %           ctxML/AP/DepthPrctile, mat, reference value of
            %           corteical ML/AP/DepthPrecitle
            % Output:   mlapd, mat, coordiantes in ML/AP/Depth
                        
            disp('Function getMLAPD on progress...')
            
            % For cell input
            cellTF = iscell(xyz);            
            if cellTF
                sz = cellfun(@(X) size(X,1),xyz);
                xyz = vertcat(xyz{:});
            end
            
            mlapd = zeros(size(xyz));
            
            % Coordinates outside the cortex
            xyz2 = round(xyz);
            ind = sub2ind(size(ctxML),xyz2(:,2),xyz2(:,1),xyz2(:,3));
            TF = ctxML(ind)~=0;
            
            if ~any(TF)
                
                if cellTF
                    mlapd = mat2cell(mlapd,sz,3);
                end
                
                return
            end
            
            % For coordinates inside the cortex ---------------------------
            xyz = xyz(TF,:);
            
            % Add padding regions for interp coordinates close to the edge
            % Padding repeat value, otherwise error message 'Insufficient
            % finite values to interpolate.'
            ctxML = prepareCtxV(ctxML);
            ctxAP = prepareCtxV(ctxAP);
            ctxDepthPrctile = prepareCtxV(ctxDepthPrctile);
            
            % Change the right hemisphere to minus ML
            leftIm = TBS.isLeftIm(ctxML);
            ctxML(~leftIm) = ctxML(~leftIm).*(-1);
            
            ML = interp3(ctxML,xyz(:,1),xyz(:,2),xyz(:,3));
            AP = interp3(ctxAP,xyz(:,1),xyz(:,2),xyz(:,3));
            detphPrctile = interp3(ctxDepthPrctile,xyz(:,1),xyz(:,2),xyz(:,3));
            
            mlapd(TF,:) = [ML,AP,detphPrctile];
            
            % Change to cell output if applicable
            if cellTF
                mlapd = mat2cell(mlapd,sz,3);
            end
            
            % Function: prepareCtxV ---------------------------------------
            % Discription: Prepare image with cortex value, for
            % interpolation (for edge effect)
            function ctxV = prepareCtxV(ctxV)
                SE = strel('sphere',2);
                SE = SE.Neighborhood;
                
                TFv = ctxV == 0;
                % nonzeroAvgFilter(im,h)
                ctxV2 = TBS.nonzeroAvgFilter(ctxV,SE);
                ctxV(TFv) = ctxV2(TFv);
            end
        end
        
        %% Function:    plotRegionRef
        % Discription:  plot flatmap region boundary and soma area
        function plotRegionRef(mlapdSoma,regionOutlineFlat)
            % Input:    mlapdSoma, mat, soma ML/AP/Depth coordinates
            %           regionOutlineFlat, region boundary ML/AP/Depth
            %           coorindates
            % Output:   h, object, handel for the plot
            
            TF = any(mlapdSoma,2);
            mlapdSoma = mlapdSoma(TF,1:2);
            
            % Injection center and radius
            injCenter = median(mlapdSoma,1);
            D = pdist2(mlapdSoma,injCenter);
            r = prctile(D',95);
            
            hold on;
            % Plot region boundaries
            scatter(regionOutlineFlat(:,1),regionOutlineFlat(:,2),1,...
                'k','filled','MarkerFaceAlpha',0.25);
            % Plot it as a disk
            rectangle('Position',[injCenter-[r r], [r r].*2],...
                'Curvature',[1 1],'FaceColor','k')
            
            g = gca; g.YDir = 'reverse';
            % Flatmap2
            g.XLim = [-9000 9000];
            daspect([1 1 1]); set(gcf,'Position',[100 100 600 300]);
        end
        
        %% Function:    flatmapSetting
        % Discription:  setting for scatter plot of cortical flatmap
        function flatmapSetting(ylim)
            % Input:    ylim, vector, min & max limit for flatmap y-axis
            
            xlabel('ML (\mum)'); ylabel('AP (\mum)');
            
            g = gca; g.YLim = [min(ylim),max(ylim)]; 
            g.XTickLabel = []; g.YTickLabel = [];
            
            set(gcf,'Position',[100 100 900 200]);
        end
        
    end
    
    methods (Static)    % Cortical analysis ===============================
        %% Function:    BCtable2cell
        % Description:  change BCtable to cell with one cell per bc
        function BCcell = BCtable2cell(BCtable,varName)
            % Input:    BCtable, table, one image per row
            %           varName, str, content of the variable into the cell
            % Output:   BCcell, cell, info from one BC per cell
            
            var = BCtable.(varName);
            id = BCtable.codeID;
            
            var = vertcat(var{:});
            id = vertcat(id{:});
            
            nBC = max(id);
            
            BCcell = cell(nBC,1);
            for i = 1:nBC
                row = id == i;
                BCcell{i} = var(row,:);
            end            
        end
        
        %% Function:    nearSomaExcl
        %  Discription: exclude rolony within a range of injection center
        function mlapdDot = nearSomaExcl(mlapdDot,mlapdSoma,p,isCtrl)
            % Input:    mlapdDot, cell, ML/AP/Depth coordinates of rolony
            %           mlapdSoma, mat, ML/AP/Depth coordinates of soma, for exclusion
            %           p, num, min percentage for exclusion
            %           isCtrl, logical, whether is a control including deleting
            %           similar area in the other side
            % Output:   mlapdDot, cell
            
            % Injection center (median) and radius
            TF = any(mlapdSoma,2);
            xy = mlapdSoma(TF,1:2);
            injCenter = median(xy,1);
            r = pdist2(xy,injCenter);
            r = prctile(r',p);
            
            % Delete dots closed to the injection site
            fh = @(X,Y) pdist2(X(:,1:2),injCenter)>= r;
            mlapdDot = cellfun(@(X) X(fh(X),:),mlapdDot,'UniformOutput',false);
            
            % LatC local-exclusion control
            if isCtrl
                fh = @(X,Y) pdist2(X(:,1:2),injCenter.*[-1 1])>= r;
                mlapdDot = cellfun(@(X) X(fh(X),:),mlapdDot,'UniformOutput',false);
            end
            
            % Delete 0
            mlapdDot = cellfun(@(X) X(any(X,2),:),mlapdDot,'UniformOutput',false);
        end
        
        %% Function:    inCtxCluster
        % Discription:  whether a rolony is within a cortical cluster
        function [clusterRng, localRng] = inCtxCluster(mlapdDot,center,...
                rng,figOutTF,mlapdSoma,regionOutlineFlat)
            % Input:    mlapdDot, cell, mlapd coordinates of rolonies
            %           center, vector, center of cluster
            %           rng, num, range to be counted within cluster
            %           figureOutTF, logical, whether have figure output
            %           mlapdSoma, mat, mlapd coordinates of soma location
            %           regionOutlineFlat, for region outline in annotation map
            % Output:   clusterRng, cell of logical, whether the rolony within the
            % cluster
            %           localRng, cell of logical, whether the rolony within the
            %           cluster or the surrounding area
            
            surroundRng = rng*sqrt(2);
            
            % Project to the cluster
            % Distance to the cluster center
            D = cellfun(@(X) pdist2(X(:,1:2),center),mlapdDot,'UniformOutput',false);
            
            % Within cluster and local (+surrounding) rolony density
            clusterRng = cellfun(@(X) X<= rng,D,'UniformOutput',false);
            localRng = cellfun(@(X) X<= surroundRng,D,'UniformOutput',false);
            
            if ~figOutTF
                return
            end
            
            % Visualize the cluster & surrounding area --------------------------------
            c = vertcat(clusterRng{:}) + vertcat(localRng{:});
            xy = vertcat(mlapdDot{:});
            TF = any(xy,2);
            
            figure; TBS.plotRegionRef(mlapdSoma,regionOutlineFlat);
            hold on; scatter(xy(TF,1),xy(TF,2),1,c(TF),'filled','MarkerFaceAlpha',0.05);
            
            % Within cluster-red; surrounding area-blue; otherwise-black
            colormap([0 0 0; 0 0 1; 1 0 0]);
            
            % Crop region of interest
            % Only include AP with data
            lim = [min(xy(TF,2)),max(xy(TF,2))];
            g = gca; g.YLim = lim;
            % Only include one hemisphere
            if center(1) > 0
                g.XLim(1) = 0;
            else
                g.XLim(2) = 0;
            end
            
            g.XTick = 0:2000:10000; 
            g.YTick = g.XTick;
            
            xlabel('ML (\mum)'); ylabel('AP (\mum)');
            TBS.axLabelSettings('Myriad Pro',12);
        end
        
        %% Function:    plotSoma
        % Description:  plot soma on ML-Depth axes
        function plotSoma(mlapdSoma,c)
            % Input:    mlapdSoma, mat, ML-AP-Depth for soma location on flatmap
            %           c, vector, index
            
            I = any(mlapdSoma,2);
            
            % Random shuffle
            I = find(I); I = TBS.shuffleRows(I);
            
            disp(['Total soma plotted: ',num2str(numel(I))]);
            
            figure; scatter(mlapdSoma(I,1),mlapdSoma(I,3),10,c(I),'filled');
            ylabel('Soma depth (%)','FontSize',15);
            g = gca; g.YDir = 'reverse'; g.YLim = [0 100]; g.XTick = [];
            set(gcf,'Position',[100 100 300 300]);
        end
        
        %% Function:    kmeansDepthHist
        % Discription:  group depth histocounts, sort index using depth
        function idx = kmeansDepthHist(X,k,C)
            % Input:    X, mat, data input
            %           k, number, output group number
            %           C, mat, start centroid
            % Output:   idx, vector, group number
            
            % X2 = cumsum(X,2);
            
            if nargin == 3
                idx = kmeans(X,k,'Start',C);
            elseif nargin == 2
                idx = kmeans(X,k,'Replicates',50);
            end
            
            % % Mean max column of the group
            [~,I] = max(X,[],2);
            
            meanI = accumarray(idx,I,[],@median);
            
            % Sort group using max column
            [~,I] = sort(meanI,'ascend');
            
            % Sorted rank
            [~,I] = sort(I,'ascend');
            
            % Use sorted rank as index
            idx = I(idx);            
        end
        
        %% Function:    depthHeatmapSetting
        % Discription:  figure settings for heatmap of projeciton histocounts
        function depthHeatmapSetting(h,X,edges)
            % Input:    h, obj, handle of the heatmap
            %           X, mat, input for heatmap
            %           edges, vector, edges for histogram
            
            h.GridVisible = 'off';
            h.CellLabelColor = 'none';
            
            % Axis setting
            h.XDisplayLabels = nan(size(X,2),1);
            % Only label the top and bottom for y-axis
            h.YDisplayLabels = nan(size(X,1),1);
            h.YDisplayLabels{1} = edges(1); h.YDisplayLabels{end} = edges(end);
            h.YLabel = 'Projection depth (%)';
            
            h.FontSize = 12;
            
            % Colormap, 0-10%
            colormap(flipud(gray)); caxis([0 10]);            
            S = struct(h); S.Colorbar.Label.String = 'Projection (%)';            
        end
        
        %% Function:    focalProjPct
        % Description:  evaludate the level of focal projection
        function D = focalProjPct(mlap,p)
            % 01262021: use mean distance of p% closest neighbors
            % Input:    mlap, mat, flatmap coordinates
            %           p, num, the clost precentage of rolony
            % Output:   D, vector, mean distance represent the focal
            % projection level of this cell
            
            mlap = mlap(:,1:2);
            
            % Number of rolony to pick this neuron
            n = size(mlap,1);
            n = round(n.*p);
            
            % Get the distance of the closet N rolony
            D = pdist2(mlap,mlap,'euclidean','Smallest',n+1);
            % Delete itself
            D = D(2:end,:);
            
            D = mean(D,1);
        end
        
    end
        
    methods (Static)    % Thalamus analysis ===============================
        %% Function:    withinStrThalFiber
        % Discription:  find whether rolony is within the str-thal fiber
        function TF = withinStrThalFiber(xyz)
            % Input:    X, mat, rolony coordinates xyzRef
            % Output:   TF, logical
            
            TF = [];
            if isempty(xyz)
                return;
            end
            
            roi = [330 360; 140 200; 270 300]';
            TF = xyz >= roi(1,:) & xyz <= roi(2,:);
            TF = all(TF,2);
        end        
        
        %% Function:    roiGroupPTCT
        % Discritpion:  grouping CTPT using fiber tract
        function idx = roiGroupPTCT(roiDot)
            % Note, group basing on fiber location, top-0, below-1
            % Input:    roiDot, cell, one row per bc; each cell is the
            % coordinates of dots within the roi
            % Output:   idx, vector, CT-0, PT-1, no dots in roi-nan;
            
            cellSize = cellfun(@(X) size(X,1),roiDot);
            
            roiDot = vertcat(roiDot{:});
            
            % Find the closest 10 roi dots as neighbors
            [~,I] = pdist2(roiDot,roiDot,'euclidean','Smallest',11);
            I = I(2:end,:); I = I';
            
            % Initial grouping
            idx = roiDot(:,2) > median(roiDot(:,2));
            
            % Grouping basing on same barcode and neigbors
            previousIdx = idx;
            for i = 1:100
                % Assign the dots with the neighbor mode idx for each roi dot
                neghborIdx = idx(I);
                idx = mode(neghborIdx,2);
                idx = mat2cell(idx,nonzeros(cellSize));
                
                % Find the mode for each BC, assign to all roi dots
                modeIdx = cellfun(@mode,idx,'Uniformoutput',false);
                % repmat2cell(matIn,cellIn)
                idx = TBS.repmat2cell(modeIdx,idx);
                idx = vertcat(idx{:});
                
                % Check convergence
                if previousIdx == idx
                    break
                else
                    previousIdx = idx;
                end
            end
            
            idx = mat2cell(idx,nonzeros(cellSize));
            idx = cellfun(@mode,idx);
            
            % Convert to the size of roi
            idx2 = nan(size(cellSize));
            idx2(cellSize > 0) = idx;
            idx = idx2;
        end
        
        %% Function:    extGroupPTCT
        % Discription:  grouping of CT and PT using neigboring rolonies
        function idxOut = extGroupPTCT(thalDot,idxIn)
            % Note 02132022: tried to use while loop to find convergence, is the same
            % Input:    thalDot, cell vector, coordinates per cell
            % Output:   idx, vector, group id for each cell
            
            cellSize = cellfun(@(X) size(X,1),thalDot);
            
            % input idx for every rolony
            idxIn = TBS.repmat2cell(idxIn,thalDot);
            idxIn = vertcat(idxIn{:});
            
            % Assign rolony idx using nearest 10 neigbors
            thalDot = vertcat(thalDot{:});
            [~,I] = pdist2(thalDot,thalDot,'euclidean','Smallest',11);
            I = I'; I = I(:,2:end);
            
            % Find the most frequent idx per rolony
            idxOut = idxIn(I);
            idxOut = mode(idxOut,2);
            
            % Find the most frequent idx per BC
            idxOut = mat2cell(idxOut,cellSize);
            idxOut = cellfun(@mode,idxOut);
        end
        
        %% Function:    dispThalGroupStat
        % Discription:  displate the grouping result of cell number
        function dispThalGroupStat(matIn)
            % Input:    matIn, logical, 2-3 groups
            
            for i = 1:size(matIn,2)
                
                n = matIn(:,i);
                n = sum(n);
                n = num2str(n);
                
                if i == 3
                    disp(['Num of non-group cells: ', n]);
                else
                    disp(['Num of group ',num2str(i),' cells: ', n]);
                end
            end
        end
                
        %% Function:    dispThalGroupCellNumStat
        function dispThalGroupCellNumStat(groupRow,inStr,inMb)
            
            fh = @(X) num2str(X);
            
            inStr = cellfun(@any,inStr);
            inMb = cellfun(@any,inMb);
                                    
            % Stats: porpotion per group
            disp(['Total cells for grouping: ',fh(size(groupRow,1))]);
            TBS.dispThalGroupStat(groupRow);
            
            % Stats: midbrain per group
            mbPos = groupRow & inMb;
            disp(['Total supCol-projecting BC: ', fh(sum(mbPos,'all'))]);
            TBS.dispThalGroupStat(mbPos)
            
            % Stats: str per group
            strPos = groupRow & inStr;
            disp(['Total str-projecting BC: ', fh(sum(strPos,'all'))]);
            TBS.dispThalGroupStat(strPos)            
        end
        
        %% Function:    thalFigSetting
        % Discription:  Setting for thalamic 3d figure
        function thalFigSetting
            alpha 0.2;
            xlabel('x'); ylabel('y'); zlabel('z');
            g = gca; g.YDir = 'reverse'; g.ZDir = 'reverse';
            g.XLim = [200 450]; g.YLim = [50 250];
            g.XTick = []; g.YTick = []; g.ZTick = [];
            daspect([1 1 1]); set(gcf, 'Position',  [100, 100, 500, 300]);
            grid off; view([-10 10]);
        end
        
        %% Function:    visThalGroup
        % Discription:  visualize the thalamic rolony of each group
        function visThalGroup(xyzCell,idx,omitnan)
            % Input:    cf, cell, rolony coordinates
            %           idx, vector,
            %           omitnan, logical, whether omit nan
            
            xyz = vertcat(xyzCell{:});
                        
            idx2 = TBS.repmat2cell(idx,xyzCell);
            idx2 = vertcat(idx2{:});
            if omitnan == 0
                idx2(isnan(idx2)) = 3;
            end
            
            % Row shuffle
            [xyz,I] = TBS.shuffleRows(xyz);
            idx2 = idx2(I);
            
            figure; scatter3(xyz(:,1),xyz(:,2),xyz(:,3),2,idx2,'filled');
            TBS.thalFigSetting;            
        end
       
        %% Function:    groupCTPT (main)
        % Discription:  group thal+ cells into CTPT
        function idx2 = groupCTPT(xyzDot,thalReg,strReg,mbReg)
            % Note 03232022: no count threshold was used in this function
            % NaN: cells with no projection to thalamus; 0: CT; 1: PT
            % Input:    xyzDot, cell, xyz coordinates in reference
            % framework
            %           thalReg/strReg/mbReg, logical stack, thal/str/mb
            %           area in registrated framework
            % Output:   idx2, vector, group number
            
            disp('Grouping thalamus+ cells....')
            
            % Settings
            
            % Colormap: corn
            % 1-CT; 2-PT; 3-NaN
            cCorn = [1 0.8 0; 0.25 0.75 0.1; 0 0 0];
                                               
            % Get CF cells ================================================
            % Dots within brain regions
            fh = @(Y) cellfun(@(X) TBS.xyz2v(X,Y),xyzDot,'UniformOutput',false);
            inThal = fh(thalReg);
            inStr = fh(strReg);
            inMb = fh(mbReg);
            
            % COI: cell of interest: thal+
            COI = cellfun(@any,inThal);      
            
            xyzCOI = xyzDot(COI);
            inThal = inThal(COI);
            inStr = inStr(COI);
            inMb = inMb(COI);
            
            % Group CF cells within roi ===================================
            % Dots within ROI
            roiDot = cellfun(@(X) X(TBS.withinStrThalFiber(X),:),...
                xyzCOI,'UniformOutput',false);
            
            % Split cells into groups using roi
            % group basing on fiber location, top-0, below-1
            idx = TBS.roiGroupPTCT(roiDot);
            
            % Stats -------------------------------------------------------
            % Stats: cell number for each group (0, 1, NaN)
            disp('Stats of grouping of all cells: ');
            row = idx == [0 1];
            row = [row,isnan(idx)];
            
            TBS.dispThalGroupCellNumStat(row,inStr,inMb)
            
            % SupFig 2. Plot the dots within roi --------------------------
            % Black: all thal dots; red: portion within roi
            
            xyzCOI2 = vertcat(xyzCOI{:});
            roiDot2 = vertcat(roiDot{:});
            
            figure; scatter3(xyzCOI2(:,1),xyzCOI2(:,2),xyzCOI2(:,3),2,'k','filled');
            hold on; scatter3(roiDot2(:,1),roiDot2(:,2),roiDot2(:,3),...
                2,'r','filled');
            TBS.thalFigSetting;
            
            % SupFig2. Visualizing grouping result ------------------------
            % visThalGroup(cf,idx,omitnan);
            TBS.visThalGroup(xyzCOI,idx,0); colormap(cCorn);
            
            % Extent the grouping to the BC with no rolony in roi =========
            
            % Thalamic rolonies
            thal = cellfun(@(X,Y) X(Y,:),xyzCOI,inThal,'UniformOutput',false);
            
            % expand idx to include cells with no rolony in roi
            % extGroupPTCT(thalDot,idxIn)
            idx2 = TBS.extGroupPTCT(thal,idx);
            
            % Report change group id
            TF = ~isnan(idx);
            disp(['ExtGroupPTCT has different grouping result: ',...
                num2str(sum(idx(TF)~=idx2(TF))),'; in ',num2str(sum(TF))]);
            
            idx = idx2;
            
            % Stats -------------------------------------------------------                      
            disp('Stats of grouping of all cells, after extension: ');
            row = idx == [0 1];
            row = [row,isnan(idx)];
            
            TBS.dispThalGroupCellNumStat(row,inStr,inMb)
            
            % Fig2. Visualizing grouping result ---------------------------
            % visThalGroup(cf,idx,omitnan);
            TBS.visThalGroup(xyzCOI,idx,0); colormap(cCorn(1:2,:));
            
            % Visualize individual group
            idx2 = idx;
            idx2(idx2 == 1) = nan;
            TBS.visThalGroup(xyzCOI,idx2,1); colormap(cCorn(1,:));
            
            idx2 = idx;
            idx2(idx2 == 0) = nan;
            TBS.visThalGroup(xyzCOI,idx2,1); colormap(cCorn(2,:));
            
            % Convert CTPT grouping to all the BC -------------------------
            idx2 = nan(size(xyzDot,1),1);
            idx2(COI) = idx;
        end
        
    end
    
    methods (Static)    % Barcoded cell reconstruction ====================
        %% Function:    pdistCluster
        % Discription:  Find the shortest distance between cluster
        function [D,I] = pdistCluster(dotXYZ,clusterID)
            % Input:    dotXYZ, mat, dot coordinates in xyz
            %           clusterID, cell, one cluster per cell
            % Output:   D, mat, min distance between cluster
            %           I, cell, pair of dots with min distance in each
            %           cluster; 1st column, current cluster (row); 2nd
            %           column, the index of other cluster
            
            D = []; I = {};
            for i = 1:size(clusterID,1)
                % Dot coordinates of the current cluster
                iCluster = clusterID{i};
                iCluster = dotXYZ(iCluster,:);
                
                % Shortest distance of dots in other cluster to every dot of the
                % current cluster-iD, dot index in the partner cluster-Ip
                [iD,Ip] = cellfun(@(X) pdist2(dotXYZ(X,:),iCluster,'euclidean',...
                    'Smallest',1),clusterID,'Uniformoutput',false);
                
                % Minimum distance in each cluster
                % and the cooresponding dot index of current cluster
                [iD,Ic] = cellfun(@(X) min(X),iD,'Uniformoutput',false);
                
                % Dot index in the partner cluster with short distance
                Ip = cellfun(@(X,Y) X(Y),Ip,Ic,'UniformOutput',false);
                
                iD = cell2mat(iD);
                % Index pair: current (Ic)-partner (Ip)
                Ic = cell2mat(Ic); Ip = cell2mat(Ip);
                
                % Find cluster with minimum distance
                iD(i) = inf;     % Exclude itself
                % iD: min distance; iI: idx in clusterID
                [iD,iI] = min(iD);
                
                % Get the cooresponding cell using iI
                Ic = Ic(iI); Ip = Ip(iI);
                
                % Convert to dot ID   
                Ic = clusterID{i}(Ic); Ip = clusterID{iI}(Ip); 
                
                D(i) = iD; I{i} = [Ic,Ip];
            end
        end
        
        %% Function:    uniqueCluster
        % Discription:  delete repeat cluster with identical ID
        function clusterID = uniqueCluster(clusterID)
            % Input/Output:    clusterID, cell, with ID for each cluster
            
            i = 1;
            % (The last one doesnt have similar below)
            while i < size(clusterID,1)
                iCluster = clusterID{i};
                
                % Find the other cluster complete match the current clusterID
                repeatCluster = cellfun(@(X) all(ismember(X,iCluster)),clusterID);
                repeatCluster = find(repeatCluster);
                
                % Delete all the other cluster
                % Skip the itself-first one
                repeatCluster = repeatCluster(2:end);
                clusterID(repeatCluster) = [];
                
                i = i+1;
            end
        end
        
        %% Function:    fuseCluster
        % Disceirption: fuse cluster with shared id
        function clusterID = fuseCluster(clusterID)
            % Input/Output:    clusterID, cell, with ID for each cluster
            
            i = 1;
            while i < size(clusterID,1)
                iCluster = clusterID{i};
                
                % Find cluster with shared id
                partnerCluster = cellfun(@(X) any(ismember(X,iCluster)),clusterID);
                partnerCluster = find(partnerCluster);
                
                % 1st is itself, move on to the next cluster
                if numel(partnerCluster) == 1
                    i = i + 1;
                    continue
                end
                
                % If there is a partner, add it to current cluster
                % (Add one per iteration)
                partnerCluster = partnerCluster(2);
                
                % Add partner cluster to the current cluster
                clusterID{i} = [clusterID{i},clusterID{partnerCluster}];
                clusterID{i} = unique(clusterID{i});
                
                % Delete the origianl position of partner cluster
                clusterID(partnerCluster) = [];
            end
        end
        
        %% Function:    getLinePair
        % Discription:  get pairs of dots for rolony-line
        function linePair = getLinePair(dotXYZ,distThreshold)
            % Input:    dotXYZ, mat, dot coordinates
            %           distThreshold, num, max length threshold
            % Output:   linePair, cell
            
            tic;
            
            n = size(dotXYZ,1);
            
            % dot ID: row number
            dotID = (1:n)';
            
            % Cluster: group of dots connected together
            % Initiation: one cluster with one rolony
            clusterID = num2cell(dotID);
            
            linePair = {};
            while size(clusterID,1) > 1
                % Shortest distance between cluster
                % I, pair of dots belongs to a cluster
                [D,I] = TBS.pdistCluster(dotXYZ,clusterID);
                
                % Whether the shortest distance pass threshold
                row = D <= distThreshold;                
                if ~any(row)
                    break
                end
                
                I = I(row);
                
                % Add the new pairs to clusterID
                clusterID = [clusterID; I'];
                
                % Delete the initial single clusterID
                row = cellfun(@(X) numel(X),clusterID) > 1;
                clusterID = clusterID(row);
                
                % Delete repeat
                clusterID = TBS.uniqueCluster(clusterID);
                
                % Line for drawing lines between dots
                linePair = [linePair; clusterID];
                % Only the new paris have two dots, fused one has more than
                % 2, exclud fused clusters (dot number > 2)
                row = cellfun(@(X) numel(X),linePair) == 2;
                linePair = linePair(row,:);
                
                clusterID = TBS.fuseCluster(clusterID);
            end
            
            disp('Done getLinePair'); toc;
        end
        
        %% Function:    plotLineIm
        % Discription:  plot line image in voxel
        function im = plotLineIm(sz,dotXYZ)
            % Input:    sz, vector, image size
            %           dotXYZ, cell, dot coordinates for each line pair
            % Output:   im, logical stack, with lines   
            
            xyz = {};            
            for i = 1:numel(dotXYZ)
                iDot = dotXYZ{i};
                
                diffDot = diff(iDot,1);
                
                % interp1 basing on the max length
                len = ceil(max(abs(diffDot)));
                
                if len > 1
                    iDot =  mat2cell(iDot,2,[1 1 1]);
                    % 02142022: change 1:len to 0:len
                    iDot = cellfun(@(X) interp1([0 len],X,0:len)',iDot,...
                        'UniformOutput',false);
                    iDot = horzcat(iDot{:});
                end
                
               xyz{i,1} = iDot; 
            end
            
            xyz = vertcat(xyz{:});
            
            xyz = round(xyz);
            xyz = unique(xyz,'row');
            
            % Delete out of range dots (dots close to the edge due to scaling)
            TF = TBS.isOutOfRngeDot(sz,xyz);
            xyz = xyz(~TF,:);
            
            im = TBS.xyzv2im(sz,xyz,[]);
        end
        
        %% Function:    getBCmodel
        % Discription:  get voxel location for individual barcode
        function bcVoxel = getBCmodel(somaXYZ,dotXYZ,distThreshold,sz,dilateR)
            % Input:    somaXYZ, mat, soma locaiton coordinates in micron
            %           dotXYZ, mat, rolony location coordinates in micron
            %           distThreshold, num, max distance for rolony-line
            %           sz, row vector, output image size
            %           dilateR, num, dilation size for soma
            % Output:   bcVoxel, mat, voxel location for output,
                                    
            % 10202021: do not include soma if no soma location
            % Include soma during connecting dots
            somaTF = any(somaXYZ);
            
            if somaTF
                dotXYZ = [dotXYZ; somaXYZ];
            end
            
            linePair = TBS.getLinePair(dotXYZ,distThreshold);
                                   
            dotXYZ2 = cellfun(@(X) dotXYZ(X,:),linePair,...
                'UniformOutput',false);
            
            % Plot lines
            im = TBS.plotLineIm(sz,dotXYZ2);
            
            % Plot somas
            if somaTF
                imSoma = TBS.xyzv2im(sz,somaXYZ,[]);
                
                % Make soma visiable
                SE = strel('sphere',dilateR);
                imSoma = imdilate(imSoma,SE);
                
                im = im | imSoma;
            end
            
            % Compress image into sparse matrix
            [y,x,z,~] = TBS.find3(im);
            
            % Coordinates for barcode voxel
            bcVoxel = [x y z];            
        end
                               
        %% Function:    BCmodel (main)
        % Discription:  get model of barcoded neurons for visualization
        function [im, outline] = BCmodel(xyzDot,xyzSoma,c,scaleFactor,annoMap,refSetting)
            % Input:    xyzDot, cell, rolony coordinates in registrated xyz-space
            %           xyzSoma, mat, soma coordinates in registrated xyz-space
            %           c, vector, color index
            %           scaleFactor, num, scale factor for output
            %           annoMap, mat stack, CCF annotation image stack, for brain
            %           outline
            %           refSetting, struct
            % Output:   im, mat, image stack for BC model (scaled)
            %           outline, mat, brain outline (scaled)
            
            % Threshold for combine into a cluster (um)
            distThreshold = 1000;
            % Convert to micron/pixel in registrated xyz-space
            distThreshold = distThreshold.*refSetting.refScale;
            
            % Voxel dilation for soma
            dilateR = 4./scaleFactor;
            
            % Scaling
            xyzSoma = xyzSoma.*scaleFactor;
            xyzDot = cellfun(@(X) X.*scaleFactor,xyzDot,'Uniformoutput',false);
            distThreshold = distThreshold.*scaleFactor;
            dilateR = dilateR.*scaleFactor; dilateR = ceil(dilateR);
            
            % Brain outline
            outline = TBS.getBrainOutline(annoMap,5);
            outline = imresize3(outline,scaleFactor);
            
            sz = size(outline);
            
            im = {};
            parfor i = 1:numel(xyzDot) % parfor
                % getBCmodel(somaXYZ,dotXYZ,distThreshold,sz,dilateR)
                im{i,1} = TBS.getBCmodel(xyzSoma(i,:),xyzDot{i},distThreshold,sz,dilateR);
            end
            
            % Color for individual barcode cell
            if isempty(c)
                % BC id as color
                c = (1:numel(xyzDot))';
            end
            c = TBS.repmat2cell(c,im);
            
            im = vertcat(im{:});
            c = vertcat(c{:});
            
            % Assemble to image
            im = TBS.xyzv2im(sz,im,c);
            
            % Thicker axon for visualization
            SE = strel('sphere',1);
            im = imdilate(im,SE);
            
        end
        
    end
    
end