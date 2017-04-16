classdef BioformatsImage
    %BioformatsImage class
    
    properties (AbortSet)   %Filename and series number
        
        filename = '';
        series = 1;
        
    end
    
    properties (Dependent)  %Image file attributes
        width
        height
        
        seriesCount
        
        sizeZ
        sizeC
        sizeT
        
        channelNames
        
        pxSize
        pxUnits
    end
    
    properties (Transient, Hidden, SetAccess = private)
        bfReader
        metadata
    end
    
    properties (Access = private, Hidden)   %Bioformats toolbox URL
        bfTbxURL = 'http://downloads.openmicroscopy.org/bio-formats/5.4.0/artifacts/bfmatlab.zip';
        
        thisVersion = '0.9.0';
    end
    
    methods %Constructor, getters and setters
        
        function obj = BioformatsImage(varargin)
            %Class constructor
            
            if nargin > 0
                obj.filename = varargin{1};
            end
            
        end
        
        function obj = set.filename(obj,filename)
            
            %Verify that filename is a string
            if ~ischar(filename)
                error('BioformatsImage:FilenameNotString',...
                    'Input filename must be a string.');
            elseif isempty(filename)
                error('BioformatsImage:NoFilename',...
                    'Please specify an image file.');
            end
            
            %Check that the file exists
            if ~exist(filename,'file')
                error('BioformatsImage:CannotFindFile',...
                    'Could not find %s. Provide the full filename or make sure the image folder is on the MATLAB path.',filename);
            end
            
            %Make sure that the full path to the file is stores
            [fPath,fName,fExt] = fileparts(filename);
            
            if ~isempty(fPath)
                %Reconstruct the full path (this will correct any
                %system-dependent file separator strings)
                obj.filename = fullfile(fPath,[fName,fExt]);
                
            else
                %Since the file must exist on the MATLAB path (we checked
                %this above), we can use which to determine the full path
                %to the file
                obj.filename = which(filename);
                
            end
            
            %Get a reader for the image
            obj = obj.getReader;
            
        end
        
        function obj = set.series(obj,iS)
            %Set series number
            
            if isempty(obj.bfReader) %#ok<MCSUP>
                obj = obj.getReader;
            end
            
            if rem(iS,1) ~= 0
                error('BioformatsImage:SeriesNumberMustBeInteger',...
                    'Expected series number to be an integer.')
            end
            
            if iS <= 0 || iS > obj.seriesCount %#ok<MCSUP>
                error('BioformatsImage:RequestedSeriesOutOfBounds',...
                    'Requested series must be > 0 and <= than %d.',obj.seriesCount); %#ok<MCSUP>
            end
            
            obj.bfReader.setSeries(iS - 1); %#ok<MCSUP>
                        
        end
        
        function seriesOut = get.series(obj)
            
            if ~obj.bfReaderExist
                seriesOut = NaN;                
            else
                seriesOut = obj.bfReader.getSeries + 1;                
            end
            
        end
        
        function height = get.height(obj)
            
            if ~isempty(obj.bfReader)
                height = obj.bfReader.getSizeY;
            end
            
        end
        
        function width = get.width(obj)
            
            if ~isempty(obj.bfReader)
                width = obj.bfReader.getSizeX;
            end
            
        end
        
        function sizeZ = get.sizeZ(obj)
            
            if ~isempty(obj.bfReader)
                sizeZ = obj.bfReader.getSizeZ;
            end
            
        end
        
        function sizeC = get.sizeC(obj)
            
            if ~isempty(obj.bfReader)
                sizeC = obj.bfReader.getSizeC;
            end
            
        end
        
        function sizeT = get.sizeT(obj)
            
            if ~isempty(obj.bfReader)
                sizeT = obj.bfReader.getSizeT;
            end
            
        end
        
        function seriesCount = get.seriesCount(obj)
           
            seriesCount = obj.bfReader.getSeriesCount;
            
        end
        
        function metadata = get.metadata(obj)
            
            metadata = obj.bfReader.getMetadataStore;
            
        end
        
        function chNames = get.channelNames(obj)
            chNames = cell(1,obj.sizeC);
            
            for iC = 1:obj.sizeC
                chNames{iC} = obj.metadata.getChannelName(0, iC - 1).char;
            end
        end
        
        function pxSize = get.pxSize(obj)
            
            pxSizeX = double(obj.metadata.getPixelsPhysicalSizeX(0).value);
            pxSizeY = double(obj.metadata.getPixelsPhysicalSizeY(0).value);
            
            pxSize = [pxSizeX, pxSizeY];
            
        end
        
        function pxUnit = get.pxUnits(obj)
            
            pxInfoStr = char(obj.metadata.getPixelsPhysicalSizeX(0).toString);
            
            %Find the 'unit[' pattern
            unitStartIdx = strfind(pxInfoStr,'unit[');
            %Offset by 5 characters to get start of unit string
            unitStartIdx = unitStartIdx + 5;    
            
            %Find the closing brace, starting from the unit string
            unitEndIdx = strfind(pxInfoStr(unitStartIdx:end),']');
            %Reduce by 1 character to not include the closing brace
            unitEndIdx = unitEndIdx - 2;
            
            %Output the pixel unit
            pxUnit = pxInfoStr(unitStartIdx: unitStartIdx + unitEndIdx);
            
        end

    end
    
    methods %Base functions
        
        function obj = getReader(obj)
            %Get a Bioformats Reader object
            %
            % You should not need to call the function again, but I'm
            % leaving it here in case it is necessary to reboot it.
            
            %Check if the Bioformats toolbox is installed. If not, install
            %it from the web.
            if ~obj.bfTbxExist
                obj.installBFtbx;
            end
            
            if ~isempty(obj.filename)
                obj.bfReader = bfGetReader(obj.filename);
            else
                error('BioformatsImage:EmptyFilename',...
                    'Set filename first.');
            end
            
            %Set current series to 1
            obj.series = 1;
            
        end
        
        function [imgOut, timestamp] = getPlane(obj,varargin)
            %Get image at specified index
            %
            %  imgOut = getImage(iZ, iC, iT)
            %
            %  imgOut = getImage(iS, iZ, iC, iT);
            %
            %  imgOut = getImage(iS, iZ, iC, iT, ROI);
            %
            %  ROI = [XMIN YMIN WIDTH HEIGHT];
            
            ip = inputParser;
            ip.addRequired('iZ',@(x) x > 0 && x <= obj.sizeZ);
            ip.addRequired('iC',@(x) x > 0 && x <= obj.sizeC);
            ip.addRequired('iT',@(x) x > 0 && x <= obj.sizeT);
            ip.addOptional('iS', [], @(x) x > 0 && x <= obj.seriesCount);
            ip.addParameter('ROI',[],@(x) numel(x) == 4 && all(x > 0));
            
            ip.parse(varargin{:});
            
            %Check that reader object already exists
            if ~bfReaderExist(obj)
                %Get a reader for the file
                obj = obj.getReader;
            end
            
            if ~isempty(ip.Results.iS)
                obj.series = ip.Results.iS;
            end            
            
            %Get actual image index
            bfIndex = obj.bfReader.getIndex(ip.Results.iZ - 1,...
                ip.Results.iC - 1, ip.Results.iT -1) + 1;
            
            %Get image
            if ~isempty(ip.Results.ROI)
                xMin = ip.Results.ROI(1);
                yMin = ip.Results.ROI(2);
                roiWidth = ip.Results.ROI(3);
                roiHeight = ip.Results.ROI(4);
                imgOut = bfGetPlane(obj.bfReader,bfIndex,xMin,yMin,roiWidth,roiHeight);
            else
                %Get full image
                imgOut = bfGetPlane(obj.bfReader,bfIndex);
            end
            
            if nargout > 1                
                %Get the timestamp
                timestamp = double(obj.metadata.getPlaneDeltaT(obj.series - 1,bfIndex).value);
            end
            
        end
        
        function [imgOut, timestamp] = getXYplane(obj, channel, XYloc, frame, varargin)
            %GETXYPLANE  Get plane from multi-XY images
            %
            %  I = O.GETXYPLANE(Channel, Location, Frame)
            %
            %  Assumes that Z = 1.
            % This is to handle the interleaving
            
            %Check that reader object already exists
            if ~bfReaderExist(obj)
                %Get a reader for the file
                obj = obj.getReader;
            end
            
            %Resolve the channel index
            iC = obj.channelname2ind(channel);
            
            %Calculate the image index
            imgIndex = XYloc + (frame - 1) * obj.seriesCount;
            
            %Calculate the actual series number and timepoint
            iS = floor(imgIndex/obj.sizeT);
            
            iT = imgIndex - (iS * obj.sizeT);
            
            %Get the image plane
            obj.series = iS;
            
            %Get the image (passing in the varargin, assumed to be ROI
            %information)
            imgOut = obj.getPlane(1, iC, iT, varargin);
            
            %Get timestamp if output is assigned            
            if nargout > 1
                %Resolve the bioformats index
                bfIndex = obj.bfReader.getIndex(0, iC - 1, iT - 1) + 1;
                
                %Get the timestamp
                timestamp = double(obj.metadata.getPlaneDeltaT(iS - 1,bfIndex).value);
            end
        end
        
        function timestamps = getTimestamps(obj,channel)
            
            iC = obj.channelname2ind(channel);
            
            timestamps = zeros(1, obj.sizeT);
            
            for iT = 1:obj.sizeT
                %Resolve the bioformats index
                bfIndex = obj.bfReader.getIndex(0, iC - 1, iT - 1) + 1;
      
                timestamps(iT) = double(obj.metadata.getPlaneDeltaT(obj.series - 1,bfIndex).value);
                
            end
        end
        
        function versionOut = version(obj)
            
            versionOut = sprintf('BioformatsImage class Version %s',obj.thisVersion);
            
        end
        
    end

    methods %Stack functions
        
        function imgOut = getZstack(obj, channel, timepoint, varargin)
            %Gets a z-stack for the specified channel and timepoint
            %
            %  I = O.GETZSTACK(Channel, Timepoint)
            %
            %  If Timepoint is not set, it defaults to 1.
            
            ip = inputParser;
            ip.addRequired('iC',@(x) isscalar(x) || ischar(x));
            ip.addOptional('iT', 1, @(x) x > 0 && x <= obj.sizeT);
            ip.addParameter('ROI',[],@(x) numel(x) == 4 && all(x > 0));
            ip.addParameter('ZRange',1:obj.sizeZ,@(x) min(x) > 0 && max(x) <= obj.sizeZ);
            
            ip.parse(channel, timepoint, varargin{:});
            
            %Get the results of the parser
            zRange = ip.Results.ZRange;
            roi = ip.Results.ROI;
            iT = ip.Results.iT;
            
            %Check that reader object already exists
            if ~bfReaderExist(obj)
                %Get a reader for the file
                obj = obj.getReader;
            end
            
            %Convert channel name to index (if it is already an index, the
            %function will just return the same number)
            iC = obj.channelname2ind(channel);
           
            %Initialize a storage matrix
            if ~isempty(ip.Results.ROI)
                imgOut = zeros(roi(4),roi(3),numel(zRange),'uint16');
            else
                imgOut = zeros(obj.height,obj.width,(zRange),'uint16');
            end
            
            for iZ = zRange
                if ~isempty(ip.Results.ROI)                
                    imgOut(:,:,iZ) = obj.getPlane(iZ, iC, iT, 'ROI', roi);                
                else
                    imgOut(:,:,iZ) = obj.getPlane(iZ, iC, iT);                
                end
            end
            
        end
        
        function imgOut = getChannel(obj, channel, varargin)
            %Gets all timepoints from the specified channel (and z-plane)
            ip = inputParser;
            ip.addRequired('Channel',@(x) isscalar(x) || ischar(x));
            ip.addOptional('iZ', 1, @(x) x > 0 && x <= obj.sizeZ);
            ip.addParameter('ROI',[],@(x) numel(x) == 4 && all(x > 0));
            
            ip.parse(channel, varargin{:});
                        
            %Check that reader object already exists
            if ~bfReaderExist(obj)
                %Get a reader for the file
                obj = obj.getReader;
            end
            
            %Convert channel name to index (if it is already an index, the
            %function will just return the same number)
            iC = obj.channelname2ind(ip.Results.Channel);
            
            roi = ip.Results.ROI;
            
            %Initialize a storage matrix
            if ~isempty(ip.Results.ROI)
                imgOut = zeros(roi(4),roi(3),obj.sizeT,'uint16');
            else
                imgOut = zeros(obj.height,obj.width,obj.sizeT,'uint16');
            end
            
            for iT = 1:obj.sizeT
                if ~isempty(ip.Results.ROI)
                    imgOut(:,:,iT) = obj.getPlane(ip.Results.iZ, iC, iT, 'ROI', roi);
                else
                    imgOut(:,:,iT) = obj.getPlane(ip.Results.iZ, iC, iT);
                end
            end
            
        end
        
        function imgOut = getTimepoint(obj, timepoint, varargin)
            %Get all channels from a specific timepoint (and z-plane)
            
            ip = inputParser;
            ip.addRequired('iT',@(x) x > 0 && x <= obj.sizeT);
            ip.addOptional('iZ', 1, @(x) x > 0 && x <= obj.sizeZ);
            ip.addParameter('ROI',[],@(x) numel(x) == 4 && all(x > 0));
            
            ip.parse(timepoint, varargin{:});
            
            %Check that reader object already exists
            if ~bfReaderExist(obj)
                %Get a reader for the file
                obj = obj.getReader;
            end
                        
            %Initialize a storage matrix
            roi = ip.Results.ROI;
          
            if ~isempty(ip.Results.ROI)
                imgOut = zeros(roi(4),roi(3),obj.sizeC,'uint16');
            else
                imgOut = zeros(obj.height,obj.width,obj.sizeC,'uint16');
            end
            
            for iC = 1:obj.sizeC
                if  ~isempty(ip.Results.ROI)
                    imgOut(:,:,iC) = obj.getPlane(ip.Results.iZ, iC, ip.Results.iT, 'ROI', roi);
                else
                    imgOut(:,:,iC) = obj.getPlane(ip.Results.iZ, iC, ip.Results.iT);
                end
            end
        end
        
        function MAPout = getZmap(obj,channel, varargin)
            %Get maximum amplitude projection of an image along the z-plane
            %
            % Default is 'maximum'
            %
            % Could also be 'mean' or 'median'
            
            mapIP = inputParser;
            mapIP.addRequired('channel',@(x) isnumeric(x) || ischar(x));
            mapIP.addOptional('iT',1, @(x) isscalar(x) && x > 0 && x <= obj.sizeT);
            mapIP.addParameter('binMode','max',@(x) ischar(x));
            mapIP.addParameter('ROI',[1 1 obj.width obj.height],@(x) numel(x) == 4);
            mapIP.addParameter('ZRange',1:obj.sizeZ,@(x) min(x) > 0 && max(x) <= obj.sizeZ);
            
            mapIP.parse(channel, varargin{:});
            
            %Check that reader object already exists
            if ~bfReaderExist(obj)
                %Get a reader for the file
                obj = obj.getReader;
            end
            
            zStack = obj.getZstack(mapIP.Results.channel, mapIP.Results.iT, 'ROI', mapIP.Results.ROI,'zRange',mapIP.Results.ZRange);
            
            switch lower(mapIP.Results.binMode)
            
                case {'max', 'maximum'}
                    MAPout = max(zStack,[],3);
                    
                case {'mean','average','avg'}
                    MAPout = mean(zStack,3);
                    
                case 'median'
                    MAPout = median(zStack,3);
                    
                case {'min', 'minimum'}
                    MAPout = min(zStack,[],3);
                    
                otherwise
                    error('BioformatsImage:getZmap:UnknownMode',...
                        'Mode should be either ''max'', ''mean'', or ''median''');
                    
            end
            
        end
        
    end
    
    methods  %Tiling functions
        
        function [tileDataOut roiOut] = getTile(obj, zcts, numTiles, tileIndex)
            
            %Parse the first input
            ip = inputParser;
            ip.addOptional('iZ',1)
            ip.addRequired('iC')
            ip.addOptional('iT',1);
            ip.addOptional('iS',obj.series);
            
            ip.addRequired('NumTiles',@(x) numel(x) == 2 && all(x > 0));
            ip.addRequired('TileRange', @(x) all(x > 0));
            
            ip.parse(zcts(:),numTiles,tileIndex);
            
            %Check that the tile range fits within the total number of
            %tiles
            if any(ip.Results.TileRange > prod(ip.Results.NumTiles))
                error('BioformatsImage:IndexExceedsTileDimensions',...
                    'Index exceeds tile dimensions. Maximum index is %d.',...
                    prod(ip.Results.NumTiles))
            end
            
            %Create a storage cell (for multiple tile data, which might have different sizes)
            tileDataOut = cell(numel(ip.Results.TileRange),1);
            
            %Create a storage for the roi
            roiOut = zeros(numel(ip.Results.TileRange),4);
            
            for iT = 1:numel(ip.Results.TileRange)
                tileIndex = ip.Results.TileRange(iT);
                roi = getTileIndices(obj, ip.Results.NumTiles, tileIndex);
                
                %Get the ROI
                tileDataOut{iT} = obj.getPlane(ip.Results.iC,...
                    ip.Results.iZ, ip.Results.iT, 'ROI', roi);
                
                roiOut(iT,:) = roiOut;
            end
            
            %If it's only a single tile, output a matrix rather than a cell
            if numel(tileDataOut) == 1
                tileDataOut = tileDataOut{:};
            end
                
        end
        
        function zMAPsOut = getTileZMap(obj, ct, numTiles, tileIndex, varargin)
            
            ip = inputParser;
            ip.addRequired('iC')
            ip.addOptional('iT',1);
%             ip.addOptional('iS',obj.series);
            
            ip.addRequired('NumTiles',@(x) numel(x) == 2 && all(x > 0));
            ip.addRequired('TileRange', @(x) all(x > 0));

            ip.addParameter('binMode','max',@(x) ischar(x));
            ip.addParameter('ZRange',1:obj.sizeZ,@(x) min(x) > 0 && max(x) <= obj.sizeZ);
            
            ip.parse(ct,numTiles,tileIndex,varargin{:});
            
            tileRange = ip.Results.TileRange;
            
            zMAPsOut = cell(numel(tileRange,1));
            for iT = 1:numel(tileRange)
                currTileIdx = tileRange(iT);
                
                %Resolve the roi
                roi = obj.getTileIndices(ip.Results.NumTiles, currTileIdx);
                
                %Get the z-stack MAP
                zMAPsOut{iT} = obj.getZmap(ip.Results.iC, ip.Results.iT, 'roi', roi,...
                        'binMode', ip.Results.binMode, 'ZRange', ip.Results.ZRange);
                
            end
            
            if numel(zMAPsOut) == 1
                zMAPsOut = zMAPsOut{:};
            end
               
        end
        
    end
    
    methods (Hidden, Static)    %Toolbox exists check
        
        function tbxExists = bfTbxExist
            %Check if the BioFormats toolbox is installed
            
            if ~exist('bfGetReader','file')
                
                %Check if there is a bfmatlab folder in the same folder
                if exist('bfmatlab','dir')
                    %Folder exists so add it to the path
                    addpath('bfmatlab');
                else
                    tbxExists = false;
                    return;
                end
            end
            
            tbxExists = true;
            
        end
        
    end
    
    methods (Access = private)
        
        function installBFtbx(obj)
            %Downloads and installs the Bioformats Toolbox
            
            %Create a temporary directory
            if ~exist('temp','dir')
                mkdir('temp')
            end
            
            %Download the toolbox from the OME website
            fprintf('Downloading Bioformats toolbox...');
            
            fn = websave('temp/bfmatlab.zip',obj.bfTbxURL);
            fprintf('Done.\n');
            
            fprintf('Installing toolbox...')
            %Unzip the toolbox
            unzip(fn);
            
            %Add toolbox to path
            addpath('bfmatlab')
            fprintf('Done.\n');
            
            %Remove temporary files and directory
            if exist('temp','dir')
                rmdir('temp','s')
            end
            
        end
       
        function bfrExists = bfReaderExist(obj)
            %Checks to see if the bfReader property contains a bfreader
            %object
            
            if ~isempty(obj.bfReader) && isa(obj.bfReader,'loci.formats.ChannelSeparator')
                bfrExists = true;
            else
                bfrExists = false;
            end
            
        end
        
        function chInd = channelname2ind(obj,channelIn)
            %Convert channel name to index
            
            if isnumeric(channelIn)
                chInd = channelIn;
                return
            end
            
            %Find the matching name
            chTF = strcmpi(channelIn, obj.channelNames);
            
            if ~any(chTF)
                error('BioformatsImage:ChannelNameNotFound',...
                    'Could not find channel %s.',channelIn);
            end
            
            chInd = find(chTF);
            
        end
        
        function roiOut = getTileIndices(obj, numTiles, tileIndex)
            %Given a tile size and index, calculate the resulting ROI
            %vector
            %ROI = [XMIN YMIN WIDTH HEIGHT];
            
            if numel(tileIndex) == 1
                %If it's a single index, resolve it into row,col
                %coordinates
                [indRow, indCol] = ind2sub(numTiles,tileIndex);
            else
                indRow = tileIndex(1);
                indCol = tileIndex(2);
            end
            
            %Calculate the tile height/width start and end indices
            %
            % Tile indices are:
            %     start floor(height/numrows) * (tile row number - 1) + 1
            %     end floor(height/numrows) * (tile row number)
            
            %Row:
            tileHeight = floor(obj.height/numTiles(1));
            
            rowStart = tileHeight * (indRow - 1) + 1;
            
            if indRow ~= numTiles(1)
                rowEnd = tileHeight * indRow;
            else
                %If it's the last row, then just use the image height
                rowEnd = obj.height;
            end
            
            %Convert to ROI height
            roiHeight = rowEnd - rowStart + 1;
            
            %Col:
            tileWidth = floor(obj.width/numTiles(2));
            
            colStart = tileWidth * (indCol - 1) + 1;
            
            if indCol ~= numTiles(2)
                colEnd = tileWidth * indCol;
            else
                %If it's the last row, then just use the image height
                colEnd = obj.height;
            end
            
            %Convert to ROI width
            roiWidth = colEnd - colStart + 1;
            
            roiOut = [colStart, rowStart, roiWidth, roiHeight];
        
        end
    end
        
end

%TODO: Tidy up ROI functions