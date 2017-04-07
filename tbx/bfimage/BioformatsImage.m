classdef BioformatsImage
    
    properties (AbortSet)   %Filename
        
        filename = '';
        
    end
    
    properties  %Current series number
        series
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
    
    methods %Public functions
        
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
            
            %Check that reader object already exists
            if ~bfReaderExist(obj)
                %Get a reader for the file
                obj = obj.getReader;
            end
            
            %Combine all inputs into a single vector
            imgIndex = cell2mat(varargin);
            
            %Validate the inputs
            if any(imgIndex <= 0)
                error('BioformatsImage:getImage:IndexIsZeroOrLess',...
                    'Values should be greater than 0.');
            end
            
            getROI = false;
            
            switch numel(imgIndex)
                
                case 3
                    %Input is [iZ, iC, iT]
                    
                    iZ = imgIndex(1);
                    iC = imgIndex(2);
                    iT = imgIndex(3);
                    
                    %Check if current series number is set. If not, default to
                    %series 1.
                    if isnan(obj.series)
                        obj.series = 1;
                    end
                    
                case 4
                    %Input is [iS, iZ, iC, iT];
                    
                    iS = imgIndex(1);
                    iZ = imgIndex(2);
                    iC = imgIndex(3);
                    iT = imgIndex(4);
                    
                    %Set series
                    obj.series = imgIndex(iS);
                    
                case 7
                    %Set ROI, no series change
                    
                    iZ = imgIndex(1);
                    iC = imgIndex(2);
                    iT = imgIndex(3);
                    
                    xMin = imgIndex(4);
                    yMin = imgIndex(5);
                    roiWidth = imgIndex(6);
                    roiHeight = imgIndex(7);
                    
                    %Check if current series number is set. If not, default to
                    %series 1.
                    if isnan(obj.series)
                        obj.series = 1;
                    end
                    
                    getROI = true;
                    
                case 8
                    
                    iS = imgIndex(1);
                    iZ = imgIndex(2);
                    iC = imgIndex(3);
                    iT = imgIndex(4);
                    
                    xMin = imgIndex(5);
                    yMin = imgIndex(6);
                    roiWidth = imgIndex(7);
                    roiHeight = imgIndex(8);
                    
                    getROI = true;
                    
                    %Set series
                    obj.series = imgIndex(iS);
                    
                otherwise                    
                    
                    error('BioformatsImage:getImage:InvalidArguments',...
                        'Invalid number of input arguments. Expected [(iS), iZ, iC, iT, (ROI)] (values in parentheses are optional).')
            end
            
            %Get actual image index
            bfIndex = obj.bfReader.getIndex(iZ - 1, iC - 1, iT -1) + 1;

            %Get image
            if getROI
                imgOut = bfGetPlane(obj.bfReader,bfIndex,xMin,yMin,roiWidth,roiHeight);
            else
                %Get full image
                imgOut = bfGetPlane(obj.bfReader,bfIndex);
            end
            
            if nargin > 1                
                %Resolve the bioformats index
                bfIndex = obj.bfReader.getIndex(0, iC - 1, iT - 1) + 1;
                
                %Get the timestamp
                timestamp = double(obj.metadata.getPlaneDeltaT(obj.series - 1,bfIndex).value);
            end
            
        end
        
        function [imgOut, timestamp] = getXYplane(obj, channel, XYloc, frame)
            %GETXYPLANE  Get plane from multi-XY images
            %
            %  I = O.GETXYPLANE(Channel, Location, Frame)
            %
            %  Assumes that Z = 1.
            
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
            imgOut = obj.getPlane([1, iC, iT]);
            
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
        
        function imgOut = getZstack(obj, channel, timepoint)
            %Gets a z-stack for the specified channel and timepoint
            %
            %  I = O.GETZSTACK(Channel, Timepoint)
            %
            %  If Timepoint is not set, it defaults to 1.
            
            %Check that reader object already exists
            if ~bfReaderExist(obj)
                %Get a reader for the file
                obj = obj.getReader;
            end
            
            %Convert channel name to index (if it is already an index, the
            %function will just return the same number)
            iC = obj.channelname2ind(channel);
            
            if ~exist('timepoint','var')
                iT = 1;
            else
                if timepoint <= 0 || timepoint > obj.sizeT
                    error('BioformatsImage:getZstack:InvalidTimepoint',...
                        'Timepoint must be >0 and <= %d.',obj.sizeT);
                end
                iT = timepoint;
            end
            
            imgOut = zeros(obj.height,obj.width,obj.sizeZ,'uint16');
            for iZ = 1:obj.sizeZ
                imgOut(:,:,iZ) = obj.getPlane([iZ, iC, iT]);                
            end
            
        end
        
        function imgOut = getChannel(obj, channel, zplane)
            %Gets all timepoints from the specified channel (and z-plane)
            
            %Check that reader object already exists
            if ~bfReaderExist(obj)
                %Get a reader for the file
                obj = obj.getReader;
            end
            
            %Convert channel name to index (if it is already an index, the
            %function will just return the same number)
            iC = obj.channelname2ind(channel);
            
            if ~exist('zplane','var')
                iZ = 1;
            else
                if zplane <= 0 || zplane > obj.sizeZ
                    error('BioformatsImage:getTstack:InvalidZplane',...
                        'Z-plane index must be > 0 and <= %d.',obj.sizeZ);
                end
                
                iZ = zplane;
            end
            
            imgOut = zeros(obj.height,obj.width,obj.sizeT,'uint16');
            for iT = 1:obj.sizeT
                imgOut(:,:,iT) = obj.getPlane([iZ, iC, iT]);
            end
            
        end
        
        function imgOut = getFrame(obj, timepoint, zplane)
            %Get all channels from a specific timepoint (and z-plane)
            
            %Check that reader object already exists
            if ~bfReaderExist(obj)
                %Get a reader for the file
                obj = obj.getReader;
            end
            
            %Validate the timepoint
            if timepoint <= 0 || timepoint > obj.sizeT
                error('BioformatsImage:getZstack:InvalidTimepoint',...
                    'Timepoint must be >0 and <= %d.',obj.sizeT);
            end
            iT = timepoint;

            %Validate the zplane
            if ~exist('zplane','var')
                iZ = 1;
            else
                if zplane <= 0 || zplane > obj.sizeZ
                    error('BioformatsImage:getTstack:InvalidZplane',...
                        'Z-plane index must be > 0 and <= %d.',obj.sizeZ);
                end
                
                iZ = zplane;
            end
            
            imgOut = zeros(obj.height,obj.width,obj.sizeC,'uint16');
            for iC = 1:obj.sizeC
                imgOut(:,:,iC) = obj.getPlane([iZ, iC, iT]);
            end
        end
        
        function MAPout = getZmap(obj,channel, timepoint, varargin)
            %Get maximum amplitude projection of an image along the z-plane
            %
            % Default is 'maximum'
            %
            % Could also be 'mean' or 'median'
            
            %Check that reader object already exists
            if ~bfReaderExist(obj)
                %Get a reader for the file
                obj = obj.getReader;
            end
            
            %Validate the timepoint input
            if exist('timepoint','var') && ~ischar(timepoint)
                if timepoint > 0 && timepoint <= obj.sizeT
                    iT = timepoint;
                else
                    error('BioformatsImage:getZmap:InvalidTimepoint',...
                        'Timepoint must be >0 and <= %d.',obj.sizeT);
                end
            else
               iT = 1;
            end
            
            %Set the operating mode
            if ischar(timepoint)
                %If the mode was specified instead of a timepoint
                mode = timepoint; 
            elseif ~isempty(varargin)
                %If the mode was specified at the end
                mode = varargin{1};
            else
                %Default mode is maximum
                mode = 'max';
            end               
            
            zStack = obj.getZstack(channel, iT);
            
            switch lower(mode)
            
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
    
    methods (Hidden, Access = private)    %Installing toolbox
        
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
        
    end
        
end














