classdef imager < detector
    %% Imaging camera
    %
    % imgr = imager(resolution) creates an imaging camera object from the
    % detector resolution
    %
    % Example:
%     tel = telescope(1,'resolution',21,'samplingTime',1);
%     imgr = imager(21);
%     src = source.*tel*imgr;
%     figure
%     imagesc(imgr.frame)

    properties
        % reference frame
        referenceFrame;
        % imaging lens
        imgLens;
        % imaging camera
        imgCcd;
        % Strehl ratio
        strehl;
        % entrapped energy
        ee;
        % entrapped energy slit width
        eeWidth;
    end
        
    properties (Access=private)
        % integration count;
        %frameCount=0;
        % telescope
        tel;
    end
    
    methods
        
        %% Constructor
        function obj = imager(in)
            if isa(in,'telescopeAbstract')
                resolution = in.resolution;
            elseif isnumeric(in)
                resolution = in;
            else
                error('oomao:imager','Inputer is either numeric or a telescope class')
            end
            obj = obj@detector(resolution);
            if isa(in,'telescopeAbstract')
                obj.tel = in;
                obj.exposureTime = in.samplingTime;
            end
            obj.imgLens = lens;
            % Frame listener
            obj.frameListener = addlistener(obj,'frameBuffer','PostSet',...
                @(src,evnt) obj.imagesc );
            obj.frameListener.Enabled = false;
        end
        
        function relay(obj,src)
            %% RELAY source propagation
            
            relay(obj.imgLens,src)
            %src.mask = tools.piston(300);
            if all([src.timeStamp]>=obj.startDelay)
                if obj.frameCount==0
                    disp(' @(detector:relay)> starting imager integration!')
                end
                obj.startDelay = -Inf;
                obj.frameCount = obj.frameCount + 1;
%                 test = size(src.mask);
%                 disp(test)
                obj.frameBuffer = obj.frameBuffer + cat(2,src.intensity);
                if src.timeStamp>=obj.exposureTime
                    src.timeStamp = 0;
                    flush(obj,src);
                end
            end
        end
        
        function flush(obj,src)
            fprintf(' @(detector:relay)> reading out and emptying buffer (%d frames)!\n',obj.frameCount)
            frameCount = obj.frameCount;
            %readOut(obj,obj.frameBuffer)
            obj.frame = obj.frameBuffer;
            obj.frameBuffer = 0;%*obj.frameBuffer;
            maskSize = obj.resolution;
            if ~isempty(obj.referenceFrame)
                obj.imgLens.fieldStopSize = obj.imgLens.fieldStopSize*2;
                src_ = source.*(obj.referenceFrame)*obj.imgLens;
                src_.mask = tools.piston(maskSize(1)*2);
                otf =  src_.amplitude;               
                src_ = source.*(obj.frame/frameCount)*obj.imgLens;
                src_.mask = tools.piston(maskSize(1)*2);
                obj.imgLens.fieldStopSize = obj.imgLens.fieldStopSize/2;
                otfAO =  src_.amplitude;
                
                %                         figure, imagesc(real(otfAO)/max(otfAO(:)))
                % strehl ratio
                obj.strehl = sum(otfAO(:))/sum(otf(:));
                % entrapped energy
                obj.tel
                a      = (obj.eeWidth/(src.wavelength/obj.tel.D*constants.radian2arcsec))/obj.tel.D;
                nOtf   = length(otfAO);
                u      = linspace(-1,1,nOtf).*obj.tel.D;
                [x,y]  = meshgrid(u);
                eeFilter ...
                    = a^2*(sin(pi.*x.*a)./(pi.*x.*a)).*...
                    (sin(pi.*y.*a)./(pi.*y.*a));
                otfAO = otfAO/max(otfAO(:));
                obj.ee = real(trapz(u,trapz(u,otfAO.*eeFilter)));
            end
            obj.frameCount = 0;
        end
        
        function imagesc(obj,varargin)
            %% IMAGESC Display the detector frame
            %
            % imagesc(obj) displays the frame of the detector object
            %
            % imagesc(obj,'PropertyName',PropertyValue) displays the frame of
            % the detector object and set the properties of the graphics object
            % imagesc
            %
            % h = imagesc(obj,...) returns the graphics handle
            %
            % See also: imagesc
            
            disp('Got it')
            if ishandle(obj.frameHandle)
                set(obj.frameHandle,'Cdata',obj.frameBuffer,varargin{:});
                %                 xAxisLim = [0,size(obj.frame,2)]+0.5;
                %                 yAxisLim = [0,size(obj.frame,1)]+0.5;
                %                 set( get(obj.frameHandle,'parent') , ...
                %                     'xlim',xAxisLim,'ylim',yAxisLim);
            else
                obj.frameHandle = image(obj.frameBuffer,...
                    'CDataMApping','Scaled',...
                    varargin{:});
                colormap(pink)
                axis xy equal tight
                colorbar('location','SouthOutside')
            end
        end

   
    end

end