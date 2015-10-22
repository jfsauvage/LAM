classdef utilities
    % UTILITIES Collection of various useful functions
    
    methods (Static)
        
        function out = piston(Npx,varargin)
            %% PISTON piston mode
            %
            % out = piston(Npx) Computes a piston on Npx pixel across the
            % diameter
            %
            % out = piston(Npx,nOut) Computes a piston on Npx pixel across
            % the diameter inside a square array of nOutXnOut pixels.
            %
            % out = piston(Npx,nOut,xOffset,yOffset) Computes a piston on
            % Npx pixel across the diameter inside a square array of
            % nOutXnOut pixels at xOffset and yOffset pixels from the
            % center.
            %
            % out = piston( ... ,'shape','square') By default the piston is
            % a disc but here it is forced to be a square
            %
            % out = piston( ... ,'shape','hex') By default the piston is
            % a disc but here it is forced to be hexagonal, nOut is equal
            % to twice the hexagonal side
            %
            % out = piston( ... ,'type','logical') By default the piston
            % values are in double but they can be casted into any types
            % supported by Matlab like logical
            
            p = inputParser;
            p.addRequired('Npx',@isnumeric);
            p.addOptional('nOut',Npx,@isnumeric);
            p.addOptional('xOffset',0,@isnumeric);
            p.addOptional('yOffset',0,@isnumeric);
            p.addParamValue('shape','disc',@ischar);
            p.addParamValue('type','double',@ischar);
            parse(p,Npx,varargin{:});
            param = p.Results;
            
            x = -(param.nOut-1)/2:(param.nOut-1)/2;
            u = x - param.xOffset;
            v = x - param.yOffset;
            
            [x,y,r,o] = utilities.cartAndPol(2.*u./Npx,2.*v./Npx);
            
            switch param.shape
                case 'disc'
                    out = double(r <= 1);
                case 'square'
                    out = double( abs(x)<=1 & abs(y)<=1 );
                case {'hex','hexagonal'}
                    out = double( abs(x)<=sqrt(3)/2 & abs(y)<=x/sqrt(3)+1 & abs(y)<=-x/sqrt(3)+1 );
                otherwise
                    error('The piston shape is either a disc or a square')
            end
            
            switch param.type
                case 'logical'
                    out = logical(out);
                case 'double'
                otherwise
                    error('The piston type is either a double or logical')
            end
        end
        
        function varargout = cartAndPol(u,varargin)
            %% CARTANDPOL Cartesian and polar coordinate arrays
            %
            % [x,y,r,o] = cartAndPol(n) Shortcuts to u = ((1-n):2:(n-1))/n;
            % [x,y] = meshgrid(u);[o,r] = cart2pol(x,y);
            %
            % [x,y,r,o] = cartAndPol(n,R) Shortcuts to u = R*((1-n):2:(n-1))/n;
            % [x,y] = meshgrid(u);[o,r] = cart2pol(x,y);
            %
            % [x,y,r,o] = cartAndPol(u) Shortcuts to [x,y] = meshgrid(u);[o,r] =
            % cart2pol(x,y);
            %
            % [x,y,r,o] = cartAndPol(u,v) Shortcuts to [x,y] = meshgrid(u,v);[o,r] =
            % cart2pol(x,y);
            %
            % [x,y,r,o] = cartAndPol(n,[],type)
            %
            % [x,y,r,o] = cartAndPol(n,R,type)
            %
            % [x,y,r,o] = cartAndPol(u,[],type
            %
            % [x,y,r,o] = cartAndPol(u,v, type) Same as above but now the type of x, y,
            % r and o is now specified: double (default) or single
            %
            % [...] = cartAndPol(...,'offset',[xOffset,yOffset]) offsets
            % the grid by xOffset and yOffset
            %
            % [r,o] = cartAndPol(...,'output','polar')
            %
            % [r] = cartAndPol(...,'output','radius')
            
            
            p = inputParser;
            p.addRequired('u',@isnumeric);
            p.addOptional('v',[],@isnumeric);
            p.addParamValue('offset',[0,0],@isnumeric);
            p.addParamValue('type','double',@ischar);
            p.addParamValue('output','all',@ischar);
            p.parse(u, varargin{:})
            u      = p.Results.u;
            v      = p.Results.v;
            offset = p.Results.offset;
            type   = p.Results.type;
            output = p.Results.output;
            
            if isempty(v)
                if numel(u)==1
                    u = 2*( -(u-1)/2:(u-1)/2 )/u;%linspace(-1,1,u);
                end
                v=u;
            elseif (numel(u)==1) && (numel(v)==1)
                u = 2*v*( -(u-1)/2:(u-1)/2 )/u;%linspace(-v,v,u);
                v = u;
            end
            
            if strcmp(type,'single')
                [x,y] = meshgrid(single(u-offset(1)),single(v-offset(2)));
            else
                [x,y] = meshgrid(u-offset(1),v-offset(2));
            end
            [o,r] = cart2pol(x,y);
            
            switch output
                case 'all'
                    varargout{1} = x;
                    varargout{2} = y;
                    varargout{3} = r;
                    varargout{4} = o;
                case 'polar'
                    varargout{1} = r;
                    varargout{2} = o;
                case 'radius'
                    varargout{1} = r;
                otherwise
                    error('oomao:utilities:cartAndPol:wrongOutput',...
                        'Valid outputs are all, polar or radius.')
            end
            
        end
        
        function frame = toggleFrame(frame,toggle)
            %% TOGGLEFRAME 2D to 3D array reshaping
            %
            % out = toggleFrame(frame) reshapes a 2D array into a 3D array or a
            % 3d array into a 2D array
            %
            % out = toggleFrame(frame,toggle) reshapes the array into a 2D
            % array if toggle is equal to 2 or into a 3D array if toggle is
            % equal to 3
            
            n    = ndims(frame);
            dims = size(frame);
            if length(dims)==2
                dims(3) = 1;
            end
            
            if nargin<2
                if n==2
                    toggle = 3;
                else
                    toggle = 2;
                end
            end
            
            
            
            if n~=toggle || toggle==2
                switch toggle
                    case 2
                        %                         fprintf(' @(toggleFrame)> 2D: [%d,%d] !\n',dims(1)*dims(2),dims(3))
                        frame = reshape(frame,dims(1)*dims(2),dims(3));
                    case 3
                        m = sqrt(dims(1));
                        %                         fprintf(' @(toggleFrame)> 3D: [%d,%d,%d] !\n',m,m,dims(2))
                        frame = reshape(frame,[m,m,dims(2)]);
                end
            end
            
        end
        
        
        function index = rearrange(sizeArray,sizeSubArray,overlap,columnMajor)
            %5 REARRANGE Array linear index scrambling
            %
            % index = rearrange(sizeArray,sizeSubArray) Rearrange the linear index of
            % an array of size sizeArray in a 2D matrix where each row contains the
            % index of a sub-array of size sizeSubArray taken from the initial array
            %
            % index = rearrange(sizeArray,sizeSubArray,overlap) Rearrange the linear
            % index of  an array of size sizeArray in a 2D matrix where each row
            % contains the index of a sub-array of size sizeSubArray taken from the
            % initial array. The sub-arrays overlap a number overlap(1) of rows and
            % overlap(2) of columns.
            %
            % index = rearrange(sizeArray,sizeSubArray,[],'column') Same as above but
            % now the sub-array browse the array along the columns not the rows as
            % before.
            %
            % index = rearrange(sizeArray,sizeSubArray,overlap,'column') Same as above
            % with overlapping
            %
            % Example:
            % >> a = reshape(1:36,6,6)
            % a =
            %      1     7    13    19    25    31
            %      2     8    14    20    26    32
            %      3     9    15    21    27    33
            %      4    10    16    22    28    34
            %      5    11    17    23    29    35
            %      6    12    18    24    30    36
            % >> index = rearrange(size(a),[3,3])
            % ans =
            %      1     4    19    22
            %      2     5    20    23
            %      3     6    21    24
            %      7    10    25    28
            %      8    11    26    29
            %      9    12    27    30
            %     13    16    31    34
            %     14    17    32    35
            %     15    18    33    36
            % >> rearrange(size(a),[3,3],[],'column')
            % ans =
            %      1    19     4    22
            %      2    20     5    23
            %      3    21     6    24
            %      7    25    10    28
            %      8    26    11    29
            %      9    27    12    30
            %     13    31    16    34
            %     14    32    17    35
            %     15    33    18    36
            % >> reshape( a(index) ,[3,3,4])
            % ans(:,:,1) =
            %      1     7    13
            %      2     8    14
            %      3     9    15
            % ans(:,:,2) =
            %      4    10    16
            %      5    11    17
            %      6    12    18
            % ans(:,:,3) =
            %
            %     19    25    31
            %     20    26    32
            %     21    27    33
            % ans(:,:,4) =
            %     22    28    34
            %     23    29    35
            %     24    30    36
            
            % $Id: rearrange.m 409 2006-07-12 16:49:24Z aoteam $
            
            if nargin<3 || isempty(overlap)
                overlap = zeros(1,2);
            end
            
            n     = sizeArray(1);
            m     = sizeArray(2);
            k     = prod(sizeArray(3:end));
            nSub  = sizeSubArray(1);
            if numel(sizeSubArray)==1
                mSub = nSub;
            else
                mSub  = sizeSubArray(2);
            end
            
            if rem(n,2)
                % Odd n
                nNSub = (n+overlap(1))/nSub;
                mMSub = (m+overlap(2))/mSub;
            else
                % Even n
                nNSub = n/nSub + overlap(1);
                mMSub = m/mSub + overlap(2);
            end
            
            % Type the index array as an unsigned integer with the coding depending on
            % the value of the largest elements of the index array
            switch find(2.^(2.^(3:6))-1 > prod(sizeArray),1)
                case 1
                    uint = @(x) uint8(x);
                case 2
                    uint = @(x) uint16(x);
                case 3
                    uint = @(x) uint32(x);
                case 4
                    uint = @(x) uint64(x);
                otherwise
                    error('Array size to big')
            end
            
            % Sub-array index in array
            [i,j] = ndgrid(uint(1:nSub),uint(1:mSub));
            index = repmat( uint(sub2ind( [n,m] , i(:) , j(:) )) , [ 1 , nNSub*mMSub*k ] );
            
            % Step index
            indexStep = ...
                repmat( uint(0:nNSub-1).'*(nSub-overlap(1))   , [  1   , mMSub*k ] ) + ...
                repmat( uint(0:mMSub*k-1)*(mSub-overlap(2))*n , [nNSub ,    1    ] );
            
            if nargin==4
                % Column major propagation of sub-array
                indexStep = indexStep.';
            end
            
            indexStep = repmat( reshape( indexStep , [1,nNSub*mMSub*k] ) , [nSub*mSub,1] );
            
            index = index + indexStep;
        end
        
        function out = sombrero(n,x)
            %% SOMBRERO Order n sombrero function
            %
            % out = sombrero(n,x) computes besselj(n,x)/x
            
            if n==0
                out = besselj(0,x)./x;
            else
                if n>1
                    out = zeros(size(x));
                else
                    out = 0.5*ones(size(x));
                end
                u = x~=0;
                x = x(u);
                out(u) = besselj(n,x)./x;
            end
        end
        
        function out = sinc(x)
            %% SINC Sinus cardinal function
            %
            % out = sinc(x) computes sin(pi*x)/(pi*x)
            
            out = ones(size(x));
            u = x~=0;
            x = x(u);
            out(u) = sin(pi*x)./(pi*x);
        end
        
        function out = fittingError(tel,atm,dm)
            %% FITTINGERROR Deformable mirror fitting error variance
            %
            % out = fittingError(telAtm,dm) computes the fitting error
            % variance of a a deformableMirror object for given telescope
            % and atmosphere objects
            
            c = (3/5)*(gamma(11/6)^2/pi^(8/3))*(24*gamma(6/5)/5)^(5/6);
            out = c*(tel.D/atm.r0)^(5/3)*...
                (dm.nValidActuator/pi + (tel.D/atm.L0)^2)^(-5/6);
        end
        
        function out = binning(frame,outRes)
            %% BINNING Frame binning
            %
            % out = binning(frame,[n,m]) bins the frame pixels into a nXm
            % array; frame can be either a single frame or a data cube
            
            [n,m,nFrame] = size(frame);
            out          = zeros(outRes(1),outRes(2),nFrame);
            n1 = n/outRes(1);
            m2 = m/outRes(2);
            if n1==1 && m2==1
                out = frame;
                return
            end
            if n1==1
                for kFrame=1:nFrame
                    out(:,:,kFrame) = ...
                        reshape( ...
                        sum( ...
                        reshape( ...
                        frame(:,:,kFrame).', m2 , [] ) ...
                        ).' , ...
                        outRes(2) , [] ).';
                end
            elseif m2==1
                for kFrame=1:nFrame
                    out(:,:,kFrame) = ...
                        reshape( ...
                        sum( ...
                        reshape( frame(:,:,kFrame) , n1 , [] ) ...
                        ) , ...
                        outRes(1) , [] );
                end
            else
                for kFrame=1:nFrame
                    out(:,:,kFrame) = ...
                        reshape( ...
                        sum( ...
                        reshape( ...
                        reshape( ...
                        sum( ...
                        reshape( frame(:,:,kFrame) , n1 , [] ) ...
                        ) , ...
                        outRes(1) , [] ).' , ...
                        m2 , [] ) ...
                        ) , ...
                        outRes(2) , [] ).';
                end
            end
        end
        
        function out = polar3(theta,rho,z,varargin)
            %% POLAR3 Polar coordinate plot with color coded markers
            %
            % polar3(theta,rho,z) makes a plot using polar coordinates of
            % the angle THETA, in radians, versus the radius RHO. The color
            % of the markers is scaled according to the values in vector z.
            %
            % polar3(theta,rho,z,style) uses the marker specified in style
            %
            % polar3(...,'zMinMax',zBound) sets the z color scale limits to
            % the zBound values
            %
            % h = polar3(...) returns a handle to the plotted object in H.
            %
            % See also polar
            
            p = inputParser;
            p.addRequired('theta',@isnumeric);
            p.addRequired('rho',@isnumeric);
            p.addRequired('z',@isnumeric);
            p.addOptional('style','.',@ischar);
            p.addParamValue('zMinMax',[],@isnumeric);
            p.parse(theta,rho,z , varargin{:});
            style   = p.Results.style;
            zMinMax = p.Results.zMinMax;
            
            n = length(theta);
            if isempty(zMinMax)
                minZ = min(z);
                maxZ = max(z);
                fprintf(' @(utilities:polar3)> Z axis minmax: [%.2f,%.2f]\n',minZ,maxZ)
            else
                minZ = zMinMax(1);
                maxZ = zMinMax(2);
            end
            
            c    = colormap;
            nc = length(c);
            zc = fix((nc-1)*(z - minZ)/(maxZ-minZ) + 1);
            
            index = find(rho==max(rho));
            h = polar(theta(index),rho(index),'.');
            delete(h)
            
            h = zeros(n,1);
            hold on
            for k=1:n
                h(k) = polar(theta(k),rho(k),style);
                set(h(k),'zData',z(k),'color',c(zc(k),:))
            end
            hold off
            
            hc = colorbar;
            set(hc,'ylim',[minZ maxZ])
            set(get(hc,'children'),'YData',[minZ maxZ])
            
            if nargout == 1
                out = h;
            end
        end
        
        
        function out = defocusDistance(a4,focalLength,diameter,wavelength,unit)
            % DEFOCUSDISTANCE Focal point deplacement for a Zernike defocus
            %
            % out = defocusDistance(a4,focalLength,diameter,wavelength)
            % Compute the focal point relative position [meter] for the
            % Zernike (Noll normalized) focus coefficients [radian], the
            % focalLength [meter], the beam diameter [meter] and the
            % wavelength [meter]
            %
            % out = defocusDistance(a4,focalLength,diameter,wavelength,unit)
            % The result is converted into the appropriate unit: 3, 0, -3,
            % -6, -9 for example correspond to km, m, mm, micron, nm,
            % respectively
            
            out = 16*sqrt(3)*a4*(focalLength/diameter)^2/...
                ( 2*pi/wavelength - 16*sqrt(3)*a4*focalLength/diameter^2 );
            
            if nargin>4
                out = out*10^-unit;
            end
        end
        
        function out = outOfFocus(delta,focalLength,diameter,wavelength,unit)
            % OUTOFFOCUS Zernike focus for a focal point deplacement
            %
            % out = outOfFocus(delta,focalLength,diameter,wavelength)
            % Compute the Zernike (Noll normalized) focus coefficients
            % [radian] for the focal point relative position [meter], the
            % focalLength [meter], the beam diameter [meter] and the
            % wavelength [meter]
            
            out = ( 2*pi*delta/wavelength ) / ...
                ( 16*sqrt(3)*( (focalLength/diameter)^2 + focalLength*delta/diameter^2 ) );
            
            if nargin>4
                out = (wavelength/(2*pi))*out*10^-unit;
            end
            
        end
        
        function out = orbitalVelocity(h,zen)
            %% ORBITALVELOCITY Orbital angular velocity
            %
            % out = orbitalVelocity(h) computes the orbital angular in
            % [rad/s] velocity at altitude h a zenith
            %
            % out = orbitalVelocity(h,zen) computes the orbital angular in
            % [rad/s] velocity at altitude h a zenith angle zen
            
            
            if nargin==1
                zen = 0;
            end
            out = sqrt(constants.G*constants.Me/(constants.Re+h)).*...
                (1-constants.Re*sin(zen)^2/(constants.Re+h))./h;
        end
        
        function out = pointAheadAngle(h,zen)
            %% POINTAHEADANGLE Point ahead angle
            %
            % out = pointAheadAngle(h) computes the orbital angular in
            % [rad] velocity at altitude h a zenith
            %
            % out = pointAheadAngle(h,zen) computes the orbital angular in
            % [rad] velocity at altitude h a zenith angle zen
            
            
            if nargin==1
                zen = 0;
            end
            out = 2*h*utilities.orbitalVelocity(h,zen)*sec(zen)/constants.c;
        end
        
        function [vertex,center] = hexagonalArray(nCycle,pitch)
            %% HEXAGONALARRAY Array of hexagonals
            %
            % [vertex,center] = hexagonalArray(nSegment,pitch) computes the
            % vertex and center coordinates of nSegment hexagonals with a
            % the given pitch arranged in a hexagonal array
            
            if nargin<2
                pitch=1;
            end
            a = pitch/sqrt(3);
            hexCoord = a*exp(1i*((0:5)*pi/3 + pi/2));
            count = 1;
            nSegment = 3*nCycle^2+3*nCycle+1;
            vertex = zeros(6,nSegment);
            vertex(:,count) = hexCoord;
            center = zeros(nSegment,1);
            for cycle=1:nCycle
                for o=1:6
                    zo = hexCoord + cycle*a*sqrt(3)*exp(1i*(o-1)*pi/3);
                    for k=1:cycle
                        zk = zo + (k-1)*a*sqrt(3)*exp(1i*((o-1)*pi/3+2*pi/3));
                        zk_center = mean(zk);
                        count = count + 1;
                        vertex(:,count) = zk;
                        center(count)   = zk_center;
                    end
                end
            end
            v = vertex(:);
            f = reshape(1:6*nSegment,6,nSegment);
            figure(nSegment)
            patch('Faces',f','Vertices',[real(v(:)),imag(v(:))],'FaceColor',[1,1,1]*0.8);
            line(real(center),imag(center),'color','r','marker','.')
            axis square
            set(gca,'ylim',get(gca,'xlim'))
            title(sprintf('%d segments',nSegment))
        end
        
        function V = gramSchmidt(V)
            %% GRAMSCHMIDT Gram-Schmidt orthonormalization process
            %
            % V = gramSchmidt(V) orthonormalize the vector set V according
            % to the Gram-Schimdt process
            
            k = size(V,2);
            for j=1:k
                v = V(:,j);
                for i=1:j-1
                    u = V(:,i);
                    v = v - u*(u'*v)*(u'*u);
                end
                V(:,j) = v/norm(v);
            end
            
        end
        function corr = correl_fcn_WienerKhinchin(one_sided_psd)
            sp_full   = [0 one_sided_psd(length(one_sided_psd):-1:2)  one_sided_psd];   % to keep symmetry. Otherwise fft(sp_full) presents imaginary components...
            corr      = fftshift(fft(fftshift(sp_full)));                               % autocorrelation thru Wiener-Khinchine theorem
            corr      = real(corr/max(corr));
        end
        
        
        function coeff = bestfitmodel(correl,t, Npts, ord, pmin, pmax, delta)
            %{
        ------------HEADER-----------------
        Objective         ::  Find the coefficients of a disturbance model that best fit the Npts first
                              steps of the (de)correlation curve correl

        Comments:             

        INPUT VARS
        correl             :: The correlation sequence with at least Npts length, normalised
        t                  :: The temporal vector of points were the correlation is evaluated
        Npts               :: Number of points to consider in the fitting
        ord                :: Order of the model to fit

        OUTPUT VARS
         coeff             :: Coefficients of the fitting model
        Created by         :: Carlos Correia
        Creation date      :: 06/05/2009
                      
        Change Record:     ::
        ------------HEADER END----------------
            %}
            % --- parse input ---
            if nargin < 5
                pmin = 0.4;
                pmax = 3;
                delta = 0.1;
            end
            coefftrial = pmin:delta:pmax;
            
            % --- build a meaningful vector of frequencies related to the temporal sampling of the correlation
            % funcion ---
            dt = t(2) - t(1);
            ffmax = 1/(2*dt);
            ffmin = 1/t(end);
            f = linspace(ffmin, ffmax, length(t)/2);
            %tab = ones(length(coefftrial))*inf;
            tab1 = ones(1,length(coefftrial))*inf;
            p2opt = zeros(length(coefftrial),1);
            if ord == 1
                %                 for c1 = 1:length(coefftrial)
                %                     p2      = coefftrial(c1);
                %                     syswind = tf(p2,[1 p2]);
                %                     [mag]   = bode(syswind,2*pi*f);
                %                     %[mag,phase,] = bode(syswind);
                %                     mag     = reshape(mag, 1,length(mag));
                %                     sp_full = [0 mag(length(mag):-1:1)  mag];
                %                     Cphi    = fftshift(fft(fftshift(sp_full.^2)));     % autocorrelation thru Wiener-Khinchine theorem
                %                     Cphi    = real(Cphi);
                %                     Cphi    = Cphi/max(Cphi);
                %                     Np      = length(Cphi);
                %                     tab(c1) = var(correl(1:Npts) - Cphi(floor(Np/2+1):floor(Np/2+1)+Npts-1));
                %                 end
                %[val,ind]   = min(tab);
                %coeff       = coefftrial(ind);
                [p2, fval] = fminbnd(@(p2) utilities.makeCphi_ord1(f,p2,correl,Npts), pmin, pmax);
                coeff = p2;
            elseif ord ==2
                parfor c1 = 1:length(coefftrial)
                    %                     for c2     = c1:length(coefftrial)
                    %                         p1     = coefftrial(c1);
                    %                         p2     = coefftrial(c2);
                    %                         syswind= tf(p1*p2,[1 p1+p2 p1*p2]);
                    %                         [mag]  = bode(syswind,2*pi*f);
                    %                         %[mag,phase,] = bode(syswind);
                    %                         mag    = reshape(mag, 1,length(mag));
                    %                         sp_full= [0 mag(length(mag):-1:1)  mag];
                    %                         Cphi   = fftshift(fft(fftshift(sp_full.^2)));     % autocorrelation thru Wiener-Khinchine theorem
                    %                         Cphi   = real(Cphi);
                    %                         Cphi   = Cphi/max(Cphi);
                    %                         Np     = length(Cphi);
                    %
                    %                         tab(c1,c2) = var(correl(1:Npts) - Cphi(floor(Np/2+1):floor(Np/2+1)+Npts-1));
                    %                     end
                    p1     = coefftrial(c1);
                    [p2,val] = fminbnd(@(p2) utilities.makeCphi_ord2(f,p1,p2,correl,Npts), pmin, pmax);
                    tab1(c1)  = val;
                    p2opt(c1) = p2;
                end
                %[val,ind] = min(tab(:));
                %[I,J]     = ind2sub([length(coefftrial) length(coefftrial)],ind);
                %coeff     = [coefftrial(I) coefftrial(J)];
                idx = find(tab1 == min(tab1(:)));
                coeff     = [coefftrial(idx) p2opt(idx)];
            elseif ord == 3
                for p3 = 1:length(coefftrial)
                parfor c1 = 1:length(coefftrial)
                    %                     for c2     = c1:length(coefftrial)
                    %                         p1     = coefftrial(c1);
                    %                         p2     = coefftrial(c2);
                    %                         syswind= tf(p1*p2,[1 p1+p2 p1*p2]);
                    %                         [mag]  = bode(syswind,2*pi*f);
                    %                         %[mag,phase,] = bode(syswind);
                    %                         mag    = reshape(mag, 1,length(mag));
                    %                         sp_full= [0 mag(length(mag):-1:1)  mag];
                    %                         Cphi   = fftshift(fft(fftshift(sp_full.^2)));     % autocorrelation thru Wiener-Khinchine theorem
                    %                         Cphi   = real(Cphi);
                    %                         Cphi   = Cphi/max(Cphi);
                    %                         Np     = length(Cphi);
                    %
                    %                         tab(c1,c2) = var(correl(1:Npts) - Cphi(floor(Np/2+1):floor(Np/2+1)+Npts-1));
                    %                     end
                    p1     = coefftrial(c1);
                    [p2,val] = fminbnd(@(p2) utilities.makeCphi_ord3(f,p1,p2,p3,correl,Npts), pmin, pmax);
                    tab1(c1,p3)  = val;
                    p2opt(c1,p3) = p2;
                end
                end
                %[val,ind] = min(tab(:));
                %[I,J]     = ind2sub([length(coefftrial) length(coefftrial)],ind);
                %coeff     = [coefftrial(I) coefftrial(J)];
                [idx idy] = find(tab1 == min(tab1(:)));
                coeff     = [coefftrial(idx) p2opt(idx) coefftrial(idy)];
            end
        end
        function out = makeCphi_ord1(f,p2,correl,Npts)
            syswind = tf(p2,[1 p2]);
            [mag]   = bode(syswind,2*pi*f);
            %[mag,phase,] = bode(syswind);
            mag     = reshape(mag, 1,length(mag));
            sp_full = [0 mag(length(mag):-1:2)  mag];
            Cphi    = fftshift(fft(fftshift(sp_full.^2)));     % autocorrelation thru Wiener-Khinchine theorem
            Cphi    = real(Cphi);
            Cphi    = Cphi/max(Cphi);
            Np      = length(Cphi);
            
            out = norm(correl(1:Npts) - Cphi(floor(Np/2+1):floor(Np/2+1)+Npts-1));
        end
        function out = makeCphi_ord2(f,p1,p2,correl,Npts)
            syswind= tf(p1*p2,[1 p1+p2 p1*p2]);
            [mag]  = bode(syswind,2*pi*f);
            %[mag,phase,] = bode(syswind);
            mag    = reshape(mag, 1,length(mag));
            sp_full= [0 mag(length(mag):-1:2)  mag];
            Cphi   = fftshift(fft(fftshift(sp_full.^2)));     % autocorrelation thru Wiener-Khinchine theorem
            Cphi   = real(Cphi);
            Cphi   = Cphi/max(Cphi);
            Np     = length(Cphi);
            
            out = norm(correl(1:Npts) - Cphi(floor(Np/2+1):floor(Np/2+1)+Npts-1));
        end
                    function out = makeCphi_ord3(f,p1,p2,p3,correl,Npts)
            syswind= zpk([],[-p1 -p2 -p3],p1*p2*p3);
            [mag]  = bode(syswind,2*pi*f);
            %[mag,phase,] = bode(syswind);
            mag    = reshape(mag, 1,length(mag));
            sp_full= [0 mag(length(mag):-1:2)  mag];
            Cphi   = fftshift(fft(fftshift(sp_full.^2)));     % autocorrelation thru Wiener-Khinchine theorem
            Cphi   = real(Cphi);
            Cphi   = Cphi/max(Cphi);
            Np     = length(Cphi);
            
            out = norm(correl(1:Npts) - Cphi(floor(Np/2+1):floor(Np/2+1)+Npts-1));
        end
    end
end