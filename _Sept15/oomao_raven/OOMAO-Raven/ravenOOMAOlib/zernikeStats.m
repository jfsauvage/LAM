classdef zernikeStats
    % Zernike statistics static class
    
    methods (Static)
         
        function out = spectrum(f,atm)
            %% SPECTRUM Phase power spectrum density
            %
            % out = phaseStats.spectrum(f,atm) computes the phase power
            % spectrum density from the spatial frequency f and an
            % atmosphere object
            %
            % See also atmosphere
            
            out = (24.*gamma(6./5)./5).^(5./6).*...
                (gamma(11./6).^2./(2.*pi.^(11./3))).*...
                atm.r0.^(-5./3);
            out = out.*(f.^2 + 1./atm.L0.^2).^(-11./6);
            out = sum([atm.layer.fractionnalR0]).*out;
        end
        
        function out = symSpectrum(symf)
            syms r0 L0
            out = (24*gamma(sym(6/5))/5)^(5./6)*...
                (gamma(sym(11/6))^2/(2*pi^(11/3)))*r0^(5/3);
            out = out*L0^(11/3)*( (symf*L0)^2 + 1 )^(-11/6);
        end
        
                
        function out = variance(zern,atm)
            %% VARIANCE Zernike coefficients variance
            %
            % out = zernikeStats.variance(modes,atmosphere) computes the
            % variance of Zernike coefficients from the modes and the
            % atmosphere object
            %
            % out = zernikeStats.variance(zernike,atmosphere) computes the
            % variance of Zernike coefficients from the Zernike polynomials
            % object and the atmosphere object
            %
            % Example:
            % atm = atmosphere(photometry.V,0.15,30);
            % tel = telescope(10);
            % modes = 1:15;
            % figure
            % semilogy(modes,zernikeStats.variance(modes,atm),'--.')
            % xlabel('Zernike modes')
            % ylabel('Variance [rd^2]')
            %
            % See also zernike, atmosphere
            
            
            if ~isa(zern,'zernike')
                zern = zernike(zern);
            end
            r0 = atm.r0;
            L0 = atm.L0;
            D  = zern.D;
            jv = zern.j;
            nv = zern.n;
            mv = zern.m;
            nv0   = nv;
            index = diff(nv)~=0;
            jv = [jv(index) jv(end)];
            mv = [mv(index) mv(end)];
            nv = [nv(index) nv(end)];
            nf  = length(nv);
            out = zeros(length(zern.j),1);
            
            for cpt = 1:nf
                
                j = jv(cpt);
                n = nv(cpt);
                m = mv(cpt);
                
                out(nv0==n,1) = zernCovCoef(r0,L0,D,j,j,n,m,n,m);
                
            end
            function out = zernCovCoef(r0,L0,D,i,j,ni,mi,nj,mj)
                if (mi==mj) && (rem(abs(i-j),2)==0 || ((mi==0) && (mj==0)))
                    if L0==Inf
                        if i==1 && j==1
                            out = Inf;
                        else
                            out = (gamma(11./6).^2.*gamma(14./3)./(2.^(8./3).*pi)).*(24.*gamma(6./5)./5).^(5./6).*...
                                (D./r0).^(5./3).*sqrt((ni+1).*(nj+1)).*(-1).^((ni+nj-mi-mj)./2).*...
                                newGamma(-5./6+(ni+nj)./2,...
                                [23./6+(ni+nj)./2 17./6+(ni-nj)./2 17./6+(nj-ni)./2]);
                        end
                    else
                        out = (4.*gamma(11./6).^2./pi.^(14./3)).*(24.*gamma(6./5)./5).^(5./6).*...
                            (L0./r0).^(5./3).*(L0./D).^2.*...
                            sqrt((ni+1).*(nj+1)).*(-1).^((ni+nj-mi-mj)./2).*...
                            UnParamEx4q2(0,ni+1,nj+1,11./6,pi.*D./L0);
                    end
                else
                    out = 0;
                end
                function out = newGamma(a,b)
                    % NEWGAMMA Computes the function defined by Eq.(1.18) in R.J. Sasiela's book :
                    % Electromagnetic Wave Propagation in Turbulence, Springer-Verlag.
                    % out = newGamma(a,b)
                    
                    out = prod(gamma(a))./prod(gamma(b));
                end
                function out = UnParamEx4q2(mu,alpha,beta,p,a)
                    % UNPARAMEX4Q2 Computes the integral given by the Eq.(2.33) of the thesis
                    % of R. Conan (Modelisation des effets de l'echelle externe de coherence
                    % spatiale du front d'onde pour l'Observation a Haute Resolution Angulaire
                    % en Astronomie, University of Nice-Sophia Antipolis, October 2000)
                    % http://www-astro.unice.fr/GSM/Bibliography.html#thesis
                    
                    a1 = [(alpha+beta+1)./2 (2+mu+alpha+beta)./2 (mu+alpha+beta)./2];
                    b1 = [1+alpha+beta 1+alpha 1+beta];
                    a2 = [(1-mu)./2+p 1+p p];
                    b2 = [1+(alpha+beta-mu)./2+p 1+(alpha-beta-mu)./2+p 1+(beta-alpha-mu)./2+p];
                    
                    out = (1./(2.*sqrt(pi).*gamma(p))).*(...
                        newGamma([a1 p-(mu+alpha+beta)./2],b1).*a.^(mu+alpha+beta).*...
                        pochammerSeries(3,5,a1,[1-p+(mu+alpha+beta)./2 b1 1],a.^2) + ...
                        newGamma([(mu+alpha+beta)./2-p a2],b2).*a.^(2.*p).*...
                        pochammerSeries(3,5,a2,[1-(mu+alpha+beta)./2+p b2 1],a.^2));
                    function out = pochammerSeries(p,q,a,b,z,tol,nmax)
                        % POCHAMMERSERIES Computes power series in Pochammer notation
                        % pochammerSeries(p,q,a,b,z)
                        % pochammerSeries(p,q,a,b,z,tol)
                        % pochammerSeries(p,q,a,b,z,[],nmax)
                        % pochammerSeries(p,q,a,b,z,tol,nmax)
                        
                        if (p==(q+1) && abs(z)<1) || (abs(z)==1 && real(sum(a)-sum(b))<0) || p<(q+1)
                            
                            if p==length(a) && q==length(b)
                                
                                switch nargin
                                    case 6
                                        nmax = 1e3;
                                    case 7
                                        if isempty(tol)
                                            tol = 1e-6;
                                        end
                                    otherwise
                                        tol = 1e-6;
                                        nmax = 1e3;
                                end
                                
                                out = zeros(size(z));
                                
                                indz = find(z==0);
                                if ~isempty(indz)
                                    out(indz) = 1;
                                end
                                
                                indnz = find(z~=0);
                                if ~isempty(indnz)
                                    z = z(indnz);
                                    ck = 1;
                                    step = Inf;
                                    k = 0;
                                    som = ck;
                                    while (k<=nmax) && (step>tol)
                                        ckp1 = prod(a+k).*z.*ck./prod(b+k);
                                        step = abs(abs(ck)-abs(ckp1));
                                        som = som + ckp1;
                                        k = k+1;
                                        ck = ckp1;
                                    end
                                    if step>tol
                                        warning('pochammerSeries','Maximum iteration reached before convergence')
                                    end
                                    out(indnz) = som;
                                end
                                
                            else
                                error('p and q must be the same length than vectors a and b, respectively')
                                
                            end
                            
                        else
                            error('This generalized hypergeometric function doesn''t converge')
                        end
                    end
                end
                
            end
        end
        
        function out = covariance(zern,atm)
            %% COVARIANCE Zernike coefficients covariance
            %
            % out = zernikeStats.covariance(modes,atmosphere)
            % computes the covariance matrix of Zernike coefficients from
            % the modes and the atmosphere object
            %
            % out = zernikeStats.covariance(zernike,atmosphere) computes the
            % covariance matrix of Zernike coefficients from the Zernike
            % polynomials object and the atmosphere object
            %
            % Example:
            % atm = atmosphere(photometry.V,0.15,30);
            % tel = telescope(10);
            % modes = 1:15;
            % figure
            % spy(zernikeStats.covariance(modes,atm))
            %
            % See also zernike, atmosphere
            
            
            if ~isa(zern,'zernike')
                zern = zernike(zern);
            end
            r0 = atm.r0;
            L0 = atm.L0;
            D  = zern.D;
            % -- covariance --
            n = zern.n;
            m = zern.m;
            [i,j]   = meshgrid(zern.j);
            [ni,nj] = meshgrid(n);
            [mi,mj] = meshgrid(m);
            index   = (mi==mj) & (rem(abs(i-j),2)==0 | ((mi==0) & (mj==0)));
            index   = logical(index - diag(diag(index)));
            i       = i(index);
            j       = j(index);
            ni      = ni(index);
            nj      = nj(index);
            mi      = mi(index);
            mj      = mj(index);
            
            nf       = sum(index(:));
            covValue = zeros(1,nf);
            for cpt = 1:nf
                covValue(cpt) = zernCovCoef(r0,L0,D,i(cpt),j(cpt),ni(cpt),mi(cpt),nj(cpt),mj(cpt));
            end
            out          = zeros(length(zern.j));
            out(index)   = covValue;
            % -- variance --
            jv = zern.j;
            nv = zern.n;
            mv = zern.m;
            nv0   = nv;
            index = diff(nv)~=0;
            jv = [jv(index) jv(end)];
            mv = [mv(index) mv(end)];
            nv = [nv(index) nv(end)];
            nf  = length(nv);
            variance = zeros(length(zern.j),1);
            
            for cpt = 1:nf
                
                j = jv(cpt);
                n = nv(cpt);
                m = mv(cpt);
                
                variance(nv0==n,1) = zernCovCoef(r0,L0,D,j,j,n,m,n,m);
                
            end
            % -- covariance + variance --
            out = out + diag(variance);
            function out = zernCovCoef(r0,L0,D,i,j,ni,mi,nj,mj)
                if (mi==mj) && (rem(abs(i-j),2)==0 || ((mi==0) && (mj==0)))
                    if L0==Inf
                        if i==1 && j==1
                            out = Inf;
                        else
                            out = (gamma(11./6).^2.*gamma(14./3)./(2.^(8./3).*pi)).*(24.*gamma(6./5)./5).^(5./6).*...
                                (D./r0).^(5./3).*sqrt((ni+1).*(nj+1)).*(-1).^((ni+nj-mi-mj)./2).*...
                                newGamma(-5./6+(ni+nj)./2,...
                                [23./6+(ni+nj)./2 17./6+(ni-nj)./2 17./6+(nj-ni)./2]);
                        end
                    else
                        out = (4.*gamma(11./6).^2./pi.^(14./3)).*(24.*gamma(6./5)./5).^(5./6).*...
                            (L0./r0).^(5./3).*(L0./D).^2.*...
                            sqrt((ni+1).*(nj+1)).*(-1).^((ni+nj-mi-mj)./2).*...
                            UnParamEx4q2(0,ni+1,nj+1,11./6,pi.*D./L0);
                    end
                else
                    out = 0;
                end
                function out = newGamma(a,b)
                    % NEWGAMMA Computes the function defined by Eq.(1.18) in R.J. Sasiela's book :
                    % Electromagnetic Wave Propagation in Turbulence, Springer-Verlag.
                    % out = newGamma(a,b)
                    
                    out = prod(gamma(a))./prod(gamma(b));
                end
                function out = UnParamEx4q2(mu,alpha,beta,p,a)
                    % UNPARAMEX4Q2 Computes the integral given by the Eq.(2.33) of the thesis
                    % of R. Conan (Modelisation des effets de l'echelle externe de coherence
                    % spatiale du front d'onde pour l'Observation a Haute Resolution Angulaire
                    % en Astronomie, University of Nice-Sophia Antipolis, October 2000)
                    % http://www-astro.unice.fr/GSM/Bibliography.html#thesis
                    
                    a1 = [(alpha+beta+1)./2 (2+mu+alpha+beta)./2 (mu+alpha+beta)./2];
                    b1 = [1+alpha+beta 1+alpha 1+beta];
                    a2 = [(1-mu)./2+p 1+p p];
                    b2 = [1+(alpha+beta-mu)./2+p 1+(alpha-beta-mu)./2+p 1+(beta-alpha-mu)./2+p];
                    
                    out = (1./(2.*sqrt(pi).*gamma(p))).*(...
                        newGamma([a1 p-(mu+alpha+beta)./2],b1).*a.^(mu+alpha+beta).*...
                        pochammerSeries(3,5,a1,[1-p+(mu+alpha+beta)./2 b1 1],a.^2) + ...
                        newGamma([(mu+alpha+beta)./2-p a2],b2).*a.^(2.*p).*...
                        pochammerSeries(3,5,a2,[1-(mu+alpha+beta)./2+p b2 1],a.^2));
                    function out = pochammerSeries(p,q,a,b,z,tol,nmax)
                        % POCHAMMERSERIES Computes power series in Pochammer notation
                        % pochammerSeries(p,q,a,b,z)
                        % pochammerSeries(p,q,a,b,z,tol)
                        % pochammerSeries(p,q,a,b,z,[],nmax)
                        % pochammerSeries(p,q,a,b,z,tol,nmax)
                        
                        if (p==(q+1) && abs(z)<1) || (abs(z)==1 && real(sum(a)-sum(b))<0) || p<(q+1)
                            
                            if p==length(a) && q==length(b)
                                
                                switch nargin
                                    case 6
                                        nmax = 1e3;
                                    case 7
                                        if isempty(tol)
                                            tol = 1e-6;
                                        end
                                    otherwise
                                        tol = 1e-6;
                                        nmax = 1e3;
                                end
                                
                                out = zeros(size(z));
                                
                                indz = find(z==0);
                                if ~isempty(indz)
                                    out(indz) = 1;
                                end
                                
                                indnz = find(z~=0);
                                if ~isempty(indnz)
                                    z = z(indnz);
                                    ck = 1;
                                    step = Inf;
                                    k = 0;
                                    som = ck;
                                    while (k<=nmax) && (step>tol)
                                        ckp1 = prod(a+k).*z.*ck./prod(b+k);
                                        step = abs(abs(ck)-abs(ckp1));
                                        som = som + ckp1;
                                        k = k+1;
                                        ck = ckp1;
                                    end
                                    if step>tol
                                        warning('pochammerSeries','Maximum iteration reached before convergence')
                                    end
                                    out(indnz) = som;
                                end
                                
                            else
                                error('p and q must be the same length than vectors a and b, respectively')
                                
                            end
                            
                        else
                            error('This generalized hypergeometric function doesn''t converge')
                        end
                    end
                end
                
            end
        end
        
        function out = residualVariance(N,atm,tel)
            %% RESIDUALVARIANCE Zernike corrected wavefront variance 
            %
            % out = zernikeResidualVariance(N,atm,tel)
            %
            % See also: atmosphere, telescope
            
            zern = zernike(1:N,tel.D,'logging',false);
            r0 = atm.r0;
            L0 = atm.L0;
            D  = tel.D;
            aiVar = zernikeStats.variance(zern,atm);
            
            if isinf(L0)
                Delta1 = -(2.*gamma(11./6).^2./pi.^1.5).*(24.*gamma(6./5)./5).^(5./6).*...
                    (D./r0).^(5./3).*newGamma([-5./6,7./3],[23./6,17./6]);
                out = Delta1 - sum(aiVar(2:end));
            else
                a = phaseStats.variance(atm);
                b = sum(aiVar);
                out = a - b;
            end
        end
        
        function aiaj = angularCovariance(zern,atm,src,optSrc)
            %% ANGULARCOVARIANCE Zernike coefficients angular covariance
            %
            % aiaj = angularCovariance(zern,atm,src) computes
            % the covariance matrix between Zernike coefficients of Zernike
            % polynomials zern corresponding to wavefront propagating from
            % two sources src(1) and src(2) through the atmosphere atm
            %
            % See also zernike, atmosphere, source
            
            nGs = numel(src);
            if nGs>2 % then its a meta-matrix
%                 disp(' @(phaseStats.zernikeAngularCovariance)> META-MATRIX:')
                if nargin<4 % a correlation meta-matrix
%                     disp(' @(phaseStats.zernikeAngularCovariance)> AUTO CORRELATION META-MATRIX:')
                    iSrc = src;
                    jSrc = src;
                    mGs = nGs;
                    aiaj = cell(nGs,mGs);
%                     for iGs = 1:nGs
%                         fprintf(' @(Data covariance)> ');
%                         gsCurrent = iSrc(iGs);
%                         parfor jGs = iGs:mGs
%                             fprintf('gs#%d/gs#%d - ',iGs,jGs)
%                             aiaj{iGs,jGs} = phaseStats.zernikeAngularCovariance(zern,atm,[gsCurrent,jSrc(jGs)]);
%                         end
%                         fprintf('\b\b\b\n')
%                     end
                    nmGs  = [nGs mGs];
                    index = 1:nGs*mGs;
                    mask  = triu(true(nGs,mGs));
                    index = index(mask);
                    buffer = cell(1,length(index));
                    for kGs = 1:length(index)
                        ijGs = index(kGs);
                        [iGs,jGs] = ind2sub(nmGs,ijGs);
%                         fprintf(' @(phaseStats.zernikeAngularCovariance)> gs#%d/gs#%d \n',iGs,jGs);
                        buffer{kGs} = zernikeStats.angularCovariance(zern,atm,[iSrc(iGs),jSrc(jGs)]);
                    end
                    aiaj(mask) = buffer;
                    index = cellfun(@isempty,aiaj);
                    aiaj(index) = cellfun(@transpose,aiaj(triu(~index,1)),'UniformOutput',false);
                else % a cross-correlation meta-matrix
%                     disp('@(phaseStats.zernikeAngularCovariance)> CROSS CORRELATION META-MATRIX:')
                    iSrc = src;
                    jSrc = optSrc;
                    mGs = numel(jSrc);
                    aiaj = cell(nGs,mGs);
%                     for iGs = 1:nGs
%                         fprintf(' @(phaseStats.zernikeAngularCovariance)> ');
%                         gsCurrent = iSrc(iGs);
%                         for jGs = 1:mGs
%                             fprintf('gs#%d/gs#%d - ',iGs,jGs);
%                             aiaj{iGs,jGs} = phaseStats.zernikeAngularCovariance(zern,atm,[gsCurrent,jSrc(jGs)]);
%                         end
%                         fprintf('\b\b\b\n')
%                     end
                    nmGs  = [nGs mGs];
                    for kGs = 1:nGs*mGs
                        [iGs,jGs] = ind2sub(nmGs,kGs);
%                         fprintf(' @(phaseStats.zernikeAngularCovariance)> gs#%d/gs#%d \n',iGs,jGs);
                        aiaj{kGs} = zernikeStats.angularCovariance(zern,atm,[iSrc(iGs),jSrc(jGs)]);
                    end
                end
%                 aiaj = cell2mat(aiaj);
            else
                if src(1)==src(2)
                    aiaj = phaseStats.zernikeCovariance(zern,atm);
                else
                    R   = zern.R;
                    zs1 = src(1).height;
                    zs2 = src(2).height;
                    xSrc = tan(src(1).zenith).*cos(src(1).azimuth) - ...
                        tan(src(2).zenith).*cos(src(2).azimuth);
                    ySrc = tan(src(1).zenith).*sin(src(1).azimuth) - ...
                        tan(src(2).zenith).*sin(src(2).azimuth);
                    rhoSrcLayer = hypot(xSrc,ySrc);
                    thetaSrcLayer = atan2(ySrc,xSrc);
                    nMode = length(zern.j);
                    znmj = mat2cell([zern.j;zern.n;zern.m],3,ones(1,nMode));
                    znmj = repmat(znmj,nMode,1);
                    znmi = znmj';
                    psdCst = (24.*gamma(6./5)./5).^(5./6).*...
                        (gamma(11./6).^2./(2.*pi.^(11./3))).*...
                        atm.r0.^(-5./3);
                    if all( isinf( [zs1 zs2] ) ) % NGS CASE
                        a1l     = R;
                        a2l     = R;
                        denom   = pi.*a1l.*a2l;
                        sl      = [atm.layer.altitude]'.*rhoSrcLayer;
                        fr0     = [atm.layer.fractionnalR0]';
                        aiajFun = @ (znmi,znmj) ...
                            quadgk(@(x) integrand(x,znmi(1),znmi(2),znmi(3),znmj(1),znmj(2),znmj(3)), ...
                            0, Inf, 'AbsTol',1e-3, 'RelTol',1e-2);
%                         n = 201;
%                         r = linspace(0,20,n);
%                         r(1) = 1e-6;
%                         aiajFun = @ (znmi,znmj) ...
%                             trapz(r,integrandNgs(r,znmi(1),znmi(2),znmi(3),znmj(1),znmj(2),znmj(3)));                       
                    else % LGS CASE (TO DO: optimize for LGS as for NGS)
                        a1l   = zeros(1,atm.nLayer);
                        a2l   = zeros(1,atm.nLayer);
                        denom = zeros(1,atm.nLayer);
                        sl    = zeros(1,atm.nLayer);
                        for lLayer=1:atm.nLayer
                            a1l(lLayer) = R.*(1 - atm.layer(lLayer).altitude./zs1);
                            a2l(lLayer) = R.*(1 - atm.layer(lLayer).altitude./zs2);
                            denom(lLayer) = pi.*a1l(lLayer).*a2l(lLayer);
                            sl(lLayer) = atm.layer(lLayer).altitude.*rhoSrcLayer;
                        end
                        aiajFun = @ (znmi,znmj) ...
                            quadgk(@(x) integrand(x,znmi(1),znmi(2),znmi(3),znmj(1),znmj(2),znmj(3)), ...
                            0, Inf, 'AbsTol',1e-3, 'RelTol',1e-2);
%                         n = 201;
%                         r = linspace(0,20,n);
%                         r(1) = 1e-6;
%                         aiajFun = @ (znmi,znmj) ...
%                             trapz(r,integrandNgs(r,znmi(1),znmi(2),znmi(3),znmj(1),znmj(2),znmj(3)));                       
                    end
                    aiaj = zeros(nMode);
                    index = triu(true(nMode));
                    %                 tic
                    aiaj(index) = cellfun(aiajFun,znmj(index),znmi(index));
                    %                 toc
                    aiaj = aiaj + triu(aiaj,1)';
                    aiaj = bsxfun(@times,aiaj,(-1).^zern.m');
                    %                     aiaj = cellfun(aiajFun,znmj,znmi);
                end
            end
            function out = integrand(x,zi,ni,mi,zj,nj,mj)
                krkr_mi = mi==0;
                krkr_mj = mj==0;
                out = 0;
                factor1 = sqrt((ni+1)*(nj+1)).*...
                    (-1).^(0.5.*(ni+nj)).*...
                    2.^(1-0.5.*(krkr_mi+krkr_mj));%.*...
                %(-1).^mj;
                factor2 = (-1).^(1.5*(mi+mj)).*...
                    cos( ...
                    (mi+mj).*thetaSrcLayer + ...
                    pi.*( (1-krkr_mi).*((-1).^zi-1) + ...
                    (1-krkr_mj).*((-1).^zj-1) )./4 );
                factor3 = (-1).^(1.5*abs(mi-mj)).*...
                    cos( ...
                    (mi-mj).*thetaSrcLayer + ...
                    pi.*( (1-krkr_mi).*((-1).^zi-1) - ...
                    (1-krkr_mj).*((-1).^zj-1) )./4 );
                for kLayer=1:atm.nLayer
%                     a1l = R.*(1 - atm.layer(kLayer).altitude./zs1);
%                     a2l = R.*(1 - atm.layer(kLayer).altitude./zs2);
%                     denom = pi.*a1l.*a2l;
%                     sl = atm.layer(kLayer).altitude.*rhoSrcLayer;
                    red1 = a1l(kLayer).*x;
                    red2 = a2l(kLayer).*x;
                    red = sl(kLayer).*x;
%                     phasePSD = phaseStats.spectrum(0.5*x/pi,atmLayers{kLayer});
                    f = 0.5*x/pi;
                    phasePSD = atm.layer(kLayer).fractionnalR0.*...
                        psdCst.*(f.^2 + 1./atm.L0.^2).^(-11./6);
%                     phasePSD = phaseStats.spectrum(0.5*x/pi,atm.slab(lLayer));
%                     besselsRadialOrder = besselj(ni+1,red1).*besselj(nj+1,red2);
%                     tripleBessel1 = besselj(mi+mj,red);
%                     tripleBessel2 = besselj(abs(mi-mj),red);
                    besselsRadialOrder = besselmx('J',ni+1,red1).*besselmx('J',nj+1,red2);
                    tripleBessel1 = besselmx('J',mi+mj,red);
                    tripleBessel2 = besselmx('J',abs(mi-mj),red);
                    out = out +  (factor1./denom(kLayer)).*(phasePSD./x).*...
                        ( factor2.*tripleBessel1 + factor3.*tripleBessel2 ).*...
                        besselsRadialOrder;
                end
                out = real(out);
            end
            function out = integrandNgs(x,zi,ni,mi,zj,nj,mj)
                mipmj = mi+mj;
                mimmj = abs(mi-mj);
                krkr_mi = mi==0;
                krkr_mj = mj==0;
                factor1 = sqrt((ni+1)*(nj+1)).*...
                    (-1).^(0.5.*(ni+nj)).*...
                    2.^(1-0.5.*(krkr_mi+krkr_mj));%.*...
                %(-1).^mj;
                factor1 = factor1./denom;
                factor2 = (-1).^(1.5*(mi+mj)).*...
                    cos( ...
                    (mi+mj).*thetaSrcLayer + ...
                    pi.*( (1-krkr_mi).*((-1).^zi-1) + ...
                    (1-krkr_mj).*((-1).^zj-1) )./4 );
                factor2 = factor2*factor1;
                factor3 = (-1).^(1.5*abs(mi-mj)).*...
                    cos( ...
                    (mi-mj).*thetaSrcLayer + ...
                    pi.*( (1-krkr_mi).*((-1).^zi-1) - ...
                    (1-krkr_mj).*((-1).^zj-1) )./4 );
                factor3 = factor3*factor1;
                red1 = a1l.*x;
                red2 = a2l.*x;
                besselsRadialOrder = besselmx('J',ni+1,red1).*besselmx('J',nj+1,red2);
                f = 0.5*x/pi;
                phasePSD = ...
                    psdCst.*(f.^2 + 1./atm.L0.^2).^(-11./6)./x;
                red = sl*x;
                tripleBessel1 = besselmx('J',mipmj,red);
                tripleBessel2 = besselmx('J',mimmj,red);
                out = (fr0*phasePSD).*...
                    ( factor2.*tripleBessel1 + factor3.*tripleBessel2 );
                out = sum( bsxfun( @times , out , besselsRadialOrder ) );
                out = real(out);
            end
        end
        
        function out = aiaj(zi,ni,mi,zj,nj,mj,r0,L0,R)
            out = 0;
            if (mi==mj) && (rem(abs(zi-zj),2)==0 || ((mi==0) && (mj==0)))
            out = (gamma(11/6)^2/pi^(14/3))*(24*gamma(6/5)/5)^(5/6)*...
                (L0/r0)^(5/3)*(L0/R)^2*...
                sqrt((ni+1)*(nj+1)).*(-1).^((ni+nj-mi-mj)/2).*...
                aiajFun(ni,nj,2*pi*R/L0);
            end
        end
        
        function varargout = residueVarianceMap(N,r,o,atm,tel)
        %% RESIDUEVARIANCEMAP
            
            persistent aiaj zMode zN zM R L0 L0r0Ratio L0RRatio c1 c2 ...
                red1 atmVar
            if isempty(aiaj)
                fprintf(' @(zernikeStats.residueVarianceMap)> Initializing aiaj ...')
                
                R         = tel.R;
                L0        = atm.L0;
                L0r0Ratio = (L0./atm.r0).^(5./3);
                L0RRatio  = (L0./R).^(2);
                
                c1 = (24.*gamma(6./5)./5).^(5./6);
                c2 = (gamma(11/6)^2/pi^(14/3))*c1*L0r0Ratio*L0RRatio;
                
                red1 = 2*pi*R/L0;
                r    = r/R;
                atmVar    = c1.*...
                    (gamma(11./6).*gamma(5./6)./(2.*pi.^(8./3))).*L0r0Ratio;
                
                zMode = 1:N;
                zern = zernike(zMode,'logging',false);
                zN    = zern.n;%radialOrder(zMode)
                zM    = zern.m;%azimuthFrequency(zMode,zN)
                aiaj = zeros(N);
                for zi = zMode
                    ni = zN(zi);
                    mi = zM(zi);
                    for zj = 1:N
                        nj = zN(zj);
                        mj = zM(zj);
                        if (mi==mj) && (rem(abs(zi-zj),2)==0 || ((mi==0) && (mj==0)))
                            fnm = sqrt((ni+1)*(nj+1)).*(-1).^((ni+nj-mi-mj)/2);
                            aiaj(zi,zj) = c2.*fnm.*aiajFun(ni,nj,red1);
                        end
                    end
                end
                fprintf('\b\b\b: done!\n')
            else
                r    = r/R;
            end
                  
            zernCov = 0;
            for zi = 1:N
                
                ni = zN(zi);
                mi = zM(zi);
                
                for zj = 1:N
                    nj = zN(zj);
                    mj = zM(zj);
                    if (mi==mj) && (rem(abs(zi-zj),2)==0 || ((mi==0) && (mj==0)))
                        zerni = polynomials(zi,ni,mi,r,o);
                        zernj = polynomials(zj,nj,mj,r,o);
                        zernCov = zernCov + aiaj(zi,zj).*zerni.*zernj;
                    end
                end
                
            end
            
            nR = numel(r);
            zernPhaseCov = zeros(size(r));
            zernPhaseCovFun = @(y) quadgk(@(x) CphiAi(x,y), 0, Inf);
            parfor kR = 1:nR
                zernPhaseCov(kR) = zernPhaseCovFun(r(kR));
            end
            varargout{1} = atmVar + zernCov - c2*zernPhaseCov;
            if nargout>1
                varargout{2} = aiaj;
            end
            
            function out1 = CphiAi(x,rr)
                out1 = 0;
                for z = 1:N
                    n = zN(z);
                    m = zM(z);
                    krkr = m==0;
                    zerni = polynomials(z,n,m,rr,o);
                    CphiAiTerm1 = sqrt(n+1)*2^(0.5*(1-krkr))*...
                        (-1)^((n-m*(1-krkr))/2)*cos(m*o+pi*(1-krkr)*(((-1)^z-1))/4);
                    integrand = (1+(x/red1).^2).^(-11/6).*besselj(n+1,x).*besselj(m,rr*x);
                    out1 = out1 + zerni.*CphiAiTerm1.*integrand;
                end
                
            end
            
        end
                
        function varargout = residueStructureFunction(N,r1,o1,r2,o2,atm,tel)
        %% RESIDUESTRUCTUREFUNCTION
            
            persistent aiaj zMode zN zM R L0 L0r0Ratio L0RRatio c1 c2 cst ...
                red1 atmVar
            if isempty(aiaj)
                fprintf(' @(zernikeStats.residueStructureFunction)> Initializing aiaj ...')
                
                R         = tel.R;
                L0        = atm.L0;
                L0r0Ratio = (L0./atm.r0).^(5./3);
                L0RRatio  = (L0./R).^(2);
                
                c1 = (24.*gamma(6./5)./5).^(5./6);
                c2 = (gamma(11/6)^2/pi^(14/3))*c1*L0r0Ratio*L0RRatio;
                
                red1 = 2*pi*R/L0;
                atmVar    = c1.*...
                    (gamma(11./6).*gamma(5./6)./(2.*pi.^(8./3))).*L0r0Ratio;
                
                cst      = c1.*...
                    (gamma(11./6)./(2.^(5./6).*pi.^(8./3))).*...
                    L0r0Ratio;
                rho = abs(r1.*exp(1i*o1)-r2.*exp(1i*o2));
                atmCov   = ones(size(rho)).*atmVar;
                index         = rho~=0;
                u             = 2.*pi.*rho(index)./L0;
                atmCov(index) = cst.*u.^(5./6).*besselk(5./6,u);
                
                r1    = r1/R;
                r2    = r2/R;
                zMode = 1:N;
                zern = zernike(zMode,'logging',false);
                zN    = zern.n;%radialOrder(zMode)
                zM    = zern.m;%azimuthFrequency(zMode,zN)
                aiaj = zeros(N);
                for zi = zMode
                    ni = zN(zi);
                    mi = zM(zi);
                    for zj = 1:N
                        nj = zN(zj);
                        mj = zM(zj);
                        if (mi==mj) && (rem(abs(zi-zj),2)==0 || ((mi==0) && (mj==0)))
                            fnm = sqrt((ni+1)*(nj+1)).*(-1).^((ni+nj-mi-mj)/2);
                            aiaj(zi,zj) = c2.*fnm.*aiajFun(ni,nj,red1);
                        end
                    end
                end
                fprintf('\b\b\b: done!\n')
            else
                rho = abs(r1.*exp(1i*o1)-r2.*exp(1i*o2));
                atmCov   = ones(size(rho)).*atmVar;
                index         = rho~=0;
                u             = 2.*pi.*rho(index)./L0;
                atmCov(index) = cst.*u.^(5./6).*besselk(5./6,u);
                r1    = r1/R;
                r2    = r2/R;
            end
                 
            atmSF = 2*(atmVar - atmCov);
                 
            zernCov = 0;
            for zi = 1:N
                
                ni = zN(zi);
                mi = zM(zi);
                
                for zj = 1:N
                    nj = zN(zj);
                    mj = zM(zj);
                    if (mi==mj) && (rem(abs(zi-zj),2)==0 || ((mi==0) && (mj==0)))
                        zerni1 = polynomials(zi,ni,mi,r1,o1);
                        zernj1 = polynomials(zj,nj,mj,r1,o1);
                        zerni2 = polynomials(zi,ni,mi,r2,o2);
                        zernj2 = polynomials(zj,nj,mj,r2,o2);
                        zernCov = zernCov + aiaj(zi,zj).*...
                            (zerni1-zerni2).*(zernj1-zernj2);
                    end
                end
                
            end
            
            nR = numel(r1);
            zernPhaseCov = zeros(size(r1));
            zernPhaseCovFun = @(r1_,o1_,r2_,o2_) quadgk(@(x) CphiAi(x,r1_,o1_,r2_,o2_), 0, Inf);
            parfor kR = 1:nR
                zernPhaseCov(kR) = zernPhaseCovFun(r1(kR),o1(kR),r2(kR),o2(kR));
            end
            zernPhaseCov = c2*zernPhaseCov;
            varargout{1} = atmSF + zernCov - zernPhaseCov;
            if nargout>1
                varargout{2} = aiaj;
            end
            
            function out1 = CphiAi(x,rr1,oo1,rr2,oo2)
                out1 = 0;
                for z = 1:N
                    n = zN(z);
                    m = zM(z);
                    krkr = m==0;
                    zerni1 = polynomials(z,n,m,rr1,oo1);
                    zerni2 = polynomials(z,n,m,rr2,oo2);
                    gnm    = sqrt(n+1)*2^(0.5*(1-krkr))*...
                        (-1)^((n-m*(1-krkr))/2);
                    CphiAiTerm1 = cos(m*oo1+pi*(1-krkr)*(((-1)^z-1))/4);
                    CphiAiTerm2 = cos(m*oo2+pi*(1-krkr)*(((-1)^z-1))/4);
                    integrand1 = (1+(x/red1).^2).^(-11/6).*besselj(n+1,x).*besselj(m,rr1*x);
                    integrand2 = (1+(x/red1).^2).^(-11/6).*besselj(n+1,x).*besselj(m,rr2*x);
                    out1 = out1 + (zerni1-zerni2).*gnm.*...
                        ( CphiAiTerm1.*integrand1 - CphiAiTerm2.*integrand2 );
                end
                
            end
            
        end
        
        function out = residueOtf(N,rho_,gamma_,atm,tel)
            %% RESIDUEOTF
            
            out    = zeros(size(rho_));
            rho_   = rho_./tel.D;
            index  = rho_<1;
            rho_   = rho_(index);
            gamma_ = gamma_(index);
            n      = sum(index(:));
%             h = waitbar(0,'Patience ...');
            for k=1:n
                rho = rho_(k);
                gamma = gamma_(k);
                z   = rho.*exp(1i.*gamma);
%                              out = quad2d(@(o_,r_) r_, 0,pi/2,0, ro )*4
                ro = @(o) 0.5*rho.*cos(gamma-o)-0.5*sqrt(1-(rho.*sin(gamma-o)).^2);
                out(k) = quad2d(@residueOtfInt, 0,pi/2,0, ro );
%                 waitbar(k/n)
            end
%             close(h)
            out = out*16/pi;
            function out1 = residueOtfInt(oo,rr)
                z_ = rr.*exp(1i.*oo);
                z1 = z_-z/2;
                z2 = z_+z/2;
                r1 = abs(z1);
                o1 = angle(z1);
                r2 = abs(z2);
                o2 = angle(z2);
                out1 = rr.*exp(-0.5*zernikeStats.residueStructureFunction(N,r1,o1,r2,o2,atm,tel));
            end
        end
        
        function out = residueStrehlRatio(N,atm,tel)
            %% RESIDUESR
            
%             rho   = linspace(0,tel.D,21);
%             gamma = zeros(size(rho));
%             otf   = residueOtf(N,rho,gamma,atm,tel);
%             out   = trapz(rho,rho.*otf)*2*pi./tel.area;
            f   = @(x) x.*zernikeStats.residueOtf(N,x,zeros(size(x)),atm,tel);
            a = 0;b=tel.D;
            out = simpson(f,a,b,5)*2*pi./tel.area;
        end
        
        function out = residueEntrappedEnergy(N,eHalfSize,atm,tel)
            %% RESIDUEENTRAPPEDENERGY
            
%             rho   = linspace(0,tel.D,21);
%             gamma = zeros(size(rho));
%             otf   = residueOtf(N,rho,gamma,atm,tel);
%             out   = trapz(rho,rho.*otf)*2*pi./tel.area;
            f   = @(x) x.*zernikeStats.residueOtf(N,x,zeros(size(x)),atm,tel).*...
                2.*utilities.sombrero(1,2*pi.*eHalfSize.*x);
            a = 0;b=tel.D;
            out = simpson(f,a,b,5)*pi*eHalfSize^2;
        end
        
        function out = zern_aiaj(zi,ni,mi,zj,nj,mj,atm,tel)
            out = 0;
            if (ni==0) && (nj==0)
                out = Inf;
            elseif rem(abs(zi-zj),2)==0 || (mi==0 && mj==0)
                red = pi*tel.D;
                out = 2/(pi*tel.R^2)*...
                    sqrt((ni+1)*(nj+1))*(-1)^((ni+nj-mi-mj)/2)*(mi==mj).*...
                    quadgk( @(x) phaseStats.spectrum(x,atm).*...
                    besselj(ni+1,red.*x).*besselj(nj+1,red.*x)./x,0,Inf);
            end
        end
        
        function out = zern_phiai(r,o,j,n,m,atm,tel)
            Inm = @(r) quadgk( @(f) phaseStats.spectrum(f,atm).*...
                besselj(n+1,pi.*f.*tel.D).*...
                besselj(m,2.*pi.*f.*r),0,Inf);
            out = arrayfun(Inm,r);
            krkr = m~=0;
            out = out.*2.*sqrt(n+1).*2.^(0.5*krkr).*(-1).^((n-m*krkr)/2).*...
                cos(m.*o+pi.*krkr.*((-1).^j-1)./4)./tel.R;
        end
        
        
        function map = zern_resid_var(N,atm,tel,r,o)
            zern = zernike(1:N);
            
            mapFun = @(rr,oo) quadgk(@(f)integrand(f,rr,oo),0,Inf);
            tic
            map = arrayfun( mapFun , r , o );
            toc
            function out = integrand(f,rr,oo)
                cum_i = 0;
                for i=1:N
                    cum_ip = 0;
                    for ip=1:N
                        flag = (zern.m(i)==zern.m(ip)) && ...
                            ( rem( abs( zern.j(i)-zern.j(ip) ),2 )==0 || ...
                            ( (zern.m(i)==0) && (zern.m(ip)==0) ) );
                        if flag
                            cum_ip = cum_ip + ...
                                sqrt(zern.n(ip)+1).*...
                                (-1).^((zern.n(ip)-zern.m(ip))/2).*...
                                %{
                                zernike.fun(zern.j(ip),zern.n(ip),zern.m(ip),rr,oo).*...
                                %}
                                besselj(zern.n(ip)+1,pi.*tel.D.*f);
                        end
                    end
                    cum_ip = cum_ip./(pi.*f.*tel.R);
                    deltaNot = 1 - zern.m(i)==0;
                    alpha = pi.*((-1).^zern.j(i)-1).*deltaNot./4;
%                     cum_ip = cum_ip - ...
%                         2.*(-1).^(zern.m(i).*(zern.m(i)==0)/2).*...
%                         2.^(deltaNot/2).*...
%                         cos(zern.m(i).*oo + alpha).*...
%                         besselj(zern.m(i),2.*pi.*f.*rr);
                    cum_i = cum_i + ...
                        sqrt(zern.n(i)+1).*...
                        (-1).^((zern.n(i)-zern.m(i))/2).*...
                        %{
                        zernike.fun(zern.j(i),zern.n(i),zern.m(i),rr,oo).*...
                        %}
                        besselj(zern.n(i)+1,pi.*tel.D.*f).*...
                        2.*cum_ip./tel.R;
                end
%                 out = phaseStats.spectrum(f,atm).*(2.*pi.*f + cum_i);
                out = phaseStats.spectrum(f,atm).*(cum_i);
            end
            
        end
        
    end
    
            
end
    
function out = polynomials(j,n,m,r,o)
%FUN Zernike polynomials mathematical expression
%
% out = zernike.fun(j,n,m,z) computes the Zernike polynomial
% defined by the mode j, the radial order n, the azimuthal
% frequency m at the locations given by the polar coordinate r
% and o; both r and o must be column vectors and r must be in the range
% [0,1]

if n==0
    out = ones(size(r));
else
    krkr = m~=0;
    out = sqrt(n+1).*R_fun1(r,n,m).*2.^(0.5.*krkr).*...
        cos(m.*o+pi.*krkr.*((-1).^j-1)./4);
end

    function R = R_fun1(r,n,m)
        R=zeros(size(r));
        for s=0:(n-m)/2
            R = R + (-1).^s.*prod(1:(n-s)).*r.^(n-2.*s)./...
                (prod(1:s).*prod(1:((n+m)/2-s)).*prod(1:((n-m)/2-s)));
        end
    end
%     function R = R_fun2(r,n,m)
%         s=0:(n-m)/2;
%         fs = (-1).^s.*gamma(n-s+1)./...
%             ( gamma(s+1).*gamma((n+m)/2-s+1).*gamma((n-m)/2-s+1) );
%         R = bsxfun(@times,fs,bsxfun(@power,r,n-2.*s));
%     end

end

function out = aiajFun(ni,nj,a)
out = ((1/2.^(ni + nj + 1)).*...
    (729*a.^(ni + nj + 2).*...
    gamma(ni/2 + nj/2 + 1).*gamma(5/6 - nj/2 - ni/2).*...
    gamma(nj/2 - ni/2 + 17/6).*gamma(ni/2 - nj/2 + 17/6).*...
    gamma(ni/2 + nj/2 + 23/6).*gamma(2/3).*...
    hypergeom([ni/2 + nj/2 + 1, ni/2 + nj/2 + 2, ni/2 + nj/2 + 3/2], ...
    [ni + 2, nj + 2, ni + nj + 3, ni/2 + nj/2 + 1/6], a^2) + ...
    275*2^(ni + nj + 2).*3^(1/2)*pi^(1/2)*a^(11/3).*...
    gamma(ni + 2).*gamma(nj + 2).*...
    gamma(ni/2 + nj/2 - 5/6).*gamma(5/6)^2.*...
    hypergeom([11/6, 7/3, 17/6], ...
    [11/6 - nj/2 - ni/2, nj/2 - ni/2 + 17/6, ni/2 - nj/2 + 17/6, ni/2 + nj/2 + 23/6], a^2)))./...
    (2430*gamma(ni + 2).*gamma(nj + 2).*...
    gamma(nj/2 - ni/2 + 17/6).*gamma(ni/2 - nj/2 + 17/6).*...
    gamma(ni/2 + nj/2 + 23/6).*gamma(2/3)*gamma(5/6));
end

function out = radialOrder(zi)
    out = ceil( (-3 + sqrt( 9 + 8*(zi-1) ))/2 );
end
function out = azimuthFrequency(zi,ni)
    out = abs((ni - 2*( zi - (ni+1).*ni/2 -1 )));
end