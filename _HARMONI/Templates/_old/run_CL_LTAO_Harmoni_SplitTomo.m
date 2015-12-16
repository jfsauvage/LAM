%% SOURCE
ngs =source;

%% ATMOSPHERE

atm = atmosphere(photometry.V,20e-2,30,...
    'fractionnalR0',[0.5,0.3,0.2],'altitude',[0e3,5e3,12e3],...
    'windSpeed',[10,5,20],'windDirection',[0,pi/2,pi]);
%% TELESCOPE
nL   = 20;
nPx  = 6;
nRes = nL*nPx;
D    = 8;
d    = D/nL; % lenslet pitch
samplingFreq = 500;

tel = telescope(D,'resolution',nRes,...
    'fieldOfViewInArcsec',120,'samplingTime',1/samplingFreq);

%%
QuadCell = shackHartmann(1,8,1);
ngs = ngs.*tel*QuadCell;

QuadCell.INIT

+QuadCell;
figure
imagesc(QuadCell.camera,'parent',subplot(3,2,[1,4]))
slopesDisplay(QuadCell,'parent',subplot(3,2,[5,6]))

%% WFS gain calibration
% <latex>
% The WFS must be calibrated such as for 1rd of tip--tilt wavefront , it will
% measured a slopes of 1rd.
% To do so, the \oop{shackHartmann}{pointingDirection} property is set on-axis
% </latex>
QuadCell.pointingDirection = zeros(2,1);
ngs.wavelength = photometry.K;

pixelScale = ngs.wavelength/...
    (2*D*QuadCell.lenslets.nyquistSampling);
tipStep = pixelScale/2;
nStep   = floor(nPx/3)*2;
sx      = zeros(1,nStep+1);
u       = 0:nStep;
QuadCell.camera.frameListener.Enabled = false;
QuadCell.slopesListener.Enabled = false;
warning('off','oomao:shackHartmann:relay')
for kStep=u
    ngs.zenith = -tipStep*kStep;
    +ngs;
    drawnow
    sx(kStep+1) = median(QuadCell.slopes(1:end/2));
end
warning('on','oomao:shackHartmann:relay')
Ox_in  = u*tipStep*constants.radian2arcsec;
Ox_out = sx*ngs.wavelength/D/2*constants.radian2arcsec;
%figure
plot(Ox_in,Ox_out)
%grid
slopesLinCoef = polyfit(Ox_in,Ox_out,1);
QuadCell.slopesUnits = 1/slopesLinCoef(1);

ngs.zenith = 0;
QuadCell.pointingDirection = [];


%% SH WFS
wfs = shackHartmann(nL,nRes,0.85);
%wfs.lenslets.nyquistSampling = 0.5;

ngs.wavelength = photometry.R;
ngs = ngs.*tel*wfs;

wfs.INIT

+wfs;
figure
imagesc(wfs.camera,'parent',subplot(3,2,[1,4]))
slopesDisplay(wfs,'parent',subplot(3,2,[5,6]))

wfs.camera.frameListener.Enabled = true;
wfs.slopesListener.Enabled = true;
%% WFS gain calibration
% <latex>
% The WFS must be calibrated such as for 1rd of tip--tilt wavefront , it will
% measured a slopes of 1rd.
% To do so, the \oop{shackHartmann}{pointingDirection} property is set on-axis
% </latex>
wfs.pointingDirection = zeros(2,1);

pixelScale = ngs.wavelength/...
    (2*d*wfs.lenslets.nyquistSampling);
tipStep = pixelScale/2;
nStep   = floor(nPx/3)*2;
sx      = zeros(1,nStep+1);
u       = 0:nStep;
wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled = false;
warning('off','oomao:shackHartmann:relay')
for kStep=u
    ngs.zenith = -tipStep*kStep;
    +ngs;
    drawnow
    sx(kStep+1) = median(wfs.slopes(1:end/2));
end
warning('on','oomao:shackHartmann:relay')
Ox_in  = u*tipStep*constants.radian2arcsec;
Ox_out = sx*ngs.wavelength/d/2*constants.radian2arcsec;
%figure
%plot(Ox_in,Ox_out)
%grid
slopesLinCoef = polyfit(Ox_in,Ox_out,1);
wfs.slopesUnits = 1/slopesLinCoef(1);
%%
% <latex>
% The source is reset on--axis and the WFS is set to always be aligned to
% the source by setting \oop{shackHartmann}{pointingDirection} to empty.
% </latex>
ngs.zenith = 0;
wfs.pointingDirection = [];


%%
dmCrossCouplingCoef = 0.35;
bifa = influenceFunction('monotonic',dmCrossCouplingCoef);
figure,show(bifa,'parent',subplot(1,2,1))
title('Monototic influence function')
bifb = influenceFunction('overshoot',dmCrossCouplingCoef);
show(bifb,'parent',subplot(1,2,2))
title('Overshooted influence function')

dm = deformableMirror(nL+1,'modes',bifa,...
    'resolution',tel.resolution,...
    'validActuator',wfs.validActuator);


wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled = false;


bifaLowRes = influenceFunction('monotonic',dmCrossCouplingCoef);
dmLowRes = deformableMirror(nL+1,'modes',bifaLowRes,'resolution',nL+1,...
    'validActuator',wfs.validActuator);


%%
% <latex>
% and setup the optical path before the DM/WFS subsystem:
% </latex>
ngs = ngs.*tel;
%%
% <latex>
% The interaction matrix is generated by pushing all the valid actuators
% sequentially and saving the corresponding WFS measurements in a matrix.
% The process is done with the \oom{deformableMirror}{calibration} method
% of the \oo{deformableMirror} class.
% The arguments of the method are the deformable mirror, the wavefront
% sensor, the calibration source and the actuator amplitude in meter.
% In \oomao, the generation of the interaction matrix is vectorized meaning
% it can be done in 1 step at the expense of requiring a lot of memory.
% Here the process is divided in as many steps as actuators accross the pupil.
% </latex>
calibDm = calibration(dm,wfs,ngs,ngs.wavelength,nL+1,'cond',1e2);

%%
tel = tel + atm;
figure
imagesc(tel)
ngs = ngs.*tel*wfs;


%%
telLowRes = telescope(tel.D,'resolution',nL+1,...
    'fieldOfViewInArcsec',120,'samplingTime',1/500);
telLowRes.pupil = wfs.validActuator;
%%
% <latex>
% The atmosphere is now bound to the new telescope, the \matcall{ngs} is
% propagated through the atmosphere to the new telescope and the piston
% free wavefront is saved.
% </latex>
telLowRes= telLowRes + atm;
ngs = ngs.*telLowRes;
phase = ngs.meanRmOpd;

%%
lgsAst = source('asterism',{[4,arcsec(30),0]},'height',90e3);
%lgsAst = source('asterism',{[1,arcsec(0),0]},'height',90e3);

% figure, imagesc(tel, [ngs,lgsAst])
lgsAst_slmmse = slopesLinearMMSE(wfs,tel,atm,lgsAst,'mmseStar',ngs,'NF',1024);
%%
% <latex>
% The LGSs are propagated through the atmosphere and the telescope to
% the wavefront sensor.
% The LMMSE wavefront estimate is obtained by multiplying the
% \oo{slopesLinearMMSE} object and the \oo{shackHartmann} object.
% \matcall{lgsAst\_ps\_e} is removed from the NGS wavefront to obtain the
% zero--piston residual wavefront as well as the residual wavefront error
% rms.
% </latex>
lgsAst = lgsAst.*tel*wfs;
lgsAst_ps_e = tools.meanSub( lgsAst_slmmse*wfs , wfs.validActuator );
ngs = ngs.*telLowRes*{wfs.validActuator,-lgsAst_ps_e*ngs.waveNumber};
lgsAst_ps_eRes = ngs.meanRmOpd;
lgsAst_ps_eResRms = ngs.opdRms*1e9;

%%
TTAst = source('asterism',{[3,arcsec(20),0]},'wavelength',photometry.H);
%TTAst = source('asterism',{[1,arcsec(0),0]},'wavelength',photometry.H);

TTAst = TTAst .* tel*QuadCell;

%%
science = source('wavelength',photometry.J);
cam = imager(tel);
%%
% <latex>
% The \oo{atmosphere} object is detached from the telescope and the
% science star is propagated through the telescope to the science camera
% producing a perfect diffraction limited image.
% The camera display can also be set to update itself when a new camera frame
% is generated.
% </latex>
tel = tel - atm;
science = science.*tel*cam;
figure(31416)
imagesc(cam,'parent',subplot(2,1,1))
%cam.frameListener.Enabled = true;
%%
% <latex>
% The diffraction limited image is set as the reference frame allowing to
% compute the image Strehl ratio on the subsequent frame captures.
% </latex>
cam.referenceFrame = cam.frame;
+science;
fprintf('Strehl ratio: %4.1f\n',cam.strehl)


tel = tel + atm;
+science;
fprintf('Strehl ratio: %4.1f\n',cam.strehl)

%%
% <latex>
% \section{Laser Tomography Adaptive Optics}
% In this section, the performance of an LTAO system on the same telescope
% is compared to the two NGS AO systems of the former section.
% \newline
% Lets reset the atmosphere to start with the same initial conditions
% that the NGS AO systems and also the DM.
% </latex>

%reset(tel)
dm.coefs = zeros(dm.nValidActuator,1);
%%
% <latex>
% The optical paths of both the LGS constellation and the science star are
% defined and the display is updated accordingly.
% Note how the WFS display is now showing the imagelets and slopes
% corresponding to the 3 LGSs.
% The WFS frame and slope listeners are turned off to speed up the
% computation.
% </latex>
science = science.*tel*dm*cam;
lgsAst = lgsAst.*tel*dm*wfs;
figure(31416)
imagesc(cam,'parent',subplot(2,1,1))
subplot(2,1,2)
h = imagesc(catMeanRmPhase(science));
axis xy equal tight
colorbar
wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled = false;

TTAst = TTAst.*tel*dm*QuadCell;
%%
% <latex>
% The logging of the wavefront variance of the science object is turned on.
% The LGS LMMSE object is set to return the wavefront estimate in a vector and
% the iterative solver will use the previous estimate as first guess.
% </latex>

flush(cam)
cam.clockRate    = 1;
exposureTime     = 100;
cam.exposureTime = exposureTime;
startDelay       = 20;

flush(cam)
cam.startDelay = startDelay;
set(science,'logging',true)
set(science,'phaseVar',[])
lgsAst_slmmse.wavefrontSize = [dm.nValidActuator,1];
lgsAst_slmmse.warmStart = true;
cam.frameListener.Enabled = false;

gain_pol = 0.7;
F = 2*bifaLowRes.modes(wfs.validActuator,:);
iF = pinv(full(F),1e-1);

gain_cl = 0.4;
u_ngs = zeros(dmLowRes.nValidActuator,1);

QuadCell.camera.photonNoise = true;
%QuadCell.camera.readOutNoise = 0;

wfs.camera.photonNoise = true;
%wfs.camera.readOutNoise = 0;


%% piston-tip-tilt removal 
t = ones(wfs.nValidLenslet,1);
TT = [t 0*t; 0*t,t];
SlopeTTRem = eye(2*wfs.nValidLenslet) - TT*pinv(TT);

zern = zernike(2:3,'resolution',nL+1, 'pupil',dmLowRes.validActuator);
TT = iF*zern.modes(dmLowRes.validActuator,:);
DMTTRem = eye(dmLowRes.nValidActuator) - TT*pinv(TT);


%%
% <latex>
% The loop is closed for one full exposure of the science camera.
% </latex>
nIteration = startDelay + exposureTime;

for k=1:cam.startDelay + cam.exposureTime
    tic
    % Objects update
    +tel;
    +lgsAst;
    +TTAst;
    +science;
    % Pseudo-open-loop controller    
    dm.coefs = (1-gain_pol)*dm.coefs + ...
        gain_pol*iF*( lgsAst_slmmse*( bsxfun( @minus, SlopeTTRem*wfs.slopes, calibDm.D*dm.coefs ) ) );


%dm.coefs = iF*( lgsAst_slmmse*( bsxfun( @minus, SlopeTTRem*wfs.slopes, calibDm.D*dm.coefs ) ) );

% Remove TT from DM
    dm.coefs = DMTTRem*dm.coefs;
    % Add TT from TT stars, % Closed-loop controller
    u_ngs = u_ngs - TTAst(1).wavelength/8*gain_cl*TT*mean(QuadCell.slopes,2); 
    
    % add ngs controls back to dm commands
    dm.coefs = dm.coefs - u_ngs;
    % Display
     set(h,'Cdata',catMeanRmPhase(science))
     drawnow
     toc
end
imagesc(cam)
set(h,'Cdata',catMeanRmPhase(science))

%%
% <latex>
% The time series of wavefront variance is read from the
% \oop{stochasticWave}{phaseVar} property of the \matcall{science}
% object. 
% The Strehl ratio is estimated from the residual phase variance using the
% Marechal approximation and it is compared to the Strehl ratio derived
% from the long exposure image.
% </latex>
var_wfe_ltao = reshape(science.phaseVar(1:nIteration*2),2,[])';
wfe_ltao = sqrt(var_wfe_ltao)/science.waveNumber*1e6;
marechalStrehl_ltao = 1e2*exp(-mean(var_wfe_ltao(startDelay:end,2)))
psfStrehl_ltao =1e2*cam.strehl

