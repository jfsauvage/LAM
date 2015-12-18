PRO fourier_ao, sim_params

; SIM_PARAMS : structure containing all simulation parameters


; ---------------------------
; Copy Structure to parameters
; ---------------------------
anglezenith = sim_params.anglezenith & sep = sim_params.sep & nb_gs = sim_params.nb_gs
LGS_H = sim_params.LGS_H & varnoise_wfs = sim_params.varnoise_wfs & lambda_seeing = sim_params.lambda_seeing
lambda_imaging = sim_params.lambda_imaging & L0 = sim_params.L0 & seeing0 = sim_params.seeing0
cn2true = sim_params.cn2true & vent = sim_params.vent & arg_v = sim_params.arg_v
h = sim_params.h & D = sim_params.D & dimpup = sim_params.dimpup & dim = sim_params.dim
nb_act = sim_params.nb_act & pitch_m = sim_params.pitch_m & fech = sim_params.fech
td = sim_params.td & h_dm = sim_params.h_dm & pitchs_dm = sim_params.pitchs_dm
alpha = sim_params.alpha & tempo = sim_params.tempo & fitting_kw = sim_params.fitting_kw
recons_h = sim_params.recons_h & cn2recons = sim_params.cn2recons & err_R0 = sim_params.err_R0
theta = sim_params.theta
beta_tab = sim_params.beta_tab & LSE = sim_params.LSE & aliasing_kw = sim_params.aliasing_kw
tomo = sim_params.tomo & file = sim_params.file ;& ngs = sim_params.ngs


; geometrical projection
zenith =  anglezenith*!pi/180.
seeing = seeing0/(cos(zenith))^(0.6)

; Frequencies array definition
local_L = D*dim/dimpup        ; Screen size [m]
local_dim = dim               ; Screen size [px]

; Frequency array : variable fx, fy, f, arg_f
IF local_dim mod 2 then $
    fx = eclat(((indgen(local_dim)-(local_dim-1)/2)/local_L)#$
               (1.+fltarr(1,local_dim))) $ 
ELSE fx = eclat(((indgen(local_dim)-local_dim/2)/local_L)#$
                (1.+fltarr(1,local_dim)))
fy = transpose(fx)

; Frequency modulus and argument
f     = sqrt(fx*fx+fy*fy)
arg_f = atan(fy,fx)
s = (size(f))(1)              ; linear size of F
;------------------------------------------------------------


;------------------------------------------------------------
; Definition of atmospherical profile
;------------------------------------------------------------
; normalized Cn2
Cn2_init =  Cn2true/total(cn2true)
h_init   = h/cos(zenith)
nl       = n_elements(h_init)       ; turbulent layers number

IF NOT keyword_set(LGS_H) THEN lgs_h =  90000. ; 90km by default
dilat    =  (lgs_h - h_init)/lgs_h
h_init   = (h_init)/dilat 

; true profile is initial profile :
h_vrai =  h_init
Cn2_vrai =  Cn2_init

; r0 computation 
rad2arcsec = 180.*3600./!pi
r0_seeing  = 0.976*lambda_seeing/seeing*rad2arcsec ; Fried's param. [m]
r0         = (r0_seeing*(lambda_imaging/lambda_seeing)^1.2) ; imaging Fried's param. [m]
hbar       = (total(cn2_init*(h)^(5/3.)))^(3./5.)


; wind profile : 
; TBC : as many wind layers as turbulent layers
; wind speed modulous given by VENT
; random repartition of winds direction
; arg_v  = ((randomu(-12345, nl)-0.5)*!pi/180.) ; wind dir.  [rad]
wind =  transpose([[vent*cos(arg_v)], [vent*sin(arg_v)]])

;------------------------------------------------------------
; AO system generation 
;------------------------------------------------------------
; wfs pitch:
pitchs_wfs = replicate(pitch_m,Nb_gs)
; WFS sampling frequencies
fech_tab = replicate(fech,Nb_gs) ; in Hz
ti = 1/fech_tab

; display field of view
sepchampvisu = max([sep/2., 1])

; number of optimisation directions
nb_dir_theta = n_elements(theta)/2.

; number of performance direction
nb_dir_perf = n_elements(beta_tab)/2.

; arrays initialisation for residual (nm^2) in each performance direction
resall = fltarr(nb_dir_perf)

; Arrays initialisation for Strehl Ratio (in %) in each performance direction
SR = fltarr(nb_dir_perf)

; arrays initialisation for res dsp (nm^2) in each performance direction
dspall = fltarr(dim, dim, nb_dir_perf)



;-----------------------------------------------------
; LGS part
;-----------------------------------------------------
; default value 90 km altitude
IF NOT keyword_set(LGS_H) THEN lgs_h =  90000.
hxxx     =   recons_h 
hxxx = hxxx/cos(zenith)

; dilatation factor due to LGS geometry
dilat =  (lgs_h - hxxx)/lgs_h
hxxx = (hxxx)/dilat
h_recons = hxxx 

Cn2_recons =  cn2recons
Cn2_recons =  Cn2_recons/total(cn2_recons)

; Normalized estimated profile
; Cn2_recons = Cn2_recons/total(Cn2_recons)
DSP_tab_recons = fltarr(n_elements(h_recons)*s, s)
cst = 0.0229                  ; normalisation constant, power spectrum of
                              ; phase aberrations after propagation through
                              ; kolmogorov turbulence 0.0229 r0^5/3 k-11/3
; if an error is introduced on Cn2, it is accounted for here
FOR i = 0, n_elements(h_recons)-1 DO $
    DSP_tab_recons(s*i:s*(i+1)-1, *)=cst*(Cn2_recons[i]^(-3./5.)*r0/err_R0)^(-5./3.)*(f^2+(1/L0)^2)^(-11./6.) ;

DSP_tab_vrai = fltarr(n_elements(h_vrai)*s, s)
FOR i = 0, n_elements(h_vrai)-1 DO $
    DSP_tab_vrai(s*i:s*(i+1)-1, *)=cst*(Cn2_vrai[i]^(-3./5.)*r0)^(-5./3.)*(f^2+(1/L0)^2)^(-11./6.) ;

;-----------------------------------------------------
;Computation of global reconstructor WMAP = POPT # WTOMO
;-----------------------------------------------------
; without aliasing => Wmap direct with frequency modulous only.
;                  => Fastest
; with    aliasing => Wmap with (x,y) Wmap, 2 x more computation
; For reconstruction error and noise,
; without aliasing, Noise and reconstruction error is the same 
 
IF aliasing_kw EQ 0 THEN BEGIN 
   Wflag = 'W2'
   
; option for (weighted) least squared error
; if LSE =1B, WLSE : (MtCb-1M)-1MtCb-1
; IF LSE THEN seuil = 1000. ELSE seuil =  1e9
; SEUIL to be adapted. 
; mostly for LSE case
; MMSE : 1e9
; LSE : 100, or less...
   condmax = 1e6
   IF keyword_set(LSE) THEN BEGIN
      print, '--------------------------'
      print, 'LSE RECONSTRUCTOR'
      condmax = 1e0           ; CONDMAX used for inverse computation in Popt (done by LA_TSVD)
      seuil = 1e2
      print,  'condmax = ', condmax
      print,  'seuil = ', seuil
      print, '--------------------------'
   ENDIF ELSE BEGIN
      print, '--------------------------'
      print, 'MMSE RECONSTRUCTOR'
      condmax = 1e0           ; CONDMAX used for inverse computation in Popt (done by LA_TSVD)
      seuil = 1e6
      print,  'condmax = ', condmax
      print,  'seuil = ', seuil
      print, '--------------------------'
   ENDELSE
      

wmap = calc_reconstructor_noalias(f, arg_f, pitchs_wfs,pitchs_DM, $
                                  Nb_gs,alpha, theta,varnoise_wfs, $
                                  DSP_tab_recons, $
                                  h_recons, h_dm, $
                                  seuil, condmax, Wflag, $
                                  fitting_kw = fitting_kw, $
                                  LSE = LSE)

ENDIF ELSE BEGIN
   Wflag = 'W2'

; to include aliasing, we need to compute reconstruction matrix including
; aliasing of measurement.
; compute Cb_alias, following a model of SH aliasing expression (B.Neichel p 143):
   
   a = max(pitchs_wfs)
   Cb_Alias = calc_cbalias(f, fx, fy, arg_f, a,Nb_gs,alpha, h_recons ,Cn2_recons, L0, r0)

;then compute W_alias including knownledge of Cb_alias and / or of noise
;Cb_alias*0. noise only
;sig*0. aliasing only
   
   IF LSE THEN seuil = 1000. ELSE seuil =  1e9
;seuil = 1e9
   print,  'seuil = ', seuil
;Seuil is for inversion of W with Cbalias only. In this case, some frequencies
;are not inversible, then tsvd is required. (Cb_alias = 0 for f_x ou f_y =0)
   IF keyword_set(LSE) THEN BEGIN
      print, '--------------------------'
      print, 'LSE RECONSTRUCTOR'
      print, '--------------------------'
   ENDIF
; LSE reconstructor does not account for Cb_alias   
;Cb_alias = fltarr(2*nb_gs*s, 2*nb_gs*s)
   
; NOISE vector definition in case of aliasing : interlace X and Y noises
   varnoise_wfs_xy = reform(transpose(rebin(varnoise_wfs, nb_gs, 2)), nb_gs * 2.)

; case aliasing_kw = 1 : No account for Cb_alias in reconstructor (KW-1. = 0.000)
; case aliasing_kw = 2 : Accounts for Cb_alias in reconstructor (KW-1. = 1.000)
   Wmap = calc_reconstructor_alias(f, fx, fy, arg_f,pitchs_wfs,pitchs_DM,Nb_gs,alpha, theta,varnoise_wfs_xy, DSP_tab_recons,$
                                   h_recons, h_dm, $
                                   seuil, cb_alias*(aliasing_kw - 1.), $
                                   condmax, $
                                   Wflag, $
                                   fitting_kw = fitting_kw, $
                                   LSE = LSE)

ENDELSE 


;-----------------------------------------------------
; Computation of tomographic residual DSP
;-----------------------------------------------------
ech = 2.                      ;surech
rt = dimpup/2.
ri = 0.                       ;obstruction centrale
polaire2, RT = rt, RHO= rho,  PHI=phi, /ENTRE4PIXELS, $
          MASQUE= pup,oc = ri, largeur= dimpup


; Loop on the performance direction wanted
FOR bbb = 0, nb_dir_perf-1 DO BEGIN   
   beta = [beta_tab(0, bbb), beta_tab(1, bbb)]
      
   IF aliasing_kw EQ 0 THEN BEGIN
      
      dsp_res = calc_dspresidual_noalias(f, arg_f, pitchs_wfs, pitchs_dm, Nb_gs,alpha,beta,varnoise_wfs, DSP_tab_vrai, h_vrai, h_dm, Wmap, $
                                         tempo = tempo, td, ti, Wind, $
                                         fitting_kw = fitting_kw,  err_recons = err_recons, err_noise = err_noise)
      
   ENDIF ELSE BEGIN
      
; NOTE : an errored Cb_alias' could be introduced here to perform the
;        scenario:
; direct model with Cb_alias, reconstruction with Cb_alias'
      dsp_res = calc_dspresidual_alias(f, fx, fy, arg_f, pitchs_wfs, pitchs_dm, Nb_gs,alpha,beta,varnoise_wfs_xy, DSP_tab_vrai, $
                                       h_vrai, h_dm, Wmap, $
                                       tempo = tempo, td, ti, Wind, Cb_alias, $
                                       fitting_kw = fitting_kw, $
                                       alias_alone = alias_alone,  err_recons = err_recons, err_noise = err_noise)
      
   ENDELSE
   
   ; size of screen [m]
   L =  D*dim/dimpup
   pixsize =  1./float(L)
   resva1 = calc_var_from_psd(dsp_res, pixsize =  pixsize, $
                              DD = D)
   
   print,'Performance direction : ', strc(bbb), $
         '  Residual [nm] : ', strc(sqrt(resva1)*lambda_imaging * 1e9/2./!pi), $
         '  Strehl Ratio [%] : ',  strc(exp(-(resva1))*100.)

   resall(bbb) =  sqrt(resva1)*lambda_imaging * 1e9/2./!pi
   SR(bbb) = exp(-(resva1))*100.
   dspall(*,*, bbb) =  dsp_res 
ENDFOR  

; RESALL is the residual [nm] for each performance direction
; DSPALL are the 2D-DSP for each performance direction

save, file = file, resall, SR, dspall, f, dsp_tab_vrai, sim_params, /compress
print, 'DATA saved in: ', file

END 

