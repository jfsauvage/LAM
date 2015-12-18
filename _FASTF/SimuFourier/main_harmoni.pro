PRO main_harmoni, varnoise_wfs = varnoise_wfs, lambda_imaging = lambda_imaging, $
                  dimpup = dimpup, anglezenith = anglezenith, adim = adim, $
                  ind_turb_case = ind_turb_case, beta_tab = beta_tab, L0 = L0

rad2arcsec = 180./!dpi*60.*60.

; -- Simulation geometry --
anglezenith = anglezenith     ; [degrees!!!] Observation angle relative to zenith  --> DO FROM 30° and 60°
sep = 0.5                      ; [arcmin] Careful > 1.6 arcmin to ensure NO unseen area
nb_gs = 6.                    ; Number of Guide Stars
LGS_H = 9e4; 9e4              ; Altitude for LGS (0 for infinity) 

; -- Simulation photometry --
lambda_wfs = 0.600e-6                                                             ; wfs wavelength [m]
lambda_seeing = 0.5e-6                                                            ; r0 wavelength  [m]
lambda_imaging =  lambda_imaging                                                  ; imaging wavelength [m] 

varnoise_wfs = varnoise_wfs * (fltarr(nb_gs)+1)*(lambda_wfs/lambda_imaging)^2.  ;analyse noise in rad^2
; [1] == Variance de bruit "pleine pupille" (?) a la longueur d'onde ASO
; Si c'est stationnaire: c'est la meme chose que variance de bord-a-bord de sous-pupille... 
; Cela dit [1] represente la somme du RON, bruit de photon etc... donc c'est pas exactement la meme chose...


; -- Simulation turbulence definition --
IF NOT KEYWORD_SET(L0) THEN BEGIN
    L0 = 25.;                  ; outer scale    [m];
ENDIF
;seeing0 = 0.80                ; zenithal seeing (arcsec)
;cn2true =  [53.6, 2.5, 4.3, 11.3, 9.6, 2.9, 5.8, 4.2, 5.8]/100.      ; Cn2 profile 
;vent = [15., 13., 13., 9., 9., 15., 25., 40., 21.]                   ; wind speed modulous [m/s]
;arg_v = [1.75,-1.67,-0.95,0.23,0.63,0.44,-0.66,-0.38,-2.29]          ; wind dir.  [rad]
;h = [47., 140., 281., 562., 1125., 2250., 4500., 9000., 18000.]      ; layers altitude [m]

turb_profile = readfits('C:\Users\jfsauvage\Desktop\HARMONI\ATC-E2E_FOURRIER\FOURIER\SimuFourier\profil_turbulent_eso.fits')
;turb_profile = readfits('D:\Simulation Codes\IDL\profil_turbulent_eso.fits')
turb_case = (['Jmedian', 'JQ1', 'JQ2', 'JQ3', 'JQ4'])[ind_turb_case]
print, '--- The turbulence case is: ', turb_case, ' ---'
h = reform(turb_profile[1, 0:34])                         ; layers altitude [m]
vent = reform(turb_profile[2, 0:34])                      ; wind speed modulous [m/s]
cn2true = reform(turb_profile[3+ind_turb_case, 0:34])     ; Cn2 profile
cn2true = cn2true/total(cn2true) 
arg_v = reform(turb_profile[8, 0:34])                     ; wind dir.  [rad] (choosen randomly but once and for all)

r0 = reform(turb_profile[3: 7,35])                        ; Fried parameter in [m]
tau0 = reform(turb_profile[3: 7,36])                      ; Coherence time in [s]
alpha = reform(turb_profile[3: 7,37])                     ; A single common wind profile is provided which has to be multiplied by the $
                                                          ; corresponding alpha coefficients to fit to the tau0 measured at the site by MASS. 
seeing0 = 0.98*lambda_seeing/r0[ind_turb_case]*rad2arcsec ; Zenithal seeing (arcsec)
vent = vent * alpha[ind_turb_case]                        ; Wind speed modulous [m/s]

;; Anisoplanatic Angle
;print, 'Anisoplanatic Angle', 0.314 * r0[ind_turb_case]/(total(h^(5./3)*cn2true)/total(cn2true))^(3./5.)*rad2arcsec
;; Greenwood time delay
;print, 'Greenwood time delay', 0.314 * r0[ind_turb_case]/(total(vent^(5./3)*cn2true)/total(cn2true))^(3./5.)*rad2arcsec

; -- Simulation parameters --
D = 37.0;                             ; Telescope diameter [m]
dimpup = dimpup;                      ; Pupil size [px]
dim = adim                            ; Screen size [px]                                      

; -- AO system definition --
nb_act = 74.                          ; number of act
pitch_m = float(D)/(float(nb_act));float(D)/(float(nb_act)-1)  ; WFS pitch [m]
fech = 500.                           ; AO temporal sampling frequency [Hz]
td = 3.0*1.e-3                        ; loop latency [s]
h_dm = [0.0]                          ; DM altitude
pitchs_dm = [pitch_m]                 ; DM pitch [m]

; -- Guide sources Geometry--
aa = 2.*!pi/Nb_gs                                       ; etoiles en cercle
alpha  = transpose([[cos(indgen(nb_gs)*aa)],  $         ; analyse directions [arcmin]
                    [sin(indgen(nb_gs)*aa)]] *sep/2.)   
ngs = [[0., 0.]]/60                                     ; NGS Geometry [arcmin]--> Not implemented yet!

; Reconstruction turbulence a priori 
recons_h = [0.0]              ; supposed altitudes for estimation 
cn2recons = [1.0]             ; supposed CN2 for estimation
; Introduction of error on global r0 in reconstruction
; if err_r0 EQ 1, then no error. for a 10% error, use err_R0 = 1.1
err_R0 = 1.

theta = [[0., 0.]]/60.      ; Optimization direction [arcmin]
IF NOT KEYWORD_SET(beta_tab) THEN BEGIN
      print, 'beta_tab not set, choosing default of [0,0]'
      beta_tab = [[0.0, 0.0]]/60.0;, $    ; Performance direction [arcmin]          <----------- Check this value
      ;            [10., 0.], $
      ;            [20., 0.], $
      ;            [30., 0.]]/60;, $
      ;            [0., 10.], $
      ;            [0., 20.], $
      ;            [0., 30.]]/60 
END  


; Reconstruction method
LSE = 0B                      ; LSE = 1B for LSE reconstruction

; TEMPORAL Error
tempo = 1B
IF keyword_set(tempo) THEN $
    print,  'TEMPORAL ERROR: YES' $
ELSE print,  'TEMPORAL ERROR: NO'


; define keywords for ALAISING
; 0 for no error in reconstructor, no error in residual phase
; 1 for no error in reconstructor,    error in residual phase
; 2 for    error in reconstructor,    error in residual phase
aliasing_kw = 0

IF aliasing_kw EQ 0 THEN $
   print, 'NO ALIASING'
IF aliasing_kw EQ 1 THEN $
   print, 'ALIASING ERROR & NO ALIASING IN RECONSTRUCTION'
IF aliasing_kw EQ 2 THEN $
   print, 'ALIASING ERROR & ALIASING IN RECONSTRUCTION'

; define keywords for FITTING
; 0 for no error in reconstructor, no error in residual phase
; 1 for    error in reconstructor,    error in residual phase
fitting_kw = 1

IF fitting_kw EQ 0 THEN $
   print, 'NO DM FITTING ERROR & PURE TOMO CASE (by def: Nb_DM=Nb_Reconstructed Turb_Layer)'
IF fitting_kw EQ 1 THEN $
   print, 'DM FITTING ERROR (Nb_DM!=Nb_Reconstructed Turb_Layer)'

; Case PURE TOMOGRAPHIC reconstruction
tomo = 0B                     ; Ideal reconstruction h_dm = h_recons
IF keyword_set(tomo) THEN begin
      h_dm = h_recons
      print, 'PURE tomographic reconstruction, h_dm = h_recons'
ENDIF


; Output
dir = 'C:\Users\jfsauvage\Desktop\HARMONI\ATC-E2E_FOURRIER\FOURIER\SimuFourier\'
;get_date, date,/time
;date = STRJOIN(STRSPLIT(date, ':', /extra), '.')
suffix = turb_case+'_'+$
;         strn(D, format='(I4.0)')+'m_'+$
;         strn(dimpup, format='(I4.0)')+'px'+strn(dim, format='(I4.0)')+'px_'+$
;         strn(nb_act, format='(I2.0)')+'x'+strn(nb_act, format='(I2.0)')+'_'+$
;         strn(fech, format='(I4.0)')+'Hz_'+$
;         strn(td*1000., format='(I4.0)')+'ms_'+$
;         strn(seeing0, format='(f5.2)')+'arcsec_'+$
         strn(anglezenith, format='(I2.0)')+'deg_'+$
         strn(L0, format='(I4.0)')+'m_'+$
         strn(varnoise_wfs[0] /(lambda_wfs/lambda_imaging)^2, format='(F7.2)')+'rad2_'+$
         strn(lambda_imaging*1e6, format='(F5.2)')+'um'
;suffix = turb_case+'_'+$
;         strn(dim, format='(I4.0)')+'px_'+$
;         strn(anglezenith, format='(I2.0)')+'deg_'+$
;         strn(varnoise_wfs /(lambda_wfs/lambda_imaging)^2, format='(F7.2)')+'rad2_'+$
;         strn(lambda_imaging*1e6, format='(F5.2)')+'um'
;bettt = strn(beta_tab[0]*60, format='(f5.2)')+'arcsec'
;file = dir+'SCAO_HARMONI_'+suffix+'_'+bettt+'.dat' ; save file string
file = dir+'PY-WFS_SCAO_HARMONI_'+suffix+'.dat' ; save file string

sim_params = create_struct('anglezenith', anglezenith, 'sep', sep, 'nb_gs', nb_gs, $
                           'LGS_H', LGS_H, 'varnoise_wfs', varnoise_wfs, 'lambda_seeing', lambda_seeing, $
                           'lambda_imaging', lambda_imaging, 'L0', L0, 'seeing0', seeing0, 'cn2true', cn2true, $
                           'vent', vent, 'arg_v', arg_v, 'h', h, 'D', D, 'dimpup', dimpup, $
                           'dim', dim, 'nb_act', nb_act, 'pitch_m', pitch_m,$
                           'fech', fech, 'td', td, 'h_dm', h_dm, 'pitchs_dm', pitchs_dm, $
                           'alpha', alpha, 'tempo', tempo, $
                           'fitting_kw', fitting_kw, $
                           'recons_h', recons_h, 'cn2recons', cn2recons, $
                           'err_R0', err_R0, 'theta', theta, $
                           'beta_tab', beta_tab, 'LSE', LSE, $
                           'aliasing_kw', aliasing_kw, $
                           'tomo', tomo, 'ngs', ngs, 'file', file)

; -- Display simulation main geometry / parameters --
display_simu, sim_params
display_geo, alpha, beta_tab, theta, nb_gs = nb_gs, window_size = [600,600]    
display_turb, Cn2true, h, Cn2recons, recons_h, vent, arg_v, h_dm, window_number = 2, window_size = [650,600], /inv_color


; launch simulation
fourier_ao, sim_params

END

