PRO main

; Simulation geometry
anglezenith = 0               ; observation angle relative to zenith 
posfov = 0.                   ; unused
deltah = 0.                   ; unused
sep = 0.01                    ; arcmin Careful > 1.6 arcmin to ensure NO
                              ; unseen area
nb_gs = 1.
LGS_H = 0; 9e4                ; altitude for LGS. 

; Simulation photometry
lambda_wfs = 0.589e-6         ; wfs wavelength [m]
lambda_seeing = 0.5e-6        ; r0 wavelength  [m]
lambda_imaging =  1.65e-6     ; imaging wavelength [m]
varnoise_wfs = 1e-3 * (fltarr(nb_gs)+1)*(lambda_wfs/lambda_imaging)^2 ;analyse noise in rad^2

; Simulation turbulence definition
L0     = 25.                  ; outer scale    [m]; TEM
seeing0 = 0.85                ; zenithal seeing (arcsec)
cn2true =  [1.0]              ; cn2 profile 
vent = [10]                   ; wind speed modulous [m/s]
arg_v  = [0]                  ; wind dir.  [rad]
h = [0]                       ; layers altitude [m]

; Simulation parameters
D = 42.                       ; Telescope diameter [m]
dimpup = 168.                 ; Pupil size [px]
dim = 512.                    ; Screen size [px] 

; AO system definition
nb_act = 84.                  ; number of act
pitch_m =  float(D)/(float(nb_act)-1)       ; WFS pitch [m]
fech = 500.                   ; AO temporal sampling frequency [Hz]
td = 3*1.e-3                  ; loop latency [s]
h_dm = [100]                  ; DM altitude
pitchs_dm = [pitch_m]         ; DM pitch [m]
aa = 2.*!pi/Nb_gs             ; etoiles en cercle
alpha  = transpose([[cos(indgen(nb_gs)*aa)],  [sin(indgen(nb_gs)*aa)]] $
                   *sep/2.)   ; analyse directions [arcmin]
; NGS Geometry [arcmin]
ngs = [[0., 20.],  $
      [-10., -35.]]/60  

tempo = 1B                    ; temporal error 

; Reconstruction turbulence a priori 
recons_h = [0]                ; supposed altitudes for estimation 
cn2recons = [1.0]             ; supposed CN2 for estimation
; Introduction of error on global r0 in reconstruction
; if err_r0 EQ 1, then no error. for a 10% error, use err_R0 = 1.1
err_R0 = 1.

theta = [[0, 0]]/60.          ; optimization direction [arcmin]
beta_tab = [[0., 0.], [30., 30.]]/60      ; performance direction [arcmin]

; Reconstruction method
LSE = 0B                      ; LSE = 1B for LSE reconstruction
tomo = 0B                     ; Ideal reconstruction h_dm = h_recons

; Output
file = 'file_save_oa_test.dat' ; save file string

IF keyword_set(tempo) THEN $
    print,  'TEMPORAL ERROR: YES' $
ELSE print,  'TEMPORAL ERROR: NO'


; define keywords for aliasing
; 0 for no error in reconstructor, no error in residual phase
; 1 for no error in reconstructor,    error in residual phase
; 2 for    error in reconstructor,    error in residual phase
aliasing_kw = 1

IF aliasing_kw EQ 0 THEN $
   print, 'NO ALIASING'
IF aliasing_kw EQ 1 THEN $
   print, 'ALIASING ERROR & NO ALIASING IN RECONSTRUCTION'
IF aliasing_kw EQ 2 THEN $
   print, 'ALIASING ERROR & ALIASING IN RECONSTRUCTION'

; define keywords for fitting
; 0 for no error in reconstructor, no error in residual phase
; 1 for    error in reconstructor,    error in residual phase
fitting_kw = 1

IF fitting_kw EQ 0 THEN $
   print, 'NO DM FITTING ERROR & PURE TOMO CASE (by def: Nb_DM=Nb_Reconstructed Turb_Layer)'
IF fitting_kw EQ 1 THEN $
   print, 'DM FITTING ERROR (Nb_DM!=Nb_Reconstructed Turb_Layer)'

; case pure tomogrpahic reconstruction
IF keyword_set(tomo) THEN begin
      h_dm = h_recons
      print, 'PURE tomographic reconstruction, h_dm = h_recons'
ENDIF

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

; display simulation main geometry / parameters
; display_simu, sim_params
; display_geo, alpha, beta_tab, theta, ngs = ngs, window_number = 1, window_size = [650,650], /inv_color       
; display_turb, Cn2true, h, Cn2recons, recons_h, vent, arg_v, h_dm, window_number = 2, window_size = [750,650], /inv_color

; launch simulation
fourier_ao, sim_params

END
