; Simulation geometry
anglezenith = 0               ; observation angle relative to zenith 
posfov = 0.                   ; unused
deltah = 0.                   ; unused
sep = 0.01                     ; arcmin Careful > 1.6 arcmin to ensure NO
                              ; unseen area
nb_gs = 1.
LGS_H = 0; 9e4                   ; altitude for LGS. 

; Simulation photometry
signoise_wfs = (fltarr(nb_gs)+1)*(0.589/1.65)^2 ;analyse noise in rad^2
lambda_seeing = 0.5e-6        ; r0 wavelength  [m]
lambda =  1.65e-6             ; imaging wavelength [m]

; Simulation turbulence definition
L0     = 25.                  ; outer scale    [m]; TEM
seeing0 = 0.85                ; zenithal seeing (arcsec)
cn2true =  [1.0]         ; cn2 profile 
vent = [10]               ; wind speed modulous [m/s]
arg_v  = [0]          ; wind dir.  [rad]
h = [0]                 ; layers altitude [m]

; Simulation parameters
D = 42.                         ; Telescope diameter [m]
dimpup = 168.                 ; Pupil size [px]
dim = 512.                    ; Screen size [px] 

; AO system definition
nb_act = 84.                 ; number of act
pitch_m =  float(D)/(float(nb_act)-1)       ; WFS pitch [m]
fech = 500.                   ; AO temporal sampling frequency [Hz]
td = 3*1.e-3                  ; loop latency [s]
h_dm = [0]                    ; DM altitude
pitchs_dm = [pitch_m]         ; DM pitch [m]
aa = 2.*!pi/Nb_gs             ; etoiles en cercle
alpha  = transpose([[cos(indgen(nb_gs)*aa)],  [sin(indgen(nb_gs)*aa)]] $
                   *sep/2.)   ; analyse directions [arcmin]
tempo = 1B                    ; temporal error 
fitting = 1B                  ; fitting error 

; Reconstruction turbulence a priori 
recons_h = [0]                ; supposed altitudes for estimation 
cn2recons = [1.0]               ; supposed CN2 for estimation

theta = [[0, 0]]/60. ; optimization direction [arcmin]
beta_tab = [[0., 0.]]/60       ; performance direction [arcmin]

; Reconstruction method
LSE = 0B                      ; LSE = 1B for LSE reconstruction
optim = 0                     ; ???
tomo = 0B                     ; Ideal reconstruction h_dm = h_recons

; Output
visu =  0
file = 'file_save_oa_test.dat' ; save file string

; define keywords for aliasing
; 0 for no error in reconstructor, no error in residual phase
; 1 for no error in reconstructor,    error in residual phase
; 2 for    error in reconstructor,    error in residual phase
aliasing_kw = 2

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
