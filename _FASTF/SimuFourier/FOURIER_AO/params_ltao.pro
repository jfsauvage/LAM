; Simulation geometry
anglezenith = 0               ; observation angle relative to zenith 
posfov      = 0               ; unused
deltah      = 0.              ; unused
sep         = 1.              ; arcmin Careful > 1.6 arcmin to ensure NO

nb_gs = 5.
LGS_H = 9e4                   ; altitude for LGS. 

; Simulation photometry
; creation of noise array (same noise value for each WFS)
signoise_wfs  = (fltarr(nb_gs)+1)*(0.589/1.65)^2 ;analyse noise in rad^2
lambda_seeing = 0.5e-6        ; r0 wavelength  [m]
lambda        =  1.65e-6             ; imaging wavelength [m]

; Simulation turbulence definition
L0      = 25.                  ; outer scale    [m]; TEM
seeing0 = 0.85                ; zenithal seeing (arcsec)
Cn2true = [52.24, 2.60, 4.44, 11.60, 9.89,  2.95,  5.98, 4.30, 6.00]/100.
vent    = [15.,  13.,  13., 9.,  9.,  15.,  25.,  40.,  21.]
arg_v   = fltarr(n_elements(vent))          ; wind dir.  [rad]
h       = [47., 140., 281., 562., 1125., 2250., 4500., 9000., 18000.]

; Simulation parameters
D      = 42                   ; Telescope diameter [m]
dimpup = 168/8.                 ; Pupil size [px]
dim    = 512/8.                 ; Screen size [px] 

; AO system definition
nb_act    = 84.               ; number of act
pitch_m   = D/(nb_act-1)     ; WFS pitch [m]
fech      = 500.              ; AO temporal sampling frequency [Hz]
td        = 3*1.e-3           ; loop latency [s]
h_dm      = [0, 5e3, 10e3]    ; DM altitude
pitchs_dm = [pitch_m, pitch_m, pitch_m] ; DM pitch [m]
aa        = 2.*!pi/Nb_gs      ; etoiles en cercle
alpha     = transpose([[cos(indgen(nb_gs)*aa)],  [sin(indgen(nb_gs)*aa)]] $
                      *sep/2.) ; analyse directions [arcmin]
tempo   = 0B                  ; temporal error 
fitting = 1B                  ; fitting error 

; Reconstruction turbulence a priori 
recons_h = h                  ; supposed altitudes for estimation 
cn2recons = cn2true           ; supposed CN2 for estimation

theta = [[0, 0]]/60. ; optimization direction [arcmin]
beta_tab = [[0., 0.], [0, 15], [-15, 0]]/60       ; performance direction [arcmin]

; Reconstruction method
LSE = 0B                      ; LSE = 1B for LSE reconstruction
alias = 'yes'                  ; Keyword for aliasing
tomo = 0B                     ; Ideal reconstruction h_dm = h_recons
gfit =  1B                    ; generalized fitting => Popt x Wtomo

; Output
visu =  1
file = 'file_save_ltao_test_lse.dat'        ; save file string

; fitting
; noise
; tempo
; aliasing 
; aliasing généralise

; test config
IF n_elements(h_dm) EQ 1 THEN BEGIN
   IF (nb_gs) EQ 1 THEN print, 'AO'
   IF n_elements(theta)/2. EQ 1. THEN print, 'LTAO'
   IF n_elements(theta)/2. GT 1. THEN print, 'GLAO'
ENDIF

IF n_elements(h_dm) GT 1 THEN print, 'MCAO'
