; Simulation geometry
anglezenith = 0               ; observation angle relative to zenith 
posfov = 0                    ; unused
deltah = 0.                   ; unused
sep = 0.4                     ; arcmin Careful > 1.6 arcmin to ensure NO
                              ; unseen area
nb_gs = 6.
LGS_H = 9e4                   ; altitude for LGS. 

; Simulation photometry
signoise = [1]*(0.589/1.65)^2 ;analyse noise [rad^2]
lambda_seeing = 0.5e-6        ; r0 wavelength  [m]
lambda =  1.65e-6             ; imaging wavelength [m]

; Simulation turbulence definition
L0     = 25.                  ; outer scale [m]
seeing0 = 0.85                ; zenithal seeing [arcsec]
Cn2true = [52.24, 2.60, 4.44, 11.60, 9.89,  2.95,  5.98, 4.30, 6.00]/100.
vent = [15.,  13.,  13., 9.,  9.,  15.,  25.,  40.,  21.]
arg_v  = fltarr(n_elements(vent))          ; wind dir.  [rad]
h = [47., 140., 281., 562., 1125., 2250., 4500., 9000., 18000.]
; cn2true =  [0.6, 0.2, 0.2]    ; cn2 profile 
; vent = [15, 10, 10]           ; wind speed modulous [m/s]
; arg_v  = [0, !pi/2., 0]          ; wind dir.  [rad]
; h = [0, 1e3, 5e3]                 ; layers altitude [m]

; Simulation parameters
D = 42                        ; Telescope diameter [m]
dimpup = 168.                 ; Pupil size [px]
dim = 512.                    ; Screen size [px] 

; AO system definition
nb_act = 84.                 ; number of act
pitch_m =  D/(nb_act-1)       ; WFS pitch [m]
fech = 500.                   ; AO temporal sampling frequency [Hz]
td = 3*1.e-3                  ; loop latency [s]
h_dm = [0]                    ; DM altitude
pitchs_dm = [pitch_m]         ; DM pitch [m]
aa = 2.*!pi/Nb_gs             ; etoiles en cercle
alpha  = transpose([[cos(indgen(nb_gs)*aa)],  [sin(indgen(nb_gs)*aa)]] $
                   *sep/2.)   ; analyse directions [arcmin]
tempo = 0B                    ; temporal error 
fitting = 1B                  ; fitting error 

; Reconstruction turbulence a priori 
recons_h = [0]                ; supposed altitudes for estimation 
cn2recons = [1]               ; supposed CN2 for estimation

theta = [[0, 0], [0, 0]]/60. ; optimization direction [arcmin]
beta_tab = [[00., 0.],  $
;             [10., 0.],  $
;             [20., 0.],  $
;             [30., 0.],  $
;             [40., 0.],  $
;             [50., 0.],  $
            [60., 0.]]/60       ; performance direction [arcmin]

; Reconstruction method
LSE = 0B                      ; LSE = 1B for LSE reconstruction
alias = 'no'                  ; Keyword for aliasing
optim = 0                     ; ???
tomo = 0B                     ; Ideal reconstruction h_dm = h_recons
gfit =  1B                    ; generalized fitting

condmax = 1e6                 ; CONDMAX used for inverse computation in Popt (done by LA_TSVD)


; fitting
; noise
; tempo
; aliasing 
; aliasing généralise



