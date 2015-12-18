; Simulation geometry
anglezenith = 0               ; observation angle relative to zenith 
posfov = 0                    ; unused
deltah = 0.                   ; unused
sep = .01                     ; arcmin Careful > 1.6 arcmin to ensure NO
                              ; unseen area
nb_gs = 6.
LGS_H = 9e4                   ; altitude for LGS. 

; Simulation photometry
signoise = [1]*(0.589/1.65)^2 ;.e-6
lambda_seeing = 0.5e-6        ; r0 wavelength  [m]
lambda =  1.65e-6             ; imaging wavelength [m]

; Simulation turbulence definition
L0     = 25.                  ; outer scale    [m]; TEM
seeing0 = 0.85                ; zenithal seeing (arcsec)
cn2true =  [0.5, 0.5]         ; cn2 profile 
vent = [10, 10]               ; wind speed modulous [m/s]
h = [0, 5e3]                 ; layers altitude [m]

; Simulation parameters
D = 42                        ; Telescope diameter [m]
dimpup = 168.                 ; Pupil size [px]
dim = 512.                    ; Screen size [px] 

; AO system definition
nb_act = 21.                 ; number of act
pitch_m =  D/(nb_act-1)       ; WFS pitch [m]
fech = 500.                   ; AO temporal sampling frequency [Hz]
td = 3*1.e-3                  ; loop latency [s]
h_dm = [0]                    ; DM altitude
pitchs_dm = [pitch_m]         ; DM pitch [m]
aa = 2.*!pi/Nb_gs             ; etoiles en cercle
alpha  = transpose([[cos(indgen(nb_gs)*aa)],  [sin(indgen(nb_gs)*aa)]] $
                   *sep/2.)   ; analyse directions [arcmin]
tempo = 1B
fitting = 1B                  ; fitting

; Reconstruction turbulence a priori 
recons_h = [0]                ; supposed altitudes for estimation 
cn2recons = [1]               ; supposed CN2 for estimation

theta = [[0, 0], [0, 0]]/60. ; optimization direction [arcmin]
beta_tab = [[0., 0.],  $
            [0., 10.],  $
            [0., 20.],  $
            [0., 30.],  $
            [0., 40.],  $
            [0., 50.],  $
            [0., 60.]]/60       ; performance direction [arcmin]

; Reconstruction method
alias = 'corr'                  ; Keyword for aliasing
optim = 0                     ; ???
tomo = 0B                     ; Ideal reconstruction h_dm = h_recons


; Output
visu =  1
file = 'file_save_oa_5km.dat'        ; save file string



