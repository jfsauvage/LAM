PRO phase2psf, phase, maskpupil, sim_params, psf = psf

; This procedure generates PSF from phase measurement

; phase                 ; phase maps for each performance direction 
;                       ; with residual turbulence, 
;                       ; TTF filtered (if asked), and static contribution (if asked) [nm]
; maskpupil             ; 2D pupil mask
; SIM_PARAMS            ; structure containing all simulation parameters
; psf                   ; psf for each performance direction 

; gets basic numbers from simulation structure 
nb_dirperf = n_elements(sim_params.beta_tab)/2.
pupil_size_px = sim_params.dimpup
array_size = sim_params.dim

; sets phase array (output)
psf = fltarr(sim_params.dim, sim_params.dim, nb_dirperf)

; sets conversion unit
nm2rad = 1.e-9 * 2. * !pi/sim_params.lambda_imaging

; loop on performance directions
FOR ind_dir = 0, nb_dirperf-1 DO $
   psf[*, *, ind_dir] = $
       norme(abs2(fftshift($
       maskpupil * exp(nm2rad * !i * phase[*, *, ind_dir]))))


END
