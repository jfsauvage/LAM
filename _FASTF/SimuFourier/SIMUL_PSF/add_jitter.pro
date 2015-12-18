FUNCTION add_jitter, jitter_mas, sim_params, maskpupil, $
                     jitter_value_mas = jitter_value_mas, $
                     phase = phase, $
                     psf = psf, $
                     verbose = verbose

;.run ['D:\Simulation Codes\IDL\ELTF2\FOURIER_AO\carre'] 

; This function adds residual jitter to a phase map (if KW PHASE present), or
; to a PSF (if KW PSF present). 
; Adding jitter to a phase is done by adding randomn tip-tilt phase screen of
; the correct amoung (for short exposure mode), 
; or convolving by a gaussian (for long exposure mode)

; JITTER_mas            ; jitter RMS value [mas]
; JITTER_value_mas      ; jitter value [mas]
; SIM_PARAMS            : structure containing all simulation parameters
; PHASE                 ; phase maps for each performance direction 
;                       ; with residual turbulence, 
;                       ; TTF filtered (if asked), and static contribution (if asked) [nm]
; psf                   ; psf for each performance direction 

rad2arcsec = 1.*180./!dpi*3600.
nb_dirperf = n_elements(sim_params.beta_tab)/2.

; diffraction size in mas
diff_mas = 1000. * rad2arcsec * sim_params.lambda_imaging / sim_params.D

; diffraction size in pix
diff_px = sim_params.dim/sim_params.dimpup

; conversion of jitter in diffraction size (focal plane)
jitter_diff = jitter_mas /diff_mas

; conversion of jitter : pixel to nm RMS of TT coefficients
; 1 lambda RMS of Tip coefficient is 4 lambda PV, and is 4 diffraction shift
; in focal plane 
jitter_nm_zern = 1e9 * 0.25 * sim_params.lambda_imaging * jitter_diff

IF verbose EQ 2 THEN begin
   print, 'DIFFRACTION SIZE : ' + strc(diff_mas) + ' [mas]'
   print, 'IMAGE SAMPLING : ' + strc(diff_px) + ' [px]'
   print, 'FOCAL PLANE JITTER : ' + strc(jitter_diff[0]) + ' [lambda]'
   IF jitter_diff[0] GT 1 THEN print, 'Poor system performance...Jitter too high'
ENDIF



; case PHASE simulation (short exposure)
IF keyword_set(phase) THEN begin
   zern = calc_mode_zernike(nbmode = 2, pupdiam = sim_params.dimpup, largeur = $
                           sim_params.dim, /silent)
   IF NOT(keyword_set(jitter_value_mas)) THEN  jitter_shortexposure = $
       randomn(seed, 2) ELSE $
           jitter_shortexposure = 1e9 * 0.25 * jitter_value_mas / diff_mas 
   phase_jitter = maskpupil * carre(zern ## ((jitter_nm_zern * jitter_shortexposure)))
   FOR ind_perf = 0, nb_dirperf-1 DO $
       phase[*, *, ind_perf] += phase_jitter
   return, phase
ENDIF

; case PSF simulation (short exposure)
IF keyword_set(psf) THEN begin
   IF keyword_set(verbose) THEN print, 'FOURIER JITTER ON LONG EXPOSURE'
   jitter_px = jitter_diff * diff_px
   
   ; generation of gaussian filter in fourier space
   gauss2d_xy, sim_params.dim, sim_params.dim, $
               (sim_params.dim)/2., (sim_params.dim)/2., $
               8. * alog(2.) / 2./ !pi / jitter_px*sim_params.dim, $
               fft_gauss_test
   
   ; computation of jittered PSF   
   FOR ind_dirperf = 0, nb_dirperf-1 DO $
       psf[*, *, ind_dirperf] = norme(real(fftshift(fftshift(psf[*, *, ind_dirperf]) * fft_gauss_test, +1)))
   
   return, psf
ENDIF



END
