PRO main_psf

; program simulating AO corrected PSF with different options
; uses the output of code FOURIER_AO

; opens FOURIER_AO output
; FOURIER_AO files contain
; RESALL       : fltarr(nb_dir_perf)  containing residual for each performance direction [nm]
; DSPALL       : fltarr(nb_f, nb_f, nb_dir_perf) containing 2D DSP for each perf direction
; F            : fltarr(nb_f, nb_f) containing frequency array
; DSP_TAB_VRAI : fltarr(nb_f, nb_f, nb_L) containing 2D DSP of turbulent layers
; SIM_PARAMS   : structure containing all simulation parameters
dir = 'C:\Users\jfsauvage\Desktop\HARMONI\ATC-E2E_FOURRIER\FOURIER\SimuFourier\FOURIER_AO\'
file = 'file_save_oa_test'
restore, file = dir + file + '.dat'

nb_dirperf = n_elements(sim_params.beta_tab)/2.
rad2arcsec = 1.*180./!dpi*3600.

; sets parameters for DSP2PHASE simulation 
seed = 1234L                  ; use a given seed for random generation of
                              ; static phase
pupil_kw = 1                  ; sets a round pupil mask in simulation (1) or
                              ; not (0), or give a user map (2), then the map
                              ; MASKPUPIL has to be provided by the user
kw_static = 1B                ; adds a static phase in simulation (1B) or not
                              ; (0B)
rms_static_phase = 300.       ; sets the RMS value of static phase [nm] to be
                              ; added 
TTF_kw = 0                    ; KW for filtering TT (1) or TTF (2) or nothing (0)
                              ; in the turbulent phase, also used for keeping
                              ; TT only (3) or TTF only (4)
short = 0                     ; KW for short (1) or long (0) exposure computation
verbose = 1B                  ; prints each step information

; sets parameters for ADD_JITTER simulation 
jitter_kw = 2                 ; additionnal jitter in long exposure images (2)
                              ; in turbulent phase (1) or not (0)
jitter_mas = [1.0, 30.0]      ; jitter value [X, Y] [mas]

; sets parameters for MAIN simulation 
number_exposure = 10000.         ; number of short exposures to be averaged

; checks pupil parameters if user map triggered
IF pupil_kw EQ 2 THEN BEGIN
   IF verbose EQ 1 THEN print, 'PUPIL USER MAP MODE REQUESTED'
   ; check array undefined
   IF NOT(keyword_set(maskpupil)) THEN BEGIN
      print, 'PUPIL USER MAP UNDEFINED --> STOP'
      return
   ENDIF
   ; check wrong array size
   IF n_elements(maskpupil) NE (sim_params.dim)^2 THEN begin
      print, 'USER PUPIL MAP UNCORRECT SIZE : ' + strc((size(maskpupil))[1]) $
             + ' // Requested ' + strc(sim_params.dim) + ' --> STOP'
      return
   ENDIF
   ; check wrong pupil size
   IF abs(total(maskpupil)/(!pi*(sim_params.dimpup/2.)^2)-1) GT 0.1 THEN begin
      print, 'USER PUPIL SEEMS UNCORRECT SIZE, PLEASE CHECK --> STOP'
      return
   ENDIF
ENDIF

; sets pupil mask
; if KW = 0, mask is 1 everywhere
IF pupil_kw EQ 0 THEN BEGIN
   IF keyword_set(verbose) THEN print, 'PUPIL : NO'
   maskpupil = fltarr(sim_params.dim, sim_params.dim)+1.
ENDIF

; if KW = 1, mask is 1 inside a round pupil
IF pupil_kw EQ 1 THEN begin
   IF keyword_set(verbose) THEN print, 'PUPIL : YES'
   polaire2, rt = sim_params.dimpup/2., largeur = sim_params.dim , masque = maskpupil, oc = $
             0, prolixe = 0
   maskpupil = float(maskpupil)
ENDIF

; loop on short exposures
IF short EQ 1 THEN begin
   FOR ind_exposure = 0, number_exposure-1 DO begin
                              ; creation of turbulent phase map
      dsp2phase, dspall, resall, f, sim_params, $
                 seed = seed, $
                 pupil_kw = pupil_kw, $
                 maskpupil = maskpupil, $
                 kw_static = kw_static, $
                 ttf_kw = ttf_kw, $
                 phase = phase, $
                 static_phase = static_phase, $
                 rms_static_phase = rms_static_phase, $
                 verbose = verbose
      
                              ; additionnal jitter on short exposure image
      IF jitter_kw EQ 1 THEN $
          phase = add_jitter(jitter_mas, sim_params, maskpupil, $
                             jitter_value = jitter_value, $
                             phase = phase, $
                             verbose = verbose)
      
      IF pupil_kw EQ 1 THEN begin
         
; creation of short exposure PSF
         phase2psf, phase, maskpupil, sim_params, psf = psf
         
; initialisation of long exposure
         IF ind_exposure EQ 0 THEN longexposure = psf * 0.
         
; integration for long exposure
         longexposure += psf
      ENDIF   
   ENDFOR
ENDIF

; computation of long exposure PSF directly
IF short EQ 0 THEN BEGIN
   IF kw_static EQ 0 THEN rms_static_phase = 0. ; [nm]
   static_phase = calc_static_phase(rms_static_phase, f, sim_params, seed, maskpupil)
   longexposure = dsp2psf(dspall, maskpupil, sim_params, $
                          static_phase = static_phase)
 ENDIF


; addition of jitter on long exposure image
IF pupil_kw EQ 1 OR pupil_kw EQ 2 THEN begin
   IF jitter_kw EQ 2 THEN $
       longexposure = add_jitter(jitter_mas, sim_params, maskpupil, $
                                 jitter_value = jitter_value, $
                                 psf = longexposure, $
                                 verbose = verbose)
   

; computation of encircled energy
   pix2rad = sim_params.lambda_imaging/sim_params.D $
             / (sim_params.dim/sim_params.dimpup )
   width = 80.*rad2arcsec * sim_params.lambda_imaging/sim_params.D
   ee = calc_ee(longexposure, sim_params, pix2rad = pix2rad, width = width, /square)
   
; save to file : same name, addition of PSF at the end of name
   save, file =  dir + file + '_PSF_LONG' + '.dat', longexposure, maskpupil, static_phase, $
         rms_static_phase, ttf_kw, number_exposure, ee
ENDIF

END
