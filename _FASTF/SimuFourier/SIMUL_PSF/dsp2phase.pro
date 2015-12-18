PRO dsp2phase, dspall,  resall, f, sim_params, $
               seed = seed, $
               pupil_kw = pupil_kw, $
               maskpupil = maskpupil, $
               kw_static = kw_static, $
               ttf_kw = ttf_kw, $
               phase = phase, $
               static_phase = static_phase, $
               rms_static_phase = rms_static_phase, $
               verbose = verbose

; This procedure performs a phase simulation from DSP 
; The output is an array containing phase maps for each performance direction
; of DSPALL
; The phase map is residual turbulent, with eventual TTF filtering, with
; eventual static contribution

; INPUT ARGUMENTS
; RESALL       : fltarr(nb_dir_perf)  containing residual for each performance direction [nm]
; DSPALL       : fltarr(nb_f, nb_f, nb_dir_perf) containing 2D DSP for each perf direction
; F            : fltarr(nb_f, nb_f) containing 2D array of frequencies
; SIM_PARAMS   : structure containing all simulation parameters

; INPUT OPTIONNAL 
; seed                  ; use a given seed for random generation
; pupil_kw              ; sets a pupil mask in simulation (1B) or not (0B)
; kw_static             ; adds a static phase in simulation (1B) or not
;                       ; (0B)
; TTF_kw                ; KW for filtering TT (1B) or TTF (2B) or nothing (0B)
                        ; in the turbulent phase

; OUTPUTS
; phase                 ; phase maps for each performance direction 
;                       ; with residual turbulence, 
;                       ; TTF filtered (if asked), and static contribution (if asked) [nm]
; static_phase          ; sole static contribution in phase [nm]
; maskpupil             ; 2D pupil mask

; gets basic numbers from simulation structure 
nb_dirperf = n_elements(sim_params.beta_tab)/2.
pupil_size_px = sim_params.dimpup
array_size = sim_params.dim

; sets default parameters
IF kw_static EQ 0 THEN rms_static_phase = 0.  ; [nm]
IF keyword_set(verbose) THEN print, 'Static phase rms [nm] = '+strc(rms_static_phase)

IF keyword_set(verbose) THEN begin
   IF ttf_kw EQ 0 THEN print,  'FILTERING TT : NO'
   IF ttf_kw EQ 1 OR ttf_kw EQ 2 THEN print,  'FILTERING TT : YES'
   IF ttf_kw EQ 2 THEN print,  'FILTERING FOCUS : YES'
ENDIF

; sets phase array (output)
phase = fltarr(sim_params.dim, sim_params.dim, nb_dirperf)

; sets pupil mask
; if KW = 0, mask is 1 everywhere
IF pupil_kw EQ 0 THEN BEGIN
   IF keyword_set(verbose) THEN print, 'PUPIL : NO'
   maskpupil = fltarr(sim_params.dim, sim_params.dim)+1.
ENDIF

; if KW = 1, mask is 1 inside a round pupil
IF pupil_kw EQ 1 THEN begin
   IF keyword_set(verbose) THEN print, 'PUPIL : YES'
   polaire2, rt = pupil_size_px/2., largeur = array_size, masque = maskpupil, oc = $
             0, prolixe = 0
   maskpupil = float(maskpupil)
ENDIF

; sets static phase map
IF keyword_set(verbose) THEN print, 'Static phase generation'
static_phase = calc_static_phase(rms_static_phase, f, sim_params, seed, maskpupil)
static_phase *= maskpupil

; loop on performance directions
FOR ind_dir = 0, nb_dirperf-1 DO begin
   IF keyword_set(verbose) THEN print, 'Generate phase for direction '+strc(ind_dir)
   phase[*, *, ind_dir] = calc_turb_phase(dspall[*, *, ind_dir], $
                                          resall[ind_dir], $
                                          sim_params, maskpupil) 
   ; sets pupil on turbulent phase
   phase[*, *, ind_dir] *= maskpupil
   
   ; filters TTF
   IF TTF_kw EQ 1 THEN phase[*, *, ind_dir] = $
       ttf_filtering(phase[*, *, ind_dir], sim_params, maskpupil, /tt)
   IF TTF_kw EQ 2 THEN phase[*, *, ind_dir] = $
       ttf_filtering(phase[*, *, ind_dir], sim_params, maskpupil, /tt, /focus)
   
   IF TTF_kw EQ 3 THEN phase[*, *, ind_dir] = $
       ttf_filtering(phase[*, *, ind_dir], sim_params, maskpupil, /tt, /inverse)
   IF TTF_kw EQ 4 THEN phase[*, *, ind_dir] = $
       ttf_filtering(phase[*, *, ind_dir], sim_params, maskpupil, /tt, /focus, $
                    /inverse)
   
   ; adds static phase (that may contain TTF)
   phase[*, *, ind_dir] += static_phase
   ; aff, [static_phase, phase[*, *, ind_dir]], /nosample
ENDFOR

END
