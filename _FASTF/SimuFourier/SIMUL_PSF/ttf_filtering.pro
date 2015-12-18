FUNCTION ttf_filtering, phase, sim_params, maskpupil, tt = tt, focus = focus, $
                        inverse = inverse

; This function performs a Tip-Tilt or Tip-Tilt-Focus filtering in the phase
; map PHASE

; phase : 2D array containing phase map to filter
; SIM_PARAMS   : structure containing all simulation parameters
; maskpupil : 2D array containing the pupil map
; TT is the KW for Tip-Tilt filtering
; F is the KW for Tip-Tilt-Focus filtering
; inverse is the KW for keeping TT only, or TTF only

IF NOT(keyword_set(tt)) THEN tt = 0
IF NOT(keyword_set(focus)) THEN focus = 0

; gets phase array dimension
size_phase = sim_params.dim

; case for TT and F filtering (removing)
IF inverse EQ 0 THEN begin
   IF tt EQ 1 THEN begin
; generates tip tilt focus
      zern = calc_mode_zernike(nbmodes = 3, largeur = size_phase, pupdiam = $
                               sim_params.dimpup, /silent)
; filtering TT
      phase -=  carre(zern[0, *] * total(zern[0,*] * ligne(phase)/total(maskpupil)) $
                      +zern[1, *] * total(zern[1,*] * ligne(phase)/total(maskpupil)))
      IF focus EQ 1 THEN BEGIN
; filtering focus      
         phase -=  carre(zern[2, *] * total(zern[2,*] * ligne(phase)/total(maskpupil)))
      ENDIF
   ENDIF
ENDIF

; case for TT and F keeping (removes the rest)
IF inverse EQ 1 THEN begin
   IF tt EQ 1 THEN begin
; generates tip tilt focus
      zern = calc_mode_zernike(nbmodes = 3, largeur = size_phase, pupdiam = $
                               sim_params.dimpup, /silent)
; keeping TT
      local_phase =  carre(zern[0, *] * total(zern[0,*] * ligne(phase)/total(maskpupil)) $
                           +zern[1, *] * total(zern[1,*] * ligne(phase)/total(maskpupil)))
      IF focus EQ 1 THEN BEGIN
; filtering focus      
         local_phase +=  carre(zern[2, *] * total(zern[2,*] * ligne(phase)/total(maskpupil)))
      ENDIF
      phase = local_phase
   ENDIF
ENDIF

return, phase

END
