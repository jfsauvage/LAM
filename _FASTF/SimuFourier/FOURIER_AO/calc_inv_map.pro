FUNCTION calc_inv_map, MAP, s, seuil, nb_gs, h_recons, h_dm, WFlag

; computes the inverse of MAP matrix
; This inversion is used to compute the term INV in calc_w1_* and calc_w2_*

; MAP
; s
; seuil
; nb_gs
; h_recons
; Wflag

; Inversion of MAP matrix frequency by frequency
IF keyword_set(verbose) THEN    print, 'W1 : Inversion of MAP matrix frequency by frequency'
inv = MAP*0.0

IF strmatch(WFLAG, 'W1') THEN size_MAP = n_elements(h_recons)
IF strmatch(WFLAG, 'W2') THEN size_MAP = nb_gs
IF strmatch(WFLAG, 'POPT') THEN size_MAP = n_elements(h_dm)

tmp = complexarr(size_map, size_map)
IF strmatch(WFLAG, 'POPT') THEN tmp = dcomplexarr(size_map, size_map)

FOR i = 0, s-1 DO BEGIN
   FOR j = 0, s-1 DO BEGIN
; extract MAP for frequency [i,j]
      FOR k = 0, size_MAP-1 DO $
          FOR l = 0, size_MAP-1 DO $
              tmp(k, l) = MAP(i+k*s, j+l*s)
; inversion of each sub matrix         
      IF total(tmp) EQ 0. THEN tmp_inv = tmp
      IF total(tmp) NE 0. THEN BEGIN
         IF size_MAP GT 1 THEN $
             la_tsvd, mat = tmp, inverse = tmp_inv, condmax = seuil, /SILENT $
         ELSE tmp_inv = invert(tmp)
; filtering of piston
         IF NOT(strmatch(WFLAG, 'POPT')) THEN IF i EQ 0. AND j EQ 0. THEN tmp_inv = tmp*0.
; fill INV matrix with result
         FOR k = 0, size_map-1 DO $
             FOR l = 0, size_map-1 DO $
                 inv(i+k*s, j+l*s) = tmp_inv(k, l)
      ENDIF
   ENDFOR
ENDFOR

return, inv

end
