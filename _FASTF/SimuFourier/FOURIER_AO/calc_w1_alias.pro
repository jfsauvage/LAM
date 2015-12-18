PRO calc_w1_alias, f, fx, fy, arg_f, $
                   pitchs_wfs,pitchs_dm, $
                   Nb_gs,alpha,theta,varnoise_wfs_xy, $
                   DSP_tab_recons, h_recons, h_dm, $
                   seuil, condmax, Wflag, $
                   LSE = LSE, verbose = verbose, $
                   wtomo = wtomo, MPalpha_xy = MPalpha_xy, MPalpha_t_xy = MPalpha_t_xy, $
                   cb_alias = cb_alias, Cb_noisetot = Cb_noisetot

; this function produces the first form of tomographic reconstructor WTOMO,
; in the presence of aliasing

; wtomo
; MPalpha
; t_MPalpha

;------------------------------------------------------------------------------------
; W1 = ((MPalphat#Cb_inv_recons#MPalpha + Cphi_inv_recons)^-1)#MPalphat#Cb_inv_recons 
;        ----------------------       |                 |    | ----------------------
;                res_tmp              |                 |    |          res_tmp
;                ----------------------                 |    |              |
;                           model_r                     |    |              |
;                           -----------------------------    |              |
;                                       MAP                  |              |
;                                       ----------------------              |
;                                                INV                        |
;                                                ----------------------------
;                                                           WTOMO
;------------------------------------------------------------------------------------

nbre = nb_gs
f_x = fx
f_y = fy
s = (size(f))(1)


;------------------------------------------------------------------------------------
; if LSE then simple Cb inverse, WFS noise only (no aliasing)
;------------------------------------------------------------------------------------
IF keyword_set(LSE) THEN BEGIN
;noise only  Cb and NOT the aliasing noise:
   
   Cb_inv_tot = fltarr(2*nb_gs*s, 2*nb_gs*s)
   FOR i = 0, 2*nb_gs-1 DO $
       Cb_inv_tot(s*i:s*(i+1)-1, s*i:s*(i+1)-1) = 1./varnoise_wfs_xy(i)
   
ENDIF ELSE BEGIN
   
   
      
;------------------------------------------------------------------------------------
;Inversion of Cb_alias + Cb_recons
;------------------------------------------------------------------------------------
;freq by freq and case by case inversion
      cb_inv_tot = fltarr(2*nb_gs*s, 2*nb_gs*s)
      tmp = fltarr(2*nb_gs, 2*nb_gs)
      FOR i = 0, s-1 DO BEGIN
         FOR j = 0, s-1 DO BEGIN
            
            FOR k = 0, 2*nb_gs-1 DO $
               FOR l = 0, 2*nb_gs-1 DO $
                  tmp(k, l) = cb_noisetot(i+k*s, j+l*s)
               
            IF total(tmp) NE 0 THEN BEGIN
               tsvd, mat = tmp, inverse = tmp_inv, condmax = seuil
               
               IF i EQ 0. AND j EQ 0. THEN tmp_inv = tmp*0.
               
               FOR k = 0, 2*nb_gs-1 DO $
                  FOR l = 0, 2*nb_gs-1 DO $
                     cb_inv_tot(i+k*s, j+l*s) = tmp_inv(k, l)
            ENDIF
         ENDFOR
      ENDFOR
   ENDELSE
      
;------------------------------------------------------------------------------------
; Computation of Cphi_inv_recons
;------------------------------------------------------------------------------------
   Cphi_inv_recons = fltarr(n_elements(h_recons)*s, n_elements(h_recons)*s)
   FOR i = 0, n_elements(h_recons)-1 DO $
            Cphi_inv_recons(s*i:s*(i+1)-1, s*i:s*(i+1)-1) = 1./DSP_tab_recons(s*i:s*(i+1)-1, *)

   
   IF keyword_set(LSE) THEN BEGIN
      Cphi_inv_recons *= 0.
   ENDIF
   
   
;----------------------------------------------------------------------------------
; computation of res_tmp term : MPalphat#Cb_inv_recons
;----------------------------------------------------------------------------------
   res_tmp = MPalpha_t_xy*0.0
   
   FOR i = 0, nb_gs*2.-1 DO $
      FOR j = 0, n_elements(h_recons)-1 DO $
         FOR k = 0, nb_gs*2.-1 DO $
            res_tmp(i*s:(i+1)*s-1, j*s:(j+1)*s-1) += MPalpha_t_xy(k*s:(k+1)*s-1, j*s:(j+1)*s-1)*Cb_inv_tot(i*s:(i+1)*s-1,k*s:(k+1)*s-1 )
      
   undefine, cb_inv_tot
   
   
   
;----------------------------------------------------------------------------------
; Computation of model_r term : MPalphat#Cb_inv#MPalpha
;----------------------------------------------------------------------------------
   model_r = complexarr(n_elements(h_recons)*s, n_elements(h_recons)*s)
   
   FOR i = 0, n_elements(h_recons)-1 DO $
      FOR j = 0, n_elements(h_recons)-1 DO $
         FOR k = 0, nb_gs*2.-1 DO $
            model_r(i*s:(i+1)*s-1, j*s:(j+1)*s-1) += res_tmp(k*s:(k+1)*s-1, j*s:(j+1)*s-1)*MPalpha_xy(i*s:(i+1)*s-1,k*s:(k+1)*s-1 )
   
;the matrix to invert:
   MAP = model_r + Cphi_inv_recons 
   undefine, model_r
   undefine, cphi_inv_recons
;---------------------------------------------------------------------
;without a priori, it is a WLS
;---------------------------------------------------------------------
   
;invert the matrix
   inv = MAP*0.0
   
   tmp = complexarr(n_elements(h_recons), n_elements(h_recons))
   FOR i = 0, s-1 DO BEGIN
      FOR j = 0, s-1 DO BEGIN
         
         FOR k = 0, n_elements(h_recons)-1 DO BEGIN
            FOR l = 0, n_elements(h_recons)-1 DO BEGIN
               tmp(k, l) = MAP(i+k*s, j+l*s)
            ENDFOR
         ENDFOR
         
         IF total(tmp) NE 0 THEN BEGIN
            
            IF n_elements(h_recons) GT 1 THEN la_tsvd, mat = tmp, inverse = tmp_inv, condmax = seuil, /SILENT $
            ELSE tmp_inv = invert(tmp)
            
            IF i EQ 0. AND j EQ 0. THEN tmp_inv = tmp*0.
            
            FOR k = 0, n_elements(h_recons)-1 DO BEGIN
               FOR l = 0, n_elements(h_recons)-1 DO BEGIN
                  inv(i+k*s, j+l*s) = tmp_inv(k, l)
               ENDFOR
            ENDFOR
         ENDIF
         
      ENDFOR
   ENDFOR
   
;W1 = inv#res_tmp
   W1 = complexarr(2*nb_gs*s, n_elements(h_recons)*s)
   
   FOR i = 0, 2*nb_gs-1 DO BEGIN
      FOR j = 0, n_elements(h_recons)-1 DO BEGIN
         FOR k = 0, n_elements(h_recons)-1 DO BEGIN
            W1(i*s:(i+1)*s-1, j*s:(j+1)*s-1) += inv(k*s:(k+1)*s-1, j*s:(j+1)*s-1)*res_tmp(i*s:(i+1)*s-1,k*s:(k+1)*s-1 )
         ENDFOR
      ENDFOR
   ENDFOR
   undefine, MAP
   undefine, inv
   undefine, res_tmp
   

;------------------------------------------------------------------------------------
; returns tomographic reconstructor with aliasing into variable WTOMO
;------------------------------------------------------------------------------------
wtomo = w1

END

