PRO calc_w2_alias, f, fx, fy, arg_f, $
                   pitchs_wfs,pitchs_dm, $
                   Nb_gs,alpha,theta,varnoise_wfs_xy, $
                   DSP_tab_recons, h_recons, h_dm, $
                   seuil, condmax, Wflag, $
                   LSE = LSE, verbose = verbose, $
                   wtomo = wtomo, MPalpha_xy = MPalpha_xy, MPalpha_t_xy = MPalpha_t_xy, $
                   cb_alias = cb_alias, Cb_noisetot = Cb_noisetot


; this function produces the second form of tomographic reconstructor WTOMO


; wtomo
; MPalpha
; t_MPalpha

;----------------------------------------------------------------------------------
; W2 = Cphi_recons#MPalphat(MPalpha#Cphi_recons#MPalphat + Cb_recons)^-1
;      -------------------- |       --------------------           |   |
;                res_tmp    |            res_tmp                   |   |
;                |          --------------------                   |   |
;                |                model_r                          |   |
;                |                ----------------------------------   |
;                |                               MAP                   |
;                |                               -----------------------  
;                |                                        INV
;                --------------------------------------------
;                                    WTOMO
;----------------------------------------------------------------------------------

;----------------------------------------------------------------------------------
; W2 allows to force noise at 0 during reconstruction. 
; Interesting to study tomo rec error only, without noise.
; But is not the same than WLSE, without a priori   
;----------------------------------------------------------------------------------

s = (size(f))(1)

;----------------------------------------------------------------------------------
; Computation of Cphi_recons
;----------------------------------------------------------------------------------
Cphi_recons = fltarr(n_elements(h_recons)*s, n_elements(h_recons)*s)
FOR i = 0, n_elements(h_recons)-1 DO $
    Cphi_recons(s*i:s*(i+1)-1, s*i:s*(i+1)-1) = DSP_tab_recons(s*i:s*(i+1)-1, *)

   
;----------------------------------------------------------------------------------
; term res_tmp : Cphi_recons#MPalpha_t first
;----------------------------------------------------------------------------------
   res_tmp2 =  bloc_product(Cphi_recons,  MPalpha_t_xy, s)
   undefine, Cphi_recons
   
   
   
;----------------------------------------------------------------------------------
; term model_r2 
;----------------------------------------------------------------------------------
   model_r2 = bloc_product(MPalpha_xy, res_tmp2, s)
   
   
   
;----------------------------------------------------------------------------------
; matrix to invert : model_r2 + Cb_noisetot
;----------------------------------------------------------------------------------
   to_inv =  model_r2 + Cb_noisetot
   undefine, model_r2
   undefine, Cb_noisetot
   
   
   
;----------------------------------------------------------------------------------
; MAP computation : invert the matrix TO_INV, freq by freq, case by case
;----------------------------------------------------------------------------------
   inv = to_inv*0.0
   
   tmp = complexarr(2*nb_gs, 2*nb_gs)
   FOR i = 0, s-1 DO BEGIN
      FOR j = 0, s-1 DO BEGIN
         
         FOR k = 0, 2*nb_gs-1 DO $
            FOR l = 0, 2*nb_gs-1 DO $
               tmp(k, l) = to_inv(i+k*s, j+l*s)
         
         IF total(tmp) NE 0 THEN BEGIN
            
; LA_TSVD : limits numerical pbs (compared to TSVD)   
            ; 2012 09 04 JFS add for case nb_gs=1
            IF n_elements(tmp) EQ 1 THEN tmp_inv = 1/tmp ELSE $
                la_tsvd, mat = tmp, inverse = tmp_inv, condmax = seuil, /SILENT
            IF i EQ 0. AND j EQ 0. THEN tmp_inv = tmp*0.
            
            FOR k = 0, 2*nb_gs-1 DO $
               FOR l = 0, 2*nb_gs-1 DO $
                  inv(i+k*s, j+l*s) = tmp_inv(k, l)
            
         ENDIF
         
      ENDFOR
   ENDFOR
   
   undefine, to_inv
   
;------------------------------------------------------------
; Reconstructor W2 = res_tmp2#inv
;------------------------------------------------------------
   wtomo = bloc_product(res_tmp2, inv, s)
   
   
   undefine, res_tmp2
   undefine, inv
      
END
