PRO calc_w1_noalias, f, arg_f, $
                     pitchs_wfs,pitchs_dm, $
                     Nb_gs,alpha,theta,sigr, $
                     DSP_tab_recons, h_recons, h_dm, $
                     seuil, condmax, Wflag, $
                     LSE = LSE, verbose = verbose, $
                     wtomo = wtomo, MPalpha = MPalpha, t_MPalpha = t_MPalpha

; this function produces the first form of tomographic reconstructor WTOMO

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

s = (size(f))(1)

; Construction of Cb_inv (a priori on noise)
IF keyword_set(verbose) THEN print, 'W1 : Cb_inv_recons'
Cb_inv_recons = complexarr(nb_gs*s, nb_gs*s)
FOR i = 0, nb_gs-1 DO Cb_inv_recons(s*i:s*(i+1)-1, s*i:s*(i+1)-1) = 1./sigr(i) 

; Cphi-1, a priori on turbulence layers, computed from DSP_tab_recons
IF keyword_set(verbose) THEN    print, 'W1 : Cphi_inv_recons'
Cphi_inv_recons = fltarr(n_elements(h_recons)*s, n_elements(h_recons)*s)
FOR i = 0, n_elements(h_recons)-1 DO $
    Cphi_inv_recons(s*i:s*(i+1)-1, s*i:s*(i+1)-1) = 1./DSP_tab_recons(s*i:s*(i+1)-1, *)

; Filtering of piston in reconstruction : 
Cphi_inv_recons(0) = 0.

IF keyword_set(LSE) THEN Cphi_inv_recons *= 0.




;------------------------------------------------------------------------------------
; compute res_tmp term : t_MPalpha#Cb_inv first
;------------------------------------------------------------------------------------
IF keyword_set(verbose) THEN    print, 'W1 : t_MPalpha#Cb_inv first'
res_tmp = bloc_product(t_MPalpha, Cb_inv_recons, s)

undefine, cb_inv_recons
;------------------------------------------------------------------------------------




;------------------------------------------------------------------------------------
; compute mode_r term : t_MPalpha#Cb_inv#MPalpha  :
;------------------------------------------------------------------------------------
IF keyword_set(verbose) THEN    print, 'W1 : t_MPalpha#Cb_inv#MPalpha'
model_r = bloc_product(res_tmp, MPalpha, s)

undefine, MPalpha
undefine, t_MPalpha
;------------------------------------------------------------------------------------




;------------------------------------------------------------------------------------
; compute MAP term
;------------------------------------------------------------------------------------
; matrix to be inversed :
; size (NB_FREQ x NB_LAYER) x (NB_FREQ x NB_LAYER)
MAP = model_r + Cphi_inv_recons 
undefine, cphi_inv_recons
undefine, model_r
;------------------------------------------------------------------------------------




;------------------------------------------------------------------------------------
; compute INV term
;------------------------------------------------------------------------------------
; Inversion of MAP matrix frequency by frequency
IF keyword_set(verbose) THEN    print, 'W1 : Inversion of MAP matrix frequency by frequency'
inv= calc_inv_map(MAP, s, seuil, nb_gs, h_recons, h_dm, WFlag)
undefine, MAP
;------------------------------------------------------------------------------------





;------------------------------------------------------------------------------------
; Last step W1 = inv#res_tmp
;------------------------------------------------------------------------------------
IF keyword_set(verbose) THEN print, 'W1 : W1 = inv#res_tmp'
wtomo = bloc_product(inv, res_tmp, s)
;------------------------------------------------------------------------------------




;------------------------------------------------------------------------------------
; delete variables
;----------------------------------------------------------------------------------
undefine, res_tmp
undefine, inv
;----------------------------------------------------------------------------------

END

