PRO calc_w2_noalias, f, arg_f, $
                     pitchs_wfs,pitchs_dm, $
                     Nb_gs,alpha,theta,signoise_wfs, $
                     DSP_tab_recons, h_recons, h_dm, $
                     seuil, condmax, Wflag, $
                     LSE = LSE, verbose = verbose, $
                     wtomo = wtomo, MPalpha = MPalpha, t_MPalpha = t_MPalpha

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
; W2 allows to force noise at 0 during reconstruction. 
; Interesting to study tomo rec error only, without noise.
; But is not the same than WLSE, without a priori   
;----------------------------------------------------------------------------------

s = (size(f))(1)

;----------------------------------------------------------------------------------
; Computation of Cb_recons
;----------------------------------------------------------------------------------
IF keyword_set(verbose) THEN  print, 'W2 : Cb_recons'
Cb_recons = complexarr(nb_gs*s, nb_gs*s)
FOR i = 0, nb_gs-1 DO Cb_recons(s*i:s*(i+1)-1, s*i:s*(i+1)-1) = signoise_wfs(i)



;----------------------------------------------------------------------------------
; Computation of Cphi_recons
;----------------------------------------------------------------------------------
IF keyword_set(verbose) THEN  print, 'W2 : Cphi_recons'
Cphi_recons = fltarr(n_elements(h_recons)*s, n_elements(h_recons)*s)
FOR i = 0, n_elements(h_recons)-1 DO $
    Cphi_recons(s*i:s*(i+1)-1, s*i:s*(i+1)-1) = DSP_tab_recons(s*i:s*(i+1)-1, *)



;----------------------------------------------------------------------------------
; removes piston :
;----------------------------------------------------------------------------------
Cphi_recons(0) = 0.



;----------------------------------------------------------------------------------
; Cphi_recons#t_MPalpha first
;----------------------------------------------------------------------------------
IF keyword_set(verbose) THEN  print, 'W2 : Cphi_recons#t_MPalpha '
res_tmp = bloc_product(Cphi_recons, t_MPalpha, s)
undefine, t_MPalpha
undefine, Cphi_recons



;----------------------------------------------------------------------------------
; MPalpha#Cphi_recons#MPalphaT then :
;----------------------------------------------------------------------------------
IF keyword_set(verbose) THEN  print, 'W2 : MPalpha#Cphi_recons#MPalphaT'
model_r = bloc_product(MPalpha, res_tmp, s)
undefine, MPalpha



;----------------------------------------------------------------------------------
; Cb_recons :
;----------------------------------------------------------------------------------
MAP =  model_r + cb_recons

undefine, Cb_recons
undefine, model_r



;----------------------------------------------------------------------------------
; Inverse MAP
;----------------------------------------------------------------------------------
IF keyword_set(verbose) THEN  print, 'W2 : Inversion frequency by frequency'
inv = calc_inv_map(MAP, s, seuil, nb_gs, h_recons, h_dm, WFlag)

undefine, MAP



;----------------------------------------------------------------------------------
;Last step W2 = res_tmp#inv
;----------------------------------------------------------------------------------
IF keyword_set(verbose) THEN  print, 'W2 : W2 = res_tmp#inv'
WTOMO = bloc_product(res_tmp, inv, s)



;----------------------------------------------------------------------------------
undefine, inv
undefine, res_tmp
;----------------------------------------------------------------------------------

END
