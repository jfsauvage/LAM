; This function computes the residual phase DSP, for all performance
; directions. 
; The DSP is computed following formulae of B. Neichel thesis
; 7.22 pp 182

FUNCTION calc_dspresidual_noalias, f, arg_f, $
                                   pitchs_wfs, pitchs_dm, $
                                   Nb_gs,alpha,beta,$
                                   sigv, DSP_tab_vrai,  h_vrai, h_dm, $
                                   Wmap, $
                                   tempo = tempo, td, ti, Wind, $
                                   fitting_kw = fitting_kw,  $
                                   err_recons = err_recons, err_noise = err_noise

;f          = spatial frequencies array
;arg_f      = freq argument
;pitchs_wfs = WFS pitches
;Nb_gs      = number of GS
;alpha      = position of GS in the field (cartesian (x,y) and arcmin)
;Beta       = position for performance evaluation (cartesian (x,y) and arcmin)
;sigv       = true noise used for each GS
;DSP_tab_vrai = turbulent layers DSP
;h_vrai     = true profile altitudes
;h_dm       = DM altitudes
;WMAP       = global  reconstructor = Popt # Wtomo
;td         = delay
;ti         = integration times of WFS
;Wind       = wind array

; conversion from angular radians to arcsec, famous 206265.
rad2arcsec = 1.*180./!dpi*3600.

nbre = nb_gs
f_x = f*cos(arg_f)
f_y = f*sin(arg_f)
s = (size(f))(1)

IF keyword_set(tempo) THEN BEGIN 
   wind = wind
   ti = ti
   td = td
ENDIF ELSE BEGIN
   wind = fltarr(2, n_elements(h_vrai))
   ti = fltarr(nb_gs)
   td = 0.
ENDELSE



;-----------------------------------------------------------
; Complete formulae, case without aliasing (only two terms)
; PSD_res =   (PbetaL - PbetaDM # WMAP # M.PalphaL) # C_phi #...
;             ...(PbetaL - PbetaDM # WMAP # M.PalphaL)^T
;           + (PbetaDM#WMAP) C_b (PbetaDM#WMAP)^T
;-----------------------------------------------------------


; 1- (PbetaL - PbetaDM#WMAP#M.PalphaL)
; WFS matrix 
; M matrix
; PbetaL matrix; 
; PbetaDM matrix
; PbetaDM#WMAP 
; Then proj_tmp#Mv
; (PbetaL - PbetaDM#WMAP#M.PalphaL)^T
; Middle matrix C_phi
; First term of DSP : reconstruction error contribution
; 2- Second term of DSP : noise propagation through tomographic reconstructor
; first term (PbetaDM # Wmap)
; Then transpose (PbetaDM # Wmap)^T
; noise covariance matrix C_b
; Err_noise = proj_noise#Cb#proj_noise_conj
;--------------------------------------------------------------------------------------




;-----------------------------------------------------------
; 1- (PbetaL - PbetaDM#WMAP#M.PalphaL)
;-----------------------------------------------------------
; WFS matrix 
wfs = complexarr(nb_gs*s, s)
FOR j = 0, nb_gs-1 DO BEGIN
   ccc =  2*!pi*(f)*sinc(pitchs_wfs(j)*f_x)*sinc(pitchs_wfs(j)*f_y)*complex(0, 1)

;   ; test sans WFS : CCC=1
;   print, 'CALC_DSPRESIDUAL CCC = 1'
;   ccc = (2.*sinc(pitchs_wfs(j)*f_x)*sinc(pitchs_wfs(j)*f_y))*complex(0, 1)

   fc =  1./(2.*pitchs_wfs(j)) 
;f_ind =  where((f NE 0) and (f GT fc),  count)
   f_ind =  where((f NE 0) and (abs(f_x) GT fc) OR (abs(f_y) GT fc),  count)
   IF (count NE 0) THEN ccc(f_ind) =  0.
   wfs(j*s:(j+1)*s-1, *) = ccc
ENDFOR

; M matrix
Mv = complexarr(n_elements(h_vrai)*s, nb_gs*s)
FOR i = 0, n_elements(h_vrai)-1 DO BEGIN
   FOR j = 0, nb_gs-1 DO BEGIN
      ff_x = f_x*alpha(0, j)*h_vrai(i)*60./rad2arcsec
      ff_y = f_y*alpha(1, j)*h_vrai(i)*60./rad2arcsec
      www = sinc(wind(0, i)*ti(j)*f_x+wind(1, i)*ti(j)*f_y)
      Mv(i*s:(i+1)*s-1, j*s:(j+1)*s-1) = www*wfs(j*s:(j+1)*s-1, *)*exp(complex(0, 1)*2*(ff_x+ff_y)*!pi)
   ENDFOR
ENDFOR
undefine, wfs

; PbetaL matrix; 
; The phase screens have moved by DeltaTxV
; shift is actually backward shift
deltaT = (max(ti)+td)(0)
proj_beta = complexarr(s*n_elements(h_vrai), s)
FOR j = 0, n_elements(h_vrai)-1 DO BEGIN
   proj_beta(j*s:(j+1)*s-1, *) = exp(complex(0, 1)*2*!pi*(h_vrai(j)*60./rad2arcsec*$
                                                          (beta[0]*f_x+beta[1]*f_y)-(wind(0, j)*deltaT*f_x+wind(1, j)*deltaT*f_y)))
ENDFOR

; PbetaDM matrix
proj_betaDM = complexarr(s*n_elements(h_dm), s)
FOR j = 0, n_elements(h_dm)-1 DO BEGIN
   proj_betaDM(j*s:(j+1)*s-1, *) = exp(complex(0, 1)*2*!pi*h_dm(j)*60./rad2arcsec*(beta[0]*f_x+beta[1]*f_y))
ENDFOR

; PbetaDM#WMAP 
proj_tmp = complexarr(nb_gs*s, s)
FOR i = 0, nb_gs-1 DO BEGIN
   FOR k = 0, n_elements(h_dm)-1 DO BEGIN
      proj_tmp(i*s:(i+1)*s-1, *) += Proj_betaDM(k*s:(k+1)*s-1, *)*WMAP(i*s:(i+1)*s-1,k*s:(k+1)*s-1 )
   ENDFOR
ENDFOR

; Then proj_tmp#Mv
proj_tmp2 = complexarr(n_elements(h_vrai)*s, s)
FOR i = 0, n_elements(h_vrai)-1 DO BEGIN
   FOR k = 0, nb_gs-1 DO BEGIN
      proj_tmp2(i*s:(i+1)*s-1, *) += Proj_tmp(k*s:(k+1)*s-1, *)*Mv(i*s:(i+1)*s-1,k*s:(k+1)*s-1 )
   ENDFOR
ENDFOR

; first part of first term is now complete
; (PbetaL - PbetaDM#WMAP#M.PalphaL)
proj = proj_beta - proj_tmp2

undefine, Mv
undefine, proj_tmp
undefine, proj_tmp2
undefine, proj_beta

; (PbetaL - PbetaDM#WMAP#M.PalphaL)^T
; transpose is called proj_conj
proj_conj = complexarr(s, n_elements(h_vrai)*s)
FOR j = 0, n_elements(h_vrai)-1 DO BEGIN
   proj_conj(*, j*s:(j+1)*s-1) = conj(proj(j*s:(j+1)*s-1, *))
ENDFOR

; Middle matrix C_phi
Cphi_vrai = fltarr(n_elements(h_vrai)*s, n_elements(h_vrai)*s)

FOR i = 0, n_elements(h_vrai)-1 DO $
         Cphi_vrai(s*i:s*(i+1)-1, s*i:s*(i+1)-1) = DSP_tab_vrai(s*i:s*(i+1)-1, *)


;--------------------------------------------------------------------------------------
; First term of DSP : reconstruction error contribution
; Err_recons =   (PbetaL - PbetaDM # WMAP # M.PalphaL) # C_phi #...
;                    ...(PbetaL - PbetaDM # WMAP # M.PalphaL)^T
;--------------------------------------------------------------------------------------
inter = complexarr(s*n_elements(h_vrai), s)

FOR i = 0, n_elements(h_vrai)-1 DO BEGIN
   FOR j = 0, n_elements(h_vrai)-1 DO BEGIN
      inter(i*s:(i+1)*s-1, *) +=  proj(j*s:(j+1)*s-1, *)*Cphi_Vrai(i*s:(i+1)*s-1, j*s:(j+1)*s-1)
   ENDFOR
ENDFOR

Err_recons = complexarr(s, s)
FOR j = 0, n_elements(h_vrai)-1 DO BEGIN
   tmp = complexarr(s, s)
   tmp = inter(j*s:(j+1)*s-1, *)*proj_conj(*, j*s:(j+1)*s-1)
   Err_recons = Err_recons+tmp
ENDFOR

err_recons(0) = 0.
err_recons = float(err_recons)

undefine, cphi_vrai
undefine, proj
undefine, proj_conj
undefine, inter




;--------------------------------------------------------------------------------------
; 2- Second term of DSP : noise propagation through tomographic reconstructor
; (PbetaDM # Wmap) # Cb # (PbetaDM # Wmap)^T
;--------------------------------------------------------------------------------------

; first term (PbetaDM # Wmap)
proj_noise = complexarr(nb_gs*s, s)
FOR i = 0, nb_gs-1 DO BEGIN
   FOR k = 0, n_elements(h_dm)-1 DO BEGIN
      proj_noise(i*s:(i+1)*s-1, *) += Proj_betaDM(k*s:(k+1)*s-1, *)*WMAP(i*s:(i+1)*s-1, k*s:(k+1)*s-1 )
   ENDFOR
ENDFOR

; Then transpose (PbetaDM # Wmap)^T
proj_noise_conj = complexarr(s, nb_gs*s)
FOR j = 0, nb_gs-1 DO BEGIN
   proj_noise_conj(*, j*s:(j+1)*s-1) = conj(proj_noise(j*s:(j+1)*s-1, *))
ENDFOR

; noise covariance matrix 
Cb_vrai = complexarr(nb_gs*s, nb_gs*s)
FOR i = 0, nb_gs-1 DO $
    Cb_vrai(s*i:s*(i+1)-1, s*i:s*(i+1)-1) = sigv(i)


;--------------------------------------------------------------------------------------
; Second term : 
;Err_noise = proj_noise#Cb#proj_noise_conj
;--------------------------------------------------------------------------------------
inter = complexarr(s*nb_gs, s)

FOR i = 0, nb_gs-1 DO BEGIN
   FOR j = 0, nb_gs-1 DO BEGIN
      inter(i*s:(i+1)*s-1, *) +=  proj_noise(j*s:(j+1)*s-1, *)*Cb_Vrai(i*s:(i+1)*s-1, j*s:(j+1)*s-1)
   ENDFOR
ENDFOR

Err_noise = complexarr(s, s)
FOR j = 0, nb_gs-1 DO BEGIN
   tmp = complexarr(s, s)
   tmp = inter(j*s:(j+1)*s-1, *)*proj_noise_conj(*, j*s:(j+1)*s-1)
   Err_noise = Err_noise+tmp
ENDFOR

err_noise(0) = 0.
err_noise = float(err_noise)



;--------------------------------------------------------------------------------------
; complete DSP : summation of Err_recon and Err_noise
;--------------------------------------------------------------------------------------
DSP_res = float(Err_recons + Err_noise)
DSP_res(0, 0) = 0.


; if fitting error is desired (FITTING_KW = 1 or 2)
IF fitting_kw EQ 1 OR fitting_kw EQ 2 THEN return, dsp_res

; else then no error is done outside of cut-off frequency
out_m = f*0.0

fc =  max(1./(2.*pitchs_dm))
f_ind =  where((f NE 0) and (abs(f_x) LE fc) AND (abs(f_y) LE fc),  count)
IF (count NE 0) THEN out_m(f_ind) =  DSP_res(f_ind)
return, out_m

END
