; This program computes the reconstruction matrix WMAP.
; accounting for all reconstruction parameters
; WMAP = Popt ## Wtomo

FUNCTION calc_reconstructor_noalias, f, arg_f, $
                                     pitchs_wfs,pitchs_dm, $
                                     Nb_gs,alpha,theta,varnoise_wfs, $
                                     DSP_tab_recons, h_recons, h_dm, $
                                     seuil, condmax, Wflag, $
                                     fitting_kw = fitting_kw, $
                                     LSE = LSE, alias = alias

; f = spatial frequencies array
; arg_f = F phase
; Nb_gs = Guide star number
; alpha = Guide stars positions
; theta = Optimisation directions
; varnoise_wfs = A priori on noise associated to each GS
; DSP_Tab_recons = A priori, DSP on estimated turbulent layers
; h_recons = Altitudes of reconstructed layers
; h_dm : DM altitude (if not pure tomo)
; condmax : Max acceptable conditionning in inversion for POPT computation
; Wflag : Choice of reconstructor W1 or W2
; Keyword Tomo : Pure tomo
; Popt : output : optimal projector for MCAO, used later for Aliasing
; LSE : LSE instead of MAP

f_x = f*cos(arg_f)
f_y = f*sin(arg_f)

out = f*0.0
s = (size(f))(1)

; WFS used is Shack.
; Each WFS has its own cut off frequency
; Construction of WFS transfert function
wfs = complexarr(nb_gs*s, s)
FOR j = 0, nb_gs-1 DO BEGIN
   ccc =  2*!pi*(f)*sinc(pitchs_wfs(j)*f_x)*sinc(pitchs_wfs(j)*f_y)*complex(0, 1)
   ; essai : calcul de la fonction de transfert pyramide
;    tf_pyramid = (2.*sinc(pitchs_wfs(j)*f_x)*sinc(pitchs_wfs(j)*f_y))*complex(0, 1)

  ; test sans WFS : CCC=1
;   print, 'CALC_RECONSTRUCTOR CCC = PYR'
;   ccc = (2.*sinc(pitchs_wfs(j)*f_x)*sinc(pitchs_wfs(j)*f_y))*complex(0, 1)

   fc =  1./(2.*pitchs_wfs(j)) 
   f_ind =  where((f NE 0) and (abs(f_x) GE fc) OR (abs(f_y) GE fc),  count)
   IF (count NE 0) THEN ccc(f_ind) =  0.
   wfs(j*s:(j+1)*s-1, *) = ccc
ENDFOR

; NB : Here CCC is the transfert function of a SH, but something else could be
;      written here. Pyramid / Curvature / direct phase sensing (CCC=1)

;-----------------------------------------------------------
; Construction of WHAP = PoptWtomo
; Size : Ngs x Ndm
;-----------------------------------------------------------
; Starting with WTomo
; 2 writings, accessible with calc_w1_alias or calc_w2_alias
;-----------------------------------------------------------

; 3 under-matrices are needed
; Palpha' (and its transposed)
; Cb   = a priori on noise for each GS
; Cphi = a priori on turbulence profile

; M.Palpha'
MPalpha = complexarr(n_elements(h_recons)*s, nb_gs*s)
FOR i = 0, n_elements(h_recons)-1 DO BEGIN
   FOR j = 0, nb_gs-1 DO BEGIN
      ff_x = f_x*alpha(0, j)*h_recons(i)*60./206265.
      ff_y = f_y*alpha(1, j)*h_recons(i)*60./206265.
      MPalpha(i*s:(i+1)*s-1, j*s:(j+1)*s-1) = wfs(j*s:(j+1)*s-1, *)*exp(complex(0, 1)*2*(ff_x+ff_y)*!pi)
   ENDFOR
ENDFOR

; suppression of WFS
undefine, wfs

; Transpose
t_MPalpha = complexarr(nb_gs*s, n_elements(h_recons)*s)
FOR i = 0, n_elements(h_recons)-1 DO $
    FOR j = 0, nb_gs-1 DO $
        t_MPalpha(j*s:(j+1)*s-1, i*s:(i+1)*s-1) = Conj(MPalpha(i*s:(i+1)*s-1, $
                                                               j*s:(j+1)*s-1))

; LSE is from the form W1
IF keyword_set(LSE) THEN Wflag = 'W1'


;----------------------------------------------------------------------------------
; WTOMO = W1
;----------------------------------------------------------------------------------
IF Wflag EQ 'W1' THEN BEGIN
   
   calc_w1_noalias, f, arg_f, $
                    pitchs_wfs,pitchs_dm, $
                    Nb_gs,alpha,theta,varnoise_wfs, $
                    DSP_tab_recons, h_recons, h_dm, $
                    seuil, condmax, Wflag, $
                    fitting_kw = fitting_kw, $
                    LSE = LSE, /verbose, $
                    wtomo = wtomo, MPalpha = MPalpha, t_MPalpha = t_MPalpha

ENDIF
;----------------------------------------------------------------------------------




;----------------------------------------------------------------------------------
; WTOMO = W2
;----------------------------------------------------------------------------------
IF Wflag EQ 'W2' THEN BEGIN
   
   calc_w2_noalias, f, arg_f, $
                    pitchs_wfs,pitchs_dm, $
                    Nb_gs,alpha,theta,varnoise_wfs, $
                    DSP_tab_recons, h_recons, h_dm, $
                    seuil, condmax, Wflag, $
                    LSE = LSE, $
                    wtomo = wtomo, MPalpha = MPalpha, t_MPalpha = t_MPalpha
   
ENDIF
;----------------------------------------------------------------------------------




;----------------------------------------------------------------------------------
; Pure tomo case : stop here
; The global reconstructor is the tomographic reconstructor 
IF fitting_kw EQ 0 THEN return, Wtomo
;----------------------------------------------------------------------------------




;----------------------------------------------------------------------------------
; computation of optimal projection 
; depends only on optimisation directions / nb and altitudes of DM
; NOTE : 0 and 0 for f_x and f_y, recomputed inside programm)
; This is different from aliasing case
Popt = calc_popt(theta, s, h_recons, h_dm,pitchs_dm, f, arg_f, 0, 0, condmax)
;----------------------------------------------------------------------------------



;------------------------------------------------------------------------------------
; ON PEUT ALORS ECRIRE WMAP = POPT#WTOMO
; BLOC PRODUCT
Wmap = bloc_product(Popt, Wtomo, s)
;------------------------------------------------------------------------------------

return, WMAP

END
