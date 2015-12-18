; This program computes the reconstruction matrix WMAPincluding aliasing errors.
; accounting for all reconstruction parameters
; WMAP = Popt ## Wtomo

FUNCTION calc_reconstructor_alias, f, fx, fy, arg_f, $
                                   pitchs_wfs, pitchs_dm,$
                                   Nb_gs,alpha, theta,$
                                   varnoise_wfs_xy, $
                                   DSP_tab_recons, h_recons, h_dm, $
                                   seuil, cb_alias, condmax, Wflag, $
                                   fitting_kw = fitting_kw, LSE = LSE

; f = spatial frequencies array
; arg_f = F phase
; Nb_gs = Guide star number
; alpha = Guide stars positions
; theta = Optimisation directions
; signoise_wfs = A priori on noise associated to each GS
; DSP_Tab_recons = A priori, DSP on estimated turbulent layers
; h_recons = Altitudes of reconstructed layers
; h_dm : DM altitude (if not pure tomo)
; condmax : Max acceptable conditionning in inversion for POPT computation
; Wflag : Choice of reconstructor W1 or W2
; Keyword Tomo : Pure tomo
; Popt : output : optimal projector for MCAO, used later for Aliasing
; LSE : LSE instead of MAP

s = (size(f))(1)


;-----------------------------------------------------------------------------
; computation of WFS matrix
;-----------------------------------------------------------------------------
wfs = complexarr(2*nb_gs*s, s)

FOR j = 0, nb_gs-1 DO BEGIN
   ccx =  2*!pi*(fx)*sinc(pitchs_wfs(j)*fx)*sinc(pitchs_wfs(j)*fy)*complex(0, 1)
   ccy =  2*!pi*(fy)*sinc(pitchs_wfs(j)*fx)*sinc(pitchs_wfs(j)*fy)*complex(0, 1)
   fc =  1./(2.*pitchs_wfs(j)) 
   f_ind =  where((f NE 0) and (abs(fx) GT fc) OR (abs(fy) GT fc),  count)
   IF (count NE 0) THEN BEGIN
      ccx(f_ind) =  0.
      ccy(f_ind) = 0.
   ENDIF
   wfs(2*j*s:(2*j+1)*s-1, *) = ccx
   wfs((2*j+1)*s:(2*j+2)*s-1, *) = ccy
ENDFOR

;XY expressions for all the model:
;M.Palpha', contains all the phasors for the reconstruction model.
;Cartesian writing
;ccx = 2*!pi*(fx)*sinc(a*fx)*sinc(a*fy)*complex(0, 1) ;Shack X
;ccy = 2*!pi*(fy)*sinc(a*fx)*sinc(a*fy)*complex(0, 1) ;Shack Y
; see B. Neichel Thesis, page 143 or so


;-----------------------------------------------------------------------------
; computation of MPalpha_xy matrix
;-----------------------------------------------------------------------------
MPalpha_xy = complexarr(n_elements(h_recons)*s, nb_gs*2*s)
FOR i = 0, n_elements(h_recons)-1 DO BEGIN
   FOR j = 0, nb_gs-1 DO BEGIN
      ff_x = fx*alpha(0, j)*h_recons(i)*60./206265.
      ff_y = fy*alpha(1, j)*h_recons(i)*60./206265.
      MPalpha_xy(i*s:(i+1)*s-1, 2*j*s:(2*j+1)*s-1) = wfs(2*j*s:(2*j+1)*s-1, *)*exp(complex(0, 1)*2*(ff_x+ff_y)*!pi)
      MPalpha_xy(i*s:(i+1)*s-1, (2*j+1)*s:(2*j+2)*s-1) = wfs((2*j+1)*s:(2*j+2)*s-1, *)*exp(complex(0, 1)*2*(ff_x+ff_y)*!pi)
   ENDFOR
ENDFOR

undefine, wfs

;-----------------------------------------------------------------------------
; transposition
;-----------------------------------------------------------------------------
MPalpha_t_xy = complexarr(2*nb_gs*s, n_elements(h_recons)*s)
FOR i = 0, n_elements(h_recons)-1 DO $
   FOR j = 0, 2*nb_gs-1 DO $
      MPalpha_t_xy(j*s:(j+1)*s-1, i*s:(i+1)*s-1) = Conj(MPalpha_xy(i*s:(i+1)*s-1, j*s:(j+1)*s-1))

Cb_recons = complexarr(2*nb_gs*s, 2*nb_gs*s)
FOR i = 0, 2*nb_gs-1 DO Cb_recons(s*i:s*(i+1)-1, s*i:s*(i+1)-1) = varnoise_wfs_xy(i)



;-----------------------------------------------------------------------------
; noise matrices
;-----------------------------------------------------------------------------
Cb_noisetot = Cb_alias + cb_recons
undefine, Cb_recons


IF keyword_set(LSE) THEN Wflag = 'W1'

IF Wflag EQ 'W1' THEN $
   calc_w1_alias, f, fx, fy, arg_f, $
                  pitchs_wfs,pitchs_dm, $
                  Nb_gs,alpha,theta,varnoise_wfs_xy, $
                  DSP_tab_recons, h_recons, h_dm, $
                  seuil, condmax, Wflag, $
                  LSE = LSE, verbose = verbose, $
                  wtomo = wtomo, MPalpha_xy = MPalpha_xy, MPalpha_t_xy = MPalpha_t_xy, $
                  cb_alias = cb_alias, Cb_noisetot = Cb_noisetot


IF Wflag EQ 'W2' THEN $
   calc_w2_alias, f, fx, fy, arg_f, $
                   pitchs_wfs,pitchs_dm, $
                   Nb_gs,alpha,theta,varnoise_wfs_xy, $
                   DSP_tab_recons, h_recons, h_dm, $
                   seuil, condmax, Wflag, $
                   LSE = LSE, verbose = verbose, $
                   wtomo = wtomo, MPalpha_xy = MPalpha_xy, MPalpha_t_xy = MPalpha_t_xy, $
                   cb_alias = cb_alias, Cb_noisetot = Cb_noisetot
   

undefine, MPalpha_xy
undefine, MPalpha_t_xy
;----------------------------------------------------------------------------------




;----------------------------------------------------------------------------------
; Pure tomo case : stop here
; The global reconstructor is the tomographic reconstructor 
IF fitting_kw EQ 0 THEN return, Wtomo
;----------------------------------------------------------------------------------




;----------------------------------------------------------------------------------
; computation of optimal projection 
; depends only on optimisation directions / nb and altitudes of DM
; NOTE : fx and fy pre-computed in aliasing case
Popt = calc_popt(theta, s, h_recons, h_dm, pitchs_dm, f, arg_f, fx, fy, condmax)
;----------------------------------------------------------------------------------



;------------------------------------------------------------------------------------
; ON PEUT ALORS ECRIRE WMAP = POPT#WTOMO
; BLOC PRODUCT
Wmap = bloc_product(Popt, Wtomo, s)
;------------------------------------------------------------------------------------


return, Wmap

       
END
