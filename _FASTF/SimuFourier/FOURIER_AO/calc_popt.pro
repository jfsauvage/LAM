FUNCTION calc_popt, theta, s, h_recons, h_dm, pitchs_dm, f, arg_f, f_x, f_y, condmax

;----------------------------------------------------------------------------------
; Now from Wtomo to Popt, need computation of optimal projector 
; OK for several DMs or not, but needs several optimisation directions
; With only 1 optimisation direction (MOAO) no more than 1 DM, or POPT can
; not be written 
;-----------------------------------------------------------

IF N_elements(f_x) EQ 1 THEN f_x = f*cos(arg_f)
IF n_elements(f_y) EQ 1 THEN f_y = f*sin(arg_f)

nnn = n_elements(size(theta))-1
nb_dir = (size(theta))(nnn)/2.

; The formulae for Popt is given in B. Neichel 7.17, pp 180
; Popt = [<(P_dm # N)T P_dm # N>_beta]^-1 # <(P_dm # N)T # P_L>_beta
; MAT1 = <(P_dm # N)T # P_L>_beta
; MAT1 = <(proj_dm # minf)T # proj_l>_beta

; MAT2 = [<(P_dm # N)T P_dm # N>_beta]^-1
; MAT2 = [<(proj_dm # minf)T proj_dm # minf>_beta]^-1
;------------------------------------------------------------------------------------


;------------------------------------------------------------------------------------
; computation of mat1
;------------------------------------------------------------------------------------
Mat1 = complexarr(s*n_elements(h_recons), s*n_elements(h_dm))

fc_dm = 1./(2*pitchs_dm)

; definition of DM filtering functions
minf = complexarr(s*n_elements(h_dm), s)
FOR j = 0, n_elements(h_dm)-1 DO BEGIN
   indi =  where((abs(f_x) LT fc_dm(j)) OR (abs(f_y) LT fc_dm(j)))
   tab =  (minf(j*s:(j+1)*s-1, *))
   tab(indi) = replicate(complex(1., 0),  n_elements(indi))
   (minf(j*s:(j+1)*s-1, *)) =  tab
ENDFOR

; loop on optimisation directions :
FOR th = 0, nb_dir-1 DO BEGIN

; Projector of DM to pup in 1 direcction
proj_dm = complexarr(s*n_elements(h_dm), s)
FOR j = 0, n_elements(h_dm)-1 DO $
    proj_dm(j*s:(j+1)*s-1, *) = exp(complex(0, 1)*2*!pi*h_dm(j)*60./206265.*(theta[0, th]*f_x+theta[1, th]*f_y))

; and transpose :
proj_conj_dm = complexarr(s, s*n_elements(h_dm))
FOR j = 0, n_elements(h_dm)-1 DO $
    proj_conj_dm(*, j*s:(j+1)*s-1) = conj(proj_dm(j*s:(j+1)*s-1, *))

; fitting :
pminf =  complexarr(s*n_elements(h_dm), s)
FOR j = 0, n_elements(h_dm)-1 DO $
  pminf(j*s:(j+1)*s-1, *) = proj_dm(j*s:(j+1)*s-1, *) * minf(j*s:(j+1)*s-1, *)

; and transpose :
pminf_conj = complexarr(s, s*n_elements(h_dm))
FOR j = 0, n_elements(h_dm)-1 DO $
    pminf_conj(*, j*s:(j+1)*s-1) = conj(pminf(j*s:(j+1)*s-1, *))

; Projector of DM to pup in same direcction
proj_L = complexarr(s*n_elements(h_recons), s)
FOR j = 0, n_elements(h_recons)-1 DO $
    proj_L(j*s:(j+1)*s-1, *) = exp(complex(0, 1)*2*!pi*h_recons(j)*60./206265.*(theta[0, th]*f_x+theta[1, th]*f_y))

;  mat1 is product of the pminf_conj and proj_l
FOR i = 0, n_elements(h_dm)-1 DO $
    FOR j = 0,  n_elements(h_recons)-1 DO $
        Mat1(s*j:(j+1)*s-1, s*i:(i+1)*s-1) += pminf_conj(*, s*i:(i+1)*s-1)*proj_L(s*j:(j+1)*s-1, *)

ENDFOR



;------------------------------------------------------------------------------------
; mat2 computation 
;------------------------------------------------------------------------------------
Mat2 = complexarr(s*n_elements(h_dm), s*n_elements(h_dm))

; fitting
minf = complexarr(s*n_elements(h_dm), s)
FOR j = 0, n_elements(h_dm)-1 DO BEGIN
   indi =  where((abs(f_x) LT fc_dm(j)) AND (abs(f_y) LT fc_dm(j)))
   tab =  (minf(j*s:(j+1)*s-1, *))
   tab(indi) = replicate(complex(1., 0),  n_elements(indi))
   (minf(j*s:(j+1)*s-1, *)) =  tab
ENDFOR

; loop on optimisation directions
FOR th = 0, nb_dir-1 DO BEGIN

; Projector of DM to pup in 1 direcction
proj_dm = complexarr(s*n_elements(h_dm), s)
FOR j = 0, n_elements(h_dm)-1 DO BEGIN
proj_dm(j*s:(j+1)*s-1, *) = exp(complex(0, 1)*2*!pi*h_dm(j)*60./206265.*(theta[0, th]*f_x+theta[1, th]*f_y))
ENDFOR

; and transpose : 
proj_conj_dm = complexarr(s, s*n_elements(h_dm))
FOR j = 0, n_elements(h_dm)-1 DO BEGIN
proj_conj_dm(*, j*s:(j+1)*s-1) = conj(proj_dm(j*s:(j+1)*s-1, *))
ENDFOR

; fitting
pminf =  complexarr(s*n_elements(h_dm), s)
FOR j = 0, n_elements(h_dm)-1 DO BEGIN
  pminf(j*s:(j+1)*s-1, *) = proj_dm(j*s:(j+1)*s-1, *) * minf(j*s:(j+1)*s-1, *)
ENDFOR

; and transpose :
pminf_conj = complexarr(s, s*n_elements(h_dm))
FOR j = 0, n_elements(h_dm)-1 DO BEGIN
pminf_conj(*, j*s:(j+1)*s-1) = conj(pminf(j*s:(j+1)*s-1, *))
ENDFOR

; mat2 is the product of minf_conj and pminf
FOR i = 0, n_elements(h_dm)-1 DO BEGIN
FOR j = 0,  n_elements(h_dm)-1 DO BEGIN

Mat2(s*j:(j+1)*s-1, s*i:(i+1)*s-1) += pminf_conj(*, s*i:(i+1)*s-1)*pminf(s*j:(j+1)*s-1, *)

ENDFOR
ENDFOR

ENDFOR


; computatin of inverse of Popt
inv_Popt = calc_inv_map(Mat2, s, condmax, nb_gs, h_recons, h_dm, 'POPT')




;------------------------------------------------------------------------------------
; adding MINF à INV_OPPT
FOR j = 0, n_elements(h_dm)-1 DO $
    FOR k = 0, n_elements(h_dm)-1 DO $
        inv_Popt(k*s:(k+1)*s-1, j*s:(j+1)*s-1) *= minf(j*s:(j+1)*s-1, *)

;------------------------------------------------------------------------------------




;------------------------------------------------------------------------------------
; computation of DM projection
Popt = bloc_product(inv_Popt, Mat1, s)
;------------------------------------------------------------------------------------


return, Popt

END
