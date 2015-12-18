FUNCTION calc_cbalias, f, fx,fy, arg_f, a,Nb_gs,alpha, h_recons ,Cn2_recons, L0, r0

; this function computes Cb_Alias :
nbre = nb_gs
f_x = fx
f_y = fy
fc =  1./(2.*a) 
out = f*0.0
s = (size(f))(1)

;--------------------------------------------------------------------------------------
;COvariance matrix of aliased measures : SST
;--------------------------------------------------------------------------------------

;construct a matrix with all pairs of measures X and Y (annex B, B. Neichel
;thesis, p313)
tab = fltarr(2*nb_gs, 2*nb_gs, 2)
FOR i = 0, nb_gs-1 DO BEGIN
   FOR j = 0, nb_gs-1 DO BEGIN

      tab(i*2, j*2, *) = [alpha(0, j)-alpha(0, i),alpha(1, j)-alpha(1, i) ]
      tab(i*2+1, j*2+1, *) = [alpha(0, j)-alpha(0, i),alpha(1, j)-alpha(1, i) ]
      tab(i*2, j*2+1, *) = [alpha(0, j)-alpha(0, i),alpha(1, j)-alpha(1, i) ]
      tab(i*2+1, j*2, *) = [alpha(0, j)-alpha(0, i),alpha(1, j)-alpha(1, i) ]
      
   ENDFOR
ENDFOR
; TAB contains angles, X and Y, for all direction combinations
; to be used as DALPHA

; aliasing : summation of shifted spectrum (shifted of 2*fc).

; summation should be done on infinite number (sum_k=1...infinity):
; but 3 is enough (only neighbour spectrum)
nmax = 3.
; (like in PAOLA)

Cov_S = complexarr(nb_gs*2.*s, nb_gs*2.*s)
;on rempli la matrice ssT
FOR i = 0, nb_gs-1 DO BEGIN
   FOR j = 0, nb_gs-1 DO BEGIN
      flag = 'xx'
      dalpha = reform(tab(i*2, j*2, *))
      tmp = make_sum3(f, f_x, f_y, arg_f, cn2_recons, a, r0, L0, flag, nmax, dalpha, h_recons)
      Cov_s(i*2*s, j*2*s) = (tmp)
      
      flag = 'yy'
      dalpha = reform(tab(i*2+1, j*2+1, *))
      tmp = make_sum3(f, f_x, f_y, arg_f, cn2_recons, a, r0, L0, flag, nmax, dalpha, h_recons)
      Cov_s((i*2+1)*s, (j*2+1)*s) = (tmp)
      
      flag = 'xy'
      dalpha = reform(tab(i*2, j*2+1, *))
      tmp = make_sum3(f, f_x, f_y, arg_f, cn2_recons, a, r0, L0, flag, nmax, dalpha, h_recons)
      Cov_s(i*2*s, (j*2+1)*s) = (tmp)
      
      flag = 'xy'
      dalpha = reform(tab(i*2+1, j*2, *))
      tmp = make_sum3(f, f_x, f_y, arg_f, cn2_recons, a, r0, L0, flag, nmax, dalpha, h_recons)
      Cov_s((i*2+1)*s, j*2*s) = (tmp)
      
   ENDFOR
ENDFOR

return, float(Cov_S)


END
