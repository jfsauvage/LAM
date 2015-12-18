PRO  tsvd,  mat =  mat, inverse =  inverse,  condmax = condmax,  condition = $
            condition,  visu =  visu,  val_propre =  val_propre



svdc,  mat,  u, v, w,  /double
sz =  (size(u))(1)
IF keyword_set(visu) THEN plot,  u(sort(u))
IF min(abs(u)) NE 0 THEN BEGIN
    condition =  max(abs(u))/min(abs(u)) 
  ;  print,  'conditionnement avant filtrage  =',  condition
ENDIF ELSE BEGIN
    condition =  0
   ; print,  'conditionnement infini , 0 par défaut'
ENDELSE 
limite =  max(abs(u))/condmax
val =  fltarr(sz)
s =  where(abs(u) LT limite, count)
;print,  'filtrage de',  count,  'modes'
val(where(abs(u) GE limite)) =  1/(u(where(abs(u GE limite))))
diag =  diag_matrix(val)
val_propre =  u
inverse =  w ## diag ## transpose(v)

END
