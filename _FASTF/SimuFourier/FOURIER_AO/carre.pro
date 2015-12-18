FUNCTION carre,  vect

; reforme un vecteur sensé être une image carrée. si elle l'est pas,
; retourne 0

taille = (size(reform(vect)))[1]
cote = floor(sqrt(taille))
IF sqrt(taille)-cote GT 1.e-12 THEN BEGIN
   print, 'Image non carrée...je sors : coté = ' + strc(sqrt(taille))
   return, 0
ENDIF
return, reform(vect, cote, cote)
END
