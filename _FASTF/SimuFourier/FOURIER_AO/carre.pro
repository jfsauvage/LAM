FUNCTION carre,  vect

; reforme un vecteur sens� �tre une image carr�e. si elle l'est pas,
; retourne 0

taille = (size(reform(vect)))[1]
cote = floor(sqrt(taille))
IF sqrt(taille)-cote GT 1.e-12 THEN BEGIN
   print, 'Image non carr�e...je sors : cot� = ' + strc(sqrt(taille))
   return, 0
ENDIF
return, reform(vect, cote, cote)
END
