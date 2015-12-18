;	$Id: crop.pro,v 1.2 2012/10/02 13:17:04 sauvage Exp $	
FUNCTION crop, image, ncrop = ncrop, $
               centre = centre, milieu = milieu, silent = $
               silent, extend = extend

;+
;NOM :
;     AFF_VIDEO - 
; proc�dure permettant d'extraire une sous-image d'une image donn�e.
; .oO NEW Oo. permet aussi d'agrandir une image en la plongeant dans une
; support plus grand si besoin est.
;   
;CATEGORIE :
;     image manipulation
;
;SYNTAXE :
;     
;
;DESCRIPTION : 
; l'image en question peut �tre une matrice ou un nom de fichier, quelle
; que soit l'extension. (on peut lui sp�cifier l'extension "spe" ou "fits" si
; le nom de fichier n'est pas explicite.
; l'extraction est centr�e sur le max de l'image, ou sur "centre", et est
; renvoy�e dans la variable de sortie "image". dans le cas d'une
; extraction centr�e sur le max qui a n�cessit� une remise � z�ro de certains
; pixels (dit pixels pisseux)
;     
;   
;   ARGUMENTS :

; IMAGE : le tableau � croper
; CROP : la taille du tableau apr�s crop
; MILIEU : MC pour croper autour du centre
; CENTRE : position du centre sinon
; SILENT : �vite de polluer un Xterm avec des lignes de comment
;          incompr�hensibles de toutes fa�ons. 
; EXTEND : MC Pour pouvoir aggrandir le support au lieu de le raccourcir
;          (seulement autour du centre pour le moment)
; ERASE : remplace la zone crop�e de l'image par la valeur moyenne au
;         bord, ou par VALUE si fournie.
; VALUE : valeur � donner aux points cropp�s.


IF NOT(keyword_set(ncrop)) THEN BEGIN
    IF NOT(keyword_set(silent)) THEN print, 'Je crope � 128 pixels !! '
    ncrop = 128
ENDIF
j = 0B
cr = 20.

; sauvegarde de l'ancienne image dans image_old
image_old = image

NP = (size(image))(1)
NP2 = (size(image))[2]
IF NOT(keyword_set(silent)) THEN print, 'Taille de l''image : ', strcompress(string(NP), /remove_all), ' pixels'

IF keyword_set(milieu) THEN centre = [NP/2., NP2/2.]

IF (NOT(keyword_set(centre)) or n_elements(centre) EQ 1) THEN BEGIN
    j = 1B
    centre = dblarr(2)
    maxi = max(image[*, *, 0], pos)
    centre(0) = pos MOD floor(NP)
    centre(1) = pos / floor(NP)
    IF NOT(keyword_set(silent)) THEN print, $
        'centre=['+strc(floor(centre[0]))+','+strc(floor(centre[1]))+']'
ENDIF

pts = [centre(0)-ncrop/2., centre(0)+ncrop/2.-1, centre(1)-ncrop/2., $
       centre(1)+ncrop/2.-1]
index = where(pts GT NP-1 OR pts LT 0) 
IF total(index) GT 0 AND NOT(keyword_set(extend)) THEN BEGIN
   print, 'le crop sort de l''image, je sors'
   return, 0
ENDIF

; juste pour le cas o� on souhaite �tendre le support au lieu de le
; raccourcir, voici les lignes correspondantes. attention, seulement pour un
; objet centr�. si la position du crop est d�centr�e et sans /extend, on est
; dans le cas du dessus.
IF keyword_set(extend) THEN BEGIN
   im = dblarr(ncrop, ncrop)
   im[ncrop/2.-NP/2.:ncrop/2.+NP/2.-1, ncrop/2.-NP/2.:ncrop/2.+NP/2.-1] = image
ENDIF

IF NOT(keyword_set(extend)) THEN im = image[centre(0)-ncrop/2.:centre(0)+ncrop/2.-1, $
                                            centre(1)-ncrop/2.:centre(1)+ncrop/2.-1, *, *, *]


return, im
END

