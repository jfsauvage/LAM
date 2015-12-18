PRO affa_3d, image, linear = linear, _REF_EXTRA = extra, sortie = sortie

; procédure d'affichage log d'une pile d'image à la meme dynamique.
; les images sont affichées juxtaposées.

	NP = (size(image))[1]
    NL = (size(image))[2]
    
    IF (size(image))[0] EQ 2 THEN BEGIN
       IF keyword_set(sortie) THEN sortie = image
       affa, image, _extra = ['z', 'title'], linear = linear 
    ENDIF ELSE BEGIN 
       Nimage = (size(image))[3]
       IF NP NE NL THEN BEGIN
          sortie = fltarr(NP * Nimage, NL)
          FOR i = 0, Nimage-1 DO sortie[NP*i:NP*(i+1)-1, *] = image[*, *, i]
       ENDIF ELSE sortie = transpose(reform(transpose(image, [1, 0, 2]), NP, NL*Nimage))
       affa, sortie, _EXTRA = ['z', 'title'], linear = linear
    ENDELSE
    
END
    
