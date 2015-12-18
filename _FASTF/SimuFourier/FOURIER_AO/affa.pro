PRO affa, image, seuil = seuil, _ref_extra = extra, linear = linear

; affichage logarithmique d'une image. seuillage à 1.e-10 si non défini.

IF NOT(keyword_set(seuil)) THEN seuil = max(image)*1.e-10

IF keyword_set(linear) THEN im = (image) ELSE im = alog(image > seuil)

tvwin, im, _extra = ['title', 'z']

END
