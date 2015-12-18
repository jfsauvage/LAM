;Programme pour prendre en compte la multi-analyse les geometries
;d'etoiles et la position de la galaxie

FUNCTION dsp_to_psf_old, dsp, pup, local_L, osamp
;+
;FUNCTION dsp_to_psf, dsp, pup, local_L, osamp
;; calcule une PSF d'apres une dsp de phase dsp et une pupille pup
;;dsp: tableau 2D contenant la DSP de la phase
;;pup: tableau 2D contenant la fonction pupille
;;local_L: taille physique du tableau dsp (en metres)
;;osamp: facteur de sur-echantillonnage dans le plan focal (2 par defaut)
;-

;;creation de la fonction de structure de la phase
Bg = fft(dsp, 1)/(local_L)^2.
;;on cree alors la fonction de structure
Dphi = real(2*(Bg[0, 0]-Bg))
Dphi = eclat(Dphi)

;on prend ce qui correspond a la pupille
dim = (size(dsp))[1]
npup = (size(pup))[1]
Dphi2 = $
    Dphi[dim/2-osamp*npup/2:dim/2+osamp*npup/2-1,dim/2-osamp*npup/2:dim/2+osamp*npup/2-1]

;creation de la FTO diff limited (autocorrelation de la pupille)
tab = fltarr(osamp*npup, osamp*npup)
tab[0,0] = pup
dlFTO = float(FFT(abs(FFT(tab,1))^2, -1))
dlFTO = eclat(abs(dlFTO)/total(pup))

;creation de la FTO AO
aoFTO = exp(-Dphi2/2)

;;calcul de la FTO resultante
sysFTO = aoFTO*dlFTO
sysFTO = eclat(sysFTO)

;;calcul de la PSF
sysPSF = real(eclat((fft(sysFTO, /INVERSE))))
sysPSF = sysPSF/total(sysPSF) ;normalisation en energie

return, sysPSF

END
