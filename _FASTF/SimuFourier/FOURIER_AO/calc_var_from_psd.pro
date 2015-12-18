FUNCTION calc_var_from_psd,  psd,  pixsize =  pixsize,  $
                 VERBOSE =  verbose, DD = DD

; Normalisation of DSP
dim =  (size(psd))[2]
psdtemp =  eclat(psd)*pixsize^2.

;Calcul de Fp 
FD =  1./float(DD)

; BOXSIZE : L/D : SCREEN_SIZE/PUP_SIZE => "sampling factor"
boxsize =  FD/pixsize
polaire2, RT =boxsize/2., largeur = dim,  /entre4, masque = ppp
maskdsp = abs(ppp-1)
psdtemp = psdtemp*maskdsp

res =   total(psdtemp)
return, res
delvar, psdtemp
END 

