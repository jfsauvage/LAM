
FUNCTION dsp_to_psf, dsp, maskpupil, sim_params, $
                     static_phase = static_phase


; This function computes the long exposure PSF directly from the DSP of
; residual phase

;; local_L: taille physique du tableau dsp (en metres)

local_L = sim_params.D*sim_params.dim/sim_params.dimpup ;screen size [m]
osamp = sim_params.dim/sim_params.dimpup ; optical sampling factor

; extraction of pupil
pup = crop(maskpupil, /m, nc = sim_params.dimpup)

; residual structure function creation
Bg = fft(dsp, 1)/(local_L)^2.
Dphi = real(2*(Bg[0, 0]-Bg))
Dphi = eclat(Dphi)

; extract pupil only content
dim = sim_params.dim
npup = sim_params.dimpup
Dphi2 = $
    Dphi[dim/2-osamp*npup/2:dim/2+osamp*npup/2-1,dim/2-osamp*npup/2:dim/2+osamp*npup/2-1]

; creation of static FTO array
dlfto = fftshift(norme(abs2(fftshift($
        maskpupil * exp(1.e-9 * 2.*!pi * !i * static_phase/sim_params.lambda_imaging)))))

; creation of turbulent FTO
aoFTO = exp(-Dphi2/2)

; total FTO computation
sysFTO = aoFTO*dlFTO
sysFTO = eclat(sysFTO)

; PSF computation
sysPSF = real(eclat((fft(sysFTO, /INVERSE))))
sysPSF = sysPSF/total(sysPSF) 

return, sysPSF

END
