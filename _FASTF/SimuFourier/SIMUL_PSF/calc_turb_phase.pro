FUNCTION calc_turb_phase, dsp, res, sim_params, pupil

; This function generates a phase map according to a given DSP
; The generation is done in a Fourier space

; RES       : residual value [nm]
; DSP       : spectrum of turbulent phase
; SIM_PARAMS   : structure containing all simulation parameters

; sets user defined spectrum for turbulent aberrations
spectrum = dsp

; sets zero frequency correctly (=No piston)
spectrum[0, 0] = 0.

; generates white random noise in direct (phase) space
noise = randomn(seed, sim_params.dim, sim_params.dim)

; computes white noise in fourier space
noise_fft = fft(noise)

; colors the noise in fourier space
pinknoise_fft = sqrt(spectrum) * real(noise_fft)$
                + !i * sqrt(spectrum) * imaginary(noise_fft)

; gets back in direct space
turb_phase = real(fft(pinknoise_fft))

; filters piston
turb_phase[where(pupil)] -= mean(turb_phase[where(pupil)])

; normalise correctly
return, turb_phase/sqrt(mean(abs2(turb_phase[where(pupil)]))) * res

END
