FUNCTION calc_static_phase, static_phase_rms, f, sim_params, seed, pupil

; This function generates a static phase, according to a 1/f^6 spectrum (may be
; modified easily)

; static_phase_rms      : value of RMS for static phase [nm]
; seed                  : use a given seed for random generation
; F                     : fltarr(nb_f, nb_f) containing 2D array of frequencies
; SIM_PARAMS            : structure containing all simulation parameters
; PUPIL                 : pupil map 

; sets local variable (avoid modifying input variables)
local_f = f
local_seed = seed

; avoid pb for null frequency
local_f[0, 0] = 1.0

; sets 1/f6 spectrum for static aberrations (very low orders)
spectrum = 1./(local_f)^6

; sets zero frequency correctly (=No piston)
spectrum[0, 0] = 0.

; generates white random noise in direct (phase) space
noise = randomn(local_seed, sim_params.dim, sim_params.dim)

; computes white noise in fourier space
noise_fft = fft(noise)

; colors the noise in fourier space
pinknoise_fft = sqrt(spectrum) * real(noise_fft)$
                + !i * sqrt(spectrum) * imaginary(noise_fft)

; gets back in direct space
static_phase = real(fft(pinknoise_fft))

; filters piston
static_phase[where(pupil)] -= mean(static_phase[where(pupil)])

; normalise correctly
return, static_phase/sqrt(mean(abs2(static_phase[where(pupil)]))) * static_phase_rms

END

