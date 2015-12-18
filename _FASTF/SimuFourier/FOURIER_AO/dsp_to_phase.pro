FUNCTION dsp_to_phase, file = file

file = 'file_save_oa_nonoise_nowind.dat'        ; save file string
restore, file = file, /verbose

; number of iterations requested
nb_iter = 10.
dim = (size(dspall))[1]
nb_dir_perf = n_elements(beta_tab)/2.

; initialisation of output array : 
; TAB_PHASE SCREEN x SCREEEN x NB_DIR_PERF
tab_phase = fltarr(dim, dim, nb_dir_perf, nb_iter)

FOR ind_iter = 0, nb_iter-1 DO BEGIN
   FOR ind_perf = 0, nb_dir_perf-1 DO begin
      wnoise = randomn(seed, dim, dim)
      wnoisefft = fft(wnoise, -1, /over) * dim * D
      pnoisefft = wnoisefft*sqrt( dspall[*, *, ind_perf] )
      tab_phase[*, *, ind_perf, ind_iter] = real(vfft(pnoisefft, 1, /overwrite)*(1./D)^2)
   ENDFOR
ENDFOR

ENDFOR
