FUNCTION  calc_ee, psf, sim_params, $
                   pix2rad = pix2rad, $
                   width = width, $
                   square = square, verbose = verbose

; This function computes the ensquared energy in a square box.

; pix2rad ; size of pixel in angular radians
; width ; size of box in arcsec
; square ; KW for square box (1) or round box (0)

nb_dirperf = n_elements(sim_params.beta_tab)/2.
rad2arcsec = 1.*180./!dpi*3600.
width_px = width /rad2arcsec /pix2rad
ee = fltarr(nb_dirperf)

; loop on performance direction
FOR ind_dirperf = 0, nb_dirperf -1 DO begin
   
; first extract a box larger than needed (2 times, even number of pixels)
   extract_px = round(width_px*2.0 < sim_params.dim)
   IF extract_px MOD 2 NE 0 THEN extract_px += 1
   psf_extract = crop(psf[*, *, ind_dirperf], /m, nc = extract_px, /silent)
   
; then congrid it
   zoom = 10.
   psf_extract_rebin = congrid(psf_extract, zoom * extract_px, zoom * extract_px, $
                               /interp, /minus)
   
; then extract the desired box for EE computation, even number of pixels
   width_px_rebin = round(zoom * width_px)
   IF width_px_rebin MOD 2 NE 0 THEN width_px_rebin += 1
   psf_extract_rebin_ee = crop(psf_extract_rebin, /m, nc = width_px_rebin, /silent)
   
; multiply by pupil : square if square, round if not
   IF (keyword_set(square)) THEN $
       masque = fltarr(width_px_rebin, width_px_rebin)+1. $
   ELSE polaire2, rt = width_px_rebin/2., largeur = width_px_rebin, masque = $
                  masque
   
   ee[ind_dirperf] = total(masque * psf_extract_rebin_ee)/zoom^2 / $
                     total(psf[*, *, ind_dirperf])
ENDFOR

return,  ee
END
