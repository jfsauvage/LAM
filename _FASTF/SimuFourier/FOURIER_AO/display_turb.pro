PRO display_turb, Cn2true, h, Cn2recons, recons_h, wind, wind_dir, h_dm, $; INPUTS
                  ; OPTIONAL INPUTS
                  window_number = window_number, window_size = window_size, $
                  thick = thick, $
                  ; OPTIONAL INPUT KEYWORD
                  inv_color = inv_color, help = help                            
;+
; NAME:
;     DISPLAY_TURB
; PURPOSE:
;     This procedure displays in a single plot window: 
;          * Turbulence parameters (i.e. true Cn2 and altitudes, reconstructed Cn2 and altitudes)
;          * Wind parameters (modulus and direction)
;          * DM altitudes.
; EXPLANATION:
;
; CALLING SEQUENCE:
;     DISPLAY_TURB, CN2TRUE, H, CN2RECONS, RECONS_H, WIND, WIND_DIR, H_DM, $
;                   [[THICK = THICK], [WINDOW_NUMBER = WINDOW_NUMBER], [WINDOW_SIZE = WINDOW_SIZE], $
;                   [/INV_COLOR], [/HELP]]
;
; INPUTS:
;     CN2TRUE   - True turbulence distribution strength used for the simulation. Values in %
;     H         - True altitude of turbulence layers used for the simulation. Values in meters.
;     CN2RECONS - Turbulence distribution used to estimated turbulence (reconstuct). Values in %
;     RECONS_H  - Altitude of turbulence layers used in the reconstruction. Values in meters.
;     WIND      - True wind modulus used in simulation. Values in meters per second. 
;     WIND_DIR  - True wind direction used in simulation. Values in radians.
;     H_DM      - Altitude of DM. Values in meters.
;
; OPTIONAL INPUTS:
;     WINDOW_NUMBER - (default: 0). Number of the display window. 
;     WINDOW_SIZE   - (default: [640, 510]). Size of the display window. 
;                     i.e.: "window, window_number, xs = window_size[0], ys = window_size[1]"
;     THICK         - (default: calculated from data). Size of the display bars for turbulence. 
;
; OPTIONAL INPUT KEYWORD:
;     INV_COLOR - By default plot white on black background (default INV_COLOR=0).
;                 By setting keyword (INV_COLOR=1), plot black and white background (for printing).
;     HELP      - Print procedure ussage and returns.
;     
; OUTPUTS:
;     None
;
; EXAMPLES:
;     Plot Cn2true, h, Cn2recons, recons_h, wind, wind_dir and h_dm with default setting
;     IDL> DISPLAY_TURB, Cn2true, h, Cn2recons, recons_h, wind, wind_dir, h_dm
;     
;     Plotting on window number 10, a 650 by 750 window and on a white background
;     IDL> DISPLAY_TURB, Cn2true, h, Cn2recons, recons_h, wind, wind_dir, h_dm,$
;                        window_number = 10, window_size = [650, 750], /Inv_color
;
; TO DO:
;     - Find a better way to print wind direction.
; 
; SEE ALSO:
;   FSC_COLOR; PLOTSYM;
;     
; REVISION HISTORY
;    01/10/2012 : Creation by Noah SCHWARTZ
;-

IF keyword_set(help) THEN BEGIN
     print,'Syntax - '
     print,'DISPLAY_TURB, CN2TRUE, H, CN2RECONS, RECONS_H, WIND, WIND_DIR, H_DM, $'
     print,'              [[WINDOW_NUMBER = WINDOW_NUMBER], [WINDOW_SIZE = WINDOW_SIZE], $'
     print,'              [THICK = THICK], [/INV_COLOR], [/HELP]]'
     RETURN
 ENDIF

; -------------------
; INPUT KEYWORD
; -------------------
; Missing keywords    
IF NOT keyword_set(window_number) THEN window_number = 1               ; Default window number
IF NOT keyword_set(window_size) THEN window_size = [640, 510]          ; Default window size
; Color convention
IF NOT keyword_set(inv_color) THEN BEGIN
    bckgrnd_col = 'black' 
    frgrnd_col = 'Snow' 
ENDIF ELSE BEGIN
    bckgrnd_col = 'Snow' 
    frgrnd_col = 'black' 
ENDELSE 


; -------------------
; DISPLAY WINDOW
; -------------------
window, window_number, xs = window_size[0], ys = window_size[1]
tek_color
!p.multi = [0,1,2,0,0]

; -- Calculate plot parameters --
minh = 200. ; minimum altitude on plot in [m]
IF NOT keyword_set(thick) THEN thick = max([h, recons_h, h_dm, minh])/180. ; Height of display bars
m2km = 1000. ; Convert meters in kilometers
xr = [0, 1.1]
yr = [-2.*thick/m2km, (max([h, recons_h, h_dm, minh])*1.1)/m2km]
nbr_Cn2True = nbr2str(n_elements(Cn2True))
nbr_Cn2recons = nbr2str(n_elements(Cn2recons))
nbr_DM = nbr2str(n_elements(h_DM))

; Plot coordinates for Turbulence
plot, Cn2true, /nodata, xr = xr, yr = yr, position = [[0.1, 0.45], [0.9, 0.95]], $
      Background = fsc_color(bckgrnd_col), color = fsc_color(frgrnd_col), $
      xtit = 'TURBULENCE STRENGTH [%]', ytit = 'LAYER ALTITUDE [km]', $
      Tit = nbr_DM +' DMs, '+nbr_Cn2True +' True Layers and '+nbr_Cn2recons +' Reconstructed Layers'


; -- True Turbulence -- 
FOR ind_layer=0, n_elements(h)-1 DO BEGIN
    Y0 = (h[ind_layer]-thick/2)/m2km
    Y1 = (h[ind_layer]+thick/2)/m2km
    X0 = 0
    X1 = Cn2true[ind_layer]
    ; Draws a box whose corners are (X0, Y0) and (X1, Y1)
    POLYFILL, [X0, X0, X1, X1], [Y0, Y1, Y1, Y0], COL = fsc_color('red')
END 

; -- Reconstructed Turbulence -- 
FOR ind_layer=0, n_elements(recons_h)-1 DO BEGIN
    Y0 = (recons_h[ind_layer]+thick/2.)/m2km
    Y1 = (recons_h[ind_layer]+thick*3./2)/m2km
    X0 = 0
    X1 = Cn2recons[ind_layer]
    ; Draws a box whose corners are (X0, Y0) and (X1, Y1)
    POLYFILL, [X0, X0, X1, X1], [Y0, Y1, Y1, Y0], COL = fsc_color('green')
END

; -- Altitude DM --
Axis, YAxis=1, yr = yr, /Save, YTitle='DM CONJUGATION ALTITUDE [km]', color = fsc_color(frgrnd_col)
FOR ind_layer=0, n_elements(h_DM)-1 DO BEGIN
    Y0 = (h_DM[ind_layer]-thick)/m2km
    Y1 = (h_DM[ind_layer]+thick)/m2km
    X0 = xr[1]
    X1 = max([max([Cn2true, Cn2recons]) + 0.05, 0.1])
    ; Draws a box whose corners are (X0, Y0) and (X1, Y1)
    POLYFILL, [X0, X0, X1, X1], [Y0, Y1, Y1, Y0], COL = fsc_color('blue') 
END


; --- DISPLAY LEGEND ---
X0 = 0.63*xr[1] & X1 = 0.95*xr[1]
Y0 = 0.90*yr[1] & Y1 = 0.92*yr[1]
PLOTSYM, 8, /FILL
oplot, [X0], [Y0], psym = 8, color = fsc_color('red')
xyouts, X0+0.02, Y0*0.89/0.9, 'True Turbulence Layer', color = fsc_color(frgrnd_col)

oplot, [X0], [Y0*0.85/0.9], psym = 8, color = fsc_color('green')
xyouts, X0+0.02, Y0*0.84/0.9, 'Reconstructed Turbulence', color = fsc_color(frgrnd_col)

oplot, [X0], [Y0*0.80/0.9], psym = 8, color = fsc_color('blue')
xyouts, X0+0.02, Y0*0.79/0.9, 'DM Conjugation', color = fsc_color(frgrnd_col)

X0 = 0.67 & X1 = 1.05
Y0 = 0.77*yr[1] & Y1 = 0.94*yr[1]
oplot, [X0, X1], [Y0, Y0], color = fsc_color(frgrnd_col)
oplot, [X0, X0], [Y0, Y1], color = fsc_color(frgrnd_col)
oplot, [X1, X1], [Y0, Y1], color = fsc_color(frgrnd_col)
oplot, [X0, X1], [Y1, Y1], color = fsc_color(frgrnd_col)


; -------------------
;      WIND SPEED
; -------------------
xr = [0, max(abs(wind))*1.1]
plot, Cn2true, /nodata, xr = xr, yr = yr, position = [[0.1, 0.1], [0.9, 0.35]], $
      xtit = 'WIND SPEED [m/s]', ytit = 'LAYER ALTITUDE [km]', $
      Background = fsc_color(bckgrnd_col), color = fsc_color(frgrnd_col)
FOR ind_layer=0, n_elements(wind)-1 DO BEGIN
    Y0 = (h[ind_layer]-thick)/m2km
    Y1 = (h[ind_layer]+thick)/m2km
    X0 = xr[0]
    X1 = (abs(wind))[ind_layer]
    ; Draws a box whose corners are (X0, Y0) and (X1, Y1)
    POLYFILL, [X0, X0, X1, X1], [Y0, Y1, Y1, Y0], color = fsc_color(frgrnd_col)
    xyouts, X1+0.02, Y1, nbr2str(wind_dir[ind_layer]*!RADEG)+'!Uo!N', color = fsc_color(frgrnd_col)
END


; -- Reset multiplot to normal value --
!p.multi = [0,0,0,0,0]
END
; END PROCEDURE





; *******************
;     Test MAIN
; *******************
; -- TURBULENCE --
; Simulation turbulence definition
Cn2true = [52.24, 2.60, 4.44, 11.60, 9.89,  2.95,  5.98, 4.30, 6.00]/100.
h = [47., 140., 281., 562., 1125., 2250., 4500., 9000., 18000.]
; Reconstruction turbulence a priori 
recons_h = h                 ; supposed altitudes for estimation 
cn2recons = Cn2true          ; supposed CN2 for estimation

; -- WIND --
wind = [15.,  13.,  13., 9.,  9.,  15.,  25.,  40.,  21.] ; wind module [m/s]
wind_dir  = randomn(seed, n_elements(wind))                       ; wind dir.  [rad]

; -- DM(s) ALTITUDE --
h_dm = [0, 500, 5000]   
;display_turb, Cn2true, h, Cn2recons, recons_h, wind, wind_dir, h_dm


display_turb, Cn2true, h, Cn2recons, recons_h, wind, wind_dir, h_dm, window_size = [750,650]
display_turb, Cn2true, h, Cn2recons, recons_h, wind, wind_dir, h_dm, window_number = 2, /INV, thick = 200
end

