PRO display_geo, $
                ; INPUTS
                alpha, beta, theta, ngs = ngs, $
                ; OPTIONAL INPUTS
                nb_gs = nb_gs, nb_ngs = nb_ngs, Nb_dir_perf = Nb_dir_perf,  nb_dir_theta = nb_dir_theta, $ 
                sepchampvisu = sepchampvisu, window_number = window_number, window_size = window_size, $
                ; OPTIONAL INPUT KEYWORD
                inv_color = inv_color, help = help                              
;+
; NAME:
;     DISPLAY_GEO
; PURPOSE:
;     Display the Simulation Geometry for LGS, NGS, Optimisation directions and Performance directions in a single window.
; EXPLANATION:
;
; CALLING SEQUENCE:
;     DISPLAY_GEO, ALPHA, BETA, THETA,[[NGS = NGS], [NB_GS = NB_GS], [NB_NGS = NB_NGS], [NB_DIR_PERF = NB_DIR_PERF], $
;                 [NB_DIR_THETA = NB_DIR_THETA], [SEPCHAMPVISU = SEPCHAMPVISU], [WINDOW_NUMBER = WINDOW_NUMBER], $
;                 [WINDOW_SIZE = WINDOW_SIZE], [/INV_COLOR], [/HELP]]
;
; INPUTS:
;     ALPHA -  Laser Guide Stars Directions in Arcmin. Example: alpha = [[0, 1], [1, 1], [-1, 0], [-1, -1]]
;     BETA  -  Performance Directions in Arcmin. Example: beta_tab = [[0., 0.],[0., 50.],[0., 80.]]/60 
;     THETA -  Optimization Directions in Arcmin. Example: theta = [[0, 0], [0, 0]]/60. 
;
; OPTIONAL INPUTS:
;     NGS           - Position of Natural Guide Stars (NGS) in the Field.
;     NB_GS         - (default: calculated from ALPHA). Number of LGS directions.
;     NB_NGS        - default: calculated from NGS). Number of NGS directions.
;     NB_DIR_PERF   - (default: calculated from BETA). Number of Performance directions.
;     NB_DIR_THETA  - (default: calculated from THETA). Number of Optimization directions.
;     SEPCHAMPVISU  - (default: 0) Visualisation FoV. The actual visual FoV is max(SEPCHAMPVISU, ALPHA, BETA, THETA).
;     WINDOW_NUMBER - (default: 0). Number of the display window. 
;     WINDOW_SIZE   - (default: [600, 600]). Size of the display window. 
;                     i.e.: "window, window_number, xs = window_size[0], ys = window_size[1]"
;
; OPTIONAL INPUT KEYWORD:
;     INV_COLOR     - By default plot white on black background (default INV_COLOR=0).
;                     By setting keyword (INV_COLOR=1), plot black and white background (for printing).
;     HELP          - Print procedure ussage and returns.
;     
; OUTPUTS:
;     None
;
; EXAMPLES:
;     Plot alpha, theta and beta with default setting
;     IDL> DISPLAY_GEO, alpha, theta, beta
;     
;     Plotting on window number 10, a 650 by 750 window and on a white background
;     IDL> DISPLAY_GEO, alpha, theta, beta, ngs = ngs, window_number = 10, window_size = [650, 750], /Inv_volor
;
; TO DO:
; 
; SEE ALSO:
;   FSC_COLOR; PLOTSYM;
;     
; REVISION HISTORY
;    01/03/2012 : Creation by Jean-Francois SAUVAGE.  
;    01/10/2012 : Modification by Noah SCHWARTZ
;                   - Addition of the 'Test MAIN'
;                   - Added comments and description
;                   - Important change in calling sequence (clarify inputs, optional inputs and keywords).
;                   - Added new keywords: 
;                           /INV_COLOR for inverting colors (white background).
;                           WINDOW_NUMBER = WINDOW_NUMBER; Specify window number (prevents erasing currently use window).
;                           WINDOW_SIZE = WINDOW_SIZE: Specify window size.
;                   - Change procedure name from DISPLAY.PRO to DISPLAY_GEO.PRO.
;-


; -------------------
; ERROR MESSASGE
; -------------------
IF (N_ELEMENTS(ALPHA) LT 1) OR (N_ELEMENTS(BETA) LT 1) OR (N_ELEMENTS(THETA) LT 1) $
   OR keyword_set(help) THEN BEGIN
     print,'Syntax - '
     print,'DISPLAY_GEO, ALPHA, BETA, THETA,[[NGS = NGS], [NB_GS = NB_GS], [NB_NGS = NB_NGS], [NB_DIR_PERF = NB_DIR_PERF], $'
     print,'             [NB_DIR_THETA = NB_DIR_THETA], [SEPCHAMPVISU = SEPCHAMPVISU], [WINDOW_NUMBER = WINDOW_NUMBER], $'
     print,'             [WINDOW_SIZE = WINDOW_SIZE], [/INV_COLOR], [/HELP]]'
     RETURN
 ENDIF
 
 
; -------------------
; INPUT KEYWORD
; -------------------
; -- Missing keywords --
; Number of Laser Guide Stars
IF NOT keyword_set(nb_gs) THEN BEGIN          
    IF n_elements(size(alpha)) EQ 5 THEN ind = 2 ELSE ind = 0   
    nb_gs = (size(alpha))[ind]
ENDIF 
; Number of Natural Guide Stars
IF NOT keyword_set(nb_ngs) AND keyword_set(ngs) THEN BEGIN          
    IF n_elements(size(ngs)) EQ 5 THEN ind = 2 ELSE ind = 0   
    nb_ngs = (size(ngs))[ind]
ENDIF      
; Number of performance Directions 
IF NOT keyword_set(nb_dir_perf) THEN BEGIN
    IF n_elements(size(beta)) EQ 5 THEN ind = 2 ELSE ind = 0   
    nb_dir_perf = (size(beta))[ind] 
ENDIF   
; Number of Optimization Directions
IF NOT keyword_set(nb_dir_theta) THEN BEGIN
    IF n_elements(size(theta)) EQ 5 THEN ind = 2 ELSE ind = 0   
    nb_dir_theta = (size(theta))[ind] 
ENDIF        
IF NOT keyword_set(sepchampvisu) THEN sepchampvisu = 0                 ; Size of dFoV     
IF NOT keyword_set(window_number) THEN window_number = 0               ; Default window number
IF NOT keyword_set(window_size) THEN window_size = [600, 600]          ; Default window size
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
; Create display window
window, window_number, xs = window_size[0], ys = window_size[1]
tek_color  ;Load a color table similar to the default Tektronix 4115 color table.

; Calculate size of Field of View (dFoV) to display 
dFoV = [-1.5, 1.5] * max([sepchampvisu, max(abs(alpha)), max(abs(theta)), max(abs(beta))])

; Plot display frame with appropriate FoV
plot, dFoV, dFoV, /nodata, position = [[0.1, 0.1], [0.95, 0.95]], $
      xtit = 'ARCMIN', ytit = 'ARCMIN', Background = fsc_color(bckgrnd_col), $
      ticklen = 1, xgridstyle = 2, ygridstyle = 2, color = fsc_color('steel blue')


; -------------------
; DISPLAY 
; LGS (ALPHA), 
; NGS (NGS), 
; Optimisation directions (THETA) 
; Performance directions (BETA)
; -------------------

; --- Display LGS analyse directions  (Gold Ellispes) --- 
FOR iii = 0, nb_gs-1 DO BEGIN
    ; Create symbols with correct elongation direction
    pos_ang =  angle(COMPLEX(alpha[0, iii], alpha[1, iii]))*!RADEG ;Position angle in degrees
    ang = pos_ang/!RADEG            ;Position angle in radians
    phi = 2*!pi*(findgen(49)/(49-1));Divide circle into 50 points
    psize = 1.5                     ;Set the ellipse size (Make it at least bigger than 0.1)
    rmin = 0.5*psize                ;Set Elongation
    rmax = 1.5*psize
    xprime = rmax*cos(phi)*cos(ang) - rmin*sin(phi)*sin(ang)    ;Rotate to desired position angle
    yprime = rmax*cos(phi)*sin(ang) + rmin*sin(phi)*cos(ang)
    usersym, xprime, yprime, FILL = 1,thick = 1 
    
    ; Plot symbols
    oplot, [alpha[0, iii]], [alpha[1, iii]], psym = 8, color =  fsc_color('gold')
ENDFOR
 
; --- Display NGS analyse directions (sky blue stars) --- ==> To be added!
PLOTSYM, 3, 1.5, /FILL
IF keyword_set(ngs) THEN BEGIN 
    FOR iii = 0, nb_ngs-1 DO BEGIN
       oplot, [ngs[0, iii]], [ngs[1, iii]], psym = 8, color =  fsc_color('sky blue'), thick = 2
    ENDFOR
ENDIF

; --- Display performance direction (White/Black circles) ---
PLOTSYM, 0, .8, /FILL; PLOTSYM: Define useful plotting symbols not in the standard !PSYM definitions
FOR iii = 0, Nb_dir_perf-1 DO BEGIN 
   oplot, [beta[0, iii]], [beta[1, iii]], psym = 8, color = fsc_color(frgrnd_col)
ENDFOR

; --- Display optimization directions (Red squares) --- 
PLOTSYM, 8, 1.
FOR iii = 0, nb_dir_theta-1 DO BEGIN
   oplot, [theta[0, iii]], [theta[1, iii]], psym = 8, color =  fsc_color('red'), thick = 2
ENDFOR



; -------------------
; DISPLAY LEGEND
; -------------------
PLOTSYM, 8
oplot, [0.9*dFoV[0]], [0.85*dFoV[1]], psym = 8, color = fsc_color('red')
xyouts, 0.85*dFoV[0], 0.84*dFoV[1], 'Optimisation directions (Theta)', color = fsc_color(frgrnd_col)

usersym, xprime, yprime, FILL = 1, thick = 1 
;PLOTSYM, 8, /FILL
oplot, [0.9*dFoV[0]], [0.8*dFoV[1]], psym = 8, color = fsc_color('gold')
xyouts, 0.85*dFoV[0], 0.79*dFoV[1], 'Laser Guide Stars (Alpha)', color = fsc_color(frgrnd_col)

PLOTSYM, 0, /FILL
oplot, [0.9*dFoV[0]], [0.75*dFoV[1]], psym = 8, color = fsc_color(frgrnd_col)
xyouts, 0.85*dFoV[0], 0.74*dFoV[1], 'PSF estimation directions (Beta)', color = fsc_color(frgrnd_col)

PLOTSYM, 3, /FILL
oplot, [0.9*dFoV[0]], [0.70*dFoV[1]], psym = 8, color = fsc_color('sky blue')
xyouts, 0.85*dFoV[0], 0.69*dFoV[1], 'Natural Guide Stars', color = fsc_color(frgrnd_col)

X0 = 0.1*dFoV[0] & X1 = 0.95*dFoV[0]
Y0 = 0.65*dFoV[1] & Y1 = 0.90*dFoV[1]
oplot, [X0, X1], [Y0, Y0], color = fsc_color(frgrnd_col)
oplot, [X0, X0], [Y0, Y1], color = fsc_color(frgrnd_col)
oplot, [X1, X1], [Y0, Y1], color = fsc_color(frgrnd_col)
oplot, [X0, X1], [Y1, Y1], color = fsc_color(frgrnd_col)

END
; END OF DISPLAY




; *******************
;     Test MAIN
; *******************
; --------------------
; Simulation geometry
; --------------------
;LGS Geometry
sep = 1.6 
nb_gs = 6.
; Analyse directions [arcmin]
aa = 2.*!pi/Nb_gs   ; etoiles en cercle
alpha  = transpose([[cos(indgen(nb_gs)*aa)],  [sin(indgen(nb_gs)*aa)]]*sep/2.)
; Optimization direction [arcmin]         
theta = [[0, 0], [0, 0]]/60. 
; Performance direction [arcmin]
beta_tab = [[0., 0.],  $
           [0., 50.],  $
           [0., 80.]]/60  
           
; NGS Geometry [arcmin]
ngs = [[0., 20.],  $
      [-10., -35.]]/60  
; --------------------
; Calculate display values from Simulation geometry
; --------------------
; Display field of view
sepchampvisu = sep/2.
; Number of performance direction
nb_dir_perf = (size(beta_tab))(2)
; Number of optimistion direction
nb_dir_theta = (size(theta))(2)

; --------------------
; DISPLAY SIMULATION GEOMETRY
; --------------------
display_geo, alpha, beta_tab, theta, $
            ; OPTIONAL INPUTS
            sepchampvisu = sepchampvisu, nb_gs = nb_gs, Nb_dir_perf = Nb_dir_perf, $
            nb_dir_theta = nb_dir_theta, window_number = window_number, $
            window_size = window_size, $ 
            ; OPTIONAL INPUT KEYWORD
            inv_color = inv_color   

alpha = [[-0.7, -0.5],  $
       [-0.7,  0.2],  $
       [0.3 ,  0.7],  $
       [0.4 , -0.7]]
display_geo, alpha, beta_tab, theta, ngs = ngs, window_number = 1, window_size = [650,650], /inv_color       

end 