;+
; NAME:
;     DISPLAY_SIMU
; PURPOSE:
;     This procedure displays in the form of a text file all simulation parameters.
; EXPLANATION:
;
; CALLING SEQUENCE:
;     DISPLAY_SIMU, SIM_PARAMS[, /HELP]
;
; INPUTS:
;     SIM_PARAMS - Structure containing all the data to be printed.
;
; OPTIONAL INPUTS:
;     None
;
; OPTIONAL INPUT KEYWORD:
;     HELP      - Print procedure ussage and returns.
;     
; OUTPUTS:
;     None
;
; EXAMPLES:
;     See 'MAIN' below
;
; TO DO:
;     - Add additionnal information such as meaning of parameters...
;     - NGS option
;     - Static Phase filename
;     - Pupil shape filename
;     - Filename of the file where simulation output is stored
; 
; SEE ALSO:
;   NBR2STR
;     
; REVISION HISTORY
;    09/10/2012 : Creation by Noah SCHWARTZ
;-



; *******************************************
;       RETURNS THE TYPE OF OA FLAVOR
; *******************************************
FUNCTION ao_flavors, sim_params

; Number of Laser Guide Stars
IF n_elements(size(sim_params.alpha)) EQ 5 THEN ind = 2 ELSE ind = 0   
nb_gs = (size(sim_params.alpha))[ind] 
; Number of performance Directions
IF n_elements(size(sim_params.beta_tab)) EQ 5 THEN ind = 2 ELSE ind = 0   
nb_dir_perf = (size(sim_params.beta_tab))[ind] 

; Number of Optimization Directions
IF n_elements(size(sim_params.theta)) EQ 5 THEN ind = 2 ELSE ind = 0   
nb_dir_theta = (size(sim_params.theta))[ind] 

flavor = '?????'

; -- Single-Conjugate Adaptive Optics case --
IF sim_params.nb_gs EQ 1 AND fix(n_elements(sim_params.h_dm)) EQ 1 AND $
   (size(sim_params.cn2recons))[1] EQ 1 AND nb_dir_theta EQ 1 THEN flavor = 'SCAO'
   
; -- Ground-Layer Adaptive Optics case --
IF sim_params.nb_gs GT 1 AND fix(n_elements(sim_params.h_dm)) EQ 1 AND $
   (size(sim_params.cn2recons))[1] EQ 1 THEN flavor = 'GLAO'

; -- Laser Tomographic Adaptive Optics case --
IF sim_params.nb_gs GT 1 AND fix(n_elements(sim_params.h_dm)) EQ 1 AND $
   (size(sim_params.cn2recons))[1] GT 1 AND nb_dir_theta EQ 1 THEN flavor = 'LTAO'

; -- Multi-Conjugate Adaptive Optics case --
IF sim_params.nb_gs GT 1 AND fix(n_elements(sim_params.h_dm)) GT 1 AND $
   (size(sim_params.cn2recons))[1] GT 1 THEN flavor = 'MCAO'


return, flavor
END 
; End ao_flavors


PRO base_destroy, ev
print, ev
END 


; *******************************************
;      EVENT: PLOT SIMULATION GEOMETRY
; *******************************************
PRO button_geo_event, ev
  WIDGET_CONTROL, ev.id, GET_UVALUE = uv
  IF ev.SELECT THEN display_geo, uv.alpha, uv.beta_tab, uv.theta, ngs = uv.ngs
END
; *******************************************
;    EVENT: PLOT TURBULENCE DISTRIBUTION
; *******************************************
PRO button_turb_event, ev
  WIDGET_CONTROL, ev.id, GET_UVALUE = uv
  IF ev.SELECT THEN BEGIN
      display_turb, uv.Cn2true, uv.h, uv.Cn2recons, $
                    uv.recons_h, uv.wind, uv.wind_dir, uv.h_dm, window_size = [750,650], thick = 100
                                  
;      ; Create widget
;      Table_DM = Widget_Base(ROW=3, Title='DM, Cn2 Values')
;      IF STRUPCASE(!VERSION.OS_FAMILY) EQ 'WINDOWS'  THEN fancyFont = 'Times*18*Italic*Bold' 
;      SCR_XSIZE = 250
;      
;      
;      ; Create DM Altitude Table
;      thisLabel = ' DM Altitude ' 
;      bulletinBoardBase = Widget_Base(Table_DM)
;      label = Widget_Label(bulletinBoardBase, VALUE=thisLabel, XOFFSET=5, FONT=fancyFont)
;      labelGeometry = Widget_Info(label, /GEOMETRY)
;      labelYSize =  labelGeometry.ysize
;      fancyBase = Widget_Base(bulletinBoardBase, /column, /FRAME, YOFFSET=labelYSize/2, YPAD=10, XPAD=10)
;      label = Widget_Label(bulletinBoardBase, VALUE=thisLabel, XOFFSET=5, FONT=fancyFont)
;      
;      IF n_elements(uv.h_dm) EQ 1 then value = [nbr2str(uv.h_dm)] ELSE value = uv.h_dm
;      Result = WIDGET_TABLE(fancyBase, /NO_ROW_header, /NO_COLUMN_header, VALUE = value)
;      
;      ; True Cn2
;      thisLabel = ' True Cn2 ' 
;      bulletinBoardBase = Widget_Base(Table_DM)
;      label = Widget_Label(bulletinBoardBase, VALUE=thisLabel, XOFFSET=5, FONT=fancyFont)
;      labelGeometry = Widget_Info(label, /GEOMETRY)
;      labelYSize =  labelGeometry.ysize
;      fancyBase = Widget_Base(bulletinBoardBase, /column, /FRAME, YOFFSET=labelYSize/2, YPAD=10, XPAD=10)
;      label = Widget_Label(bulletinBoardBase, VALUE=thisLabel, XOFFSET=5, FONT=fancyFont)
;      
;      IF n_elements(uv.Cn2true) EQ 1 then value = [[nbr2str(uv.Cn2true)], [nbr2str(uv.h)], [nbr2str(uv.wind)], [nbr2str(uv.wind_dir*!RADEG)]] $
;      ELSE value = [[uv.Cn2true], [nbr2str(uv.h)],[uv.wind], [nbr2str(uv.wind_dir*!RADEG)]]
;      Result = WIDGET_TABLE(fancyBase, /NO_COLUMN_header, ROW_labels=['Cn2 [%]','Altitude [m]', 'Wind Strength [m/s]','Wind Direction [deg]'], VALUE = value)
;      
;      ; Reconstructed Cn2 
;      thisLabel = ' Reconstructed Cn2 ' 
;      bulletinBoardBase = Widget_Base(Table_DM)
;      label = Widget_Label(bulletinBoardBase, VALUE=thisLabel, XOFFSET=5, FONT=fancyFont)
;      labelGeometry = Widget_Info(label, /GEOMETRY)
;      labelYSize =  labelGeometry.ysize
;      fancyBase = Widget_Base(bulletinBoardBase, /column, /FRAME, YOFFSET=labelYSize/2, YPAD=10, XPAD=10)
;      label = Widget_Label(bulletinBoardBase, VALUE=thisLabel, XOFFSET=5, FONT=fancyFont)
;      
;      IF n_elements(uv.Cn2recons) EQ 1 then value = [[nbr2str(uv.Cn2recons)], [nbr2str(uv.recons_h)]] ELSE value = [[uv.Cn2recons], [nbr2str(uv.recons_h)]]
;      Result = WIDGET_TABLE(fancyBase, /NO_COLUMN_header, ROW_labels=['Cn2 [%]','Altitude [m]'], VALUE = value)
;      
;      ; Realise widget
;      Widget_Control, Table_DM, /REALIZE
  ENDIF 
END 
; End button_geo_event



; ===================================================
;  MAIN PROCEDURE: DISPLAY SIMULATION PARAMS
; ===================================================
PRO display_simu, sim_params, help = help

IF keyword_set(help) THEN BEGIN
     print,'Syntax - '
     print,'DISPLAY_SIMU, SIM_PARAMS[, /HELP]
     print, ' '
     print,'=============================== To create the SIM_PARAMS structure  ==============================='
     print, "im_params = create_struct('anglezenith', anglezenith, 'sep', sep, 'nb_gs', nb_gs, $"
     print, "                          'LGS_H', LGS_H, 'varnoise_wfs', varnoise_wfs, 'lambda_seeing', lambda_seeing, $"
     print, "                          'lambda', lambda, 'L0', L0, 'seeing0', seeing0, 'cn2true', cn2true, $"
     print, "                          'vent', vent, 'arg_v', arg_v, 'h', h, 'D', D, 'dimpup', dimpup, $"
     print, "                          'dim', dim, 'nb_act', nb_act, 'pitch_m', pitch_m,$"
     print, "                          'fech', fech, 'td', td, 'h_dm', h_dm, 'pitchs_dm', pitchs_dm, $"
     print, "                          'alpha', alpha, 'tempo', tempo, 'fitting', fitting, $"
     print, "                          'recons_h', recons_h, 'cn2recons', cn2recons, 'theta', theta, $"
     print, "                          'beta_tab', beta_tab, 'LSE', LSE, 'alias', alias, $"
     print, "                          'optim', optim, 'tomo', tomo)"
     print, "==================================================================================================="
     RETURN
 ENDIF


; ---------------------------
; Copy Structure to parameters
; ---------------------------
anglezenith = sim_params.anglezenith & sep = sim_params.sep & nb_gs = sim_params.nb_gs
LGS_H = sim_params.LGS_H & varnoise_wfs = sim_params.varnoise_wfs & lambda_seeing = sim_params.lambda_seeing
lambda_imaging = sim_params.lambda_imaging & L0 = sim_params.L0 & seeing0 = sim_params.seeing0
cn2true = sim_params.cn2true & vent = sim_params.vent & arg_v = sim_params.arg_v
h = sim_params.h & D = sim_params.D & dimpup = sim_params.dimpup & dim = sim_params.dim
nb_act = sim_params.nb_act & pitch_m = sim_params.pitch_m & fech = sim_params.fech
td = sim_params.td & h_dm = sim_params.h_dm & pitchs_dm = sim_params.pitchs_dm
alpha = sim_params.alpha & tempo = sim_params.tempo & fitting_kw = sim_params.fitting_kw
recons_h = sim_params.recons_h & cn2recons = sim_params.cn2recons & theta = sim_params.theta
beta_tab = sim_params.beta_tab & LSE = sim_params.LSE & aliasing_kw = sim_params.aliasing_kw
tomo = sim_params.tomo & ngs = sim_params.ngs


; ---------------------------
; Create Base Widget
; ---------------------------
base = Widget_Base(COLUMN=3, Title='SIMULATION PARAMETERS SUMMARY')
; Fancy Label
IF STRUPCASE(!VERSION.OS_FAMILY) EQ 'WINDOWS'   THEN fancyFont = 'Times*18*Italic*Bold' 
IF STRUPCASE(!VERSION.OS_FAMILY) EQ 'WINDOWS'   THEN TextFont = 'Arial*14*Bold' 
IF STRUPCASE(!VERSION.OS_FAMILY) EQ 'WINDOWS'   THEN ErrorText = 'Times*14*Italic*Bold'


; ---------------------------
; --- GEOMETRY --
; ---------------------------
; Create Base for turbulence
thisLabel = ' Geometry ' 
bulletinBoardBase = Widget_Base(base)
label = Widget_Label(bulletinBoardBase, VALUE=thisLabel, XOFFSET=5, FONT=fancyFont)
labelGeometry = Widget_Info(label, /GEOMETRY)
labelYSize =  labelGeometry.ysize
fancyBase = Widget_Base(bulletinBoardBase, /column, /FRAME, YOFFSET=labelYSize/2, YPAD=10, XPAD=10)
label = Widget_Label(bulletinBoardBase, VALUE=thisLabel, XOFFSET=5, FONT=fancyFont)

; Write information
SCR_XSIZE1 = 105
SCR_XSIZE2 = 35
SCR_XSIZE3 = 40
; - Analysis Direction, LGS, alpha -
row0 = widget_base(fancyBase, /ROW)
text = Widget_Label(row0, Value='Analysis (LGS) Dir.', FONT=TextFont, SCR_XSIZE=SCR_XSIZE1)
IF n_elements(size(alpha)) EQ 5 THEN ind = 2 ELSE ind = 0  
text = Widget_Text(row0, Value=nbr2str((size(alpha))[ind]), SCR_XSIZE=SCR_XSIZE2*2)
; - Altitude LGS -
row0 = widget_base(fancyBase, /ROW)
label = widget_label(row0, VALUE = 'Altitude LGS', FONT=TextFont, SCR_XSIZE=SCR_XSIZE1) 
text = Widget_Text(row0, Value=nbr2str(LGS_H/1000.), SCR_XSIZE=SCR_XSIZE2*2)
text = widget_label(row0, Value='[km]', SCR_XSIZE=SCR_XSIZE3)
; - Separation LGS -
row0 = widget_base(fancyBase, /ROW)
label = widget_label(row0, VALUE = 'LGS Separation', FONT=TextFont, SCR_XSIZE=SCR_XSIZE1) 
text = Widget_Text(row0, Value=nbr2str(sep), SCR_XSIZE=SCR_XSIZE2*2)
text = widget_label(row0, Value='[arcmin]', SCR_XSIZE=SCR_XSIZE3)
; - NGS -
row0 = widget_base(fancyBase, /ROW)
label = widget_label(row0, VALUE = 'Nb NGS', FONT=TextFont, SCR_XSIZE=SCR_XSIZE1)
IF n_elements(size(ngs)) EQ 5 THEN ind = 2 ELSE ind = 0
text = Widget_Text(row0, Value=nbr2str((size(ngs))[ind]), SCR_XSIZE=SCR_XSIZE2*2)
; - Optimisation Direction, theta -
row0 = widget_base(fancyBase, /ROW)
text = Widget_Label(row0, Value='Optimisation Dir.', FONT=TextFont, SCR_XSIZE=SCR_XSIZE1)
IF n_elements(size(theta)) EQ 5 THEN ind = 2 ELSE ind = 0  
text = Widget_Text(row0, Value=nbr2str((size(theta))[ind]), SCR_XSIZE=SCR_XSIZE2*2)
; - Performance Direction, beta -
row0 = widget_base(fancyBase, /ROW) 
label = widget_label(row0, VALUE = 'Performance Dir.', FONT=TextFont, SCR_XSIZE=SCR_XSIZE1) 
IF n_elements(size(beta)) EQ 5 THEN ind = 2 ELSE ind = 0  
text = Widget_Text(row0, Value=nbr2str((size(beta))[ind]), SCR_XSIZE=SCR_XSIZE2*2)
; - Observation angle relative to zenith -
row0 = widget_base(fancyBase, /ROW) 
label = widget_label(row0, VALUE = 'Observation Angle', FONT=TextFont, SCR_XSIZE=SCR_XSIZE1) 
text = Widget_Text(row0, Value=nbr2str(anglezenith), SCR_XSIZE=SCR_XSIZE2*2)
text = widget_label(row0, Value='[arcmin]', SCR_XSIZE=SCR_XSIZE3)
; - AO Flavours -
row0 = widget_base(fancyBase, /ROW) 
label = widget_label(row0, VALUE = 'AO Flavour', FONT=TextFont, SCR_XSIZE=SCR_XSIZE1) 
text = Widget_Text(row0, Value=ao_flavors(sim_params), SCR_XSIZE=SCR_XSIZE2*2)
; - Plot Simulation Geometry -
row0 = widget_base(fancyBase, /ROW)  
uvalue = create_struct('alpha', alpha, 'beta_tab', beta_tab, 'theta', theta, 'ngs', ngs)
button_geo = widget_button(row0, VALUE = 'Plot Simulation Geometry', SCR_XSIZE=SCR_XSIZE1+SCR_XSIZE2+SCR_XSIZE3*2, uvalue = uvalue) 
XMANAGER, 'button_geo', button_geo



; -----------------------------
; --- WARNING --
; -----------------------------
bulletinBoardBase = Widget_Base(base)
thisLabel = ' Warning ' 
label = Widget_Label(bulletinBoardBase, VALUE=thisLabel, XOFFSET=5, FONT=fancyFont)
labelGeometry = Widget_Info(label, /GEOMETRY)
labelYSize =  labelGeometry.ysize
fancyBase = Widget_Base(bulletinBoardBase, COLUMN=1, /FRAME, YOFFSET=labelYSize/2, YPAD=10, XPAD=11)
label = Widget_Label(bulletinBoardBase, VALUE=thisLabel, XOFFSET=5, FONT=fancyFont)

row0 = widget_base(fancyBase, /ROW) 
label = widget_label(row0, VALUE = 'LGS Separation', FONT=TextFont, SCR_XSIZE=SCR_XSIZE1) 
IF sep GT tan(D/LGS_H)*!RADEG*60 THEN val='Ok' ELSE val='Increase LGS Separation'
text = Widget_Text(row0, Value=Val, SCR_XSIZE=SCR_XSIZE3+SCR_XSIZE2*2, FONT=ErrorText)

row0 = widget_base(fancyBase, /ROW) 
label = widget_label(row0, VALUE = 'AO Geometry Error', FONT=TextFont, SCR_XSIZE=SCR_XSIZE1)
IF STRMATCH(ao_flavors(sim_params), '?????', /FOLD_CASE) EQ 1 THEN val = 'Check AO Geometry' ELSE val = 'Ok'
text = Widget_Text(row0, Value=Val, SCR_XSIZE=SCR_XSIZE3+SCR_XSIZE2*2, FONT=ErrorText)

row0 = widget_base(fancyBase, /ROW) 
label = widget_label(row0, VALUE = 'Other Errors', FONT=TextFont, SCR_XSIZE=SCR_XSIZE1)
text = Widget_Text(row0, Value='TBD', SCR_XSIZE=SCR_XSIZE3+SCR_XSIZE2*2, FONT=ErrorText)


; -----------------------------
;   --- Adaptive Optics ---
; -----------------------------
bulletinBoardBase = Widget_Base(base)
thisLabel = ' Adaptive Optics ' 
label = Widget_Label(bulletinBoardBase, VALUE=thisLabel, XOFFSET=5, FONT=fancyFont)
labelGeometry = Widget_Info(label, /GEOMETRY)
labelYSize =  labelGeometry.ysize
fancyBase = Widget_Base(bulletinBoardBase, COLUMN=1, /FRAME, YOFFSET=labelYSize/2, YPAD=14, XPAD=10)
label = Widget_Label(bulletinBoardBase, VALUE=thisLabel, XOFFSET=5, FONT=fancyFont)


; - Number of DM -
row0 = widget_base(fancyBase, /ROW) 
label = widget_label(row0, VALUE = 'Nb DM', FONT=TextFont, SCR_XSIZE=SCR_XSIZE1)
text = Widget_Text(row0, Value=nbr2str(fix(n_elements(h_dm))), SCR_XSIZE=SCR_XSIZE2*2)

; - Number of actuators per DM -
row0 = widget_base(fancyBase, /ROW) 
label = widget_label(row0, VALUE = 'Nb of Actuators', FONT=TextFont, SCR_XSIZE=SCR_XSIZE1)  
nb_act_str =''
FOR ind=0, n_elements(nb_act)-1 DO BEGIN
    IF ind NE n_elements(nb_act)-1 THEN nb_act_str += nbr2str(nb_act[ind])+'; ' $
    ELSE nb_act_str += nbr2str(nb_act[ind]) 
ENDFOR
text = Widget_Text(row0, Value=nb_act_str, SCR_XSIZE=SCR_XSIZE2*2)

; - DM altitudes [m] -
row0 = widget_base(fancyBase, /ROW) 
label = widget_label(row0, VALUE = 'DM Altitudes', FONT=TextFont, SCR_XSIZE=SCR_XSIZE1) 
h_DM_str =''
FOR ind=0, n_elements(h_DM)-1 DO BEGIN
    IF ind NE n_elements(h_DM)-1 THEN BEGIN
        tmp = nbr2str(h_DM[ind]/1000.)
        IF STRMATCH(tmp, '*.00', /FOLD_CASE) EQ 0 THEN h_DM_str += tmp+'; ' $
        ELSE h_DM_str += nbr2str(h_DM[ind]/1000., FORMAT='(I)')+'; '
    ENDIF ELSE BEGIN
        h_DM_str += nbr2str(h_DM[ind]/1000.)
    ENDELSE 
ENDFOR
text = Widget_Text(row0, Value=h_DM_str, SCR_XSIZE=SCR_XSIZE2*2)
text = widget_label(row0, Value='[km]', SCR_XSIZE=SCR_XSIZE3/2)

; - DM pitch [m] -
row0 = widget_base(fancyBase, /ROW) 
label = widget_label(row0, VALUE = 'DM Pitch', FONT=TextFont, SCR_XSIZE=SCR_XSIZE1) 
pitchs_dm_str =''
FOR ind=0, n_elements(pitchs_dm)-1 DO BEGIN
    IF ind NE n_elements(pitchs_dm)-1 THEN pitchs_dm_str += nbr2str(fix(pitchs_dm[ind]))+'; ' $
    ELSE pitchs_dm_str += nbr2str(pitchs_dm[ind]) 
ENDFOR
text = Widget_Text(row0, Value=pitchs_dm_str, SCR_XSIZE=SCR_XSIZE2*2)
text = widget_label(row0, Value='[m]', SCR_XSIZE=SCR_XSIZE3)

; - AO temporal sampling frequency [Hz] -
row0 = widget_base(fancyBase, /ROW) 
label = widget_label(row0, VALUE = 'AO Sampling Freq.', FONT=TextFont, SCR_XSIZE=SCR_XSIZE1) 
text = Widget_Text(row0, Value=nbr2str(fech, FORMAT='(F25.0)'), SCR_XSIZE=SCR_XSIZE2*2)
text = widget_label(row0, Value='[Hz]', SCR_XSIZE=SCR_XSIZE3)
; - Loop latency [s] -
row0 = widget_base(fancyBase, /ROW) 
label = widget_label(row0, VALUE = 'Loop Latency', FONT=TextFont, SCR_XSIZE=SCR_XSIZE1) 
text = Widget_Text(row0, Value=nbr2str(td*1000.), SCR_XSIZE=SCR_XSIZE2*2)
text = widget_label(row0, Value='[ms]', SCR_XSIZE=SCR_XSIZE3)

; - WFS pitchs [m] -
row0 = widget_base(fancyBase, /ROW) 
label = widget_label(row0, VALUE = 'WFS Pitch', FONT=TextFont, SCR_XSIZE=SCR_XSIZE1) 
pitch_m_str =''
FOR ind=0, n_elements(pitch_m)-1 DO BEGIN
    IF ind NE n_elements(pitch_m)-1 THEN pitch_m_str += nbr2str(fix(pitch_m[ind]))+'-' $
    ELSE pitch_m_str += nbr2str(pitch_m[ind]) 
ENDFOR
text = Widget_Text(row0, Value=pitch_m_str, SCR_XSIZE=SCR_XSIZE2*2)
text = widget_label(row0, Value='[m]', SCR_XSIZE=SCR_XSIZE3)

; - analyse noise in [rad^2] -
row0 = widget_base(fancyBase, /ROW) 
label = widget_label(row0, VALUE = 'WFS Noise', FONT=TextFont, SCR_XSIZE=SCR_XSIZE1) 
text = Widget_Text(row0, Value=nbr2str(varnoise_wfs), SCR_XSIZE=SCR_XSIZE2*2)
text = widget_label(row0, Value='[rad^2]', SCR_XSIZE=SCR_XSIZE3)

; - Telescope Diametre -
row0 = widget_base(fancyBase, /ROW) 
label = widget_label(row0, VALUE = '_________________', FONT=TextFont, SCR_XSIZE=SCR_XSIZE1) 
row0 = widget_base(fancyBase, /ROW) 
label = widget_label(row0, VALUE = 'Telescope Diam.', FONT=TextFont, SCR_XSIZE=SCR_XSIZE1) 
text = Widget_Text(row0, Value=nbr2str(D), SCR_XSIZE=SCR_XSIZE2*2)
text = widget_label(row0, Value='[m]', SCR_XSIZE=SCR_XSIZE3)




; -----------------------------
;       -- Photometry --
; -----------------------------
bulletinBoardBase = Widget_Base(base)
thisLabel = ' Photometry ' 
label = Widget_Label(bulletinBoardBase, VALUE=thisLabel, XOFFSET=5, FONT=fancyFont)
labelGeometry = Widget_Info(label, /GEOMETRY)
labelYSize =  labelGeometry.ysize
fancyBase = Widget_Base(bulletinBoardBase, COLUMN=1, /FRAME, YOFFSET=labelYSize/2, YPAD=10, XPAD=10)
label = Widget_Label(bulletinBoardBase, VALUE=thisLabel, XOFFSET=5, FONT=fancyFont)

; - lambda - 
row0 = widget_base(fancyBase, /ROW) 
label = widget_label(row0, VALUE = 'Imaging Lambda', FONT=TextFont, SCR_XSIZE=SCR_XSIZE1) 
text = Widget_Text(row0, Value=nbr2str(lambda_imaging*1.e6), SCR_XSIZE=SCR_XSIZE2*2)
text = widget_label(row0, Value='[um]', SCR_XSIZE=SCR_XSIZE3)
; - r0 wavelength- 
row0 = widget_base(fancyBase, /ROW) 
label = widget_label(row0, VALUE = 'r0 wavelength', FONT=TextFont, SCR_XSIZE=SCR_XSIZE1) 
text = Widget_Text(row0, Value=nbr2str(lambda_seeing*1.e6), SCR_XSIZE=SCR_XSIZE2*2)
text = widget_label(row0, Value='[um]', SCR_XSIZE=SCR_XSIZE3)


; -----------------------------
;         -- TURBULENCE --
; -----------------------------
bulletinBoardBase = Widget_Base(base)
thisLabel = ' Turbulence ' 
label = Widget_Label(bulletinBoardBase, VALUE=thisLabel, XOFFSET=5, FONT=fancyFont)
labelGeometry = Widget_Info(label, /GEOMETRY)
labelYSize =  labelGeometry.ysize
fancyBase = Widget_Base(bulletinBoardBase, COLUMN=1, /FRAME, YOFFSET=labelYSize/2, YPAD=10, XPAD=10)
label = Widget_Label(bulletinBoardBase, VALUE=thisLabel, XOFFSET=5, FONT=fancyFont)

; - L0 - 
row0 = widget_base(fancyBase, /ROW) 
label = widget_label(row0, VALUE = 'L0', FONT=TextFont, SCR_XSIZE=SCR_XSIZE1) 
text = Widget_Text(row0, Value=nbr2str(L0), SCR_XSIZE=SCR_XSIZE2*2)
text = widget_label(row0, Value='[m]', SCR_XSIZE=SCR_XSIZE3)
; - seeing0 - 
row0 = widget_base(fancyBase, /ROW) 
label = widget_label(row0, VALUE = 'Seing at Zenith', FONT=TextFont, SCR_XSIZE=SCR_XSIZE1) 
text = Widget_Text(row0, Value=nbr2str(seeing0), SCR_XSIZE=SCR_XSIZE2*2)
text = widget_label(row0, Value='[arcsec]', SCR_XSIZE=SCR_XSIZE3)
; - Number of True layers -
row0 = widget_base(fancyBase, /ROW) 
label = widget_label(row0, VALUE = 'Nb True Layers', FONT=TextFont, SCR_XSIZE=SCR_XSIZE1) 
text = Widget_Text(row0, Value=nbr2str((size(cn2true))[1]), SCR_XSIZE=SCR_XSIZE2*2)
; - Number of Reconstructed Layers -
row0 = widget_base(fancyBase, /ROW) 
label = widget_label(row0, VALUE = 'Nb Recons Layers', FONT=TextFont, SCR_XSIZE=SCR_XSIZE1) 
text = Widget_Text(row0, Value=nbr2str((size(cn2recons))[1]), SCR_XSIZE=SCR_XSIZE2*2)
;  
; - Plot Turbulence Distribution -
row0 = widget_base(fancyBase, /ROW)  
uvalue = create_struct('Cn2true', Cn2true, 'h', h, 'Cn2recons', Cn2recons, $
                       'recons_h', recons_h, 'wind', vent, 'wind_dir', arg_v, 'h_dm', h_dm)
button_turb = widget_button(row0, VALUE = 'Plot Turbulence Distribution', SCR_XSIZE=SCR_XSIZE1+SCR_XSIZE2+SCR_XSIZE3*2, uvalue = uvalue) 
XMANAGER, 'button_turb', button_turb


; -----------------------------
;         -- SIMULATION --
; -----------------------------
bulletinBoardBase = Widget_Base(base)
thisLabel = ' Simulation Parameters ' 
label = Widget_Label(bulletinBoardBase, VALUE=thisLabel, XOFFSET=5, FONT=fancyFont)
labelGeometry = Widget_Info(label, /GEOMETRY)
labelYSize =  labelGeometry.ysize
fancyBase = Widget_Base(bulletinBoardBase, COLUMN=1, /FRAME, YOFFSET=labelYSize/2, YPAD=10, XPAD=10)
label = Widget_Label(bulletinBoardBase, VALUE=thisLabel, XOFFSET=5, FONT=fancyFont)

; - Dimpup - 
row0 = widget_base(fancyBase, /ROW) 
label = widget_label(row0, VALUE = 'DIMPUP', FONT=TextFont, SCR_XSIZE=SCR_XSIZE1) 
text = Widget_Text(row0, Value=nbr2str(DIMPUP), SCR_XSIZE=SCR_XSIZE2*2)
text = widget_label(row0, Value='[pix]', SCR_XSIZE=SCR_XSIZE3)
; - Dim [pix] -
row0 = widget_base(fancyBase, /ROW) 
label = widget_label(row0, VALUE = 'DIM', FONT=TextFont, SCR_XSIZE=SCR_XSIZE1) 
text = Widget_Text(row0, Value=nbr2str(Dim), SCR_XSIZE=SCR_XSIZE2*2)
text = widget_label(row0, Value='[pix]', SCR_XSIZE=SCR_XSIZE3)
; - Kwrd tempo - 
row0 = widget_base(fancyBase, /ROW) 
label = widget_label(row0, VALUE = 'TEMPO Kwrd', FONT=TextFont, SCR_XSIZE=SCR_XSIZE1) 
text = Widget_Text(row0, Value=nbr2str(tempo, FORMAT='(B2.0)'), SCR_XSIZE=SCR_XSIZE2*2)
; - Kwrd fitting -
row0 = widget_base(fancyBase, /ROW) 
label = widget_label(row0, VALUE = 'FITTING Kwrd', FONT=TextFont, SCR_XSIZE=SCR_XSIZE1) 
text = Widget_Text(row0, Value=nbr2str(Fitting_kw, FORMAT='(B2.0)'), SCR_XSIZE=SCR_XSIZE2*2)
; - Kwrd LSE -
row0 = widget_base(fancyBase, /ROW) 
label = widget_label(row0, VALUE = 'LSE Kwrd', FONT=TextFont, SCR_XSIZE=SCR_XSIZE1) 
text = Widget_Text(row0, Value=nbr2str(LSE, FORMAT='(B2.0)'), SCR_XSIZE=SCR_XSIZE2*2)
; - Kwrd Alias - 
row0 = widget_base(fancyBase, /ROW) 
label = widget_label(row0, VALUE = 'ALIAS Kwrd', FONT=TextFont, SCR_XSIZE=SCR_XSIZE1) 
text = Widget_Text(row0, Value=nbr2str(aliasing_kw), SCR_XSIZE=SCR_XSIZE2*2)
; - Kwrd tomo - 
row0 = widget_base(fancyBase, /ROW) 
label = widget_label(row0, VALUE = 'TOMO Kwrd', FONT=TextFont, SCR_XSIZE=SCR_XSIZE1) 
text = Widget_Text(row0, Value=nbr2str(tomo, FORMAT='(B2.0)'), SCR_XSIZE=SCR_XSIZE2*2)



; ---------------------------
; Realise Base
; ---------------------------
Widget_Control, base, /REALIZE, /KILL_NOTIFY
widget_control, base, redraw = 1
xmanager, 'Widget_base', base, /no_block



; AO configuration (avec le code de Jeff)
; 
; -- ALIASING --
; No -> No really useful
; Yes -> Long to simulate
; Corr -> Very long to simulate and doesn't bring much more than Yes
; Simplified -> Fast, less accurate than Corr (to be implemented)
;  
; WARNINGS (separation des etoiles laser)

END
; END PROCEDURE




; *******************
;     Test MAIN
; *******************
; Simulation geometry
anglezenith = 0               ; observation angle relative to zenith 
sep = 1.65                     ; arcmin Careful > 1.6 arcmin to ensure NO
                              ; unseen area
nb_gs = 6.
LGS_H = 9e4                   ; altitude for LGS. 

; Simulation photometry
lambda_wfs = 0.589e-6         ; wfs wavelength [m]
lambda_seeing = 0.5e-6        ; r0 wavelength  [m]
lambda_imaging =  1.65e-6     ; imaging wavelength [m]
varnoise_wfs = (fltarr(nb_gs)+1)*(lambda_wfs/lambda_imaging)^2 ;analyse noise in rad^2

; Simulation turbulence definition
L0     = 25.                  ; outer scale    [m]; TEM
seeing0 = 0.85                ; zenithal seeing (arcsec)
cn2true =  [0.5, 0.5]         ; cn2 profile 
vent = [10, 10]               ; wind speed modulous [m/s]
arg_v  = [0, !pi/2.]          ; wind dir.  [rad]
h = [0., 10.e3]                 ; layers altitude [m]

; Simulation parameters
D = 42.                        ; Telescope diameter [m]
dimpup = 168.                 ; Pupil size [px]
dim = 512.                    ; Screen size [px] 

; AO system definition
nb_act = 84                  ; number of act
pitch_m =  D/(nb_act-1)       ; WFS pitch [m]
fech = 500.                   ; AO temporal sampling frequency [Hz]
td = 3*1.e-3                  ; loop latency [s]
h_dm = [1000]                    ; DM altitude
pitchs_dm = [pitch_m]         ; DM pitch [m]
aa = 2.*!pi/Nb_gs             ; etoiles en cercle
alpha  = transpose([[cos(indgen(nb_gs)*aa)],  [sin(indgen(nb_gs)*aa)]] $
                   *sep/2.)   ; analyse directions [arcmin]
tempo = 0B                    ; temporal error 
fitting = 1B                  ; fitting error 

; Reconstruction turbulence a priori 
recons_h = [0]                ; supposed altitudes for estimation 
cn2recons = [1]               ; supposed CN2 for estimation

theta = [[0, 0], [50, 50]]/60. ; optimization direction [arcmin]
beta_tab = [[0., 0.],  $
;             [0., 10.],  $
;             [0., 20.],  $
;             [0., 30.],  $
;             [0., 40.],  $
             [0., 50.],  $
            [0., 60.]]/60       ; performance direction [arcmin]

; Reconstruction method
LSE = 0B                      ; LSE = 1B for LSE reconstruction
alias = 'no'                ; Keyword for aliasing
optim = 0                     ; ???
tomo = 0B                     ; Ideal reconstruction h_dm = h_recons


; Output
visu =  1
file = 'D:\Simulation Codes\IDL\ELTF\Data\file_save_oa_8m.dat'        ; save file string

; fitting
; noise
; tempo
; aliasing 
; aliasing g�n�ralis

; NGS Geometry [arcmin]
ngs = [[0., 20.],  $
      [-10., -35.]]/60  
      
sim_params = create_struct('anglezenith', anglezenith, 'sep', sep, 'nb_gs', nb_gs, $
                           'LGS_H', LGS_H, 'varnoise_wfs', varnoise_wfs, 'lambda_seeing', lambda_seeing, $
                           'lambda_imaging', lambda_imaging, 'L0', L0, 'seeing0', seeing0, 'cn2true', cn2true, $
                           'vent', vent, 'arg_v', arg_v, 'h', h, 'D', D, 'dimpup', dimpup, $
                           'dim', dim, 'nb_act', nb_act, 'pitch_m', pitch_m,$
                           'fech', fech, 'td', td, 'h_dm', h_dm, 'pitchs_dm', pitchs_dm, $
                           'alpha', alpha, 'tempo', tempo, 'fitting', fitting, $
                           'recons_h', recons_h, 'cn2recons', cn2recons, 'theta', theta, $
                           'beta_tab', beta_tab, 'LSE', LSE, 'alias', alias, $
                           'optim', optim, 'tomo', tomo, 'ngs', ngs)


display_simu, sim_params
end
  
            
            
