PRO main_psf_harmoni, filename = filename, SR = SR, EE = EE, savefile = savefile, $
                      sim_params = sim_params, longexposure = longexposure, new_lambda = new_lambda

defsysv, '!I', complex(0,1), 1

; Program simulating AO corrected PSF with different options
; uses the output of code FOURIER_AO

; opens FOURIER_AO output
; FOURIER_AO files contain
; RESALL       : fltarr(nb_dir_perf)  containing residual for each performance direction [nm]
; DSPALL       : fltarr(nb_f, nb_f, nb_dir_perf) containing 2D DSP for each perf direction
; F            : fltarr(nb_f, nb_f) containing frequency array
; DSP_TAB_VRAI : fltarr(nb_f, nb_f, nb_L) containing 2D DSP of turbulent layers
; SIM_PARAMS   : structure containing all simulation parameters
; NEW_LAMBDA   : Wavelength where PSF is evaluated (DSP SCALING DSP2 = DSP1 *(lambda1/lambda2)^2)


; --- Get File name and restore data ---
;dir = 'D:\Simulation Codes\IDL\ELTF2\FOURIER_AO\'
;file = 'file_save_oa_harmoni'
;restore, file = dir + file + '.dat'
IF NOT KEYWORD_SET(savefile) THEN savefile = 0
IF NOT KEYWORD_SET(filename) THEN BEGIN
    PATH = 'E:\HARMONI Simulations\ONERA-UKATC Simulation Verification\'
    PATH = 'C:\Users\ukatcuser\Desktop\IDL\HarmoniSimulations\'
    filename = DIALOG_PICKFILE(GET_PATH = dir, PATH = PATH, DEFAULT_EXTENSION='dat', /multi)
    IF filename[0] EQ '' THEN BEGIN
        print, 'No file selected. Returning'
        return
    ENDIF
ENDIF

;-------------------------------------
Nnew = 512.;1040.
psf_total = fltarr(Nnew,Nnew,N_ELEMENTS(filename)) ; add the 8th of sept
;airy_total = fltarr(Nnew,Nnew,N_ELEMENTS(filename)) ; add the 8th of sept
;-------------------------------------

FOR ii=0, N_ELEMENTS(filename)-1 DO BEGIN 
    restore, file = filename[ii]
    file = FILE_BASENAME(filename[ii])
    file = Strmid(file, 0, strlen(file)-4)
        
        
    ;****************  NEW 03 April 2014 ***************************
    DSP_ALL_OLD = dspall
    lambda_imaging_tab = 2.2e-6; [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5] * 1.0e-6
    psf_final = fltarr(Nnew,Nnew, N_elements(lambda_imaging_tab))
    airy_final = psf_final
    FOR iiii = 0, N_elements(lambda_imaging_tab)-1 DO BEGIN
    new_lambda = lambda_imaging_tab[iiii]
    print, new_lambda

    IF NOT KEYWORD_SET(new_lambda) THEN new_lambda = sim_params.lambda_imaging
    dspall  = DSP_ALL_OLD * (sim_params.lambda_imaging/new_lambda)^2.0
    ;***************************************************************
    
    nb_dirperf = n_elements(sim_params.beta_tab)/2.
    rad2arcsec = 1.*180./!dpi*3600.
    
    
    ; --- Sets parameters for DSP2PHASE simulation --- 
    seed = 1234L                  ; Use a given seed for random generation of
                                  ; static phase
    pupil_kw = 2;1;               ; Sets a round pupil mask in simulation (1) or
                                  ; not (0), or give a user map (2), then the map
                                  ; MASKPUPIL has to be provided by the user
    kw_static = 2B;0B;            ; Adds a static phase in simulation (1B) or not
                                  ; (0B). (2B) for User defined static phase
    rms_static_phase = 0.0        ; Sets the RMS value of static phase [nm] to be
                                  ; added 
    TTF_kw = 0                    ; KW for filtering TT (1) or TTF (2) or nothing (0)
                                  ; in the turbulent phase, also used for keeping
                                  ; TT only (3) or TTF only (4)
    short = 0                     ; KW for short (1) or long (0) exposure computation
    verbose = 1B                  ; prints each step information
    
    ; Sets parameters for ADD_JITTER simulation 
    jitter_kw = 2;0;              ; additionnal jitter in long exposure images (2)
                                  ; in turbulent phase (1) or not (0)
    jitter_mas = [3.0, 3.0]/2.;[1.0, 1.65]      ; jitter value [X, Y] [mas]
    
    ; sets parameters for MAIN simulation 
    number_exposure = 3000.         ; number of short exposures to be averaged
      
    ; --- Start of IDL code ---
    print, '****************************************************'
    
    ; USER-DEFINED PUPIL MASK - Values in [m rms]
    IF (kw_static EQ 2) OR (pupil_kw EQ 2) THEN BEGIN
        ;;;im = readfits('D:\Harmoni-Atlas\Simulations\CombinedError_WavefrontAfter.fits')
        im = readfits('C:\Users\ukatcuser\Desktop\IDL\ELTF2\CombinedError_WavefrontAfter.fits')
        ; Get Not-A-Number and Number indexes
        index_nan = where(im NE im)
        index_n   = where(im EQ im)
        ; Remove NAN and replace by 0
        im[index_nan] = 0
            
        ; --- User defined static aberrations ---
        IF (kw_static EQ 2) THEN BEGIN
            ; Extract useful area
            user_static_phase = im[8:1040-9, 8:1040-9]
            ; Extend support to dimpup
            user_static_phase = CONGRID(user_static_phase, sim_params.dimpup, sim_params.dimpup)
            ; Zero padding for sampling, convert from [m] to [nm]
            user_static_phase = crop(user_static_phase, ncrop=sim_params.dim, /centre, /extend, /silent)* 1.0e9
        ENDIF
        
        ; --- User defined static aberrations ---
        IF (pupil_kw EQ 2) THEN BEGIN
            ; Binarize
            im[index_n] = 1
            ; Extract useful area
            maskpupil = im[8:1040-9, 8:1040-9]
            ; Extend support to dimpup
            maskpupil = CONGRID(maskpupil, sim_params.dimpup, sim_params.dimpup)
            maskpupil = crop(maskpupil, ncrop=sim_params.dim, /centre, /extend, /silent)
        END 
        undefine, im
    END 

              
    ; sets pupil mask
    ; if KW = 0, mask is 1 everywhere
    IF pupil_kw EQ 0 THEN BEGIN
       IF keyword_set(verbose) THEN print, '    PUPIL: (NO)'
       maskpupil = fltarr(sim_params.dim, sim_params.dim)+1.
    ENDIF
    
    ; if KW = 1, mask is 1 inside a round pupil
    IF pupil_kw EQ 1 THEN begin
       IF keyword_set(verbose) THEN print, '    PUPIL: (YES) CIRCULAR PUPIL'
       polaire2, rt = sim_params.dimpup/2., largeur = sim_params.dim , masque = maskpupil, $
                 oc = 0, prolixe = 0
       maskpupil = float(maskpupil)
    ENDIF
    
    ; checks pupil parameters if user map triggered
    ; if KW = 2
    IF pupil_kw EQ 2 THEN BEGIN
        IF verbose EQ 1 THEN print, '   PUPIL: (YES) USER DEFINED'
;        ; Create User defined pupil ==> Other option for Caculation
;        polaire2, rt = sim_params.dimpup/2., largeur = sim_params.dim , masque = maskpupil, oc = 11.208/36.903, prolixe = 0
;        size_spider = sim_params.dimpup/36.903*0.5 ; Size of Spider in pixels (0.5m thickness) 
;        tmp = maskpupil*0.0+1.0
;        tmp[*, sim_params.dim/2 - round(size_spider/2.):sim_params.dim/2 + round(size_spider/2.)] = 0.0 
;        maskpupil = tmp * ROT(tmp, 60) * ROT(tmp, -60) * maskpupil
;        undefine, tmp
       
       ; check array undefined
       IF NOT(keyword_set(maskpupil)) THEN BEGIN
          print, 'PUPIL USER MAP UNDEFINED --> STOP'
          return
       ENDIF
       ; check wrong array size
       IF n_elements(maskpupil) NE (sim_params.dim)^2. THEN begin
          print, 'USER PUPIL MAP UNCORRECT SIZE : ' + strc((size(maskpupil))[1]) $
                 + ' // Requested ' + strc(sim_params.dim) + ' --> STOP'
                 stop
          return
       ENDIF
       ; check wrong pupil size
       IF abs(total(maskpupil)/(!pi*(sim_params.dimpup/2.)^2)-1) GT 0.3 THEN begin
          print, 'USER PUPIL SEEMS UNCORRECT SIZE, PLEASE CHECK --> STOP'
          return
       ENDIF
    ENDIF
    
    ; loop on short exposures
    IF short EQ 1 THEN BEGIN
        IF verbose EQ 1 THEN print, '   PSF: ADDITION of SHORT EXPOSURES'
        
        FOR ind_exposure = 0, number_exposure-1 DO begin
            ; creation of turbulent phase map
              dsp2phase, dspall, resall, f, sim_params, $
                          seed = seed, $
                          pupil_kw = pupil_kw, $
                          maskpupil = maskpupil, $
                          kw_static = kw_static, $
                          ttf_kw = ttf_kw, $
                          phase = phase, $
                          static_phase = static_phase, $
                          rms_static_phase = rms_static_phase, $
                          verbose = verbose
          
            ; additionnal jitter on short exposure image
            IF jitter_kw EQ 1 THEN phase = add_jitter(jitter_mas, sim_params, maskpupil, $
                                                      jitter_value = jitter_value, $
                                                      phase = phase, verbose = verbose)
          
            IF pupil_kw EQ 1 THEN begin
                ; creation of short exposure PSF
                phase2psf, phase, maskpupil, sim_params, psf = psf
             
                ; initialisation of long exposure
                IF ind_exposure EQ 0 THEN longexposure = psf * 0.
             
                ; integration for long exposure
                longexposure += psf
            ENDIF   
        ENDFOR
    ENDIF
    
    ; computation of long exposure PSF directly
    IF short EQ 0 THEN BEGIN
        IF verbose EQ 1 THEN print, '   PSF: ON LONG EXPOSURE'
       
        ; Static Phase Creation
        IF kw_static EQ 0 THEN BEGIN
            IF verbose EQ 1 THEN print, '   STATIC ABERRATIONS: NONE'
            rms_static_phase = 0. ; [nm]
            static_phase = calc_static_phase(rms_static_phase, f, sim_params, seed, maskpupil)
        ENDIF 
        IF kw_static EQ 1 THEN BEGIN
            IF verbose EQ 1 THEN print, '   STATIC ABERRATIONS: 1/f6 SPECTRUM'
            static_phase = calc_static_phase(rms_static_phase, f, sim_params, seed, maskpupil)
        ENDIF        
        IF kw_static EQ 2 THEN BEGIN
            IF verbose EQ 1 THEN print, '   STATIC ABERRATIONS: USER DEFINED'
            static_phase = user_static_phase ; User defined static_phase in nm
            ; normalise correctly
;            static_phase = static_phase/sqrt(mean(abs2(static_phase[where(maskpupil)]))) * static_phase_rms
        ENDIF
       
        ; Compute Long Exposure PSF
        longexposure = dsp2psf(dspall, maskpupil, sim_params, static_phase = static_phase)
     ENDIF
    
    ; addition of jitter on long exposure image
    IF pupil_kw EQ 1 OR pupil_kw EQ 2 THEN BEGIN
        IF jitter_kw EQ 2 THEN BEGIN
            IF verbose EQ 1 THEN print, '   JITTER: (YES)', jitter_mas
            longexposure = add_jitter(jitter_mas, sim_params, maskpupil, $
                                      jitter_value = jitter_value, $
                                      psf = longexposure, verbose = verbose)
        ENDIF ELSE BEGIN  
            IF verbose EQ 1 THEN print, '   JITTER: (NO)'
        ENDELSE
        
;        atv, [maskpupil/max(maskpupil), static_phase/max(static_phase)]
    
    ; computation of encircled energy
       pix2rad = sim_params.lambda_imaging/sim_params.D / (sim_params.dim/sim_params.dimpup)
       width = 80.*rad2arcsec * sim_params.lambda_imaging/sim_params.D
       EE = CALC_EE(longexposure, sim_params, pix2rad = pix2rad, width = width, /square, verbose = verbose)
              
       
     ; Computation of Strehl Ratio
       SR = fltarr((size(dspall))[3])
       IF pupil_kw EQ 1 OR pupil_kw EQ 2 THEN BEGIN
          airy = norme(abs2(fftshift(float(maskpupil))))
          Indices = where(airy EQ max(airy))
          WhereToMulti, airy, Indices, Col, Row
          IF (size(dspall))[0] EQ 2 THEN nb_dspall=1 ELSE  nb_dspall=(size(dspall))[3]
          FOR ind_sr=0, nb_dspall-1 DO BEGIN
              SR[ind_sr] = (longexposure(col, row, ind_sr))/(airy(col, row))*100.
          ENDFOR 
       ENDIF
          
    ;   ; save to file : same name, addition of PSF at the end of name
    ;   save, file =  dir + file + '_PSF_LONG' + '.dat', longexposure, maskpupil, static_phase, $
    ;         rms_static_phase, ttf_kw, number_exposure, ee, width, SR
    ENDIF
    
    
    ;Print SR & EE
    L =  sim_params.D*sim_params.dim/sim_params.dimpup
    pixsize =  1./float(L)
    FOR bbb = 0, n_elements(sim_params.beta_tab)/2-1 DO BEGIN  
        resva1 = calc_var_from_psd(dspall(*,*, bbb), pixsize =  pixsize, DD = sim_params.D)
        print,'Performance direction : ', strc(bbb), $
              '  Residual [nm] : ', strc(sqrt(resva1)*sim_params.lambda_imaging * 1e9/2./!pi), '-->', strc(round(sqrt(-alog(SR[bbb]/100))*sim_params.lambda_imaging * 1e9/2./!pi)), $
              '  Strehl Ratio [%] : ',  strc(exp(-(resva1))*100.), '-->', SR[bbb]
    ENDFOR  
    print, 'Reference Wavelength', sim_params.lambda_imaging;, ' Wavelength', new_lambda
    print, 'EE 10mas', calc_ee(longexposure, sim_params, pix2rad = pix2rad, width = 10./1000, /square)
    print, 'EE 50mas', calc_ee(longexposure, sim_params, pix2rad = pix2rad, width = 50./1000, /square)     
    print, '****************************************************'     
             stop
;    ;------ 04/09/2014 ----------
;    over = oversamp_psf(longexposure[*,*,0], sim_params.D, sim_params.dimpup, sim_params.lambda_imaging, newpixsize = 1.0)
;    over = oversamp_psf(longexposure[*,*,0], sim_params.D, sim_params.dimpup, new_lambda, newpixsize = 1.0)
;    psf_final[*,*,iiii] = crop(over, ncrop=Nnew)  
    psf_final[*,*,iiii] = PSF_RESCALE(longexposure[*,*,0], sim_params.lambda_imaging, new_lambda, /crop) 
    psf_final[*,*,iiii] = abs(psf_final[*,*,iiii]/total(psf_final[*,*,iiii], /double))
    shift = [0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.5, $
             1.5, 1.5, 1.5, 1.5, 2.0, 2.5, 2.0, 2.5, 2.5, 2.5, $
             2.5, 2.5]
    psf_final[*,*,iiii] = shiftbyfft(psf_final[*,*,iiii], -shift[iiii], -shift[iiii])
    
;;    over = oversamp_psf(airy, sim_params.D, sim_params.dimpup, sim_params.lambda_imaging, newpixsize = 1.0)   
;    airy_final[*,*,iiii] = PSF_RESCALE(airy, sim_params.lambda_imaging, new_lambda, /crop)
;    airy_final[*,*,iiii] = abs(airy_final[*,*,iiii]/total(airy_final[*,*,iiii], /double))
;    airy_final[*,*,iiii] = shiftbyfft(airy_final[*,*,iiii], -shift[iiii], -shift[iiii])

ENDFOR
    

    ;----------------------------
             
    ; SAVE DATA
    r0 = sim_params.lambda_imaging/(sim_params.seeing0/rad2arcsec)
    Hmean = (total(sim_params.h^(5./3)*sim_params.cn2true)/total(sim_params.cn2true))^(3./5.)  ; Average turbulence height 
    Vmean = (total(sim_params.vent^(5./3)*sim_params.cn2true)/total(sim_params.cn2true))^(3./5.)  ; Average wind speed
    TAU0 = 0.314 * r0/Vmean *1000.
    THETA0 = 0.314*(r0*cos(sim_params.anglezenith/(180./!dpi*60.)))/(Hmean)*rad2arcsec
    res_jitter = 0.
    IF jitter_kw NE 0 THEN res_jitter = 'X='+strc(jitter_mas[0])+'Y='+strc(jitter_mas[1])
    beta = '['+strc(sim_params.beta_tab[0])+','+strc(sim_params.beta_tab[1])+']'
    
    ; -- Write the FITS Files
    IF KEYWORD_SET(savefile) THEN BEGIN
        IF jitter_kw NE 0 THEN jitter = '_Jitter' ELSE  jitter = '_NoJitter'
        IF keyword_set(CoPhasing_error) OR kw_static EQ 2 THEN coPhase = '_M1Error' ELSE coPhase = '_NoM1Error'
        file = dir + file + jitter + coPhase +'.fits'
        mkhdr, h, psf_final
        ;sxaddpar, h, 'CDELT1','???','What is this value?', before='COMMENT'
        ;sxaddpar, h, 'CDELT2','???','What is this value?', before='COMMENT'
        ;sxaddpar, h, 'CTYPE1','???','What is this value?', before='COMMENT'
        ;sxaddpar, h, 'CTYPE2','???','What is this value?', before='COMMENT'
        sxaddpar, h, 'r0',strc(sim_params.lambda_imaging/(sim_params.seeing0/rad2arcsec)),'Fried Parameter @ Zenith (calculated)', before='COMMENT'
        sxaddpar, h, 'WAVELENG',strc(sim_params.lambda_imaging),'imaging wavelength in metres ', before='COMMENT'
        sxaddpar, h, 'TELESCOP',strc(sim_params.D),'in metres (clear IR aperture)', before='COMMENT'
        sxaddpar, h, 'PIX_SIZE',strc(pix2rad*rad2arcsec*1000.),'Size of pixel in mas', before='COMMENT'
        sxaddpar, h, 'TEL_ABER','None','No co-phasing error', before='COMMENT'
        sxaddpar, h, 'RES_JITT',res_jitter,'X and Y in mas rms', before='COMMENT'
        sxaddpar, h, 'SEEING',strc(sim_params.seeing0),'arcsec @ Zenith', before='COMMENT'
        sxaddpar, h, 'L0',strc(sim_params.L0),'in metres', before='COMMENT'
        sxaddpar, h, 'ZENITH_A',strc(sim_params.anglezenith/60.),'degrees from Zenith', before='COMMENT'
        sxaddpar, h, 'TAU0', strc(TAU0),'ms (calculated)', before='COMMENT'
        sxaddpar, h, 'THETA0',strc(THETA0),'arcsec @0.5micron (calculated)', before='COMMENT'
        sxaddpar, h, 'AO_MODE', 'SCAO-on axis optim', 'AO Flavour'  , before='COMMENT'
        sxaddpar, h, 'LGS_NUMB', '0', 'No LGS as this is SCAO' , before='COMMENT'
        sxaddpar, h, 'LGS_SEPA', '0', 'arcmin. No LGS as this is SCAO'  , before='COMMENT'
        sxaddpar, h, 'NGS_NUMB', '1', 'Only one GS' , before='COMMENT'
        sxaddpar, h, 'NGS_SEPA', '0', 'arcmin. None since SCAO', before='COMMENT'
        sxaddpar, h, 'NGS_MAGN', 'TBD', 'Jovian Moons' , before='COMMENT'
        sxaddpar, h, 'PERF_DIR', beta, 'Direction of Performance from NGS in arcmin' , before='COMMENT'
        sxaddpar, h, 'SAMP_FRE', strc(sim_params.fech), 'Hz', before='COMMENT'
        sxaddpar, h, 'LOOP_DEL', strc(sim_params.td*1000), 'ms' , before='COMMENT'
        sxaddpar, h, 'PITCH', strc(sim_params.pitch_m), 'm'  , before='COMMENT'
        sxaddpar, h, 'SR_TURB', strc(strc(exp(-(resva1[0]))*100.)), 'Turbulent effects only', before='COMMENT'
        sxaddpar, h, 'SR_TUR_T', strc(SR[0]), 'Include M1 residuals and Jitter (approximate)'  , before='COMMENT'
        sxaddpar, h, 'IMAGE_1', 'PSF with perfect M1', '', before='COMMENT'
        sxaddpar, h, 'VARNOISE', strc(sim_params.varnoise_wfs[0]), 'WFS Analysis noise in rad^2', before='COMMENT'
     
        writefits, file, abs(psf_final), h
;        file = 'E:\HARMONI Simulations\SCAO\Simon-4Sept2014\2048\Airy.fits
;        writefits, file, abs(airy_final), h
    ENDIF 
stop

ENDFOR 
stop

;****************  NEW 03 April 2014 ***************************
dspall = DSP_ALL_OLD
;***************************************************************
END

