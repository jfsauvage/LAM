
pro scao_harmoni


anglezenith_tab = [30.0]; [0.0, 15., 30., 45., 60.];                        ; in [degrees]
;lambda_imaging_tab = 2.2e-6;                                               ; in [m]
lambda_imaging_tab = 2.2e-6;[0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5] * 1.0e-6
varnoise_wfs_tab = [0.1];                                                  ; in [rad2]
ind_turb_case = 0;'Jmedian'==0; 'JQ1'==1; 'JQ2'==2; 'JQ3'==3; 'JQ4'==4;    ; ESO Turbulence Case


; Performance direction [arcmin]                                           
beta_tabxxx = [[0.0, 0.0]]/60.0
;beta_tabxxx = [[0.0, 0.0], [0.5, 0.0], [1.0, 0.0], [1.5, 0.0], [2.0, 0.0], $           ; One ligne --> 17 positions
;               [2.5, 0.0], [3.0, 0.0], [4.0, 0.0], [5.0, 0.0], [6.0, 0.0], $
;               [7.0, 0.0], [8.0, 0.0], [9.0, 0.0], [10., 0.0], [12., 0.0], $
;               [15., 0.0], [20, 0.0]]/60.0

; Outscale L0 [m]
L0_tab = [25.];[10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 100.0]

;; A whole square --> 633 positions
;tmp = [[0.0, 0.0]]
;deb = tmp
;FOR ind_tmp=1,12 DO tmp = [[tmp], [deb[0,*]+ind_tmp, deb[1,*]], [deb[0,*]-ind_tmp, deb[1,*]]]
;deb = tmp
;FOR ind_tmp=1,12 DO tmp = [[tmp], [deb[0,*], deb[1,*]+ind_tmp], [deb[0,*], deb[1,*]-ind_tmp]]
;tmp = [[tmp], [0.5, 0.0], [-0.5, 0.0], [0.5, 0.5], [0.5, -0.5], [-0.5, 0.5], [-0.5, -0.5], [0.0, 0.5], [0.0, -0.5]]
;beta_tab = tmp
; NB for r0=0.8m and h(bar)=5km: theta0=10"@2.2um and theta0=0.75"@.5um !!!

; Calculate dimpup
mas_per_pix = 1.0;3.0         ; mas per pixel
dim = 512. & adim = dim
D = 37.0 ;(need to change in main_harmoni.pro!)
rad2arcsec = 180./!dpi*60.*60.
; pix2rad = lambda/D * dimpup/dim & pix2mas = lambda/D * dimpup/dim * rad2arcsec * 1000.
dimpup_tab = 128.; ceil(D/lambda_imaging_tab *dim*mas_per_pix /(rad2arcsec*1000.)); 512.; 
pix2mas = lambda_imaging_tab/(D) * dimpup_tab/dim * rad2arcsec * 1000.
print, '---------------------------------------'
print, 'Wavelengths (um)'
print, lambda_imaging_tab*1.0e6
print, 'dimpup x dim'
print, dimpup_tab, dim
print, 'mas per pixel'
print, pix2mas
print, 'Number of Fourier Simulations (i.e. output files)'
Nb_simul = N_ELEMENTS(anglezenith_tab) * N_ELEMENTS(lambda_imaging_tab) * N_ELEMENTS(varnoise_wfs_tab) * N_ELEMENTS(L0_tab) * 5
print, Nb_simul
print, 'Number of obs. angles (beta_tab)'
print, N_ELEMENTS(beta_tab)/2 
print, '---------------------------------------'


; Loop over Zenith Angles
counter = 1
FOR ind_angle = 0, n_elements(anglezenith_tab)-1 DO BEGIN
    anglezenith = anglezenith_tab[ind_angle]
;    print, '*********************************************'
;    print, 'Angle Zenith (index, angle): ', ind_angle, anglezenith_tab[ind_angle]

    ; Loop over Outer Scale
    FOR ind_L0 = 0, n_elements(L0_tab)-1 DO BEGIN
        L0 = L0_tab[ind_L0]
;        print, 'L0 (index, L0): ', ind_L0, L0_tab[ind_L0]
;        print, 'Simulation ', strn(counter), ' of ', strn(Nb_simul)
;        print, '*********************************************'
        
;        FOR ind_beta_tab = 2,3 DO BEGIN;, n_elements(beta_tabxxx)-1 DO BEGIN
;            beta_tab = beta_tabxxx[*, ind_beta_tab]
;            print, beta_tab*60.
;            print, 'beta_tab (index, beta_tab): ', ind_beta_tab, beta_tabxxx[*, ind_beta_tab]

         ; Loop over imaging wavelength
         FOR ind_lamdba = 0, n_elements(lambda_imaging_tab)-1  DO BEGIN
          
            ; Loop over turbulence conditions
            FOR ind_turb_case = 0,0 DO BEGIN
              
                beta_tab = beta_tabxxx
                varnoise_wfs = varnoise_wfs_tab[0]
                dimpup = dimpup_tab[ind_lamdba]
                lambda_imaging = lambda_imaging_tab[ind_lamdba]
                
                print, '*********************************************'
                print, 'VARNOISE: ', varnoise_wfs
                print, 'Wavelength (index, angle): ', ind_lamdba, lambda_imaging
                print, 'DIMPUP (dimpup, dim): ', dimpup, adim
                print, 'Angle Zenith (index, angle): ', ind_angle, anglezenith_tab[ind_angle]
                print, 'TURB CASE: ', ind_turb_case
                print, 'beta_tab: ', beta_tab
                print, 'L0: ', L0_tab[ind_L0]
                print, 'Simulation ', strn(counter), ' of ', strn(Nb_simul)
                print, '*********************************************'
               
                main_harmoni, varnoise_wfs = varnoise_wfs, lambda_imaging = lambda_imaging, $
                              dimpup = dimpup, adim = adim, anglezenith = anglezenith, $
                              ind_turb_case = ind_turb_case, beta_tab = beta_tab, L0 = L0 
                counter = counter +1
            ENDFOR
        ENDFOR
    ENDFOR
ENDFOR


;; Loop over Zenith Angles
;FOR ind_angle = 0, n_elements(anglezenith_tab)-1 DO BEGIN
;    anglezenith = anglezenith_tab[ind_angle]
;    print, 'Angle Zenith (index, angle): ', ind_angle, anglezenith_tab[ind_angle]
;    
;    ; Loop over imaging wavelength
;    FOR ind_lamdba = 0, n_elements(lambda_imaging_tab)-1  DO BEGIN
;        dimpup = dimpup_tab[ind_lamdba] 
;        lambda_imaging = lambda_imaging_tab[ind_lamdba]
;        print, 'Wavelength (index, angle): ', ind_lamdba, lambda_imaging
;        
;        ; Loop over WFS noise levels
;        FOR ind_varnoise=0, n_elements(varnoise_wfs_tab)-1 DO BEGIN
;            varnoise_wfs = varnoise_wfs_tab[ind_varnoise]
;            print, 'Noise (index, angle): ', ind_varnoise, varnoise_wfs
;            
;            main_harmoni, varnoise_wfs = varnoise_wfs, lambda_imaging = lambda_imaging, $
;                          dimpup = dimpup, adim = adim, anglezenith = anglezenith, $
;                          ind_turb_case = ind_turb_case, beta_tab = beta_tab
;        ENDFOR
;    ENDFOR
;ENDFOR 


END 


;SR15 = [44.9728,44.7008,43.5906,41.7482,39.3298,33.4582,27.2022,21.3712,$
;        16.3670,9.11940,4.86532,2.52263,1.28097,0.642017,0.320116,0.159703,0.0799528,0.0402067]
;SR22 = [77.0353,76.8673,76.1743,74.9991,73.4049,69.2537,64.2801,58.9323,$
;        53.5354,43.3702,34.5898,27.3049,21.3920,16.6813,12.9837,10.1075,7.87840,6.15165]
;plot,beta_tab[0,*]*60, SR22, psym=-1, yr=[0,100], thick = 4, xtit='Angle [Arcsec]'
;oplot,beta_tab[0,*]*60, SR15, psym=-1, col=fsc_color('red'), thick = 4



;; lambda_wfs = 0.675e-6         ; wfs wavelength [m]
;main_harmoni, varnoise_wfs = 1, lambda_imaging = 1.5e-6, dimpup = 490., anglezenith = 0
;; ==> Performance direction : 0  Residual [nm] : 186.187  Strehl Ratio [%] : 54.4308



