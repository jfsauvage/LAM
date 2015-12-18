FUNCTION sr_residu, f, fx, fy, dim, dimpup, D, arg_f,pitchs_wfs,pitchs_DM,Nb_gs,alpha, $
                           beta_tab, theta, sig, $
                           DSP_tab_recons, $
                           h_recons, DSP_tab_vrai, h_vrai, h_dm, $
                           seuil, condmax, Wflag, $
                           LSE = LSE, $
                           servo = servo, td, ti, Wind, $
                           fitting = fitting, Tomo = tomo, Gfit = Gfit,$
                           alias = alias, Cn2_recons, L0, r0, err_recons = err_recons, err_noise = err_noise, $
                           turb = turb,  DSP_tab = DSP_tab, PSF_tab = PSF_tab, Var_tab = var_tab



;-----------------------------------------------------------
;CALCUL DE LA MATRICE DE RECONSTRUCTION TOMOGRAPHIQUE
;-----------------------------------------------------------
local_L = D*dim/dimpup  ;taille de l'ecran en m.
ech = 2. ;surech
rt = dimpup/2.
ri = 0. ;obstruction centrale
polaire2, RT = rt, RHO= rho,  PHI=phi, /ENTRE4PIXELS, $
            MASQUE= pup,oc = ri;, largeur= dpup


IF keyword_set(turb) THEN BEGIN 
print, 'CAS PUREMENT TURBULENT'
s = (size(f))(1)
dsp_res = fltarr(s, s)
FOR iii = 0, n_elements(h_vrai)-1 DO BEGIN
dsp_res += dsp_tab_vrai(s*iii:s*(iii+1)-1, *)
ENDFOR
GOTO, fin
ENDIF

;Remarque préliminaire:
;Dans le cas on veut pas d'aliasing, on calcul un Wmap direct en modul de f : le calcul est plus rapide
;Dans le cas où on veut de l'aliasing, on doit calculer un Wmap en "(x,y)" donc deux
;fois plus grand, plus de calculs = c'est plus long.
;Donc, je laisse les 2 options, avec ou sans Aliasing.

;Pour l'erreur de recons et de noise (et si le W n'inclut pas l'erreur d'aliasing) tout le monde donne la mm chose


IF keyword_set(alias) THEN BEGIN
print, '--------------------------'
print, 'WITH ALIASING'
print, '--------------------------'
ENDIF


IF NOT keyword_set(alias) THEN BEGIN 
;on a le choix entre les deux formes de reconstructeurs,
;La forme W1 est avec les (Cphi)^-1 (Cb)^1
;La forme W2 est plus robuste au cas de fort SNR
Wflag = wflag

;condmax est utile lors du calcul de l'inversion dans Popt (faite par la_tsvd)
;condmax = 1e6

;on ajoute une option pour sauter la projection optimale :
;tomo = 1B
;si tomo est ok alors il n'y a pas de projection sur les DMs
;C'est pour le cas tomographique pure : e.g. LTAO/MOAO
;Ca permet d'accelerer un peu les calculs
;Dans ce cas, le pitch DM est le meme pour toutes les couches

;Une autre option pour si on veut faire du moindre carré (pondéré):
;LSE = 0B
IF keyword_set(LSE) THEN BEGIN
print, '--------------------------'
print, 'LSE RECONSTRUCTOR'
print, '--------------------------'
ENDIF
;Si LSE =1B, on fait du WLSE : (MtCb-1M)-1MtCb-1

;seuil = 1e9
;pour l'inverse dans le W
;et surtout pour le cas LSE !!
;genre en MMSE, on peut mettre 1e9
;mais en LSE, c'est plutot 100, voir moins...

wmap = calc_mat_rec_finale(f, arg_f, pitchs_wfs,pitchs_DM, Nb_gs,alpha, theta,sig, DSP_tab_recons, $
                           h_recons, h_dm, $
                           seuil, condmax, Wflag, $
                           gfit = gfit, $
                           LSE = LSE)


ENDIF ELSE BEGIN
  
;Si on veut inclure l'aliasing, faut utiliser une matrice de reconstruction qui inclut la mesure selon x et y.

;En premier, on doit calculer le Cb_alias:
;C'est un Cb_alias model, le model de bruit d'aliasing qu'on utilise (ou pas) dans le reconstructeur

a = max(pitchs_wfs)
Cb_Alias = calc_CbAlias(f, fx,fy, arg_f, a,Nb_gs,alpha, h_recons ,Cn2_recons, L0, r0)

;Remarque1 : Pour l'instant y'a pas multipitch, en fait ca va etre assez chiant à inclure...
;En premiere approximation, pessimiste, on choisit le plus grand pitch
;Remarque2 : Le calcul de Cb_alias est assez long, mais tant que le pitch 
;et la GS constellation n'est pas modifié, 
;c'est toujours le meme, on peut donc l'ecrire sur le disque et venir le lire pour gagner du temps.

;writefits, 'Cb_alias.fits',Cb_alias
;Cb_alias = readfits('Cb_alias.fits', /SILENT)

;Puis on calcul un W_alias qui inclut la connaissance de Cb_alias et/ou du bruit
;pour un W_alias avec uniquement du bruit classique il faut faire Cb_alias*0.
;pour un W_alias avec uniquement la connaissance de l'aliasing il faut faire sig*0.

;seuil = 1e9
;Seuil c'est parceque qd on calcul le W avec que Cb_alias comme bruit à priori, y'a des fréquences
;qu'on ne peut inverser (Cb_alias = 0 pour f_x ou f_y =0), donc je fais une tsvd
;Tomo = 0B
;LSE = 0B
IF keyword_set(LSE) THEN BEGIN
print, '--------------------------'
print, 'LSE RECONSTRUCTOR'
print, '--------------------------'
ENDIF
;Le LSE ne prend pas en compte le bruit d'aliasing


;Cb_alias = fltarr(2*nb_gs*s, 2*nb_gs*s)

;le vecteur de bruit doit etre deux fois plus grand, car on a le bruit selon x et selon y:
;on duplique chaque bruit/GS, mais on pourrait mettre un bruit different selon x/y, je sais ce que ca veut 
;dire physiquement(Des spots tous allongés pareils), mais on peut...
sig_xy = fltarr(2*nb_gs)
FOR nmn = 0, nb_gs-1 DO BEGIN
sig_xy(nmn*2.) = sig(nmn)
sig_xy(nmn*2.+1.) = sig(nmn)
ENDFOR

;seuil = 1e12

Wmap = calc_w_alias(f, fx, fy, arg_f,pitchs_wfs,pitchs_DM,Nb_gs,alpha, theta,sig_xy, DSP_tab_recons,$
                    h_recons, h_dm, $
                    seuil, cb_alias, $
                    condmax, $
                    Wflag, $
                    gfit = gfit, $
                    LSE = LSE)

ENDELSE


;stop
;-----------------------------------------------------
;CALCUL DE LA DSP RESIDUELLE TOMO 
;-----------------------------------------------------


nb_dir_perf = round(n_elements(beta_tab)/2.)
sr_tab = fltarr(nb_dir_perf)

IF keyword_set(DSP_tab) THEN dsp_tabtab = fltarr(nb_dir_perf, dim, dim)
IF keyword_set(PSF_tab) THEN PSF_tabtab = fltarr(nb_dir_perf, dim, dim)
IF keyword_set(var_tab) THEN var_tabtab = fltarr(nb_dir_perf)

FOR bbb = 0, nb_dir_perf-1 DO BEGIN

beta = [beta_tab(0, bbb), beta_tab(1, bbb)]

;Calcul des dsp residuelle Tomographique dans chaque direction:

IF NOT keyword_set(gfit) THEN  BEGIN 
h_dm = h_recons
print, 'TOMO PURE RECONSTRUCTION'
ENDIF

IF NOT keyword_set(alias) THEN BEGIN

;tempo = 0B ;pour inclure ou pas l'erreur de servo-lag
;autre option pour dire si on veut que le fitting soit inclut ou pas
;fitting = 1B
;par défaut, le fitting s'inclut tout seul, mais si on le veut pas, on met fitting à OB.

dsp_res = calc_dsp_res_finale(f, arg_f, pitchs_wfs, Nb_gs,alpha,beta,sig, DSP_tab_vrai, h_vrai, h_dm, Wmap, $
                              tempo = servo, td, ti, Wind, $
                              fitting = fitting,  err_recons = err_recons, err_noise = err_noise)

ENDIF ELSE BEGIN
  
;Si on fait une erreur de model sur la statistique de bruit de l'aliasing, on peut recalculer un cb_alias vrai:
;Cb_Alias = calc_CbAlias(f, fx,fy, arg_f, a,Nb_gs,alpha, h_vrai ,Cn2_vrai, L0, r0)
;En toute rigeur, si on a un temps d'integration fini / WFS (i.e. ti NE 0) il faut alors recalculer
;Cb_alias en incluant ti, car le signal qui se replie c'est une signal moyenné sur ti.
;J'ai pas encore inclut ca, mais je pense que c'est peanuts (sauf si on integre des heures...)
;Donc pour l'instant, le Cb_alias d'avant suffit bien
;Surout que la calcul de Cb_alias est assez long...

;tempo = 1B ;pour inclure ou pas l'erreur de servo-lag
;autre option pour dire si on veut que le fitting soit inclut ou pas
;fitting = 1B

dsp_res = calc_dsp_res_alias(f, fx, fy, arg_f, pitchs_wfs, Nb_gs,alpha,beta,sig_xy, DSP_tab_vrai, $
                             h_vrai, h_dm, Wmap, $
                              tempo = servo, td, ti, Wind, Cb_alias, $
                              fitting = fitting, $
                              alias_alone = alias_alone,  err_recons = err_recons, err_noise = err_noise)

alias = eclat(alias_alone)
;print, total(alias_alone)*f(1)^2./(a/r0)^(5./3.)

;ce dsp_res inclut l'aliasing + propa du bruit classique + Tomo et tout.
;et alias_alone contient la dsp de l'aliasing seule, err_recons que l'err de recons et err_noise que l'erreur de noise
ENDELSE

;pour récuperer les DSP d'erreur individuelles:
IF keyword_set(fitting) THEN BEGIN
s = (size(f))(1)
turb = fltarr(s, s)
FOR iii = 0, n_elements(h_vrai)-1 DO BEGIN
turb += dsp_tab_vrai(s*iii:s*(iii+1)-1, *)
ENDFOR
fc = max(1/(2*pitchs_wfs))
ind = where(abs(fx) GE fc OR abs(fy) GE fc)
fitting = fltarr(dim, dim)
fitting(ind) = turb(ind)
fitting = eclat(fitting)
ENDIF
fc = max(1/(2*pitchs_wfs))
ind = where(abs(fx) GE fc OR abs(fy) GE fc)
err_recons(ind) = 0.
err_recons = eclat(err_recons)
err_noise = eclat(err_noise)

;-----------------------------------------------------
;Pour faire des Psfs :
;-----------------------------------------------------
;dsp_tab(bbb) = dsp_res
IF keyword_set(DSP_tab) THEN dsp_tabtab(bbb, *, *) = eclat(dsp_res)
IF keyword_set(PSF_tab) THEN psf_tabtab(bbb, *, *) = dsp_to_psf(dsp_res, pup, local_L, ech)
;tab_psf(bbb, *, *) = 
airy = dsp_to_psf(dsp_res*0., pup, local_L, ech)
psf_map = dsp_to_psf(dsp_res, pup, local_L, ech)
sr_tab(bbb) = max(psf_map)/max(airy)*100.
IF keyword_set(var_tab) THEN var_tabtab(bbb) = total(dsp_res)*f(1)^2.
;
;SR_tab(bbb) = max(psf_map)/max(airy)*100.
;var_tab(bbb) = total(dsp_res)*f(1)^2.
;print, 'Il reste : ', nb_dir_perf-bbb-1
;up

ENDFOR


IF nb_dir_perf EQ 1 AND keyword_set(DSP_tab) THEN dsp_tabtab = reform(dsp_tabtab)
IF nb_dir_perf EQ 1 AND keyword_set(PSF_tab) THEN psf_tabtab = reform(psf_tabtab)

IF keyword_set(PSF_tab) THEN PSF_tab = psf_tabtab
IF keyword_set(var_tab) THEN var_tab = var_tabtab
IF keyword_set(dsp_tab) THEN dsp_tab = dsp_tabtab

fin:

IF keyword_set(turb) THEN BEGIN
dsp_tab = eclat(dsp_res)
var_tab = total(dsp_tab)*f(1)^2.
psf_tab =  dsp_to_psf(dsp_res, pup, local_L, ech)
airy = dsp_to_psf(dsp_res*0., pup, local_L, ech)
sr_tab = max(psf_tab)/max(airy)*100.
ENDIF

return, sr_tab


END
