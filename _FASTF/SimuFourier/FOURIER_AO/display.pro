PRO display, sepchampvisu = sepchampvisu, Nb_dir_perf = Nb_dir_perf, nb_gs = $
              nb_gs, alpha = alpha, nb_dir_theta = nb_dir_theta, theta = $
             theta, beta_tab = beta_tab

window, 0, xs =  800, ys = 800
tek_color
champ = [-1.5, 1.5] * max([sepchampvisu, max(abs(alpha)), max(abs(theta)), max(abs(beta_tab))])

plot, champ, champ, /nod, position = [[0.1, 0.1], [0.95, 0.95]], $
      xtit = 'Arcmin', ytit = 'Arcmin'
PLOTSYM, 0, .5
FOR iii = 0, Nb_dir_perf-1 DO BEGIN 
   oplot, [beta_tab[0, iii]], [beta_tab[1, iii]], psym = 8
ENDFOR
PLOTSYM, 3, 2., /FILL
FOR iii = 0, nb_gs-1 DO BEGIN
   oplot, [alpha[0, iii]], [alpha[1, iii]], psym = 8
ENDFOR
PLOTSYM, 8, 1.
FOR iii = 0, nb_dir_theta-1 DO BEGIN
   oplot, [theta[0, iii]], [theta[1, iii]], psym = 8, color =  2, thick = 2
ENDFOR

PLOTSYM, 8
oplot, [0.9*champ[0], 0.9*champ[0]], [0.85*champ[1], 0.85*champ[1]], psym = 8, color = 2
xyouts, 0.8*champ[0], 0.83*champ[1], 'Optimisation directions'
PLOTSYM, 3, /FILL
oplot, [0.9*champ[0], 0.9*champ[0]], [0.8*champ[1], 0.8*champ[1]], psym = 8
xyouts, 0.8*champ[0], 0.78*champ[1], 'Guide Stars'
PLOTSYM, 0
oplot, [0.9*champ[0], 0.9*champ[0]], [0.75*champ[1], 0.75*champ[1]], psym = 8
xyouts, 0.8*champ[0], 0.73*champ[1], 'PSF estimation directions'


END
