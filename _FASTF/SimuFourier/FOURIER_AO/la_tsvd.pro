; cette proc�dure calcul l'inverse g�n�ralis�e d'une matrice (carr�e ou rectangulaire, complexe ou r�elle) en r�alisant sa d�composition en valeur singuli�re via la_svd (biblioth�que Lapack). 

; ENTREES
; mat = matrice � inverser
; condmax = condition maximal souhait� pour filtrer le seuil de troncature. 
; nbr_mode_filtre = nombre de modes � filtrer. cette m�thode de filtrage est exclusive de la pr�c�dente, laquelle est prioritaire
; visu = bool�en pour le trac� des valeurs singulieres
; silent = bool�en pour supprimer les sorties texte
; double = pour des calculs r�alis�s en double pr�cision. si vrai, alors quelque soit le format des donn�es d'entr�e, la sortie est en double pr�cision. si l'entr�e est d�j� en double, par d�faut les calculs sont en double pr�cision



; SORTIES
; condition = conditionnement initial avant filtrage
; val_propre = vecteur contenan les valeurs singulieres
; inverse = inverse g�n�ralis�
;  condfiltre = conditionnement apres filtrage

; utilisation 
; la_svd, mat = mat, inverse = inverse, [condmax = condmax], [condition =  condition], [nbr_mode_filtre =  nbr_mode_filtre],  [visu =  visu],  [val_sing =  val_sing],  [silent =  silent],  [double =  double]

PRO  la_tsvd,  mat =  mat, inverse =  inverse,  condmax = condmax,  condition = $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $
               condition, nbr_mode_filtre =  nbr_mode_filtre,  visu =  visu,  val_sing =  val_sing,  silent =  silent,  double =  double,  condfiltre =  condition_filtre






; definition du parametre status de la_svd pour connaitre le resultat de la decomposition
status =  0
; appel de la_svd
la_svd,  mat,  u, v, w,  double =  double,  status =  status

; si la_svd � converg�
IF status EQ  0 THEN BEGIN
   
   sz =  (size(u))(1)
   ; plot des valeurs singulieres ordonn�es
   IF keyword_set(visu) THEN plot,  u(sort(u))
   ; les valeurs singulieres sont reelles non negatives meme en complexe
   ; je calcule le conditionnement initial
   IF min(u) NE 0 THEN BEGIN
      condition =  max(u)/min(u) 
      IF NOT keyword_set(silent) THEN print,  'conditionnement avant filtrage  =',  condition
   ENDIF ELSE BEGIN
      condition =  0
      IF NOT keyword_set(silent) THEN print,  'conditionnement infini , 0 par d�faut'
   ENDELSE 
   val =  u*0.
   
     ; cas 1 : on tronque en fixant un conditionnement maximal, qui fixe aussi un seuil sur les valeurs singulieres. Les valeurs singulieres inf�rieures ne sont pas invers�es : on fixe leur inverse � 0  
   
   IF keyword_set(condmax) THEN BEGIN    
      limite =  max(u)/condmax    
      s =  where(u LE limite, count)
;IF total(mat) EQ 0. THEN stop
      IF NOT keyword_set(silent) THEN print,  'filtrage de',  count,  'modes'
      val(where(u GE limite)) =  1./(u(where(u GE limite)))
      condition_filtre =   max(u)/min((u(where(u GE limite)))) 
      IF NOT keyword_set(silent) THEN print, 'conditionnement apr�s filtrage =', condition_filtre
   ENDIF ELSE BEGIN
      ; cas 2 : on fixe un nombre de valeurs singulieres � tronquer
      
      IF keyword_set(nbr_mode_filtre) THEN BEGIN 
         IF nbr_mode_filtre LT sz THEN  BEGIN
            IF NOT keyword_set(silent) THEN print,'filtrage de  ',  nbr_mode_filtre,  'modes'
            val(nbr_mode_filtre:*) =  1./(u(sort(u)))(nbr_mode_filtre:*) 
            condition_filtre =  max(u)/(u(sort(u)))(nbr_mode_filtre)
            IF NOT keyword_set(silent) THEN print, 'conditionnement apr�s filtrage =', condition_filtre
         ENDIF  ELSE BEGIN $ $ $ $ $ $ $ $ $
             IF NOT keyword_set(silent) THEN print, 'nbr_de mode � filtrer trop elev� pas de filtrage'
            
            val =  1/u
            condition_filtre = condition 
         ENDELSE  
      ENDIF   ELSE  BEGIN   ; cas ou on ne tronque rien 
         IF NOT keyword_set(silent) THEN print,'pas de filtrage'
         val =  1/u         
         condition_filtre = condition 
      ENDELSE
   ENDELSE
   
    
; enfin on calcule l'inverse g�n�ralis� en distinguant le cas r�el du cas complexe
diag =  diag_matrix(val)
val_sing =  u
IF (size(v))(3) GT 5 THEN inverse =  w ## diag ## transpose(conj(v)) ELSE  inverse =  w ## diag ## transpose(v)
$ $ $ $
    ENDIF   ELSE    print,  'la svd n a  pas converg�' ; cas d'une absence de convergence de la_svd

END 
