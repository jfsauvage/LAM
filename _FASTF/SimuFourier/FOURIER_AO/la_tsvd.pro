; cette procédure calcul l'inverse généralisée d'une matrice (carrée ou rectangulaire, complexe ou réelle) en réalisant sa décomposition en valeur singulière via la_svd (bibliothèque Lapack). 

; ENTREES
; mat = matrice à inverser
; condmax = condition maximal souhaité pour filtrer le seuil de troncature. 
; nbr_mode_filtre = nombre de modes à filtrer. cette méthode de filtrage est exclusive de la précédente, laquelle est prioritaire
; visu = booléen pour le tracé des valeurs singulieres
; silent = booléen pour supprimer les sorties texte
; double = pour des calculs réalisés en double précision. si vrai, alors quelque soit le format des données d'entrée, la sortie est en double précision. si l'entrée est déjà en double, par défaut les calculs sont en double précision



; SORTIES
; condition = conditionnement initial avant filtrage
; val_propre = vecteur contenan les valeurs singulieres
; inverse = inverse généralisé
;  condfiltre = conditionnement apres filtrage

; utilisation 
; la_svd, mat = mat, inverse = inverse, [condmax = condmax], [condition =  condition], [nbr_mode_filtre =  nbr_mode_filtre],  [visu =  visu],  [val_sing =  val_sing],  [silent =  silent],  [double =  double]

PRO  la_tsvd,  mat =  mat, inverse =  inverse,  condmax = condmax,  condition = $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $
               condition, nbr_mode_filtre =  nbr_mode_filtre,  visu =  visu,  val_sing =  val_sing,  silent =  silent,  double =  double,  condfiltre =  condition_filtre






; definition du parametre status de la_svd pour connaitre le resultat de la decomposition
status =  0
; appel de la_svd
la_svd,  mat,  u, v, w,  double =  double,  status =  status

; si la_svd à convergé
IF status EQ  0 THEN BEGIN
   
   sz =  (size(u))(1)
   ; plot des valeurs singulieres ordonnées
   IF keyword_set(visu) THEN plot,  u(sort(u))
   ; les valeurs singulieres sont reelles non negatives meme en complexe
   ; je calcule le conditionnement initial
   IF min(u) NE 0 THEN BEGIN
      condition =  max(u)/min(u) 
      IF NOT keyword_set(silent) THEN print,  'conditionnement avant filtrage  =',  condition
   ENDIF ELSE BEGIN
      condition =  0
      IF NOT keyword_set(silent) THEN print,  'conditionnement infini , 0 par défaut'
   ENDELSE 
   val =  u*0.
   
     ; cas 1 : on tronque en fixant un conditionnement maximal, qui fixe aussi un seuil sur les valeurs singulieres. Les valeurs singulieres inférieures ne sont pas inversées : on fixe leur inverse à 0  
   
   IF keyword_set(condmax) THEN BEGIN    
      limite =  max(u)/condmax    
      s =  where(u LE limite, count)
;IF total(mat) EQ 0. THEN stop
      IF NOT keyword_set(silent) THEN print,  'filtrage de',  count,  'modes'
      val(where(u GE limite)) =  1./(u(where(u GE limite)))
      condition_filtre =   max(u)/min((u(where(u GE limite)))) 
      IF NOT keyword_set(silent) THEN print, 'conditionnement après filtrage =', condition_filtre
   ENDIF ELSE BEGIN
      ; cas 2 : on fixe un nombre de valeurs singulieres à tronquer
      
      IF keyword_set(nbr_mode_filtre) THEN BEGIN 
         IF nbr_mode_filtre LT sz THEN  BEGIN
            IF NOT keyword_set(silent) THEN print,'filtrage de  ',  nbr_mode_filtre,  'modes'
            val(nbr_mode_filtre:*) =  1./(u(sort(u)))(nbr_mode_filtre:*) 
            condition_filtre =  max(u)/(u(sort(u)))(nbr_mode_filtre)
            IF NOT keyword_set(silent) THEN print, 'conditionnement après filtrage =', condition_filtre
         ENDIF  ELSE BEGIN $ $ $ $ $ $ $ $ $
             IF NOT keyword_set(silent) THEN print, 'nbr_de mode à filtrer trop elevé pas de filtrage'
            
            val =  1/u
            condition_filtre = condition 
         ENDELSE  
      ENDIF   ELSE  BEGIN   ; cas ou on ne tronque rien 
         IF NOT keyword_set(silent) THEN print,'pas de filtrage'
         val =  1/u         
         condition_filtre = condition 
      ENDELSE
   ENDELSE
   
    
; enfin on calcule l'inverse généralisé en distinguant le cas réel du cas complexe
diag =  diag_matrix(val)
val_sing =  u
IF (size(v))(3) GT 5 THEN inverse =  w ## diag ## transpose(conj(v)) ELSE  inverse =  w ## diag ## transpose(v)
$ $ $ $
    ENDIF   ELSE    print,  'la svd n a  pas convergé' ; cas d'une absence de convergence de la_svd

END 
