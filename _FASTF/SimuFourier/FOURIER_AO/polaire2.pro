; $Id: polaire2.pro,v 1.6 2012/10/03 08:41:24 lmugnier Exp $
PRO POLAIRE2, RT=rt, RHO=rho, PHI=phi, MASQUE=masque,  $
              POINTS_UTILES=points_utiles, OC=oc, ENTRE4PIXELS=entre4pixels, $
              LEQ=leq, LARGEUR=largeur, CENTRE_X = centre_x, CENTRE_Y = $
              centre_y, DOUBLE=double, PROLIXE = prolixe, $   
              HELP=help, VERSION=version 


;+ 
;NOM : 
; POLAIRE2 - Calcul des coordonnees polaires et de la pupille d'un telescope.
; 
;CATEGORIE : 
; Optics Routines 
; 
;SYNTAXE : 
; POLAIRE2, RT=rt [, RHO=rho, PHI=phi, MASQUE=masque, /LEQ, OC=oc,
; POINTS_UTILES=points_utiles, /ENTRE4PIXELS, LARGEUR=largeur,
; CENTRE_X = centre_x, CENTRE_Y = centre_y, /DOUBLE, /PROLIXE, /HELP, /VERSION]
;
;DESCRIPTION : 
;  
; POLAIRE2 genere les coordonnees  polaires RHO et PHI et le MASQUE 
; d'intensite (pupille) d'un telescope de rayon RT et de taux d'occultation
; centrale lineaire OC.
; 
; En general (i.e., si on ne force pas les variables LARGEUR, CENTRE_X ou
; CENTRE_Y, cf. plus bas) les tableaux sont de largeur paire (= 2*round(RT))
; si ENTRE4PIXELS est present et non nul, de largeur impaire (= 2*fix(rt)+1)
; sinon.
; 
;PARAMETRES : 
; 
;  RT           : [entree] Rayon du Telescope (en nombre de pixels, entier ou
;                 flottant). 
;
;  RHO          : [sortie] tableau des coordonnees radiales des points du tableau.
; 
;  PHI          : [sortie]     "    "      "   angulaires des points du tableau.
; 
; MASQUE        : [sortie] tableau (byte) qui vaut  1 sur le  telescope et 0
;                 ailleurs. 
; 
; /LEQ          : [entree optionnelle] sert uniquement a definir le MASQUE ;
;                 si ce drapeau est mis, le masque est defini par les points a
;                 distance <= RT (et >= OC*RT). Par defaut il est defini par
;                 les points a distance < RT (et >= OC*RT).
;  
; OC            : [entree optionnelle] sert uniquement a definir le MASQUE ;
;                 Occultation Centrale (reel dans [0,1[) ; 0 si absent.
; 
; POINTS_UTILES : [entree & sortie optionnelle] si ce mot-cle est present et
;                 non nul, on ne renvoie dans RHO et PHI que les points utiles
;                 (i.e., la où le MASQUE est non nul) RHO et PHI sont alors
;                 des tableaux 1-D et POINTS_UTILES contient en sortie les
;                 indices des points. Sinon (i.e., par defaut) RHO et PHI
;                 contiennent tous les points du tableau.
; 
; /ENTRE4PIXELS : [entree optionnelle] si ce drapeau est mis, les tableaux ont
;                 un nombre pair de points et le centre du telescope est pris
;                 entre les 4 pixels centraux du masque. Sinon, les tableaux
;                 ont un nombre de points impair et le centre est sur un pixel
;                 (defaut). Ce mot-cle est ignore si on impose la taille du
;                 masque par le mot-cle LARGEUR ou le centrage par CENTRE_X ou
;                 CENTRE_Y.
; 
; LARGEUR       : [entree ou sortie optionnelle, entier]
;                 Si ce mot-cle est present et contient une variable definie
;                 et > 0, alors on impose la taille (largeur) des tableaux
;                 rho, phi et masque. Sinon (i.e., mot-cle absent ou <=0),
;                 largeur est calculee automatiquement, cf. ENTRE4PIXELS.
;                 Si ce mot-cle est present mais contient une variable non
;                 definie ou <= 0, alors largeur _reçoit_ la valeur de la
;                 taille des tableaux rho, phi et masque.
;
; CENTRE_X ou _Y: [entree ou sortie optionnelle, reel] 
;                 si ce mot-cle est present et contient une variable definie
;                 et > 0, alors on impose la position (centre_x, centre_y) des
;                 tableaux rho, phi et masque. Si ce mot-cle est present mais
;                 contient une variable non definie ou <= 0, alors elle
;                 _reçoit_ la valeur du centre des tableaux rho, phi et
;                 masque.
; 
; /DOUBLE       : [entree] RHO et PHI en flottant double precision.  
; 
; /PROLIXE      : [entree] = verbose en VF.  Messages d'info a l'execution  ;
;
; /VERSION      : [entree] affichage de la version avant l'execution. 
; 
; /HELP         : [entree] affichage de la syntaxe et sortie du programme. 
; 
; ATTENTION : si  on impose "CENTRE_[X,Y]" ou "LARGEUR" il faut s'assurer que  
; la  taille "LARGEUR" est  assez  grande ! 
;
;  
;DIAGNOSTIC D'ERREUR : 
; 
;AUTEUR :
;   $Author: lmugnier $
;
;HISTORIQUE :  
;   $Log: polaire2.pro,v $
;   Revision 1.6  2012/10/03 08:41:24  lmugnier
;   Ajout d'un test pour verifier que LARGEUR, s'il est passe', est bien entier.
;
;   Revision 1.5  1999/03/04 18:18:43  mugnier
;   Mot-cle PROLIXE ajoute (pour ne plus avoir de message par defaut).
;
;   Revision 1.4  1998-02-05 11:40:42+01  mugnier
; V1.4, 05/02/98 - Laurent Mugnier. Doc amelioree + fichier mis sous RCS. 
; V1.3, 15/10/96 - Laurent Mugnier. Mots-cle CENTRE_X et CENTRE_Y. 
;                  Masque calcule de toute façon pour pouvoir passer le
;                  mot-cle avec une variable non definie.
; V1.2, 31/07/96 - Laurent  Mugnier. Mot-cle VERSION.  
; V1.1, 15/12/95 - Laurent Mugnier. D'apres le POLAIRE de Ludovic Meynadier
; (mais avec uniquement des mots-cles, plus d'options, de documentation, et
;  moins de bidouilles de tailles).
;
;-

on_error,2
IF keyword_set(version) THEN $
    message, "$Revision: 1.6 $, $Date: 2012/10/03 08:41:24 $", /INFO

IF (n_params() NE 0) OR NOT keyword_set(rt) OR keyword_set(help) THEN message, $
    "Usage : POLAIRE2, RT=rt [, RHO=rho, PHI=phi, MASQUE=masque, /LEQ, OC=oc]"+$
    "[, POINTS_UTILES=points_utiles] " + $
    "[, /ENTRE4PIXELS, LARGEUR=largeur, CENTRE_X = centre_x, " + $
    "CENTRE_Y = centre_y, /DOUBLE, /PROLIXE, /HELP, /VERSION])"

IF NOT keyword_set(oc) THEN oc = 0.0      

; on impose "largeur" par un mot-cle contenant une variable definie et > 0
set_largeur = byte(0)
IF (n_elements(largeur) NE 0) THEN BEGIN
   IF (largeur NE fix(largeur)) THEN $
      message, 'LARGEUR doit etre entier !' $
   ELSE $
      ; si largeur LT 0 la var *recoit* la valeur de LARGEUR en sortie :
      IF (largeur GT 0) THEN set_largeur = byte(1) 
ENDIF 

IF NOT set_largeur THEN BEGIN ; largeur "naturelle", i.e. minimale requise :
    IF keyword_set(entre4pixels) THEN BEGIN 
        largeur = 2 * round(rt) ; largeur paire
    ENDIF ELSE BEGIN 
        largeur = 2*fix(rt) + 1 ; largeur impaire; fix=floor en int pour nb>0
    ENDELSE    
ENDIF

; on impose "centre_[x,y]" par un mot-cle contenant une variable definie et > 0
set_centre_x = byte(0)
IF (n_elements(centre_x) NE 0) THEN  $
    IF (centre_x GT 0) THEN set_centre_x = byte(1)

IF NOT set_centre_x THEN BEGIN ; centre "naturel" :
    centre_x = float((largeur-1.0)/2)
ENDIF

set_centre_y = byte(0)
IF (n_elements(centre_y) NE 0) THEN  $
    IF (centre_y GT 0) THEN set_centre_y = byte(1)

IF NOT set_centre_y THEN BEGIN ; centre "naturel" :
    centre_y = float((largeur-1.0)/2)
ENDIF

IF keyword_set(prolixe) THEN BEGIN
    printf, -2, "% LARGEUR  vaut"+string(largeur)
    printf, -2, "% CENTRE_X vaut"+string(centre_x)
    printf, -2, "% CENTRE_Y vaut"+string(centre_y)
ENDIF

;calcul de x et y pour tous les points de la pupille
x = float(dindgen(largeur,largeur) mod largeur)
y = transpose(x)
x(*) = x - centre_x
y(*) = y - centre_y

IF keyword_set(double) THEN BEGIN 
    rho = double(sqrt(x^2+y^2)) / double(RT)
    phi = double(atan(y, x + (rho EQ 0.))) ; phi doit etre defini quand x=y=0
ENDIF ELSE BEGIN
    rho = float(sqrt(x^2+y^2)) / float(RT)
    phi = float(atan(y, x + (rho EQ 0.)))
ENDELSE

;IF (keyword_set(masque)) OR (keyword_set(points_utiles))THEN BEGIN
; Calculer Masque de toute façon pour pouvoir passer le mot-cle avec une
; variable non definie.

IF keyword_set(LEQ) THEN BEGIN
        masque =  ((rho LE 1.) AND (rho GE oc))
ENDIF ELSE BEGIN
        masque =  ((rho LT 1.) AND (rho GE oc))
ENDELSE 


IF keyword_set(points_utiles) THEN BEGIN
    IF keyword_set(prolixe) THEN  $
        printf, -2, "% RHO et PHI ne conservent que les points de la pupille"
    count = 0L
    points_utiles = where(masque NE byte(0), count)
    IF (count NE 0l) THEN BEGIN
        rho = rho(points_utiles)
        phi = phi(points_utiles)
        ENDIF ELSE message, 'Erreur : POINTS_UTILES contient 0 element !'
END

END
