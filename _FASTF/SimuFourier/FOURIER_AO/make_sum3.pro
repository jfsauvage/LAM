;C'est un sous programme qui doit faire la somme sur n,m et L 
;pour ecrire les ssT etc...

;-------------------------------------------------------
;Y'a 3 cas a traiter : le cas SxSx - le cas SxSy - le cas SySy 
;et x2 parceque apres y'a les cas qd on regarde S1S2, mais
;ils sont inclus car l'exponentiel fera 1.
;-------------------------------------------------------
;on met un flag selon le cas a traiter, util pour dans la boucle.

;ICI on prend des frequences de coupures carr�es


FUNCTION make_sum3, f, fx, fy, arg_f, cn2, a, r0, L0, flag, nmax, dalpha, h_vrai

fc =  1./(2.*a) 

; what's that ? OK normalisation for slopes cov computation
; CF line 186 out(f_ind) = float(cst*r0^(-5./3.)*sum(f_ind))
cst = 0.0229

s = (size(f))(2)
out = complexarr(s, s)


;Deja, on veux travailler uniquement sur les frequences LT fc
f_ind =  where(abs(fx) LE fc AND abs(fy) LE fc,  count)
;et on a une double somme sur f_x et f_y, donc on a besoin de f_x et f_y    
;f_x =  f(f_ind)*cos(arg_f(f_ind))
;f_y =  f(f_ind)*sin(arg_f(f_ind))

;ce coeff il pondere toutes les produits SxSy etc... sauf cas limites ou (m,n)=0
coeff = (2*!pi*fx*fy*sinc(a*fx)*sinc(a*fy))^2.

;c'est lui qui recoit le resultat de la somme sur toutes les contributions de l'Aliasing
sum = complexarr(s, s)
;et somme sur les couches.

;on doit sommer sur les couches
;----------------------------------------------------------
FOR nel = 0, n_elements(Cn2)-1 DO BEGIN
summ1 = complexarr(s, s)


;Cas ou on va calculer la correlation entre des mesures selon X et selon X
;-------------------------------------------------------
IF flag EQ 'xx' THEN BEGIN 
;rint, 'xx'
FOR n=-nmax, nmax DO BEGIN
FOR m=-nmax, nmax DO BEGIN                

w =  where(abs(fx) LE fc AND abs(fy) LE fc AND fy NE 0.,  count)
;w =  where(f LT fc,  count)
;f_x =  f(f_ind)*cos(arg_f(f_ind))
;f_y =  f(f_ind)*sin(arg_f(f_ind))

IF n NE 0 OR m NE 0 THEN BEGIN  
;print, 'n', round(n), '  m', round(m)        
f_nx = fx - 2.*n*fc ;
f_my = fy - 2.*m*fc ;
f_nm_sq = f_nx*f_nx + f_my*f_my ;
                
spectrum = (f_nm_sq+(1/L0)^2)^(-11./6.)
;ca c'est la DSP(fm) pour l'altiude L.

;si y'a de l'erreur temporelle :
;temp = sinc(ti*v(cpt)*sqrt(f_nm_sq)*cos(atan(f_my,f_nx)-arg_v(cpt)))
;sinon pour l'instant on considere que non :
temp = 1.
;exponentielle complexe qui arrive qd on regarde 2 directions differentes :
ff_x = f_nx*dalpha(0)*h_vrai(nel)*60./206265.
ff_y = f_my*dalpha(1)*h_vrai(nel)*60./206265.
angle = exp(complex(0, 1)*2*(ff_x+ff_y)*!pi)
;info, angle
;C'est ici que ca change selon le FLAG :
summ1[w] = summ1[w] + spectrum[w]*temp*temp*(1./f_my[w])^2.*coeff[w]*angle[w] ;     <====== Modif!!!!      

;On traite le cas particuler :
w =  where(abs(fx) LE fc AND abs(fy) LE fc AND fy EQ 0.,  count)
coeff2 = (2*!pi*fx*sinc(a*fx)*sinc(a*f_my))^2.
;coeff2 = (2*!pi*fx*sinc(a*fx))^2.
summ1[w] = summ1[w] + spectrum[w]*temp*temp*coeff2[w]*angle[w] ;       ;     <====== Modif!!!!     

ENDIF     
     
ENDFOR
ENDFOR



ENDIF


;Cas ou on va calculer la correlation entre des mesures selon X et selon Y
;-------------------------------------------------------
IF flag EQ 'xy' THEN BEGIN 
;print, 'xy'
FOR n=-nmax, nmax DO BEGIN
FOR m=-nmax, nmax DO BEGIN                

w =  where(abs(fx) LE fc AND abs(fy) LE fc AND fx NE 0. AND fy NE 0.,  count)

IF n NE 0 OR m NE 0 THEN BEGIN  
                
f_nx = fx - 2.*n*fc ;
f_my = fy - 2.*m*fc ;
f_nm_sq = f_nx*f_nx + f_my*f_my ;
;exponentielle complexe qui arrive qd on regarde 2 directions differentes :
ff_x = f_nx*dalpha(0)*h_vrai(nel)*60./206265.
ff_y = f_my*dalpha(1)*h_vrai(nel)*60./206265.
angle = exp(complex(0, 1)*2*(ff_x+ff_y)*!pi)                
spectrum = (f_nm_sq+(1/L0)^2)^(-11./6.)
;ca c'est la DSP(fm)

;si y'a de l'erreur temporelle :
;temp = sinc(ti*v(cpt)*sqrt(f_nm_sq)*cos(atan(f_my,f_nx)-arg_v(cpt)))
;sinon
temp = 1.

;on peut l'ecrire aussi comme :
summ1[w] = summ1[w] + spectrum[w]*temp*temp*(1./(f_my[w]*f_nx[w]))*coeff[w]*angle[w] ;    ;     <====== Modif!!!!        

;Et le cas particulier ?
;C'est en (0,0) on s'en fout...
; w =  where(abs(fx) LT fc AND abs(fy) LT fc AND fx EQ 0. AND fy EQ 0.,  count)
; coeff2 = (2*!pi*sinc(a*f_my)*sinc(a*f_nx))^2.*f_nx*f_my
; summ1(w) = summ1(w) + spectrum(w)*temp*temp*coeff2(w)*angle(w)
                                
ENDIF

ENDFOR
ENDFOR


ENDIF

;Cas ou on va calculer la correlation entre des mesures selon Y et selon Y
;-------------------------------------------------------
IF flag EQ 'yy' THEN BEGIN 
;print, 'yy'

FOR n=-nmax, nmax DO BEGIN
FOR m=-nmax, nmax DO BEGIN                

w =  where(abs(fx) LE fc AND abs(fy) LE fc AND fx NE 0.,  count)

IF n NE 0 OR m NE 0 THEN BEGIN  
                
f_nx = fx - 2.*n*fc ;
f_my = fy - 2.*m*fc ;
f_nm_sq = f_nx*f_nx + f_my*f_my ;
                
spectrum = (f_nm_sq+(1/L0)^2)^(-11./6.)
;ca c'est la DSP(fm)
;exponentielle complexe qui arrive qd on regarde 2 directions differentes :
ff_x = f_nx*dalpha(0)*h_vrai(nel)*60./206265.
ff_y = f_my*dalpha(1)*h_vrai(nel)*60./206265.
angle = exp(complex(0, 1)*2*(ff_x+ff_y)*!pi)
;si y'a de l'erreur temporelle :
;temp = sinc(ti*v(cpt)*sqrt(f_nm_sq)*cos(atan(f_my,f_nx)-arg_v(cpt)))
;sinon
temp = 1.

;on peut l'ecrire aussi comme :
summ1[w] = summ1[w] + spectrum[w]*temp*temp*(1./f_nx[w])^2*coeff[w]*angle[w] ;     ;     <====== Modif!!!!       
;On traite le cas particuler :
w =  where(abs(fx) LE fc AND abs(fy) LE fc AND fx EQ 0.,  count)
coeff2 = (2*!pi*fy*sinc(a*f_nx)*sinc(a*fy))^2.
;coeff2 = (2*!pi*fy*sinc(a*fy))^2.
summ1[w] = summ1[w] + spectrum[w]*temp*temp*coeff2[w]*angle[w] ;       ;     <====== Modif!!!!          
ENDIF
         
ENDFOR
ENDFOR


ENDIF


;et ici on vient faire la somme sur les couches :  
sum = sum + Cn2(nel)*summ1 ;
ENDFOR

;AU FINAL, SUM CONTIENT LA CORRELATION SSt QUIL FAUT APPLIQUER DANS LE RECONSTRUCTEUR.
;h2 = 1.
out(f_ind) = float(cst*r0^(-5./3.)*sum(f_ind))
out(0) = 0.
;stop
return, out

END
