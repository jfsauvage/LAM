FUNCTION norme, x, max = max, centre = centre, energy = energy, silent = silent
; NORM retourne x normalis� par son total, sauf si nul.
; si max, retourne x normalis� par son max...
; si centre, retourne x normalis� par x[centre]...attention avec NP pair ou
; impair
; si energy retourne x normalis� par son energie (�nergie unitaire)

IF keyword_set(max) + keyword_set(centre) +keyword_set(energy) GT 1 THEN BEGIN
    print, 'Options incompatibles...'
    return, x
ENDIF

IF (NOT(keyword_set(energy)) AND NOT(keyword_set(max)) AND NOT(keyword_set(centre)) AND total(x) EQ 0) THEN BEGIN
    IF NOT keyword_set(silent) then    print, 'normalisation total(.) = 1.'
    print, 'total nul, je ne normalise pas'
    return, x
ENDIF

IF keyword_set(max) THEN begin
    IF NOT keyword_set(silent) then    print, 'normalisation max(.) = 1.'
    IF max(x) EQ 0 THEN BEGIN
       print, 'maximum nul, je ne normalise pas'
       return, x
    ENDIF
ENDIF

IF keyword_set(centre) THEN begin
    IF NOT keyword_set(silent) then    print, 'normalisation centre(.) = 1.'
    centros = x[(size(x))[1]/2., (size(x))[1]/2.]
    IF centros EQ 0 THEN BEGIN
    print, 'valeur au centre nulle, je ne normalise pas'
        return, x
    ENDIF
ENDIF

IF keyword_set(energy) THEN BEGIN
    IF NOT keyword_set(silent) then    print, 'normalisation energy(.) = 1.'
    norme = total(x^2)
    IF norme EQ 0 THEN BEGIN
        print, 'Energie nulle, je ne normalise pas'
        return, x
    ENDIF
ENDIF

; IF (total(x) LT 1.e-20 OR max(x) LT 1.e-20 OR centros LT 1.e-20 OR norme LT 1.e-20) $
;     THEN print, 'Attention, norme(x) = ', strc(total(x)), ', ' + $
;     'risque d''explosion ;o)'

IF keyword_set(max) THEN out = x/max(double(x))
IF keyword_set(centre) THEN out = x/centros
IF keyword_set(energy) THEN out = x/sqrt(total(x^2, /double))
IF keyword_set(centre) + keyword_set(max) + keyword_set(energy) EQ 0 THEN out = x/total(x, /double)

return, out
END
