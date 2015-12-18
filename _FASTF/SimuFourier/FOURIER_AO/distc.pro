; $Id: distc.pro,v 1.2 2005/05/25 12:50:10 mugnier Exp $

FUNCTION DistC, n, m, Cx=cx, Cy=cy

;Return a rectangular array in which each pixel = euclidian
;distance from the corner or from a given origin (Center),
;e.g., ((n-1)/2,(m-1)/2).
;+
; NAME:
;	DISTC - Version ameliorée de dist avec Centrage (origine) variable.
;
; PURPOSE:
;	Create a rectangular array in which each element is proportional
;	to its distance.  This array may be used for a variety
;	of purposes, including frequency-domain filtering and
;	making pretty pictures.
;
; CATEGORY:
;	Mathematics Routines.
;
; CALLING SEQUENCE:
;	Result = DISTC(N [, M] [,Cx=cx] [, Cy=cy])
;
; INPUTS:
;	N = number of columns in result.
;	M = number of rows in result.  If omitted, N is used to return
;		a square matrix.
;
; OUTPUTS:
;	Returns an (N,M) floating array in which:
;
;	R(i,j) = SQRT(F(i)^2 + G(j)^2)   where:
;   if neither Cx nor Cy is set:
;		 F(i) = i  IF 0 <= i <= n/2
;		      = n-i  IF i > n/2
;		 G(i) = i  IF 0 <= i <= m/2
;		      = m-i  IF i > m/2
;   if Cx or Cy is set:
;		 F(i) = i - Cx (n/2 by default)
;		 G(i) = i - Cy (m/2 by default)
;
; SIDE EFFECTS:
;	None.
;
; RESTRICTIONS:
;	None.
;
; PROCEDURE:
;	Straightforward.  The computation is done a row at a time.
;
; MODIFICATION HISTORY:
;   $Log: distc.pro,v $
;   Revision 1.2  2005/05/25 12:50:10  mugnier
;   Le test de présence sur Cx et Cy passe de keyword_set() à n_elements() pour
;   pouvoir passer Cx=0 ou Cy=0.
;
;   Revision 1.1  2005/01/10 16:35:18  mugnier
;   Initial revision
;
;               Jan. 2005 - put under version control.
;	            August 95 - Check absence of parameters.
;	            June 1995 - Changed name to DistC.
;	L. MUGNIER, July 1994 - Added Cx, Cy keywords to center array on point
;	                        (Cx,Cy)
;	DMS, July, 1992.      - Added M parameter to make non-square arrays.
; 	SMR, March 27, 1991   - Added the NOZERO keyword to increase efficiency.
;				            (Recommended by Wayne Landsman)
;	Very Old.
;-
on_error,2                      ;Return to caller if an error occurs

IF (n_params() EQ 0) THEN  $
    message, "Usage : Result = DISTC(N [, M] [,Cx=cx] [, Cy=cy])"

x=findgen(n)		             ;Make a row
IF n_elements(m) LE 0 THEN m = n ;define m if not given by user

IF (n_elements(cx) EQ 0L) AND (n_elements(cy) EQ 0L) THEN $
    x = (x < (n-x)) ^ 2 $        ;column squares
ELSE BEGIN 
    IF (n_elements(cx) EQ 0L) THEN cx = n/2.0
    IF (n_elements(cy) EQ 0L) THEN cy = m/2.0
    x = (x - cx)^2
ENDELSE 

    
a = FLTARR(n,m,/NOZERO)	;Make array

IF NOT (keyword_set(cx) OR keyword_set(cy))  THEN BEGIN 
    FOR i=0L, m/2 DO BEGIN           ;Row loop
        y = sqrt(x + i^2)            ;Euclidian distance
        a[0,i] = y	;Insert the row
        IF i ne 0 THEN a[0, m-i] = y ;Symmetrical
    ENDFOR
ENDIF ELSE BEGIN
    FOR i = 0L, m-1 DO BEGIN 
        y = sqrt(x + (i - cy)^2) 
        a[0, i] = y                  ;Insert the row 
    ENDFOR
ENDELSE

return, a
END 


