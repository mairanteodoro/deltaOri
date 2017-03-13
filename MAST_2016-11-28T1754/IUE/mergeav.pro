;
; -> EXTRACTED FROM
;     https://archive.stsci.edu/pub/iue/software/iuedac/procedures/
;
;     https://archive.stsci.edu/iue/manual/dacguide/node76.html
;
;******************************************************************************
;+
;*NAME:
;
;  	MERGEAV  (formerly IUEMERGE_AVERGE)
;
;*CLASS:
;
;	Spectral Data Reduction
;
;*CATEGORY:
;
;*PURPOSE:
;
;  	Combines two wave-flux-flag (epsilon or nuflags) groups.
;
;*CALLING SEQUENCE:
;
;   	MERGEAV,WAVE1,FLUX1,FLAG1,WAVE2,FLUX2,FLAG2,SPNT,WAVE3,FLUX3,   $
;          FLAG3,/newsips
;
;*PARAMETERS:
;
;    	WAVE1	(REQ) (I) (1) (F D)
;		First wavelength vector, monotonically increasing.
;
;    	FLUX1	(REQ) (I) (1) (F D)
;		First flux vector, with the fluxes for wave1.
;
;	FLAG1	(REQ) (I) (1) (I F D)
;		First epsilon (IUESIPS) or nu flag (NEWSIPS) vector.
;
;    	WAVE2	(REQ) (I) (1) (F D)
;		Second wavelength vector.
;
;    	FLUX2	(REQ) (I) (1) (F D)
;		Fluxes for WAVE2.
;
;	FLAG2	(REQ) (I) (1) (I F D)
;		Epsilons (IUESIPS) or nu flags (NEWSIPS) for WAVE2.
;
;    	SPNT	(REQ) (O) (1) (I)
;		Splice points (intarr(2)).
;
;    	WAVE3	(REQ) (O) (1) (F D)
;		Combined wavelength vector.
;
;    	FLUX3	(REQ) (O) (1) (F D)
;		Combined and averaged flux vector.
;
;	FLAG3	(REQ) (O) (1) (I F D)
;		Combined epsilons (IUESIPS) or nu flags (NEWSIPS).
;
;	NEWSIPS	(KEY) (I) (0) (I)
;		Keyword when set indicates that the FLAG1, FLAG2, and FLAG3
;               parameters should be treated as NEWSIPS nu flags
;
;*EXAMPLES:
;
;*SYSTEM VARIABLES USED:
;
;	None.
;
;*INTERACTIVE INPUT:
;
;	None.
;
;*SUBROUTINES CALLED:
;
;	TABINV
;	PARCHECK
;
;*FILES USED:
;
;	none
;
;*SIDE EFFECTS:
;
;*RESTRICTIONS:
;
;    	The calling routine is responsible for error checking.
;
;*PROCEDURE:
;
;   	Combines two wave-flux-flag groups as follows:  in the regions where
;	wave1 and wave2 do not overlap, use all values of wave1-flux1-flag1 and
;	wave2-flux2-flag2.  In the region of overlap, resample flux2 and flag2
;	by linear interpolation to the same wavelength values as wave1 and then
;	average flux1 with flux2 at those wavelengths.  For IUESIPS data (the
;	NEWSIPS flag is not set or equals 0) the greatest value of epsilons in
;	that region (flux points with negative quality flags having been
;	rejected), or the sum of the two quality flags, whichever is greater.
;	For NEWSIPS data (the NEWSIPS flag is set), the values of the nu flags
;	are combined.
;
;*I_HELP nn:
;
;*NOTES:
;
;*MODIFICATION HISTORY:
;
;	 7-17-89 RWT mods for SUN IDL, remove all EXTRAC and LOOKUP commands
;	22 Jul 91 LLT clean up, add parcheck, update prolog, fix assignment of
;		      eps2 and spnt, add many inline comments. tested on VAX
;	23 Jul 91 PJL tested on SUN; updated prolog
;	18 Sep 91 PJL corrected typo
;	27 Dec 91 RWT allow eps arrays to stay integers (i.e., to reduce memory)
;	13 Sep 93  PJL  added NEWSIPS keyword; changed EPS documentation
;			references to FLAG; added NEWSIPS nu flag combination
;       20 Jun 94  PJL  renamed from IUEMERGE_AVERGE to MERGEAV
;
;-
;******************************************************************************
 pro mergeav,wave1,flux1,eps1,wave2,flux2,eps2,spnt,wave3,flux3,eps3,   $
        newsips=newsips
;
 npar = n_params(0)
 if (npar eq 0) then begin
    print,'MERGEAV,WAVE1,FLUX1,EPS1,WAVE2,FLUX2,EPS2,SPNT,WAVE3,' +   $
       'FLUX3,EPS3,/newsips'
    retall
 endif  ; npar eq 0
 parcheck,npar,10,'MERGEAV'
;
 if (keyword_set(newsips)) then newsips = 1 else newsips = 0
;
 s1   = n_elements(wave1)
 s2   = n_elements(wave2)
 wmin = wave1(0) > wave2(0)          ; Min. wavelength in overlap region
 wmax = wave1(s1-1) < wave2(s2-1)    ; Max. wavelength in overlap region
 nw1  = total(wave1 lt wmin)                    ; # points shortward of wmin
 nw2  = total((wave2 le wmax)*(wave2 ge wmin))  ; # points in overlap region
 nw3  = total(wave2 lt wmin)                    ; # points shortward of wmin
 nw4  = total((wave1 le wmax)*(wave1 ge wmin))  ; # points in overlap region
;
;  create output vectors
;
 eps3  = intarr(s1+s2-nw2)
 wave3 = eps3 * wave1(0)               ; Output vectors are length of input
 flux3 = wave3                         ; wavelength vectors combined, less the
 eps3  = eps3 * eps1(0)                ; # of pts in wave2 in overlap region.
;
;  insert data into vectors
;  note flux is weighted by flag vector
;
 if (nw3 gt 0) then begin            ; If wave2 is shortward of wave1
    wave3(0) = wave2(0:nw3-1)        ; Insert wave2 in beginning of wave3
    tmp = flux2*(eps2>0)             ; Screen out fluxes with neg. qual. flags
    flux3(0) = tmp(0:nw3-1)          ; Insert fluxes into flux3
    eps3(0)  =  eps2(0:nw3-1)  ; this line fixed by LLT on 22 July 1991
 endif ; nw3 gt 0
;
 if ((nw3+nw2) lt s2) then begin  ; If #points in wave2 NOT LONGWARD of overlap
                                  ; region are less than n_elements(wave2)
    wave3(nw3+s1) = wave2(nw3+nw2:s2-1) ; Insert wave2 longward of overlap
                                    ; region in wave3 longward of overlap region
    tmp = flux2*(eps2>0)        ; Screen out fluxes with negative quality flags
    flux3(nw3+s1) = tmp(nw3+nw2:s2-1)   ; Insert fluxes into flux3
    eps3(nw3+s1)  =  eps2(nw3+nw2:s2-1) ; Inert quality flags into eps3
 endif  ; (nw3+nw2) lt s2
;
 wave3(nw3) = wave1       ; Insert wave1 in beginning of wave3 (or at wmin
                          ; if wave2 is shortward of wave1)
 flux3(nw3) = flux1*(eps1>0)      ; Insert flux1 into flux3
 eps3(nw3)  = eps1                ; Inert eps1 into eps3
;
;  interpolate and average in overlap region
;
 if (nw2 gt 0) then begin  ; have overlap region to average
    tmp = nw1>nw3          ; Use the vectors that have data shortward of wmin
    tabinv,wave2,wave3(tmp:tmp+nw4-1),ro   ; Interpolate indices - prepare to
                                           ; interpolate fluxes
    io = fix(ro)                      ; Truncate indices from tabinv
    ro = ro-io                     ; Fractional part of indices from tabinv
    tmp2 = flux2*(eps2>0)   ; Screen out fluxes (second set) with
                            ; negative quality flags
    fo = tmp2(io) * (1-ro) + tmp2(io+1) * ro  ; Interpolate fluxes (second set)
    eo = eps2(io) * (1-ro) + eps2(io+1) * ro  ; Interpolate quality flags
                                              ; (second set)
    flux3(tmp) = fo + flux3(tmp:tmp+nw4-1)  ; Add interpolated fluxes to fluxes
                                            ; (first set) already in flux3
    ro = eps3(tmp:tmp+nw4-1)      ; Quality flags (first set) already in eps3
;
;  NEWSIPS - combine eps1 and interpolated eps  (second set)
;  IUESIPS - Use greatest of eps1, interpolated eps (second set), or sum of
;            the two.
;
    if (newsips) then eps3(tmp) = -( (abs(eo)) or (abs(ro)) ) else   $
       eps3(tmp) = (eo+ro)>eo>ro
 endif  ; averge overlap region
 flux3 = flux3 / (eps3>0 + (eps3 le 0))          ; convert back to pure flux
;
;  return splice points
;
 spnt    = intarr(2)
 spnt(0) = nw1>nw3        ; This was done by LLT on 22 July 1991 - this is
                          ; the index of wmin
 spnt(1) = spnt(0)+nw4-1    ; This was done by LLT on 22 July 1991 - this is
                            ; the index of wmax
;
 return
 end  ; mergeav
