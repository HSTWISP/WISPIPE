pro new_seflag,JHN,det_H,det_J,det_comb,field,path0
;
; Set "det_J" and "det_comb" to 0 if J image is not present.
; Set "det_H" and "det_comb" to 0 if H image is not present.

;pro seflag,filter,det,path0

; This program flags problematic objects in the SExtracted catalogs
; 1. objects close to the edge
; 2. objects that contains > limitr (%) of weight=0 pixels
; 3. objects that are possibly due to star diffractions or persistence
;
;---------------------------------------------------
; Changes with respect to old version (seflag.pro):
;---------------------------------------------------
; The new version of seflag is written to flag sources extracted in
; Dual Imag mode in "new_match_cat2.pro". The difference with
; seflag.pro is that in this new version the J and H catalogs can't be
; considered separately and if a source is removed in one catalog it
; must be removed from the other too, in order to keep the catalogs in
; agreement each other.
;
; The case in which only the H image and catalog are present is
; automatically managed by the program
;
; Another important change with respect to seflag.pro is that the
; dimension of the images has changed because of a new external frame of
; zero values (made to have J-H images with identical sizes as a
; result of tweakreg and astrodrizzle). In this new version, the
; borders of the internal image are automatically identified (so that
; if they change the pipeline doesn't enter in an esistential
; crisis) and not defined once forever, as in the old pipeline
;
; Sources are flagged as:
;-  9 if close to the border;
;-  90 if they are in areas where >50% of the pixels have wheight=0;
;-  900 if theyr shape is very elongated (spikes).
; The flags can be summed up, so that, for example, 999 indicates a
; source with all the previous problems at the same time, while 909
; indicates that a source is in the border and with an elongated shape.
;
; In old seflag program, the flags were 9, 99, and 999. If one source
; was already flagged, it was not flagged anymore and the three flags
; were actually summed each other only for the first and second flag
; (was it a bug?)
;
; At the end of the flagging procedure, two catalogs are created for
; each filter and for each extraction strategy (Dual or Single image
; mode): one "full" and one "cleaned" catalog. While the full catalog
; is a copy of the input catalog, in the cleaned one, sources flagged
; are removed. 
;
; OUTPUTS:
; 
; - F14(6)0W_SIM_full.cat     (always)
; - F14(6)0W_SIM_cleaned.cat  (always)
; - F110W_SIM_full.cat        (if a J image is present)
; - F110W_SIM_cleaned.cat     (if a J image is present)
; - F14(6)0W_DIM_full.cat     (if a J image is present)
; - F14(6)0W_DIM_cleaned.cat  (if a J image is present)
; - F110W_DIM_full.cat        (if a J image is present)
; - F110W_DIM_cleaned.cat     (if a J image is present)
;
; By Ivano Baronhelli July 2016
;--------------------------------------------------
;ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ

TS='0' ; normal run
;TS='1' ; Test phase


print, "----------------------------"
print, "       new_seflag.pro       "
print, "----------------------------"

;NORMAL RUN
if TS eq '0' then begin
pathse = path0+"/aXe/"+field+'/SEX/'
endif

; TEST 
if TS eq '1' then begin
pathse='../images/SEX/'
endif


; J
img_Jn=pathse+'F110W_sci.fits'
cat_J_SIMn=pathse+'F110W_'+det_J+'.cat'
cat_J_DIMn=pathse+'F110W_DIM_'+det_comb+'.cat'

; H
; Single image mode catalog
cat_H_SIMn=pathse+'F1'+JHN+'0W_'+det_H+'.cat'
; Dual image mode catalog
cat_H_DIMn=pathse+'F1'+JHN+'0W_DIM_'+det_comb+'.cat'
; Image sci
img_Hn=pathse+'F1'+JHN+'0W_sci.fits'

; ----------------------------------------
; HOW MANY FILTERS ARE WE CONSIDERING?
; --> Check J-H image exixtence
ISTHEREJ=FILE_TEST(img_Jn)
ISTHEREH=FILE_TEST(img_Hn)
; ----------------------------------------


IF ISTHEREH gt 0 THEN BEGIN
 ; read H Single Image mode catalog
 readcol,cat_H_SIMn,ID_HSIM,X_IM_HSIM,Y_IM_HSIM,A_IM_HSIM,B_IM_HSIM,THETA_IM_HSIM,X_WO_HSIM,Y_WO_HSIM,A_WO_HSIM,B_WO_HSIM,THETA_WO_HSIM,MAG_HSIM, MAGERR_HSIM ,CLASS_STAR_HSIM ,FLAGS_HSIM,format='i,f,f,f,f,f,f,f,f,f,f,f,f,f,i',skipline=15,/silent
 ; Set magnitude uncertainties higher than 99 to 99 (to prevent some
 ; bugs ex. Par 1)
 MAGERR_HSIM[where(MAGERR_HSIM gt 99)]=99.
 ; Define two vectors used for rejection
 ; FLAGS
 FLAG2_HSIM=fltarr(n_elements(ID_HSIM))
 ; RATES
 RATE0_HSIM=fltarr(n_elements(ID_HSIM))
ENDIF


IF ISTHEREJ gt 0 THEN BEGIN
 ; read J Single Image mode catalog
 readcol,cat_J_SIMn,ID_JSIM,X_IM_JSIM,Y_IM_JSIM,A_IM_JSIM,B_IM_JSIM,THETA_IM_JSIM,X_WO_JSIM,Y_WO_JSIM,A_WO_JSIM,B_WO_JSIM,THETA_WO_JSIM,MAG_JSIM, MAGERR_JSIM ,CLASS_STAR_JSIM ,FLAGS_JSIM,format='i,f,f,f,f,f,f,f,f,f,f,f,f,f,i',skipline=15,/silent
 ; Set magnitude uncertainties higher than 99 to 99 (to prevent some
 ; bugs ex. Par 1)
 MAGERR_JSIM[where(MAGERR_JSIM gt 99.)]=99.
 ; Define two vectors used for rejection
 ; FLAGS
 FLAG2_JSIM=fltarr(n_elements(ID_JSIM))
 ; RATES
 RATE0_JSIM=fltarr(n_elements(ID_JSIM))
ENDIF


IF ISTHEREJ gt 0 and ISTHEREH gt 0 then begin
 ; Dual image mode
 ;@@@@@@@@@@@@@@@@@@@@
 NIMG=2
 ; CATALOG J Dual Image mode
 readcol,cat_J_DIMn,ID_JDIM,X_IM_JDIM,Y_IM_JDIM,A_IM_JDIM,B_IM_JDIM,THETA_IM_JDIM,X_WO_JDIM,Y_WO_JDIM,A_WO_JDIM,B_WO_JDIM,THETA_WO_JDIM,MAG_JDIM, MAGERR_JDIM ,CLASS_STAR_JDIM ,FLAGS_JDIM,format='i,f,f,f,f,f,f,f,f,f,f,f,f,f,i',skipline=15,/silent
 ; CATALOG H Dual Image mode
 readcol,cat_H_DIMn,ID_HDIM,X_IM_HDIM,Y_IM_HDIM,A_IM_HDIM,B_IM_HDIM,THETA_IM_HDIM,X_WO_HDIM,Y_WO_HDIM,A_WO_HDIM,B_WO_HDIM,THETA_WO_HDIM,MAG_HDIM, MAGERR_HDIM ,CLASS_STAR_HDIM ,FLAGS_HDIM,format='i,f,f,f,f,f,f,f,f,f,f,f,f,f,i',skipline=15,/silent
 ; Set magnitude uncertainties higher than 99 to 99 (to prevent some
 ; bugs ex. Par 1)
 MAGERR_JDIM[where(MAGERR_JDIM gt 99.)]=99.
 MAGERR_HDIM[where(MAGERR_HDIM gt 99.)]=99.
 ; FLAGS
 FLAG2_DIM=fltarr(n_elements(ID_JDIM))
 ; RATES
 RATE0_DIM=fltarr(n_elements(ID_JDIM))

ENDIF ELSE BEGIN
 ; Single image mode
 ;@@@@@@@@@@@@@@@@@@@@
 NIMG=1
ENDELSE

; ----------------------------------------


; xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
; DETECT IMAGES SIZE (external frame excluded)
; xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
; conservative values are used for the borders
; (the area common to both images
; is considered to cut the catalog)
;-----------------------
XMIN_vect=lonarr(NIMG)
XMAX_vect=lonarr(NIMG)
YMIN_vect=lonarr(NIMG)
YMAX_vect=lonarr(NIMG)
;-----------------------
IMN=0
PP=0
WHILE IMN LT 2 DO BEGIN
 IF IMN eq 0 THEN imagename=img_Hn
 IF IMN eq 1 THEN imagename=img_Jn
 
 DOIT=FILE_TEST(imagename)

 if DOIT eq '1' then begin
  imaget=mrdfits(imagename,0,HD_sci)
  IF n_elements(imaget) le 1 then imaget=mrdfits(imagename,1,HD_sci)
  NX=n_elements(imaget[*,0]) ; X size
  NY=n_elements(imaget[0,*]) ; Y size
  REFX=NX/2 ; Approximated central pixel (x)
  REFY=NY/2 ; Approximated central pixel (y)
  ; Get real image size excluding frame
  FX_min=0L
  FX_max=long(NX-1)
  FY_min=0L
  FY_max=long(NY-1)
  tst=700 ; test strip will be 2*tst

  ;------------------------------------------------------------
  internal_border=10 ; pixels to be excluded [pixels]
  ;------------------------------------------------------------

  ; DETECT X FRAME LIMITS FX_min,FX_max
  MINOK=0
  MAXOK=0
  NN=0L
  while NN lt REFX do begin
   TESTSTRIP1=imaget[NN,REFY-tst:REFY+tst]
   TESTSTRIP2=imaget[(NX-1)-NN,REFY-tst:REFY+tst]
   Check_IDX1=where(TESTSTRIP1 eq 0.)
   Check_IDX2=where(TESTSTRIP2 eq 0.)
   if n_elements(Check_IDX1) lt 0.2*float(n_elements(TESTSTRIP1)) and MINOK lt 1 then begin
    FX_min=NN+internal_border
    MINOK=1
   endif
   if n_elements(Check_IDX2) lt 0.2*float(n_elements(TESTSTRIP2)) and MAXOK lt 1 then begin
    FX_max=(NX-1)-NN-internal_border
    MAXOK=1
   endif
   NN=NN+1
  endwhile

  ; DETECT Y FRAME LIMITS FY_min,FY_max
  MINOK=0
  MAXOK=0
  NN=0L
  while NN lt REFY do begin
   TESTSTRIP1=imaget[REFX-tst:REFX+tst,NN]
   TESTSTRIP2=imaget[REFX-tst:REFX+tst,(NY-1)-NN]
   Check_IDX1=where(TESTSTRIP1 eq 0.)
   Check_IDX2=where(TESTSTRIP2 eq 0.)
   if n_elements(Check_IDX1) lt 0.2*float(n_elements(TESTSTRIP1)) and MINOK lt 1  then begin
    FY_min=NN+internal_border
    MINOK=1
   endif
   if n_elements(Check_IDX2) lt 0.2*float(n_elements(TESTSTRIP2)) and MAXOK lt 1  then begin
    FY_max=(NY-1)-NN-internal_border
    MAXOK=1
   endif 
   NN=NN+1 
  endwhile

  XMIN_vect[PP]=FX_min
  XMAX_vect[PP]=FX_max
  YMIN_vect[PP]=FY_min
  YMAX_vect[PP]=FY_max

  PP=PP+1
 endif ; if DOIT eq '1' then begin
 IMN=IMN+1
ENDWHILE

; xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

;--------------------------------------------------------------
; FLAG SOURCES CLOSE TO THE EDGES IN ALL THE AVAILABLE CATALOGS
;--------------------------------------------------------------
; Flag objects at the edge as '9'

; FLAG FOR H-SINGLE IMAGE MODE
IF ISTHEREH gt 0 then begin
 FLAG_A_HSIM=where(X_IM_HSIM le max(XMIN_vect) or X_IM_HSIM ge min(XMAX_vect) or Y_IM_HSIM le max(YMIN_vect) or Y_IM_HSIM ge min(YMAX_vect))
 IF FLAG_A_HSIM[0] NE -1 THEN FLAG2_HSIM[FLAG_A_HSIM]=9
ENDIF

; FLAG FOR J-SINGLE IMAGE MODE
IF ISTHEREJ gt 0 then begin
 FLAG_A_JSIM=where(X_IM_JSIM le max(XMIN_vect) or X_IM_JSIM ge min(XMAX_vect) or Y_IM_JSIM le max(YMIN_vect) or Y_IM_JSIM ge min(YMAX_vect))
 IF FLAG_A_JSIM[0] NE -1 THEN FLAG2_JSIM[FLAG_A_JSIM]=9
ENDIF

; FLAG FOR DUAL IMAGE MODE
IF ISTHEREJ gt 0 and ISTHEREH gt 0 then begin
 FLAG_A_DIM=where(X_IM_JDIM le max(XMIN_vect) or X_IM_JDIM ge min(XMAX_vect) or Y_IM_JDIM le max(YMIN_vect) or Y_IM_JDIM ge min(YMAX_vect))
 IF FLAG_A_DIM[0] NE -1 THEN FLAG2_DIM[FLAG_A_DIM]=9

 ; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 ; If a source is flagged in one catalog, it must be flagged in all the
 ; other catalogs too, otherwise a source would be set to magnitude
 ; limit while it is actually detected.
 dmax=1.5d-4 ; maximum distance for matching (degrees)
 SEARCH_RAD=3600.*dmax
 Y1=0
 Y2=0
 Y3=0
 Y4=0
 count=0L
 while Y1[0] ne -1 or Y2[0] ne -1 or Y3[0] ne -1 or Y4[0] ne -1 do begin
  ; DIM - JSIM
  cccpro,X_WO_JDIM[FLAG_A_DIM],Y_WO_JDIM[FLAG_A_DIM],X_WO_JSIM,Y_WO_JSIM,culo,JtoDIM_IDX,dt=SEARCH_RAD
             Y1=where(FLAG2_JSIM[JtoDIM_IDX] ne 9)
  IF Y1[0] ne -1 THEN FLAG2_JSIM[JtoDIM_IDX[Y1]]=9
  ; DIM - HSIM
  cccpro,X_WO_JDIM[FLAG_A_DIM],Y_WO_JDIM[FLAG_A_DIM],X_WO_HSIM,Y_WO_HSIM,culo,HtoDIM_IDX,dt=SEARCH_RAD
             Y2=where(FLAG2_HSIM[HtoDIM_IDX] ne 9)
  IF Y2[0] ne -1 THEN FLAG2_HSIM[HtoDIM_IDX[Y2]]=9
  ; JSIM - DIM
  cccpro,X_WO_JSIM[FLAG_A_JSIM],Y_WO_JSIM[FLAG_A_JSIM],X_WO_JDIM,Y_WO_JDIM,culo,DIMtoJ_IDX2,dt=SEARCH_RAD
             Y3=where(FLAG2_DIM[DIMtoJ_IDX2] ne 9)
  IF Y3[0] ne -1 THEN FLAG2_DIM[DIMtoJ_IDX2[Y3]]=9
  ; HSIM - DIM
  cccpro,X_WO_HSIM[FLAG_A_HSIM],Y_WO_HSIM[FLAG_A_HSIM],X_WO_JDIM,Y_WO_JDIM,culo,DIMtoH_IDX2,dt=SEARCH_RAD
             Y4=where(FLAG2_DIM[DIMtoH_IDX2] ne 9)
  IF Y4[0] ne -1 THEN FLAG2_DIM[DIMtoH_IDX2[Y4]]=9

  print, 'Cleaning (-1 at the end of the process)'
  print, "--------"
  print, "Y1",Y1
  print, "--------"
  print, "Y2",Y2
  print, "--------" 
  print, "Y3",Y3
  print, "--------"
  print, "Y4",Y4
  print, "--------"
  IF COUNT eq 10 then BEGIN
   PRINT, 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
   PRINT, ' WARNING (new_seflag): clean the sources at the border.  '
   PRINT, " Some sources couldn't be cleaned in all the catalogs"
   PRINT, " More than 10 attemps were made"
   PRINT, 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
  ENDIF
  ; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 endwhile
ENDIF




;--------------------------------------------------------------
; FLAG OBJECTS CONTAINING >50% OF WEIGHT=0 PIXELS
;--------------------------------------------------------------
; flag objects as 90 + OLD FLAG (90+0 OR 90+9)
; NOTE: 
; J and H weights images are used for the J and H single
; extractions. Instead, the combined wht image is used for the J and H
; Dual image mode extractions
; a source is excluded if more than 50% of underlying pixels in the wht
; image used are null values.


; Read weight images

; H IMAGE PRESENT xxxxxxxxxx
IF ISTHEREH gt 0 then imwht_HSIM = mrdfits(pathse+'F1'+JHN+'0W_wht.fits')
; J IMAGE PRESENT xxxxxxxxxx
IF ISTHEREJ gt 0 then imwht_JSIM = mrdfits(pathse+'F110W_wht.fits')
; J+H IMAGES PRESENT xxxxxxxxxxxxxxx
IF ISTHEREJ gt 0 and ISTHEREH gt 0 then imwht_DIM = mrdfits(pathse+'JH_combined_wht.fits')

limitr=0.5

; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; H IMAGE PRESENT xxxxxxxxxx
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IF ISTHEREH gt 0 then begin
 i=0L
 while i lt n_elements(ID_HSIM) do begin
  ; In old pipeline, flag was 99 and already flagged were skipped. 
  ; Here instead we add the two flag (9+90=99 or 0+90=90)
  ;if FLAGS_HSIM[i] lt 4 and FLAG2_HSIM[i] eq 0 then begin
  XSIZE=n_elements(imwht_HSIM[*,0])-1
  YSIZE=n_elements(imwht_HSIM[0,*])-1
  XMINS=(X_IM_HSIM[i]-2*A_IM_HSIM[i])
  XMAXS=(X_IM_HSIM[i]+2*A_IM_HSIM[i])
  YMINS=(Y_IM_HSIM[i]-2*A_IM_HSIM[i])
  YMAXS=(Y_IM_HSIM[i]+2*A_IM_HSIM[i])
  IF XMINS le 0 then XMINS=0
  IF XMAXS ge XSIZE then XMAXS=XSIZE
  IF YMINS le 0 then YMINS=0
  IF YMAXS ge YSIZE then YMAXS=YSIZE

  stamp=imwht_HSIM[XMINS:XMAXS,YMINS:YMAXS]
  ; stamp=imwht_HSIM[(X_IM_HSIM[i]-2*A_IM_HSIM[i]):(X_IM_HSIM[i]+2*A_IM_HSIM[i]),(Y_IM_HSIM[i]-2*A_IM_HSIM[i]):(Y_IM_HSIM[i]+2*A_IM_HSIM[i])]

  RATE0_HSIM[i]=float(n_elements(where(stamp eq 0)))/float(n_elements(stamp))
  if RATE0_HSIM[i] gt limitr then FLAG2_HSIM[i] = FLAG2_HSIM[i] + 90
  i=i+1
 endwhile
ENDIF

; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; J IMAGE PRESENT xxxxxxxxxx
IF ISTHEREJ gt 0 then begin
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 ; J SIM
 i=0L
 while i lt n_elements(ID_JSIM) do begin
  ; In old pipeline, flag was 99 and already flagged were skipped. 
  ; Here instead we add the two flag (9+90=99 or 0+90=90)
  ; if FLAGS_JSIM[i] lt 4 and FLAG2_JSIM[i] eq 0 then begin
  XSIZE=n_elements(imwht_JSIM[*,0])-1
  YSIZE=n_elements(imwht_JSIM[0,*])-1
  XMINS=(X_IM_JSIM[i]-2*A_IM_JSIM[i])
  XMAXS=(X_IM_JSIM[i]+2*A_IM_JSIM[i])
  YMINS=(Y_IM_JSIM[i]-2*A_IM_JSIM[i])
  YMAXS=(Y_IM_JSIM[i]+2*A_IM_JSIM[i])
  IF XMINS le 0 then XMINS=0
  IF XMAXS ge XSIZE then XMAXS=XSIZE
  IF YMINS le 0 then YMINS=0
  IF YMAXS ge YSIZE then YMAXS=YSIZE

  stamp=imwht_JSIM[XMINS:XMAXS,YMINS:YMAXS]
  ; stamp=imwht_JSIM[(X_IM_JSIM[i]-2*A_IM_JSIM[i]):(X_IM_JSIM[i]+2*A_IM_JSIM[i]),(Y_IM_JSIM[i]-2*A_IM_JSIM[i]):(Y_IM_JSIM[i]+2*A_IM_JSIM[i])]
  RATE0_JSIM[i]=float(n_elements(where(stamp eq 0)))/float(n_elements(stamp))
  if RATE0_JSIM[i] gt limitr then FLAG2_JSIM[i] = FLAG2_JSIM[i] + 90
  i=i+1
 endwhile
ENDIF

; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; J+H IMAGES PRESENT xxxxxxxxxx
IF ISTHEREJ gt 0 then begin
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 ; J and H DIM
 i=0L
 while i lt n_elements(ID_JDIM) do begin
  ; In old pipeline, flag was 99 and already flagged were skipped. 
  ; Here instead we add the two flag (9+90=99 or 0+90=90)
  ; if FLAGS_JDIM[i] lt 4 and FLAG2_JDIM[i] eq 0 then begin
  XSIZE=n_elements(imwht_DIM[*,0])-1
  YSIZE=n_elements(imwht_DIM[0,*])-1
  XMINS=(X_IM_JDIM[i]-2*A_IM_JDIM[i])
  XMAXS=(X_IM_JDIM[i]+2*A_IM_JDIM[i])
  YMINS=(Y_IM_JDIM[i]-2*A_IM_JDIM[i])
  YMAXS=(Y_IM_JDIM[i]+2*A_IM_JDIM[i])
  IF XMINS le 0 then XMINS=0
  IF XMAXS ge XSIZE then XMAXS=XSIZE
  IF YMINS le 0 then YMINS=0
  IF YMAXS ge YSIZE then YMAXS=YSIZE

  stamp=imwht_DIM[XMINS:XMAXS,YMINS:YMAXS]
  ; stamp=imwht_DIM[(X_IM_JDIM[i]-2*A_IM_JDIM[i]):(X_IM_JDIM[i]+2*A_IM_JDIM[i]),(Y_IM_JDIM[i]-2*A_IM_JDIM[i]):(Y_IM_JDIM[i]+2*A_IM_JDIM[i])]
  RATE0_DIM[i]=float(n_elements(where(stamp eq 0)))/float(n_elements(stamp))
  if RATE0_DIM[i] gt limitr then FLAG2_DIM[i] = FLAG2_DIM[i] + 90
  i=i+1
 endwhile
ENDIF
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


;--------------------------------------------------------------
; FLAG OBJECTS WITH TOO ELONGATED GEOMETRIES
;--------------------------------------------------------------
; flag objects as 900 + OLD FLAG (examples: 900+0+0 or 900+90+9 or 900+0+9)

ECC_LIM=0.98

; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; H IMAGE PRESENT xxxxxxxxxx
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IF ISTHEREH gt 0 then begin
 ECC_HSIM=sqrt(1-((B_IM_HSIM^2)/(A_IM_HSIM^2)))
 TOO_ECC_HSIM=where(ECC_HSIM gt ECC_LIM)
 if TOO_ECC_HSIM[0] ne -1 then FLAG2_HSIM[TOO_ECC_HSIM] = FLAG2_HSIM[TOO_ECC_HSIM]+900
ENDIF

IF ISTHEREJ gt 0 then begin
 ECC_JSIM=sqrt(1-((B_IM_JSIM^2)/(A_IM_JSIM^2)))
 TOO_ECC_JSIM=where(ECC_JSIM gt ECC_LIM)
 if TOO_ECC_JSIM[0] ne -1 then FLAG2_JSIM[TOO_ECC_JSIM] = FLAG2_JSIM[TOO_ECC_JSIM]+900
ENDIF

IF ISTHEREJ gt 0 and ISTHEREH gt 0 then begin
 ECC_DIM=sqrt(1-((B_IM_JDIM^2)/(A_IM_JDIM^2)))
 TOO_ECC_DIM=where(ECC_DIM gt ECC_LIM)
 ;------------------------------------
 if TOO_ECC_DIM[0] ne -1 then FLAG2_DIM[TOO_ECC_DIM] = FLAG2_DIM[TOO_ECC_DIM]+900
ENDIF


 ;--------------------------------------------------------------
 ; GENERATES THE UPDATED CATALOGS
 ;--------------------------------------------------------------

 ; H-SIM
 ; DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
 ; H SINGLE IMAGE MODE CATALOG
 ; DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
IF ISTHEREH gt 0 then begin
 OPENR,lun,cat_H_SIMn,/get_lun
 header=strarr(17)
 readf,lun,header
 header(15)="#  16 BADPIXLE_RATIO         N_badpixel/N_total   "
 header(16)="#  17 FLAG2                  WISPs flags: 9 = edge, 90 = badpix, 900 = high eccentricity   "
 CLOSE,lun
 FREE_LUN,lun

 openw,u1,pathse+'F1'+JHN+'0W_SIM_full.cat',/get_lun
 openw,u2,pathse+'F1'+JHN+'0W_SIM_clean.cat',/get_lun
 printf,u1,header
 printf,u2,header[0:14]

 for i=0, n_elements(ID_HSIM)-1 do begin
  printf,u1,ID_HSIM[i],X_IM_HSIM[i],Y_IM_HSIM[i],A_IM_HSIM[i],B_IM_HSIM[i],THETA_IM_HSIM[i],X_WO_HSIM[i],Y_WO_HSIM[i],A_WO_HSIM[i],B_WO_HSIM[i],THETA_WO_HSIM[i],MAG_HSIM[i], MAGERR_HSIM[i],CLASS_STAR_HSIM[i],FLAGS_HSIM[i],RATE0_HSIM[i],FLAG2_HSIM[i],FORMAT='(I5, 4(F10.3), F10.1, 4(F20.10), F10.1, 2F10.4, F10.2, I10, F10.4, I10)'
  if FLAG2_HSIM[i] eq 0 then begin
   printf,u2,ID_HSIM[i],X_IM_HSIM[i],Y_IM_HSIM[i],A_IM_HSIM[i],B_IM_HSIM[i],THETA_IM_HSIM[i],X_WO_HSIM[i],Y_WO_HSIM[i],A_WO_HSIM[i],B_WO_HSIM[i],THETA_WO_HSIM[i],MAG_HSIM[i], MAGERR_HSIM[i],CLASS_STAR_HSIM[i],FLAGS_HSIM[i],FORMAT='(I5, 4(F10.3), F10.1, 4(F20.10), F10.1, 2F10.4, F10.2, I10)'
  endif
 endfor

 CLOSE,u1,u2
 FREE_LUN,u1,u2
ENDIF


; J (if present)
; DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
; J- SINGLE IMAGE MODE CATALOG
; DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
IF ISTHEREJ gt 0 then begin
 OPENR,lun,cat_J_SIMn,/get_lun
 header=strarr(17)
 readf,lun,header
 header(15)="#  16 BADPIXLE_RATIO         N_badpixel/N_total   "
 header(16)="#  17 FLAG2                  WISPs flags: 9 = edge, 90 = badpix, 900 = high eccentricity   "
 CLOSE,lun
 FREE_LUN,lun

 openw,u1,pathse+'F110W_SIM_full.cat',/get_lun
 openw,u2,pathse+'F110W_SIM_clean.cat',/get_lun
 printf,u1,header
 printf,u2,header[0:14]

 for i=0, n_elements(ID_JSIM)-1 do begin
  printf,u1,ID_JSIM[i],X_IM_JSIM[i],Y_IM_JSIM[i],A_IM_JSIM[i],B_IM_JSIM[i],THETA_IM_JSIM[i],X_WO_JSIM[i],Y_WO_JSIM[i],A_WO_JSIM[i],B_WO_JSIM[i],THETA_WO_JSIM[i],MAG_JSIM[i], MAGERR_JSIM[i],CLASS_STAR_JSIM[i],FLAGS_JSIM[i],RATE0_JSIM[i],FLAG2_JSIM[i],FORMAT='(I5, 4(F10.3), F10.1, 4(F20.10), F10.1, 2F10.4, F10.2, I10, F10.4, I10)'
  if FLAG2_JSIM[i] eq 0 then begin
  printf,u2,ID_JSIM[i],X_IM_JSIM[i],Y_IM_JSIM[i],A_IM_JSIM[i],B_IM_JSIM[i],THETA_IM_JSIM[i],X_WO_JSIM[i],Y_WO_JSIM[i],A_WO_JSIM[i],B_WO_JSIM[i],THETA_WO_JSIM[i],MAG_JSIM[i], MAGERR_JSIM[i],CLASS_STAR_JSIM[i],FLAGS_JSIM[i],FORMAT='(I5, 4(F10.3), F10.1, 4(F20.10), F10.1, 2F10.4, F10.2, I10)'
  endif
 endfor

 CLOSE,u1,u2
 FREE_LUN,u1,u2
ENDIF



;Dual image mode
IF ISTHEREJ gt 0 and ISTHEREH gt 0 then begin

; H (if J image is present)
; DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
; H-DUAL IMAGE MODE CATALOG
; DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
OPENR,lun,cat_H_DIMn,/get_lun
header=strarr(17)
readf,lun,header
header(15)="#  16 BADPIXLE_RATIO         N_badpixel/N_total   "
header(16)="#  17 FLAG2                  WISPs flags: 9 = edge, 90 = badpix, 900 = high eccentricity   "
CLOSE,lun
FREE_LUN,lun
;
openw,u1,pathse+'F1'+JHN+'0W_DIM_full.cat',/get_lun
openw,u2,pathse+'F1'+JHN+'0W_DIM_clean.cat',/get_lun
printf,u1,header
printf,u2,header[0:14]

for i=0, n_elements(ID_HDIM)-1 do begin
   printf,u1,ID_HDIM[i],X_IM_HDIM[i],Y_IM_HDIM[i],A_IM_HDIM[i],B_IM_HDIM[i],THETA_IM_HDIM[i],X_WO_HDIM[i],Y_WO_HDIM[i],A_WO_HDIM[i],B_WO_HDIM[i],THETA_WO_HDIM[i],MAG_HDIM[i], MAGERR_HDIM[i],CLASS_STAR_HDIM[i],FLAGS_HDIM[i],RATE0_DIM[i],FLAG2_DIM[i],FORMAT='(I5, 4(F10.3), F10.1, 4(F20.10), F10.1, 2F10.4, F10.2, I10, F10.4, I10)'
   if FLAG2_DIM[i] eq 0 then begin
      printf,u2,ID_HDIM[i],X_IM_HDIM[i],Y_IM_HDIM[i],A_IM_HDIM[i],B_IM_HDIM[i],THETA_IM_HDIM[i],X_WO_HDIM[i],Y_WO_HDIM[i],A_WO_HDIM[i],B_WO_HDIM[i],THETA_WO_HDIM[i],MAG_HDIM[i], MAGERR_HDIM[i],CLASS_STAR_HDIM[i],FLAGS_HDIM[i],FORMAT='(I5, 4(F10.3), F10.1, 4(F20.10), F10.1, 2F10.4, F10.2, I10)'
   endif
endfor

CLOSE,u1,u2
FREE_LUN,u1,u2


; J (if present)
; DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
; J- DUAL IMAGE MODE CATALOG
; DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
OPENR,lun,cat_J_DIMn,/get_lun
header=strarr(17)
readf,lun,header
header(15)="#  16 BADPIXLE_RATIO         N_badpixel/N_total   "
header(16)="#  17 FLAG2                  WISPs flags: 9 = edge, 90 = badpix, 900 = high eccentricity   "
CLOSE,lun
FREE_LUN,lun
;
openw,u1,pathse+'F110W_DIM_full.cat',/get_lun
openw,u2,pathse+'F110W_DIM_clean.cat',/get_lun
printf,u1,header
printf,u2,header[0:14]

for i=0, n_elements(ID_JDIM)-1 do begin
   printf,u1,ID_JDIM[i],X_IM_JDIM[i],Y_IM_JDIM[i],A_IM_JDIM[i],B_IM_JDIM[i],THETA_IM_JDIM[i],X_WO_JDIM[i],Y_WO_JDIM[i],A_WO_JDIM[i],B_WO_JDIM[i],THETA_WO_JDIM[i],MAG_JDIM[i], MAGERR_JDIM[i],CLASS_STAR_JDIM[i],FLAGS_JDIM[i],RATE0_DIM[i],FLAG2_DIM[i],FORMAT='(I5, 4(F10.3), F10.1, 4(F20.10), F10.1, 2F10.4, F10.2, I10, F10.4, I10)'
   if FLAG2_DIM[i] eq 0 then begin
      printf,u2,ID_JDIM[i],X_IM_JDIM[i],Y_IM_JDIM[i],A_IM_JDIM[i],B_IM_JDIM[i],THETA_IM_JDIM[i],X_WO_JDIM[i],Y_WO_JDIM[i],A_WO_JDIM[i],B_WO_JDIM[i],THETA_WO_JDIM[i],MAG_JDIM[i], MAGERR_JDIM[i],CLASS_STAR_JDIM[i],FLAGS_JDIM[i],FORMAT='(I5, 4(F10.3), F10.1, 4(F20.10), F10.1, 2F10.4, F10.2, I10)'
   endif
endfor

CLOSE,u1,u2
FREE_LUN,u1,u2

ENDIF ; IF ISTHEREJ gt 0 

end
