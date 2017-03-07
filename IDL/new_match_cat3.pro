pro new_match_cat3,field, path0

; IDL program for wisp pipeline
; Ivano baronchelli July 2016
;
; Creates and matches two catalogs of sources , one for the J (110),
; if present, and one for the H (140 or 160) bands.
;
; This IDL program sobstitutes the old following programs:
; - match_cat.pro
; - match match_cat_f140.pro
; - match match_cat_f160.pro
; - match_cat_g141 (from version 2)
;
; NOTE
; The following notes refer to the the case when both the H and J
; images are present. However, the case when only the H image is
; present is also treated in an automatic way. We refer to H as "F160"
; filter, but the H="F140" case is also automatically considered
; and treated by this program.
;
; 1) SEXTRACTION and intermediate catalogs
; This IDL program call SExtractror and creates 4 intermediate
; catalogs. 
; DIM <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
; The first two are J and H catalogs extracted in Dual Image
; Mode (DIM), using a combined J-H image as a reference (extraction)
; image. The flux computation, instead, is performed on the J and H
; images. The two catalogs created in this way are:
; - F110W_DIM_1.9.cat (or F110W_DIM_2.3.cat depending on the combined
;   exposure time) 
; - F160W_DIM_1.9.cat (or F160W_DIM_2.3.cat depending on the combined
;   exposure time) 
; During the DIM extraction, a segmentation map is created:
; - JH_combined_seg.fits'. This segmentation map is used in the
;   following steps to identify sources identified in single image
;   mode but not in dual image mode.
; 
; SIM <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
; Two additional catalogs are extracted in Single Image Mode (SIM) from the
; J and H images:
; - F110W_1.9.cat (or F110W_2.3.cat depending on the total expt)
; - F160W_2.0.cat
;
;
; 2) FLAGGING-CLEANING
; The 4 catalogs previously created are passed through new_seflag.pro
; This program flags the sources in the image border (FLAG2=9), objects
; with 50% null pixels in the underlying weight map (FLAG2=90), and
; objects with extreme eccentricities (FLAG2=900). For each input
; catalog (the 2 DIM + the 2 SIM catalogs), the sources are
; saved in 2 additional intermediate catalogs:
; - F110W_SIM_clean.cat (where flagged sources are eliminated)
; - F110W_SIM_full.cat (where all the sources are present)
; The "full" catalogs are not used in the following steps, while the
; "cleaned" catalogs are considered in their place. 
; 
; 3) CREATING FINAL OUTPUT CATALOGS
; Two final catalogs are created: one for the J and one for the H
; filter:
; - cat_F110.cat 
; - cat_F160.cat 
; CATALOGS STRUCTURE <<<<<<<<<<<<<<<<<
; Differently form previous versions of the wisp pipeline, the main
; part of the catalog is made by sources extracted in dual image
; mode. This bring to source positions and sizes identical in J and H.
; - The first part of the two catalog is made by sources detected in DIM
;   In the catalogs, only the measured magnitude and the associated
;   uncertainty differ in the two filters. 
; - The second part of the catalogs is made by sources detected in
;   single image mode (SIN), in the J band (only), but not in dual
;   image mode. Their IDs starts from 10000. The magnitude in the J
;   catalog is the measured one, while in the H filter it is set to
;   the H magnitude limit, with mag uncertainty =-99 and
;   flag=extraction flag+90.
; - The third part is made by sources detected in single image mode
;   (SIN), in the H band (only), but not in dual image mode. Their IDs
;   starts from 20000. The magnitude in the H catalog is the measured
;   one, while in the J catalog it is set to the J magnitude limit,
;   with magnitude uncertainty set to -99 and flag=extraction flag+90. 
; - The fourth part is made by sources detected in both J and H band,
;   but not in dual image mode. Their IDs starts from 30000.
; Sources identified in SIM in the J and/or H band are included in the
; catalog when they complay with the following criteria:
; A) The SIM source position overlap to a void area (pixel value = 0)
;   in the segmentation map obtained in the DIM extraction phase.
; B) The SIM source position DOES NOT correspond to a void area (pixel
;   value = VAL /= 0) in the segmentation map, BUT it lies outside a
;   certain distance from the correspondent DIM source. This threshold
;   distance correspondent to the sum of the A_IMAGE values
;   (semi-major axis) of the source extracted in SIM and the one
;   extrated in DIM (whose ID is indicated by the value of the underlying
;   segmentation map).
; CATALOG CHARACTERISTICS <<<<<<<<<<<<<<<<<<<<<<<<<
; - The magnitude of sources with measured
; magnitudes higher than the magnitude limit (measured through
; new_depth.pro - measuremaglim_ir) are set to minus magnitude limit.
; - Sources are ordered as a function of J magnitude (lower ID=
;   bright). If J magnitude is not present (J=mag limit J), they are
;   ordered as a function of H magnitude.

; Last edit:
; Ivano Baronchelli July 2016
;--------------------------------------------------

; Versions history
; 
; Version 3
; The magnitude limit for the undetected sources is set to the
; magnitude limit itself (and the extraction flag is increased of 90)
;
; Version 2
; The absence of J coverage is automatically considered. 
; Now this program also 
; sobstitutes match_cat_g141

; Version 1
;  version 0 is used only for the test phase and is not complete

;ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ



TS='0' ; normal run
;TS='1' ; Test phase


MAG_ZEROPOINT_110=26.83
MAG_ZEROPOINT_140=26.46
MAG_ZEROPOINT_160=25.96

;MAG_ZEROPOINT_110_str=strmid(strcompress(string(MAG_ZEROPOINT_110),/remove_all),0,5)
;MAG_ZEROPOINT_140_str=strmid(strcompress(string(MAG_ZEROPOINT_140),/remove_all),0,5)
;MAG_ZEROPOINT_160_str=strmid(strcompress(string(MAG_ZEROPOINT_160),/remove_all),0,5)


;NORMAL RUN
if TS eq '0' then begin
   path=expand_path(path0)+'/aXe/'+field+'/DATA/DIRECT_GRISM/'
   path1=expand_path(path0)+'/aXe/'+field+'/'
endif

; TEST 
if TS eq '1' then begin
   path0='../images/'
   path='../images/'
   path1='../images/'
   filed='Par364'
endif


;images 110 and (140-or-160)
img_110n=path+'F110W_sci.fits'
img_140n=path+'F140W_sci.fits'
img_160n=path+'F160W_sci.fits'

;vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
TEST_J=file_test(img_110n)  ; 1 if exists
TEST_H1=file_test(img_140n) ; 1 if exists 
TEST_H2=file_test(img_160n) ; 1 if exists 
;vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

img_JN=img_110n
ZPM_J=MAG_ZEROPOINT_110
ZPM_J_str=strmid(strcompress(string(MAG_ZEROPOINT_110),/remove_all),0,5)

JHN='none' ; neither F140 nor F160 present 

IF TEST_H1 gt 0 or TEST_H2 gt 0 THEN BEGIN

 switch_160_140=FILE_TEST(img_140n)
 IF switch_160_140 eq 0 then begin
  img_Hn=img_160n
  ZPM_H=MAG_ZEROPOINT_160
  ZPM_H_str=strmid(strcompress(string(MAG_ZEROPOINT_160),/remove_all),0,5)
  JHN='6' ; letter to make "F160"
 ENDIF
 IF switch_160_140 eq 1 then begin
  img_Hn=img_140n
  ZPM_H=MAG_ZEROPOINT_140 
  ZPM_H_str=strmid(strcompress(string(MAG_ZEROPOINT_140),/remove_all),0,5)
  JHN='4' ; letter to make "F140"
 ENDIF

ENDIF

; ----------------------------------------
; HOW MANY FILTERS ARE WE CONSIDERING?
; --> Check J image exixtence
ISTHEREH=TEST_H1 + TEST_H2
ISTHEREJ=TEST_J;FILE_TEST(img_Jn)
; ----------------------------------------



;HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
; Calculate the depth of the IR images
;HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

IF ISTHEREH gt 0 then limit_H=measuremaglim_ir(img_HN,ZPM_H, 2., 0.01, 5, 1) ; plot parameter added to the end, default = 0 so no plot is shown
IF ISTHEREJ gt 0 then limit_J=measuremaglim_ir(img_JN,ZPM_J, 2., 0.01, 5, 1) ; Same as H, if J is present


;HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
; Match objects in the two grisms
;HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

dmax=1.5d-4 ; maximum distance for matching (degrees)



;Read image exposure times and compute weigths
;*********************************************
; copy both F140 (or 160) and combined image to the SEX folder
IF ISTHEREJ gt 0 and ISTHEREH gt 0 then spawn,'cp '+path+'JH_combined_*.fits '+path1+'SEX/'
spawn,'cp '+path+'F1*W_*fits '+path1+'SEX/'

IF ISTHEREH gt 0 then begin
 ; Exposure time H (F140 or F160)
 h_H=headfits(path1+'SEX/F1'+JHN+'0W_sci.fits')
 exptime_H=sxpar(h_H,'EXPTIME')
 expt_H_str=strcompress(string(exptime_H),/remove_all)
ENDIF

IF ISTHEREJ gt 0 then begin
 ; Exposure time J (F110)
 h_J=headfits(path1+'SEX/F110W_sci.fits')
 exptime_J=sxpar(h_J,'EXPTIME')
 expt_J_str=strcompress(string(exptime_J),/remove_all)
ENDIF


IF ISTHEREJ gt 0 and ISTHEREH gt 0 then begin
; read weights used to combine the images
File_wheights=path+'combine_weights.dat'
readcol,File_wheights,WHEIGHTS,format='a'
W_J=double(WHEIGHTS[0])
W_H=double(WHEIGHTS[1])

; GET COMBINED EXPOSURE TIME - 
; Both J and H exposure times have to be considered

; combined exposure times and wheights
;expt_comb=(exptime_J*W_J + exptime_H*W_H)
; IMPORTANT: The H image is normalized to the J standard.
;            the exposure time of the H image has no meaning here.
;            The exptime of the J image has to be considered instead,
;            with a wheight that keeps into account the different S/N
;    THE COMBINED IMAGE IS LIKE A J IMAGE WITH A LONGER EXPOSURE TIME.
exptime_comb=(exptime_J*W_J + exptime_J*W_H)/W_J ; NOTE: W_J=1
expt_comb_str=strcompress(string(exptime_comb),/remove_all)

;HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
; DUAL IMAGE SEXtraction
;HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

; NOTE: 
; smooth_and_combine.pro combines the J and H images making a final
; image that is similar to a J smoothed image (with longer exposure
; time). For this reason we can keep the same reference thresholds
; used when the extraction was performed on the J image (old
; match_cat.pro versions) 

if exptime_comb gt 1041 then det_comb='1.9'
if exptime_comb le 1041 then det_comb='2.3'

; DUAL IMAGE SEXtraction for combined J-H on J image
; output catalog name will be:
; SEX/F110W_DIM_det.cat (DIM=Dual Image Mode, det =1.9 or 2.3)
; NOTE: for the GAIN use the exposure tie of each speific image (not
; the "combined" one

; Segmentation map name (for both J and H, the same)
;segmap_name=path+'SEX/DIM_'+det_comb+'_seg.fits'
segmap_name=path1+'SEX/JH_combined_seg.fits'

 spawn,'sex '+path1+"SEX/JH_combined_sci.fits,"+path1+'SEX/F110W_sci.fits -c '+path1+'SEX/config.sex -catalog_name '+path1+$
      'SEX/F110W_DIM_'+det_comb+'.cat -mag_zeropoint '+ZPM_J_str+' -WEIGHT_TYPE MAP_WEIGHT,MAP_RMS -weight_image '+$
      path1+'SEX/JH_combined_wht.fits,'+path1+'SEX/JH_combined_rms.fits -parameters_name '+path1+'SEX/config.param '+$
      '-FILTER Y -filter_name '+path1+'SEX/gauss_2.0_5x5.conv -detect_minarea 12 -detect_thresh '+det_comb+$
      ' -ANALYSIS_THRESH 2.0 -DEBLEND_NTHRESH 64 -DEBLEND_MINCONT 0.005 -GAIN '+expt_J_str+' -STARNNW_NAME '+path1+$
      'SEX/default.nnw -CHECKIMAGE_NAME '+segmap_name


; DUAL IMAGE SEXtraction for combined J-H on H image
; output catalog name will be:
; SEX/F160W_DIM_det.cat (DIM=Dual Image Mode, det =1.9 or 2.3)  OR: 
; SEX/F140W_DIM_det.cat (DIM=Dual Image Mode, det =1.9 or 2.3)
; Segmentation map overwritten (they are the identical)
 spawn,'sex '+path1+"SEX/JH_combined_sci.fits,"+path1+'SEX/F1'+JHN+'0W_sci.fits -c '+path1+'SEX/config.sex -catalog_name '+path1+$
      'SEX/F1'+JHN+'0W_DIM_'+det_comb+'.cat -mag_zeropoint '+ZPM_H_str+' -WEIGHT_TYPE MAP_WEIGHT,MAP_RMS -weight_image '+$
      path1+'SEX/JH_combined_wht.fits,'+path1+'SEX/JH_combined_rms.fits -parameters_name '+path1+'SEX/config.param '+$
      '-FILTER Y -filter_name '+path1+'SEX/gauss_2.0_5x5.conv -detect_minarea 12 -detect_thresh '+det_comb+$
      ' -ANALYSIS_THRESH 2.0 -DEBLEND_NTHRESH 64 -DEBLEND_MINCONT 0.005 -GAIN '+expt_H_str+' -STARNNW_NAME '+path1+$
      'SEX/default.nnw -CHECKIMAGE_NAME '+segmap_name
 
ENDIF ;IF ISTHEREJ gt 0 and ISTHEREH gt 0 



;HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
; SEXtraction on single image J
;HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
IF ISTHEREJ gt 0 then begin

; SEXtraction on J image (SINGLE IMAGE MODE)

;2015.02.04 updated with new parameters
if exptime_J gt 1041 then det_J='1.9'
if exptime_J le 1041 then det_J='2.3'

spawn,'sex '+path1+'SEX/F110W_sci.fits -c '+path1+'SEX/config.sex -catalog_name '+path1+$
      'SEX/F110W_'+det_J+'.cat -mag_zeropoint '+ZPM_J_str+' -WEIGHT_TYPE MAP_WEIGHT,MAP_RMS -weight_image '$
      +path1+'SEX/F110W_wht.fits,'+path1+'SEX/F110W_rms.fits -parameters_name '+path1+$
      'SEX/config.param -FILTER Y -filter_name '+path1+'SEX/gauss_2.0_5x5.conv -detect_minarea 12 -detect_thresh '+det_J+$
      ' -ANALYSIS_THRESH 2.0 -DEBLEND_NTHRESH 64 -DEBLEND_MINCONT 0.005 -GAIN '+expt_J_str+' -STARNNW_NAME '+path1+$
      'SEX/default.nnw -CHECKIMAGE_NAME '+path1+'SEX/F110W_'+det_J+'_seg.fits '

ENDIF ;IF ISTHEREJ gt 0 


;HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
; SEXtraction on single image H
;HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
IF ISTHEREH gt 0 then begin

; SEXtraction on H image (SINGLE IMAGE MODE)

;F160-140W
det_H='2.0'
spawn,'sex '+path1+'SEX/F1'+JHN+'0W_sci.fits -c '+path1+'SEX/config.sex -catalog_name '+path1+$
      'SEX/F1'+JHN+'0W_'+det_H+'.cat -mag_zeropoint '+ZPM_H_str+' -WEIGHT_TYPE MAP_WEIGHT,MAP_RMS -weight_image '$
      +path1+'SEX/F1'+JHN+'0W_wht.fits,'+path1+'SEX/F1'+JHN+'0W_rms.fits -parameters_name '+path1+$
      'SEX/config.param -FILTER Y -filter_name '+path1+'SEX/gauss_2.0_5x5.conv -detect_minarea 12 -detect_thresh '+det_H+$
      ' -ANALYSIS_THRESH 2.0 -DEBLEND_NTHRESH 64 -DEBLEND_MINCONT 0.005 -GAIN '+expt_H_str+' -STARNNW_NAME '+path1+$
      'SEX/default.nnw -CHECKIMAGE_NAME '+path1+'SEX/F1'+JHN+'0W_'+det_H+'_seg.fits '

ENDIF ;IF ISTHEREH gt 0 


; NEW SEFLAG: eliminates sources in the borders, seto to 0 sources not
; covered in one of the two filters

; Flag for dual image mode and for J an H single image mode
IF ISTHEREJ gt 0 AND ISTHEREH gt 0 then new_seflag,JHN,det_H,det_J,det_comb,field,expand_path(path0)+'/'
; FLAG for H single image mode only 
IF ISTHEREJ eq 0 then new_seflag,JHN,det_H,'0','0',field,expand_path(path0)+'/'
; FLAG for J single image mode only 
IF ISTHEREH eq 0 then new_seflag,JHN,'0',det_J,'0',field,expand_path(path0)+'/'

;HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
; READ SExtracted catalogs (cleaned)
;HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
; READ CLEANED CATALOGS (cleaned by new_seflag.pro)
;-------------------------
; SIM = Single Image Mode
; DIM = Dual Image Mode
;-------------------------


IF ISTHEREH gt 0 then BEGIN
 cat_H_SIMn=path1+'SEX/'+'F1'+JHN+'0W_SIM_clean.cat' 
 readcol,cat_H_SIMn,ID_HSIM,X_IM_HSIM,Y_IM_HSIM,A_IM_HSIM,B_IM_HSIM,THETA_IM_HSIM,X_WO_HSIM,Y_WO_HSIM,A_WO_HSIM,B_WO_HSIM,THETA_WO_HSIM,MAG_HSIM, MAGERR_HSIM ,CLASS_STAR_HSIM ,FLAGS_HSIM   ,format='i,f,f,f,f,f,f,f,f,f,f,f,f,f,i',skipline=15,/silent
ENDIF

IF ISTHEREJ gt 0 then BEGIN
 cat_J_SIMn=path1+'SEX/'+'F110W_SIM_clean.cat'
 readcol,cat_J_SIMn,ID_JSIM,X_IM_JSIM,Y_IM_JSIM,A_IM_JSIM,B_IM_JSIM,THETA_IM_JSIM,X_WO_JSIM,Y_WO_JSIM,A_WO_JSIM,B_WO_JSIM,THETA_WO_JSIM,MAG_JSIM, MAGERR_JSIM ,CLASS_STAR_JSIM ,FLAGS_JSIM   ,format='i,f,f,f,f,f,f,f,f,f,f,f,f,f,i',skipline=15,/silent
ENDIF

IF ISTHEREJ gt 0 AND ISTHEREH gt 0 then BEGIN
 cat_J_DIMn=path1+'SEX/'+'F110W_DIM_clean.cat'
 readcol,cat_J_DIMn,ID_JDIM,X_IM_JDIM,Y_IM_JDIM,A_IM_JDIM,B_IM_JDIM,THETA_IM_JDIM,X_WO_JDIM,Y_WO_JDIM,A_WO_JDIM,B_WO_JDIM,THETA_WO_JDIM,MAG_JDIM, MAGERR_JDIM ,CLASS_STAR_JDIM ,FLAGS_JDIM   ,format='i,f,f,f,f,f,f,f,f,f,f,f,f,f,i',skipline=15,/silent
 cat_H_DIMn=path1+'SEX/'+'F1'+JHN+'0W_DIM_clean.cat'
 readcol,cat_H_DIMn,ID_HDIM,X_IM_HDIM,Y_IM_HDIM,A_IM_HDIM,B_IM_HDIM,THETA_IM_HDIM,X_WO_HDIM,Y_WO_HDIM,A_WO_HDIM,B_WO_HDIM,THETA_WO_HDIM,MAG_HDIM, MAGERR_HDIM ,CLASS_STAR_HDIM ,FLAGS_HDIM   ,format='i,f,f,f,f,f,f,f,f,f,f,f,f,f,i',skipline=15,/silent

 ; CHECKTEST 1
 IF n_elements(ID_JDIM) ne n_elements(ID_HDIM) then begin
  PRINT, 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
  PRINT, ' MAIN ERROR: The dimension of the two catalogs '
  PRINT, ' extracted in dual image mode and cleaned using'
  PRINT, " new_seflag.pro DOESN'T MATCH. The reduction will"
  PRINT, ' stop here until the bug is fixed'
  PRINT, 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
  stop
 ENDIF
 ; CHECKTEST 1 - END

 ; CHECKTEST 2
 nel=0
 while nel lt n_elements(ID_JDIM) do begin
  IF X_IM_JDIM[nel] ne X_IM_HDIM[nel] then begin
   PRINT, 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
   PRINT, ' MAIN ERROR: The source positions in the two catalogs '
   PRINT, " extracted in dual image mode DOESN'T MATCH. "
   PRINT, " The reduction will stop here until the bug is fixed"
   PRINT, 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
   stop
  ENDIF
  ; CHECKTEST 1 - END
  nel=nel+1
 endwhile
ENDIF

;HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
; SET MAGNITUDE LIMIT 
;HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
; All the sources with magnitude higher than the magnitude limit are
; set equal to the magnitude limit in the specfic band considered

; The negative value indicates that the magnitude is a magnitude limit.

; IMPORTANT NOTE: These are not the only sources that must be set to
; these values. When writing the catalog, further sources are flagged
; and set to mag_lim (for each band, the ones detected in SIM in the other band)

IF ISTHEREJ gt 0 and ISTHEREH gt 0 then begin
 IDX_BAD_JDIM=where(MAG_JDIM gt limit_J)
 if IDX_BAD_JDIM[0] ne -1 then begin
  MAG_JDIM[IDX_BAD_JDIM]=limit_J
  MAGERR_JDIM[IDX_BAD_JDIM]=-99
  FLAGS_JDIM[IDX_BAD_JDIM]=FLAGS_JDIM[IDX_BAD_JDIM]+90
 endif
 IDX_BAD_HDIM=where(MAG_HDIM gt limit_H)
 if IDX_BAD_HDIM[0] ne -1 then begin
  MAG_HDIM[IDX_BAD_HDIM]=limit_H
  MAGERR_HDIM[IDX_BAD_HDIM]=-99
  FLAGS_HDIM[IDX_BAD_HDIM]=FLAGS_HDIM[IDX_BAD_HDIM]+90
 endif
ENDIF

IF ISTHEREJ gt 0 then begin
 IDX_BAD_JSIM=where(MAG_JSIM gt limit_J)
 if IDX_BAD_JSIM[0] ne -1 then begin
  MAG_JSIM[IDX_BAD_JSIM]=limit_J
  MAGERR_JSIM[IDX_BAD_JSIM]=-99
  FLAGS_JSIM[IDX_BAD_JSIM]=FLAGS_JSIM[IDX_BAD_JSIM]+90
 endif
ENDIF

IF ISTHEREH gt 0 then begin
 IDX_BAD_HSIM=where(MAG_HSIM gt limit_H)
 if IDX_BAD_HSIM[0] ne -1 then begin
  MAG_HSIM[IDX_BAD_HSIM]=limit_H
  MAGERR_HSIM[IDX_BAD_HSIM]=-99
  FLAGS_HSIM[IDX_BAD_HSIM]=FLAGS_HSIM[IDX_BAD_HSIM]+90
 endif
ENDIF
; TEST TEST TEST TEST TEST TEST TEST TEST
;rad_MID=7
;mk_regionfile,X_WO_JDIM,Y_WO_JDIM,rad_MID,file="../regions/JH_DIM_cleaned.reg",color='cyan'
; TEST TEST TEST TEST TEST TEST TEST TEST


;HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
;OPEN FILES TO WRITE CATALOGS
;HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
; Two catalogs are made:
; F110 catalog  for F110 magnitudes
; F140 (160) for F160 magnitudes
; Since these catalogs were extracted in Dual Image Mode, they are
; already matched and no-deblending problems are present.

; JDIM JSIM HDIM ----------------------------------------

; 1 -- Copy the matched Dual image mode catalog in the new catalog
IF ISTHEREJ gt 0 and ISTHEREH gt 0 then begin
 ; Read header in J DIM
 OPENR,lun,cat_J_DIMn,/get_lun
 header_JDIM=strarr(15)
 readf,lun,header_JDIM
 CLOSE,lun
 FREE_LUN,lun
 ; Read header in H DIM
 OPENR,lun,cat_H_DIMn,/get_lun
 header_HDIM=strarr(15)
 readf,lun,header_HDIM
 CLOSE,lun
 FREE_LUN,lun

 ; MODIFY HEADERS
 ;*******************************
 header_JDIM[11]="#  12 MAG_F1153W             Kron-like elliptical aperture magnitude-AUTO               [mag]"
 IF switch_160_140 eq 1 then header_HDIM[11]="#  12 MAG_F1392W             Kron-like elliptical aperture magnitude-AUTO               [mag]"
 IF switch_160_140 eq 0 then header_HDIM[11]="#  12 MAG_F1537W             Kron-like elliptical aperture magnitude-AUTO               [mag]"

 openw,u1,path1+'SEX/cat_F110.cat',/get_lun
 openw,u2,path1+'SEX/cat_F1'+JHN+'0.cat',/get_lun

 ; WRITE MODIFIED HEADERS
 printf, u1,header_JDIM
 printf, u2,header_HDIM
 ;SORT DIM CATALOGS BY MAGNITUDE (J)
 ;******************************
 ;  <>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>
 ; OLD sorting method 
 ; Problem: when J=Jlim, there is no sort
 ;;; SORT_IDXJ=sort(abs(MAG_JDIM))   ; JDIM and HDIM ARE ALREADY MATCHED
 ; The "abs", when sorting is needed since limit magnitude are negative
 ; and they would end up at the top of the list
 ;  <>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>
 ; NEW sorting method
 ; When MAG J is =mag lim, sort as a function of H magnitude
 ; Since limit magnitude are negative, when sorting, these
 ; would end up at the top of the list --> use "abs"
 MAGSORT=abs(MAG_JDIM)
 IDX_EQJ=where(float(abs(MAG_JDIM)) eq float(limit_J)) ; NOTE: use same float type (not double)
 if IDX_EQJ[0] ne -1 then MAGSORT[IDX_EQJ]=MAGSORT[IDX_EQJ]+abs(MAG_HDIM[IDX_EQJ])
 SORT_IDXJ=sort(MAGSORT)
 ;  <>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>


 DD=0L
 WHILE DD lt n_elements(ID_JDIM) do begin
  ; JDIM and HDIM ARE ALREADY MATCHED
  ; WRITING IN J CATALOG
  printf, u1, DD+1,X_IM_JDIM[SORT_IDXJ[DD]],Y_IM_JDIM[SORT_IDXJ[DD]],A_IM_JDIM[SORT_IDXJ[DD]],B_IM_JDIM[SORT_IDXJ[DD]],THETA_IM_JDIM[SORT_IDXJ[DD]],X_WO_JDIM[SORT_IDXJ[DD]],Y_WO_JDIM[SORT_IDXJ[DD]],A_WO_JDIM[SORT_IDXJ[DD]],B_WO_JDIM[SORT_IDXJ[DD]],THETA_WO_JDIM[SORT_IDXJ[DD]],MAG_JDIM[SORT_IDXJ[DD]], MAGERR_JDIM[SORT_IDXJ[DD]],CLASS_STAR_JDIM[SORT_IDXJ[DD]],FLAGS_JDIM[SORT_IDXJ[DD]],format='((I5, 4(F10.3), F10.1, 4(F20.10), F10.1, 2F10.4, F10.2, I10))'
  ; WRITING IN H CATALOG
  printf, u2, DD+1,X_IM_HDIM[SORT_IDXJ[DD]],Y_IM_HDIM[SORT_IDXJ[DD]],A_IM_HDIM[SORT_IDXJ[DD]],B_IM_HDIM[SORT_IDXJ[DD]],THETA_IM_HDIM[SORT_IDXJ[DD]],X_WO_HDIM[SORT_IDXJ[DD]],Y_WO_HDIM[SORT_IDXJ[DD]],A_WO_HDIM[SORT_IDXJ[DD]],B_WO_HDIM[SORT_IDXJ[DD]],THETA_WO_HDIM[SORT_IDXJ[DD]],MAG_HDIM[SORT_IDXJ[DD]], MAGERR_HDIM[SORT_IDXJ[DD]],CLASS_STAR_HDIM[SORT_IDXJ[DD]],FLAGS_HDIM[SORT_IDXJ[DD]],format='((I5, 4(F10.3), F10.1, 4(F20.10), F10.1, 2F10.4, F10.2, I10))'
  DD=DD+1
 ENDWHILE


 ; OPEN SEGMENTATION MAP
 SEGMAP=mrdfits(segmap_name,0,HDSEG)

 ;----------------------------------------------------------------
 ; FIND J AND H SIM SOURCES WITHOUT COUNTERPART IN DIM 
 ;----------------------------------------------------------------
 ; --> IDs: 10000+ identification in J SIM only
 ; --> IDs: 20000+ identification in H SIM only
 ; --> IDs: 30000+ identification in both J and H SIM 

 ; Match position and check in the segmentation map if they are
 ; superimposed to a detected source.
 ; - If they are not (underlying segmentation = 0), then they are
 ; isolated  and they have to be included in the catalog;
 ; - If they are (underlying segmentation != 0), check if they are
 ; inside a 2 A_IMAGE distance from that source, otherwise
 ; consider them as isoleted (and include them in the catalog)
 ; - Then check for possible J-H SIM counterparts among the isolated
 ;   sources. In this case remove both and put a new 30000+ J-H source

 ; FIRST:
 ; All the sources with underlying segmentation map =0 are
 ; automatically considered, The ones with non null underlying seg map
 ; are included only if they are distant more than 2
 ; A_IMAGE from the source indicated in the seg. map a.
 
 ; ----------------- ADDITIONAL J SIM --------------------------
 
 SEGVECT_J=intarr(n_elements(ID_JSIM))
 SEGVECT_J[*]=0
 
 NJ=0L
 while NJ lt n_elements(ID_JSIM) do begin
 
 ; All images, after tweakreged, have same x,y sizes as the
 ; segmentation map, we can use their position in the image and in the
 ; segmentation map indifferently
 xp=round(X_IM_JSIM[NJ])-1 ; (because IDL starts from 0)
 yp=round(Y_IM_JSIM[NJ])-1 ; (because IDL starts from 0)
 SEGSTAMP=SEGMAP[xp-2:xp+2,yp-2:yp+2]
 NONULL=where(SEGSTAMP gt 0)

 IF NONULL[0] eq -1 then begin
  ; underlying seg map is null (and the source has to be included)
  SEGVECT_J[NJ]=1
 ENDIF ELSE BEGIN
  ; underlying seg map is not null
  SEGVAL=SEGSTAMP[uniq(SEGSTAMP)]
  ; Which sources are underlying?
  SEGVAL=SEGVAL[where(SEGVAL gt 0)]
  CC=0L
  while CC lt n_elements(SEGVAL) do begin
   ; Check if the source is outside two A_image from the source position
   SIDX=where(ID_JDIM eq SEGVAL[CC]); Index of the underlying source (from seg. map)
   DIST_THRESH=2*(A_IM_JDIM[SIDX]+A_IM_JSIM[NJ])
   ; JSIM-JHDIM Radial distance
   DIST_RAD=sqrt( ((X_IM_JSIM[NJ]-X_IM_JDIM[SIDX])^2) + ((Y_IM_JSIM[NJ]-Y_IM_JDIM[SIDX])^2) )
   IF DIST_RAD gt DIST_THRESH THEN BEGIN
    SEGVECT_J[NJ]=1 ; outside the search radius, but continue the search
   ENDIF ELSE BEGIN
    ; There is a DIM identified source in which the JSIN
    ; identification can be included. Exit the cycle
    SEGVECT_J[NJ]=0
    cc=n_elements(SEGVAL)
   ENDELSE 
   cc=cc+1
  endwhile
 ENDELSE

 NJ=NJ+1
 endwhile

 ; ----------------- ADDITIONAL H SIM --------------------------
 
 SEGVECT_H=intarr(n_elements(ID_HSIM))
 SEGVECT_H[*]=0

 NH=0L
 while NH lt n_elements(ID_HSIM) do begin

  ; All images, after tweakreged, have same x,y sizes as the
  ; segmentation map, we can use their position in the image and in the
  ; segmentation map indifferently
  xp=round(X_IM_HSIM[NH])-1 ; (because IDL starts from 0)
  yp=round(Y_IM_HSIM[NH])-1 ; (because IDL starts from 0)
  SEGSTAMP=SEGMAP[xp-2:xp+2,yp-2:yp+2]
  NONULL=where(SEGSTAMP gt 0)

  IF NONULL[0] eq -1 then begin
   ; underlying seg map is null (and the source has to be included)
   SEGVECT_H[NH]=1
  ENDIF ELSE BEGIN 
   ; underlying seg map is not null
   SEGVAL=SEGSTAMP[uniq(SEGSTAMP)]
   ; Which sources are underlying?
   SEGVAL=SEGVAL[where(SEGVAL gt 0)]
   CC=0L
   while CC lt n_elements(SEGVAL) do begin
    ; Check if the source is outside two A_image from the source position
    SIDX=where(ID_JDIM eq SEGVAL[CC]); Index of the underlying source (from seg. map)
    DIST_THRESH=2*(A_IM_JDIM[SIDX]+A_IM_HSIM[NH])
    ; HSIM-JHDIM Radial distance
    DIST_RAD=sqrt( ((X_IM_HSIM[NH]-X_IM_JDIM[SIDX])^2) + ((Y_IM_HSIM[NH]-Y_IM_JDIM[SIDX])^2) )
    IF DIST_RAD gt DIST_THRESH THEN BEGIN
     SEGVECT_H[NH]=1 ; outside the search radius, but continue the search
    ENDIF ELSE BEGIN
     ; There is a DIM identified source in which the JSIN
     ; identification can be included. Exit the cycle
     SEGVECT_H[NH]=0
     cc=n_elements(SEGVAL)
    ENDELSE 
    cc=cc+1
   endwhile
  ENDELSE
  NH=NH+1
 endwhile

 ; SECOND: 
 ; Check for J-H possible counterparts among the isolated ones
 ; identified before. These will have a 30000+ ID

 J_H_IDX=-1
 H_J_IDX=-1
 ISOL_J=where(SEGVECT_J gt 0)
 ISOL_H=where(SEGVECT_H gt 0)
 IF ISOL_J[0] ne -1 and ISOL_H[0] ne -1 then begin 
  SEARCH_RAD=3600.*dmax ; Search radius in arcseconds
  cccpro,X_WO_JSIM[ISOL_J],Y_WO_JSIM[ISOL_J],X_WO_HSIM[ISOL_H],Y_WO_HSIM[ISOL_H],J_H_IDX,H_J_IDX,dt=SEARCH_RAD
  PRINT, 'NOTE: Zero associations here is not an error!'
  print, ' It just means that there are 0 sources with J and H identifications without a countrpart in the J+H combined image. '
  IF J_H_IDX[0] ne -1 then begin 
   SEGVECT_J[ISOL_J[J_H_IDX]]=0 ; reset for 30000+ sources
   SEGVECT_H[ISOL_H[H_J_IDX]]=0 ; reset for 30000+ sources
  ENDIF
 ENDIF

 ; Indexes Isolated sources
 ;-------------------------------
 IDX_ISOL_J=where(SEGVECT_J gt 0)
 IDX_ISOL_H=where(SEGVECT_H gt 0)
 ;-------------------------------

 ; THIRD ADD SOURCES TO THE CATALOG

 ; 10000+ ; J ONLY SOURCES

 IF IDX_ISOL_J[0] ne -1 then begin
  ; SORT sources as a function of J magnitude
  SORT_JSIN=sort(abs(MAG_JSIM[IDX_ISOL_J]))

  JJ=0L
  while JJ lt n_elements(IDX_ISOL_J) do begin
   ; WRITING IN J CATALOG (u1)
   printf, u1, JJ+10000,X_IM_JSIM[IDX_ISOL_J[SORT_JSIN[JJ]]],Y_IM_JSIM[IDX_ISOL_J[SORT_JSIN[JJ]]],A_IM_JSIM[IDX_ISOL_J[SORT_JSIN[JJ]]],B_IM_JSIM[IDX_ISOL_J[SORT_JSIN[JJ]]],THETA_IM_JSIM[IDX_ISOL_J[SORT_JSIN[JJ]]],X_WO_JSIM[IDX_ISOL_J[SORT_JSIN[JJ]]],Y_WO_JSIM[IDX_ISOL_J[SORT_JSIN[JJ]]],A_WO_JSIM[IDX_ISOL_J[SORT_JSIN[JJ]]],B_WO_JSIM[IDX_ISOL_J[SORT_JSIN[JJ]]],THETA_WO_JSIM[IDX_ISOL_J[SORT_JSIN[JJ]]],MAG_JSIM[IDX_ISOL_J[SORT_JSIN[JJ]]], MAGERR_JSIM[IDX_ISOL_J[SORT_JSIN[JJ]]],CLASS_STAR_JSIM[IDX_ISOL_J[SORT_JSIN[JJ]]],FLAGS_JSIM[IDX_ISOL_J[SORT_JSIN[JJ]]],format='((I5, 4(F10.3), F10.1, 4(F20.10), F10.1, 2F10.4, F10.2, I10))'
   ; WRITING IN H CATALOG (u2)(in this case, mag=mag_lim_H)
   printf, u2, JJ+10000,X_IM_JSIM[IDX_ISOL_J[SORT_JSIN[JJ]]],Y_IM_JSIM[IDX_ISOL_J[SORT_JSIN[JJ]]],A_IM_JSIM[IDX_ISOL_J[SORT_JSIN[JJ]]],B_IM_JSIM[IDX_ISOL_J[SORT_JSIN[JJ]]],THETA_IM_JSIM[IDX_ISOL_J[SORT_JSIN[JJ]]],X_WO_JSIM[IDX_ISOL_J[SORT_JSIN[JJ]]],Y_WO_JSIM[IDX_ISOL_J[SORT_JSIN[JJ]]],A_WO_JSIM[IDX_ISOL_J[SORT_JSIN[JJ]]],B_WO_JSIM[IDX_ISOL_J[SORT_JSIN[JJ]]],THETA_WO_JSIM[IDX_ISOL_J[SORT_JSIN[JJ]]],limit_H,-99,CLASS_STAR_JSIM[IDX_ISOL_J[SORT_JSIN[JJ]]],FLAGS_JSIM[IDX_ISOL_J[SORT_JSIN[JJ]]]+90,format='((I5, 4(F10.3), F10.1, 4(F20.10), F10.1, 2F10.4, F10.2, I10))'
   JJ=JJ+1
  endwhile
 ENDIF

 ; 20000+ ; H ONLY SOURCES

 IF IDX_ISOL_H[0] ne -1 then begin
  ; SORT sources as a function of H magnitude
  SORT_HSIN=sort(abs(MAG_HSIM[IDX_ISOL_H]))

  HH=0L
  while HH lt n_elements(IDX_ISOL_H) do begin
   ; WRITING IN H CATALOG (u2)
   printf, u2, HH+20000,X_IM_HSIM[IDX_ISOL_H[SORT_HSIN[HH]]],Y_IM_HSIM[IDX_ISOL_H[SORT_HSIN[HH]]],A_IM_HSIM[IDX_ISOL_H[SORT_HSIN[HH]]],B_IM_HSIM[IDX_ISOL_H[SORT_HSIN[HH]]],THETA_IM_HSIM[IDX_ISOL_H[SORT_HSIN[HH]]],X_WO_HSIM[IDX_ISOL_H[SORT_HSIN[HH]]],Y_WO_HSIM[IDX_ISOL_H[SORT_HSIN[HH]]],A_WO_HSIM[IDX_ISOL_H[SORT_HSIN[HH]]],B_WO_HSIM[IDX_ISOL_H[SORT_HSIN[HH]]],THETA_WO_HSIM[IDX_ISOL_H[SORT_HSIN[HH]]],MAG_HSIM[IDX_ISOL_H[SORT_HSIN[HH]]], MAGERR_HSIM[IDX_ISOL_H[SORT_HSIN[HH]]],CLASS_STAR_HSIM[IDX_ISOL_H[SORT_HSIN[HH]]],FLAGS_HSIM[IDX_ISOL_H[SORT_HSIN[HH]]],format='((I5, 4(F10.3), F10.1, 4(F20.10), F10.1, 2F10.4, F10.2, I10))'
   ; WRITING IN J CATALOG (u1)(in this case, mag=mag_lim_J)
   printf, u1, HH+20000,X_IM_HSIM[IDX_ISOL_H[SORT_HSIN[HH]]],Y_IM_HSIM[IDX_ISOL_H[SORT_HSIN[HH]]],A_IM_HSIM[IDX_ISOL_H[SORT_HSIN[HH]]],B_IM_HSIM[IDX_ISOL_H[SORT_HSIN[HH]]],THETA_IM_HSIM[IDX_ISOL_H[SORT_HSIN[HH]]],X_WO_HSIM[IDX_ISOL_H[SORT_HSIN[HH]]],Y_WO_HSIM[IDX_ISOL_H[SORT_HSIN[HH]]],A_WO_HSIM[IDX_ISOL_H[SORT_HSIN[HH]]],B_WO_HSIM[IDX_ISOL_H[SORT_HSIN[HH]]],THETA_WO_HSIM[IDX_ISOL_H[SORT_HSIN[HH]]],limit_J,-99,CLASS_STAR_HSIM[IDX_ISOL_H[SORT_HSIN[HH]]],FLAGS_HSIM[IDX_ISOL_H[SORT_HSIN[HH]]]+90,format='((I5, 4(F10.3), F10.1, 4(F20.10), F10.1, 2F10.4, F10.2, I10))'
   HH=HH+1
  endwhile
 ENDIF

 ; 30000+ ; J and H SOURCES, not detected by DIM

 IF J_H_IDX[0] ne -1 then begin
  ; SORT sources as a function of J magnitude
  SORT_JHSIN=sort(abs(MAG_JSIM[ISOL_J[J_H_IDX]]))

  ;******************************
  ;  <>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>
  ; OLD sorting method 
  ; Problem: when J=Jlim, there is no other sort specified
  ;;;  SORT_JHSIN=sort(abs(MAG_JSIM[ISOL_J[J_H_IDX]]))
  ; The "abs", when sorting is needed since limit magnitude are negative
  ; and they would end up at the top of the list
  ;  <>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>
  ; NEW sorting method
  ; When MAG J is =mag lim, sort as a function of H magnitude
  ; Since limit magnitude are negative, when sorting, these
  ; would end up at the top of the list --> use "abs"
  MAGSORT2=abs(MAG_JSIM[ISOL_J[J_H_IDX]])
  IDX_EQJ2=where(float(abs(MAG_JSIM[ISOL_J[J_H_IDX]])) eq float(limit_J))
  if IDX_EQJ2[0] ne -1 then  MAGSORT2[IDX_EQJ2]=MAGSORT2[IDX_EQJ2]+abs(MAG_HSIM[ISOL_H[H_J_IDX[IDX_EQJ2]]])
  SORT_JHSIN=sort(MAGSORT2)
  ;  <>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>


  JH=0L
  while JH lt n_elements(J_H_IDX) do begin
   ; WRITING IN J CATALOG (u1)
   printf, u1, JH+30000,X_IM_JSIM[ISOL_J[J_H_IDX[SORT_JHSIN[JH]]]],Y_IM_JSIM[ISOL_J[J_H_IDX[SORT_JHSIN[JH]]]],A_IM_JSIM[ISOL_J[J_H_IDX[SORT_JHSIN[JH]]]],B_IM_JSIM[ISOL_J[J_H_IDX[SORT_JHSIN[JH]]]],THETA_IM_JSIM[ISOL_J[J_H_IDX[SORT_JHSIN[JH]]]],X_WO_JSIM[ISOL_J[J_H_IDX[SORT_JHSIN[JH]]]],Y_WO_JSIM[ISOL_J[J_H_IDX[SORT_JHSIN[JH]]]],A_WO_JSIM[ISOL_J[J_H_IDX[SORT_JHSIN[JH]]]],B_WO_JSIM[ISOL_J[J_H_IDX[SORT_JHSIN[JH]]]],THETA_WO_JSIM[ISOL_J[J_H_IDX[SORT_JHSIN[JH]]]],MAG_JSIM[ISOL_J[J_H_IDX[SORT_JHSIN[JH]]]], MAGERR_JSIM[ISOL_J[J_H_IDX[SORT_JHSIN[JH]]]],CLASS_STAR_JSIM[ISOL_J[J_H_IDX[SORT_JHSIN[JH]]]],FLAGS_JSIM[ISOL_J[J_H_IDX[SORT_JHSIN[JH]]]],format='((I5, 4(F10.3), F10.1, 4(F20.10), F10.1, 2F10.4, F10.2, I10))'
   ; WRITING IN H CATALOG (u2)
   printf, u2, JH+30000,X_IM_HSIM[ISOL_H[H_J_IDX[SORT_JHSIN[JH]]]],Y_IM_HSIM[ISOL_H[H_J_IDX[SORT_JHSIN[JH]]]],A_IM_HSIM[ISOL_H[H_J_IDX[SORT_JHSIN[JH]]]],B_IM_HSIM[ISOL_H[H_J_IDX[SORT_JHSIN[JH]]]],THETA_IM_HSIM[ISOL_H[H_J_IDX[SORT_JHSIN[JH]]]],X_WO_HSIM[ISOL_H[H_J_IDX[SORT_JHSIN[JH]]]],Y_WO_HSIM[ISOL_H[H_J_IDX[SORT_JHSIN[JH]]]],A_WO_HSIM[ISOL_H[H_J_IDX[SORT_JHSIN[JH]]]],B_WO_HSIM[ISOL_H[H_J_IDX[SORT_JHSIN[JH]]]],THETA_WO_HSIM[ISOL_H[H_J_IDX[SORT_JHSIN[JH]]]],MAG_HSIM[ISOL_H[H_J_IDX[SORT_JHSIN[JH]]]], MAGERR_HSIM[ISOL_H[H_J_IDX[SORT_JHSIN[JH]]]],CLASS_STAR_HSIM[ISOL_H[H_J_IDX[SORT_JHSIN[JH]]]],FLAGS_HSIM[ISOL_H[H_J_IDX[SORT_JHSIN[JH]]]],format='((I5, 4(F10.3), F10.1, 4(F20.10), F10.1, 2F10.4, F10.2, I10))'
   JH=JH+1
  endwhile
 ENDIF


 CLOSE,u1,u2
 FREE_LUN,u1,u2

ENDIF ; IF ISTHEREJ gt 0 and ISTHEREH gt 0 then begin



;XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
; HSIM IN CASE ONLY H IMAGE IS PRESENT --------------
;XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
IF ISTHEREH gt 0 AND ISTHEREJ le 0 then begin

 ; Read header in H SIM
 OPENR,lun,cat_H_SIMn,/get_lun
 header_HSIM=strarr(15)
 readf,lun,header_HSIM
 CLOSE,lun
 FREE_LUN,lun

 openw,u3,path1+'SEX/cat_F1'+JHN+'0.cat',/get_lun

 ; MODIFY HEADER
 ;*******************************
 IF switch_160_140 eq 1 then header_HSIM[11]="#  12 MAG_F1392W             Kron-like elliptical aperture magnitude-AUTO               [mag]"
 IF switch_160_140 eq 0 then header_HSIM[11]="#  12 MAG_F1537W             Kron-like elliptical aperture magnitude-AUTO               [mag]"

 ; WRITE MODIFIED HEADER
 printf, u3,header_HSIM
 
 SORT_IDXH=sort(abs(MAG_HSIM)) ; HSIM MAGNITUDE ORDER

 DD=0L
 WHILE DD lt n_elements(ID_HSIM) do begin
  ; HSIM
  printf, u3, DD+1,X_IM_HSIM[SORT_IDXH[DD]],Y_IM_HSIM[SORT_IDXH[DD]],A_IM_HSIM[SORT_IDXH[DD]],B_IM_HSIM[SORT_IDXH[DD]],THETA_IM_HSIM[SORT_IDXH[DD]],X_WO_HSIM[SORT_IDXH[DD]],Y_WO_HSIM[SORT_IDXH[DD]],A_WO_HSIM[SORT_IDXH[DD]],B_WO_HSIM[SORT_IDXH[DD]],THETA_WO_HSIM[SORT_IDXH[DD]],MAG_HSIM[SORT_IDXH[DD]], MAGERR_HSIM[SORT_IDXH[DD]],CLASS_STAR_HSIM[SORT_IDXH[DD]],FLAGS_HSIM[SORT_IDXH[DD]],format='((I5, 4(F10.3), F10.1, 4(F20.10), F10.1, 2F10.4, F10.2, I10))'
  DD=DD+1
 ENDWHILE
 CLOSE,u3
 FREE_LUN,u3

ENDIF



;XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
; JSIM IN CASE ONLY J IMAGE IS PRESENT --------------
;XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
IF ISTHEREJ gt 0 AND ISTHEREH le 0 then begin

 ; Read header in J SIM
 OPENR,lun,cat_J_SIMn,/get_lun
 header_JSIM=strarr(15)
 readf,lun,header_JSIM
 CLOSE,lun
 FREE_LUN,lun

 openw,u3,path1+'SEX/cat_F110.cat',/get_lun

 ; MODIFY HEADER
 ;*******************************
 header_JSIM[11]="#  12 MAG_F1153W             Kron-like elliptical aperture magnitude-AUTO               [mag]"
 
 ; WRITE MODIFIED HEADER
 printf, u3,header_JSIM
 
 SORT_IDXJ=sort(abs(MAG_JSIM)) ; JSIM MAGNITUDE ORDER

 DD=0L
 WHILE DD lt n_elements(ID_JSIM) do begin
  ; JSIM
  printf, u3, DD+1,X_IM_JSIM[SORT_IDXJ[DD]],Y_IM_JSIM[SORT_IDXJ[DD]],A_IM_JSIM[SORT_IDXJ[DD]],B_IM_JSIM[SORT_IDXJ[DD]],THETA_IM_JSIM[SORT_IDXJ[DD]],X_WO_JSIM[SORT_IDXJ[DD]],Y_WO_JSIM[SORT_IDXJ[DD]],A_WO_JSIM[SORT_IDXJ[DD]],B_WO_JSIM[SORT_IDXJ[DD]],THETA_WO_JSIM[SORT_IDXJ[DD]],MAG_JSIM[SORT_IDXJ[DD]], MAGERR_JSIM[SORT_IDXJ[DD]],CLASS_STAR_JSIM[SORT_IDXJ[DD]],FLAGS_JSIM[SORT_IDXJ[DD]],format='((I5, 4(F10.3), F10.1, 4(F20.10), F10.1, 2F10.4, F10.2, I10))'
  DD=DD+1
 ENDWHILE
 CLOSE,u3
 FREE_LUN,u3

ENDIF






if TS eq '0' then begin
 spawn,'cp '+path1+'SEX/cat_F1'+JHN+'0.cat '+path+'cat_F1'+JHN+'0.cat'
 IF ISTHEREJ gt 0 then spawn,'cp '+path1+'SEX/cat_F110.cat '+path+'cat_F110.cat'
endif


print, '------------------------------------------------------------------------------------'
print, 'match_cat.pro: NOTE:'
print, 'The possible presence of floating overflow and underflow error messages '
print, 'error messages right below this line, does not have to be considered as an error.'
print, 'If you read this sentece, the program completed its computations.'
print, '------------------------------------------------------------------------------------'


if TS eq '1' then stop
end
