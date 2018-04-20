;##############################################################
;# WISPIPE
;# Reduction Pipeline for the WISP program (Atek et al. 2010)
;# Hakim Atek 2009
;###############################################################
;# Purpose: to create the 1D spectra in .dat files
;# Input: field 
;#        trim
;#        cat_F110.cat (or 140 or 160)
;#
;# Output: Spectra/files.dat containing the spectra (flux, error,
;# contamination) 
;# Completely rewritten by Ivano Baronchelli November 2016
;###############################################################


pro new_wisp_extract_IB1,field,trim, path0

bin=1.

path = expand_path(path0)+'/aXe/'+field+"/"
path4=path+"DATA/DIRECT/"
path5=path+"DATA/GRISM/"
 
G102_path=path+"/G102_DRIZZLE/"
G141_path=path+"/G141_DRIZZLE/"
g102_path_out=path+'G102_OUTPUT/'
g141_path_out=path+'G141_OUTPUT/'

;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; OBS CHECK 
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
J_OBS='NO'    ; NO or YES observations in the F110 band
H_OBS='NO'    ; NO, F140 or F160
G102_OBS='NO' ; NO or YES observations in the G102 grism
G141_OBS='NO' ; NO or YES observations in the G141 grism
f110_list='none'
f140_list='none'
f160_list='none'
g102_list='none'
g141_list='none'
;------------------------------------------------
TEST_J=file_test(path4+'F110_clean.list',/zero_length)    ; 1 if exists but no content
TEST_JB=file_test(path4+'F110_clean.list')                ; 1 if exists
TEST_H1=file_test(path4+'F140_clean.list',/zero_length)   ; 1 if exists but no content
TEST_H1B=file_test(path4+'F140_clean.list')               ; 1 if exists but
TEST_H2=file_test(path4+'F160_clean.list',/zero_length)   ; 1 if exists but no content
TEST_H2B=file_test(path4+'F160_clean.list')               ; 1 if exists
TEST_G102=file_test(path5+'G102_clean.list',/zero_length) ; 1 if exists but no content
TEST_G102B=file_test(path5+'G102_clean.list')             ; 1 if exists
TEST_G141=file_test(path5+'G141_clean.list',/zero_length) ; 1 if exists but no content
TEST_G141B=file_test(path5+'G141_clean.list')             ; 1 if exists

;     J      +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
IF TEST_J eq 0 and TEST_JB eq 1 then begin
   readcol,path4+'F110_clean.list',f110_list,format=('A')
   if strlowcase(f110_list[0]) ne 'none' and n_elements(f110_list[0]) gt 0 then J_OBS='YES'
ENDIF
;     H      F140 / F160 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
IF TEST_H1 eq 0 and TEST_H1B eq 1 then begin
   readcol,path4+'F140_clean.list',f140_list,format=('A')
   if strlowcase(f140_list[0]) ne 'none' and n_elements(f140_list[0]) gt 0 then H_OBS='F140'
ENDIF
IF TEST_H2 eq 0 and TEST_H2B eq 1 then begin
   readcol,path4+'F160_clean.list',f160_list,format=('A')
   if strlowcase(f160_list[0]) ne 'none' and n_elements(f160_list[0]) gt 0 then H_OBS='F160'
ENDIF
;    G102    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
IF TEST_G102 eq 0 and TEST_G102B eq 1 then begin
   readcol,path5+'G102_clean.list',g102_list,format=('A')
   if strlowcase(g102_list[0]) ne 'none' and n_elements(g102_list[0]) gt 0 then G102_OBS='YES'
ENDIF
;    G141    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
IF TEST_G141 eq 0 and TEST_G141B eq 1 then begin
   readcol,path5+'G141_clean.list',g141_list,format=('A')
   if strlowcase(g141_list[0]) ne 'none' and n_elements(g141_list[0]) gt 0 then G141_OBS='YES'
ENDIF
; ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

print, 'wisp_extract'
print, 'HHHHHHHHHHHHHHHHH'
print, TEST_J,TEST_JB
print, TEST_H1,TEST_H1B
print, TEST_H2,TEST_H2B
print, TEST_G102, TEST_G102B
print, TEST_G141, TEST_G141B
print, 'J_OBS = '+J_OBS
print, 'H_OBS = '+H_OBS
print, 'G102_OBS = '+G102_OBS
print, 'G141_OBS = '+G141_OBS
print, 'HHHHHHHHHHHHHHHHH'

num110 = n_elements(f110_list)
num140 = n_elements(f140_list)
num160 = n_elements(f160_list)
numG102 = n_elements(G102_list)
numG141 = n_elements(G141_list)

IF strlowcase(f110_list[0]) eq 'none' then num110  =0
IF strlowcase(f140_list[0]) eq 'none' then num140  =0
IF strlowcase(f160_list[0]) eq 'none' then num160  =0
IF strlowcase(g102_list[0]) eq 'none' then numG102 =0
IF strlowcase(g141_list[0]) eq 'none' then numG141 =0
; ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


wave=0.
flux=0.
ferr=0.
contam=0.
zeroth=0.




IF H_OBS eq 'F160' then catname=path+'DATA/DIRECT_GRISM/cat_F160.cat'
IF H_OBS eq 'F140' then catname=path+'DATA/DIRECT_GRISM/cat_F140.cat'
IF J_OBS eq 'YES' then catname=path+'DATA/DIRECT_GRISM/cat_F110.cat'

readcol,catname,object,format=('I')


;=========================================================================
;                  READ EXTRACTED MULTI-SPECTRA FILES                    =
;=========================================================================
IF J_OBS EQ 'YES'  THEN  readcol,path+'DATA/DIRECT_GRISM/cat_F110.cat',magj,format=('x,x,x,x,x,x,x,x,x,x,x,f'),/silent
IF H_OBS EQ 'F140' THEN readcol,path+'DATA/DIRECT_GRISM/cat_F140.cat',magh,format=('x,x,x,x,x,x,x,x,x,x,x,f'),/silent
IF H_OBS EQ 'F160' THEN readcol,path+'DATA/DIRECT_GRISM/cat_F160.cat',magh,format=('x,x,x,x,x,x,x,x,x,x,x,f'),/silent

IF G102_OBS EQ "YES" THEN BEGIN
 G102_SPEC=G102_path+'aXeWFC3_G102_2_opt.SPC.fits'
 G102_STP=G102_path+'aXeWFC3_G102_2_opt.STP.fits'
 readcol,path+'DATA/DIRECT_GRISM/G102_0th.txt',g102_zx,g102_zy,g102_z_beam,g102_zmag,format=('f,f,i,f'),/silent 
 readcol,path+'DATA/DIRECT_GRISM/G102_1st.txt',g102_x,g102_y,g102_1_beam,format=('f,f,i'),/silent
 ; from (old) beam_pet.pro
 readcol,path+'G102_axeprep.lis',g102_axeprep_list,format=('A'),/silent
 tmp_102=strsplit(g102_axeprep_list[0],'.',/extract)
 root_102=tmp_102[0]
 g102_pet=g102_path_out+root_102+'_2.PET.fits'
 ; fits_open,g102_pet, fcb    ;read beam names in SPEC          
 ; NEXTEN_pet=fcb.nextend
 ; pet_IDs=fcb.extname
 ; fits_close,fcb
 fits_open,g102_pet, fcb_g102    ;read beam names in SPEC          
 NEXTEN_pet_g102=fcb_g102.nextend
 pet_IDs_g102=fcb_g102.extname
 fits_close,fcb_g102
ENDIF

IF G141_OBS EQ "YES" THEN BEGIN
 G141_SPEC=G141_path+'aXeWFC3_G141_2_opt.SPC.fits'
 G141_STP=G141_path+'aXeWFC3_G141_2_opt.STP.fits'
 readcol,path+'DATA/DIRECT_GRISM/G141_0th.txt',g141_zx,g141_zy,g141_z_beam,g141_zmag,format=('f,f,i,f'),/silent 
 readcol,path+'DATA/DIRECT_GRISM/G141_1st.txt',g141_x,g141_y,g141_1_beam,format=('f,f,i'),/silent
 ; from (old) beam_pet.pro
 readcol,path+'G141_axeprep.lis',g141_axeprep_list,format=('A'),/silent
 tmp_141=strsplit(g141_axeprep_list[0],'.',/extract)
 root_141=tmp_141[0]
 g141_pet=g141_path_out+root_141+'_2.PET.fits'
 ; fits_open,g141_pet, fcb    ;read beam names in SPEC          
 ; NEXTEN_pet=fcb.nextend
 ; pet_IDs=fcb.extname
 ; fits_close,fcb
 fits_open,g141_pet, fcb_g141    ;read beam names in SPEC          
 NEXTEN_pet_g141=fcb_g141.nextend
 pet_IDs_g141=fcb_g141.extname
 fits_close,fcb_g141
ENDIF





n=0L ; Object number
;n=646L ; Object number (TEST)
while n lt n_elements(object) do begin
 beam=object[n]
 Print,"Extraction BEAM ...............",beam
  
 ;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 ; The part below, inside the while cycle sobstitutes this function: 
 ;;;;;; axe_sing_IB1,field,beam,n,shift,trim,/save,expand_path(path0),/noplot
 ;;;; And it is used as if it was:
 ;;;;;; new_axe_sing_IB1,field,beam,n,expand_path(path0)
 ;;;; Because there is no need to specify a shift since it is not
 ;;;; plotting anymore (shift=0), files are lways saved (/save keyword
 ;;;; always set and /noplot is not needed anymore (the spectra are
 ;;;; not plotted anymore)
 ;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

 ;########################################################
 ; From beam_pet.pro
 ;########################################################
 beam_id=STRING(beam,Format='(I0.3)')+'A'
 pet_a=0
 pet_b=0

 ; G102 
 if G102_OBS EQ 'YES' then begin 
  beam_ind1=where(pet_IDs_g102 EQ beam_id) 
  ;extract PET fields
  ;***********************************************
  ftab_ext,g102_pet,'N,P_X,P_Y,X,Y,DIST,XS,YS,LAMBDA',npixa,p_xa,p_ya,xa,ya,dista,xsa,ysa,p_lambdaa,EXTEN_NO=beam_ind1[0]
  pet_a=[[p_xa],[p_ya],[p_lambdaa]]
 endif 

 ; G141
 if G141_OBS EQ 'YES'  then begin
  beam_ind2=where(pet_IDs_g141 EQ beam_id)
  ;extract PET fields
  ;***********************************************
  ftab_ext,g141_pet,'N,P_X,P_Y,LAMBDA',npixb,p_xb,p_yb,p_lambdab,EXTEN_NO=beam_ind2[0]
  pet_b=[[p_xb],[p_yb],[p_lambdab]]
 endif
 ;########################################################
 ; From beam_pet.pro (end)
 ;########################################################

 
;convert pet_a & pet_b to the larger scale, added by SD 2015.8
;--------------------------
pet_a = 1.60318*pet_a
pet_b = 1.60318*pet_b
;--------------------------

; !p.font=0
 match_count=0

 SPC_IND=[-1]
 SPC_IND2=[-1]



 
 ;################################  GRISM G102   ##################################

 IF G102_OBS EQ "YES" THEN  BEGIN

  ; READ STAMP FILES
  ;***********************************************************
  STP_ID='BEAM_'+STRING(beam,Format='(I0.3)')+'A' ;beams(k)
  print,"READING G102 ....... STAMP .......",STP_ID
  FITS_READ,G102_STP,G102_STAMP,header,EXTNAME=STP_ID,/NO_ABORT

  ; MATCH SPECTRUM & STAMP BEAMS
  ;**********************************************************
  fits_open,G102_SPEC, fcb    ;read beam names in SPEC          
  NEXTEN_SPEC=fcb.nextend
  SPEC_IDs=fcb.extname
  fits_close,fcb
  SPC_IND=where(SPEC_IDs EQ STP_ID) 

  ; READ SPEC FILES, array size= 255
  ;**********************************************************
  IF (SPC_IND[0] NE -1) THEN BEGIN
   match_count=match_count+1
   ftab_ext,G102_SPEC,'ID,lambda,flux,ferror,weight,contam',$
   SPEC_ID,G102_wave,G102_flux,G102_ferr,weight,G102_contam,EXTEN_NO=SPC_IND[0]
   print,"READING G102 ........ SPEC .......",SPEC_ID

   ; TRIM EDGES  array size 255-->134
   ;****************
   IF (MIN(G102_wave) LE 8250 and MAX(G102_wave) GE 11500) THEN BEGIN
    min=where(G102_wave GE 8250) & max=where(G102_wave GE 11500)
    min=min[0]+trim+3 & max=max[0]-trim
    G102_wave_trim=G102_wave[min:max]
    G102_flux_trim=G102_FLUX[min:max]
    G102_ferr_trim=G102_ferr[min:max]
    G102_contam_trim=G102_contam[min:max]
   ENDIF ELSE BEGIN
    max=n_elements(g102_wave)-1    ;added by SD for plot correction
    min=0
    G102_wave_trim=G102_wave
    G102_flux_trim=G102_FLUX
    G102_ferr_trim=G102_ferr
    G102_contam_trim=G102_contam
   ENDELSE

   ;REBIN  SPECTRA  array size 134-->133
   ;*******************
   dim_in=max-min+1
   dim_out=dim_in/bin
   rebin=indgen(dim_in-1/bin)*bin 
   G102_flux_rebin=congrid(G102_flux_trim,dim_out)
   G102_ferr_rebin=congrid(G102_ferr_trim,dim_out)
   G102_contam_rebin=congrid(G102_contam_trim,dim_out)
   G102_wave_rebin=G102_wave_trim[rebin]

   ; flagging zeroth order contamination and edge truncation
   ; problem targeted by CS of wrong contamination flag in G102 on 2015.8.7
   ; updated by SD to correct the contamination flag errors: add the
   ; g102_1x,g102_1y parameters
   ;*********************************************************

   ;G102 ---------
   g102_zo=G102_flux*0.

   for l=0,n_elements(g102_zx)-1 do begin
    ind_a=where(round(pet_a[*,0]) eq round(g102_zx(l)) and round(pet_a[*,1]) eq round(g102_zy(l)) and g102_zmag(l) lt 22.5 ); and g102_zx(l) gt (g102_1x-92) and  g102_zx(l) lt (g102_1x+92) and g102_zy(l) gt (g102_1y -4) and g102_zy(l) lt (g102_1y +4) )
    if ind_a[0] ne -1 then begin             
     g102_lam_zo=reform(pet_a[ind_a,2])
     for p=0,n_elements(ind_a)-1 do begin 
      g102_zo( where( abs(G102_wave-g102_lam_zo[p]) eq min( abs(G102_wave-g102_lam_zo[p]) )))=1
      print,".....G102..BEAM..CONTAMINATED..by..0th..order..of..object...",g102_z_beam(l)
     endfor   
    endif 
   endfor
  
   g102_trun=G102_flux*0.
   maxla=max(pet_a[*,2])    
   ind_truna=where(G102_wave gt maxla)
   if ind_truna[0] ne -1 then begin 
    g102_trun(ind_truna)=2
    g102_zo=g102_zo+g102_trun
   endif

  ENDIF   ; corresponds to 'IF (SPC_IND[0] NE -1)'
 ENDIF    ; corresponds to 'IF G102_OBS EQ "YES" '




 ;################################  GRISM G141   ##################################

 IF G141_OBS EQ "YES" THEN  BEGIN
   
  ; READ STAMP FILES
  ;***********************************************************
  STP_ID='BEAM_'+STRING(beam,Format='(I0.3)')+'A' ;beams(k)
  print,"READING G141 ....... STAMP .......",STP_ID
  FITS_READ,G141_STP,G141_STAMP,header,EXTNAME=STP_ID,/NO_ABORT

  ; MATCH SPECTRUM & STAMP BEAMS
  ;**********************************************************
  fits_open,G141_SPEC, fcb     ;read beam names in SPEC          
  NEXTEN_SPEC=fcb.nextend
  SPEC_IDs=fcb.extname
  fits_close,fcb
  SPC_IND2=where(SPEC_IDs EQ STP_ID) 

  ; READ SPEC FILES
  ;**********************************************************
  IF (SPC_IND2[0] NE -1) THEN BEGIN
   match_count=match_count+1
   ftab_ext,G141_SPEC,'ID,lambda,flux,ferror,weight,contam',$
   SPEC_ID,G141_wave,G141_flux,G141_ferr,weight,g141_contam,EXTEN_NO=SPC_IND2[0]
   print,"READING G141 ........ SPEC .......",SPEC_ID
  
   ;;;; MR edit - if SPEC_IND[0] eq -1 then it doesn't read the
   ;;;;           spectrum. But it still expects to see the spectrum
   ;;;;           below. Hence comment out the endif here, and add it
   ;;;;           below so it doesn't crash. What was logic here before?

   ; TRIM EDGES
   ;****************
   IF (MIN(G141_wave) LT 10800 and MAX(G141_wave) GT 16900) THEN BEGIN
    min=where(G141_wave GE 10800) & max=where(G141_wave GE 16900)
    min=min[0]+trim+11 & max=max[0]-trim
    G141_wave_trim=G141_wave[min:max]
    G141_flux_trim=G141_FLUX[min:max]
    G141_ferr_trim=G141_ferr[min:max]
    G141_contam_trim=G141_contam[min:max]
   ENDIF ELSE BEGIN
    max=n_elements(g141_wave)-1    ;added by SD for plot correction
    min=0                          
    G141_wave_trim=G141_wave
    G141_flux_trim=G141_FLUX
    G141_ferr_trim=G141_ferr
    G141_contam_trim=G141_contam
   ENDELSE

   ;REBIN  SPECTRA
   ;*******************
   dim_in=max-min+1
   dim_out=dim_in/bin
   rebin=indgen(dim_in-1/bin)*bin 
   G141_flux_rebin=congrid(G141_flux_trim,dim_out)
   G141_ferr_rebin=congrid(G141_ferr_trim,dim_out)
   G141_wave_rebin=G141_wave_trim[rebin]
   G141_contam_rebin=G141_contam_trim[rebin]


   ;G141 --zeroth--flag------
   ; updated by SD to correct the contamination flag errors:
   ;add the g141_1x,g141_1y parameters
   g141_zo=G141_flux*0.
  
   for m=0,n_elements(g141_zx)-1 do begin
    ind_b=where(round(pet_b[*,0]) eq round(g141_zx(m)) and round(pet_b[*,1]) eq round(g141_zy(m)) and g141_zmag(m) lt 22.5 ); and g141_zx(m) gt (g141_1x-92) and  g141_zx(m) lt (g141_1x+92) and g141_zy(m) gt (g141_1y -4) and g141_zy(m) lt (g141_1y +4))
    if ind_b[0] ne -1 then begin
     g141_lam_zo=pet_b[ind_b,2]
     for p=0,n_elements(ind_b)-1 do begin 
      g141_zo( where( abs(G141_wave-g141_lam_zo[p]) eq min( abs(G141_wave-g141_lam_zo[p]) )))=1
      print,".....G141..BEAM..CONTAMINATED..by..0th..order..of..object...",g141_z_beam(m)
     endfor   
    endif  
   endfor

   g141_trun=G141_flux*0.
   maxlb=max(pet_b[*,2])    
   ind_trunb=where(G141_wave gt maxlb)
   if ind_trunb[0] ne -1 then begin
    g141_trun(ind_trunb)=2
    g141_zo=g141_zo+g141_trun
   endif
  ENDIF   ;correspond to  'IF (SPC_IND2[0] NE -1)' 
 ENDIF   ; corresponds to 'IF G141_OBS EQ "YES" '








  ; SAVE SPECRA INSIDE file.dat files
  ;**********************************************************

 IF (SPC_IND[0] NE -1) then begin
  FORPRINT,G102_wave,G102_FLUX,G102_ferr,G102_contam,G102_zo $
  ,COMMENT='#     wave      flux            error           contam           zeroth' $
  ,TEXTOUT=path+'Spectra/'+field+'_G102_'+STP_ID+'.dat',/silent    
  WRITEFITS,path+'Stamps/'+field+'_G102_'+STP_ID+'.fits',G102_STAMP
 ENDIF

 IF (SPC_IND2[0] NE -1) THEN BEGIN
  FORPRINT,G141_wave,G141_FLUX,G141_ferr,G141_contam,G141_zo $
  ,COMMENT='#     wave      flux            error           contam           zeroth' $
  ,TEXTOUT=path+'Spectra/'+field+'_G141_'+STP_ID+'.dat',/SILENT
  WRITEFITS,path+'Stamps/'+field+'_G141_'+STP_ID+'.fits',G141_STAMP
 ENDIF
      
 cut=0
 dim1=1
 dim2=1
 IF SPC_IND[0] NE -1  THEN dim1=n_elements(G102_wave_trim)
 IF SPC_IND2[0] NE -1 THEN dim2=n_elements(G141_wave_trim)

 IF (SPC_IND[0] NE -1)  and (SPC_IND2[0] NE -1) THEN BEGIN
  wave=[G102_wave_trim[0:dim1-(cut+1)],G141_wave_trim[cut:dim2-1]]
  flux=[G102_flux_trim[0:dim1-(cut+1)],G141_flux_trim[cut:dim2-1]]
  ferr=[G102_ferr_trim[0:dim1-(cut+1)],G141_ferr_trim[cut:dim2-1]]
  contam=[G102_contam_trim[0:dim1-(cut+1)],G141_contam_trim[cut:dim2-1]] 
  zeroth=[g102_zo[0:dim1-(cut+1)],g141_zo[cut:dim2-1]]
 ENDIF 
 IF (SPC_IND[0] NE -1) and (SPC_IND2[0] eq -1) THEN BEGIN
  ; NOTE: As suggested by Hugh Dickinson,
  ; inside this IF statement, "dim2" must be 
  ; sobstituted by "dim1". When only G102 
  ; is present, the size of the G102 spectrum
  ; (i.e. dim1) must be considered.
  wave=G102_wave_trim[cut:dim1-1]
  flux=G102_flux_trim[cut:dim1-1]
  ferr=G102_ferr_trim[cut:dim1-1]
  contam=G102_contam_trim[cut:dim1-1]
  zeroth=g102_zo[cut:dim1-1]
 endif
 IF (SPC_IND[0] eq -1) and (SPC_IND2[0] NE -1) then begin
  wave=[G141_wave_trim[cut:dim2-1]]
  flux=[G141_flux_trim[cut:dim2-1]]
  ferr=[G141_ferr_trim[cut:dim2-1]]
  contam=[G141_contam_trim[cut:dim2-1]] 
  zeroth=[g141_zo[cut:dim2-1]]
 endif

 FORPRINT,wave,flux,ferr,contam,zeroth $
  ,COMMENT='#     wave      flux            error           contam           zeroth' $
  ,TEXTOUT=path+'Spectra/'+field+'_'+STP_ID+'.dat',/silent  





 n=n+1
endwhile

end
