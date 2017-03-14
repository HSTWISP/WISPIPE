
FUNCTION chname, In_name, g102=g102, g141=g141, oned=oned,out_name

; This function Adds digits to create uniform filenames
; with 5 digit numberings
;------------------------------------------------------
; In_name --> example: "Par439_BEAM_7A"
; /g102 --> set this keyword if g102 is present
; /g141 --> set this keyword if g141 is present
; Out_name --> output string of this function

N_DIGITS=5 ; Select here the number of digits

out_name='ciao'
part1="hlsp_wisp_hst_wfc3_par"
IDX1=strpos(In_name,'_BEAM_')
NL1=IDX1-3 ;(3 letters for "Par) ; number of letters for the Par#
part2=strmid(In_name,3,NL1)+'-'

IDX2=IDX1+6 ;(6 letters for "_BEAM_)
NL2=strlen(In_name)-1-IDX2 ; number of letters for the beam
BEAM=strmid(In_name,IDX2,NL2)
N_ZEROs=N_DIGITS-strlen(BEAM)

part3=BEAM
NZ=0L
while NZ lt N_ZEROs do begin
part3='0'+part3
NZ=NZ+1
endwhile
part3=part3+'a_'

part4='g'
if keyword_set(g102) eq 1 and keyword_set(g141) ne 1 then part4=part4+'102'
if keyword_set(g141) eq 1 and keyword_set(g102) ne 1  then part4=part4+'141'
if keyword_set(g102) eq 1 and keyword_set(g141) eq 1 then part4=part4+'102-g141'

if keyword_set(oned) then part5="_v6.2_spec1d"
if not keyword_set(oned) then part5="_v6.2_stamp2d"

out_name=part1+part2+part3+part4+part5

end



;===============================
;to generate MAST delivery for WISPs
;including drzzled fits files, catalogs, 1dspectra, and 2dspectra
;input:  folder named ParXXX
;output: folder named parXXX for MAST delivery
;updated:2017.3.6 for DR2 on V6.2
; Requirement: predefined $WISPMASTDR2
;updated:2017.3.10 for DR2 on V6.2
; No keywords are required anymore.
 
;pro mastdr2_IB,g141=g141,f140=f140, id=id
pro mastdr2_IB, id=id
; id=[xyz,nlm,...] --> List of id's to process. Example: id=[12,432,23]
  
;if n_elements(id) lt 1 then begin
if not keyword_set(id) then begin
 print, 'Please provide a list of IDs to process, e.g. mastdr2, id=[333,334,335]'
endif
  
;mastpath = '$WISPDATA/../MAST/MAST-DR2/'
mastpath = '$WISPMASTDR2/'

for j = 0,n_elements(id)-1 do begin
   i = strtrim(string(id[j]),2)

print, "working on Par"+i
; =============================================
; AUTOMATICALLY DETECT THE KIND OF OUTPUT DATA 
; =============================================
  
path='$WISPDATA/aXe/'
path2=path+'Par'+i+'/'
path3=path2+'DATA/DIRECT_GRISM/'
path4=path2+'DATA/DIRECT/'
path5=path2+'DATA/GRISM/'
path6=path2+'DATA/UVIS/'


J_OBS='NO'    ; NO or YES observations in the F110 band
H_OBS='NO'    ; NO, F140 or F160
G102_OBS='NO' ; NO or YES observations in the G102 grism
G141_OBS='NO' ; NO or YES observations in the G141 grism
UV1_OBS='NO'  ; NO or YES observations in UVIS1
UV2_OBS='NO'  ; NO or YES observations in UVIS2
UV_OBS='NO'   ; NO or YES observations in one of the UVIS
f110_list='none'
f140_list='none'
f160_list='none'
g102_list='none'
g141_list='none'
UV1_list ='none'
UV2_list ='none'
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
; UVIS
TEST_UV1=file_test(path6+'UVIS1.list',/zero_length) ; 1 if exists but no content
TEST_UV1B=file_test(path6+'UVIS1.list')             ; 1 if exists
TEST_UV2=file_test(path6+'UVIS2.list',/zero_length) ; 1 if exists but no content
TEST_UV2B=file_test(path6+'UVIS2.list')             ; 1 if exists
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
;   UVIS1    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
IF TEST_UV1 eq 0 and TEST_UV1B eq 1 then begin
   readcol,path6+'UVIS1.list',UV1_list,format=('A')
   if strlowcase(UV1_list[0]) ne 'none' and n_elements(UV1_list[0]) gt 0 then UV1_OBS='YES'
ENDIF
;   UVIS2    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
IF TEST_UV2 eq 0 and TEST_UV2B eq 1 then begin
   readcol,path6+'UVIS2.list',UV2_list,format=('A')
   if strlowcase(UV2_list[0]) ne 'none' and n_elements(UV2_list[0]) gt 0 then UV2_OBS='YES'
ENDIF
; ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

print, 'HHHHHHHHHHHHHHHHH'
;print, TEST_J,TEST_JB
;print, TEST_H1,TEST_H1B
;print, TEST_H2,TEST_H2B
;print, TEST_UV1,TEST_UV1B
;print, TEST_UV2,TEST_UV2B
;print, TEST_G102, TEST_G102B
;print, TEST_G141, TEST_G141B
print, 'J_OBS = '+J_OBS
print, 'H_OBS = '+H_OBS
print, 'UV1_OBS = '+UV1_OBS
print, 'UV2_OBS = '+UV2_OBS
print, 'G102_OBS = '+G102_OBS
print, 'G141_OBS = '+G141_OBS
print, 'UV1_OBS = '+UV1_OBS
print, 'UV2_OBS = '+UV2_OBS
print, 'HHHHHHHHHHHHHHHHH'
IF UV1_OBS eq 'YES' or UV2_OBS eq 'YES' then UV_OBS='YES'
; =============================================


spawn,'mkdir '+mastpath+'par'+strtrim(i,2)
spawn,'mkdir '+mastpath+'par'+strtrim(i,2)+'/1dspectra'
spawn,'mkdir '+mastpath+'par'+strtrim(i,2)+'/2dstamp'

;////////////////////////////////////////////////
LINE_UVIS1=8
LINE_UVIS2=9
;////////////////////////////////////////////////
uvisfilter=[0]
IF UV_OBS eq 'YES' then readcol,'$WISPDATA/axe/par'+strtrim(i,2)+'/DATA/UVIS/NOTES_ON_UVIS_FILTERS.txt',skipline=LINE_UVIS1,uvisfilter,f='x,a'

; Check
;;;;   IF UV_OBS eq 'YES' and uvisfilter[0] ne 'UVIS' then begin
;;;;      print, 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
;;;;      print, 'ERROR:'
;;;;      print,'uvis observations are present but this program is not able to recognize the filters used.'
;;;;      print,'This problem is probably due to the fact that the uvis readme file'
;;;;      print,'(/DATA/UVIS/NOTES_ON_UVIS_FILTERS.txt) has changed since this program was written.'
;;;;      print, 'In this case, reset in the program the line at which the UVIS filters are defined.'
;;;;      print, 'Search for these definitions:'
;;;;      print, 'LINE_UVIS1='+LINE_UVIS1
;;;;      print, 'LINE_UVIS2='+LINE_UVIS2
;;;;      print, 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
;;;;   ENDIF

;-------Image files------- UVIS --------
;if uvisfilter[0] ne 'UVIS' then begin
if UV_OBS eq 'YES' then begin
 IF UV1_OBS eq 'YES'THEN uvis1 = strlowcase(uvisfilter[0])
 IF UV1_OBS eq 'YES'THEN uvis2 = strlowcase(uvisfilter[1])
 if UV1_OBS eq 'YES' then begin
  spawn, 'rsync -a $WISPDATA/aXe/Par'+strtrim(i,2)+'/DATA/UVIS/IRtoUVIS/UVIS1_UVIS_drz.fits  '+mastpath+'par'+strtrim(i,2)+'/hlsp_wisp_hst_wfc3_par'+strtrim(i,2)+'-40mas_'+strtrim(uvis1,2)+'_v6.2_drz.fits'
  spawn, 'rsync -a $WISPDATA/aXe/Par'+strtrim(i,2)+'/DATA/UVIS/UVIStoIR/UVIS1_IR_drz.fits  '+mastpath+'par'+strtrim(i,2)+'/hlsp_wisp_hst_wfc3_par'+strtrim(i,2)+'-130mas_'+strtrim(uvis1,2)+'_v6.2_drz.fits'
 endif
 if UV2_OBS eq 'YES' then begin
  spawn, 'rsync -a $WISPDATA/aXe/Par'+strtrim(i,2)+'/DATA/UVIS/IRtoUVIS/UVIS2_UVIS_drz.fits  '+mastpath+'par'+strtrim(i,2)+'/hlsp_wisp_hst_wfc3_par'+strtrim(i,2)+'-40mas_'+strtrim(uvis2,2)+'_v6.2_drz.fits'
  spawn, 'rsync -a $WISPDATA/aXe/Par'+strtrim(i,2)+'/DATA/UVIS/UVIStoIR/UVIS2_IR_drz.fits  '+mastpath+'par'+strtrim(i,2)+'/hlsp_wisp_hst_wfc3_par'+strtrim(i,2)+'-130mas_'+strtrim(uvis2,2)+'_v6.2_drz.fits'
 endif
 if J_OBS eq 'YES' then begin 
  spawn, 'rsync -a $WISPDATA/aXe/Par'+strtrim(i,2)+'/DATA/UVIS/UVIStoIR/F110W_IR_drz.fits  '+mastpath+'par'+strtrim(i,2)+'/hlsp_wisp_hst_wfc3_par'+strtrim(i,2)+'-130mas_f110w_v6.2_drz.fits'
  spawn, 'rsync -a $WISPDATA/aXe/Par'+strtrim(i,2)+'/DATA/UVIS/IRtoUVIS/F110W_UVIS_drz.fits  '+mastpath+'par'+strtrim(i,2)+'/hlsp_wisp_hst_wfc3_par'+strtrim(i,2)+'-40mas_f110w_v6.2_drz.fits'
 endif
 if H_OBS eq 'F160' then begin
  spawn, 'rsync -a $WISPDATA/aXe/Par'+strtrim(i,2)+'/DATA/UVIS/UVIStoIR/F160W_IR_drz.fits  '+mastpath+'par'+strtrim(i,2)+'/hlsp_wisp_hst_wfc3_par'+strtrim(i,2)+'-130mas_f160w_v6.2_drz.fits'
  spawn, 'rsync -a $WISPDATA/aXe/Par'+strtrim(i,2)+'/DATA/UVIS/IRtoUVIS/F160W_UVIS_drz.fits  '+mastpath+'par'+strtrim(i,2)+'/hlsp_wisp_hst_wfc3_par'+strtrim(i,2)+'-40mas_f160w_v6.2_drz.fits'
 endif
 if H_OBS eq 'F140' then begin
  spawn, 'rsync -a $WISPDATA/aXe/Par'+strtrim(i,2)+'/DATA/UVIS/UVIStoIR/F140W_IR_drz.fits  '+mastpath+'par'+strtrim(i,2)+'/hlsp_wisp_hst_wfc3_par'+strtrim(i,2)+'-130mas_f140w_v6.2_drz.fits'
  spawn, 'rsync -a $WISPDATA/aXe/Par'+strtrim(i,2)+'/DATA/UVIS/IRtoUVIS/F140W_UVIS_drz.fits  '+mastpath+'par'+strtrim(i,2)+'/hlsp_wisp_hst_wfc3_par'+strtrim(i,2)+'-40mas_f140w_v6.2_drz.fits'
 endif
endif

;-------Image files------- IR --------
IF J_OBS eq 'YES' then    spawn, 'rsync -a $WISPDATA/aXe/Par'+strtrim(i,2)+'/DATA/DIRECT_GRISM/F110W_drz.fits  '+mastpath+'/par'+strtrim(i,2)+'/hlsp_wisp_hst_wfc3_par'+strtrim(i,2)+'-80mas_f110w_v6.2_drz.fits'
IF G102_OBS eq 'YES' then spawn, 'rsync -a $WISPDATA/aXe/Par'+strtrim(i,2)+'/DATA/DIRECT_GRISM/G102_drz.fits  '+mastpath+'par'+strtrim(i,2)+'/hlsp_wisp_hst_wfc3_par'+strtrim(i,2)+'-80mas_g102_v6.2_drz.fits'
IF H_OBS eq 'F140' then   spawn, 'rsync -a $WISPDATA/aXe/Par'+strtrim(i,2)+'/DATA/DIRECT_GRISM/F140W_drz.fits  '+mastpath+'par'+strtrim(i,2)+'/hlsp_wisp_hst_wfc3_par'+strtrim(i,2)+'-80mas_f140w_v6.2_drz.fits'
IF H_OBS eq 'F160' then   spawn, 'rsync -a $WISPDATA/aXe/Par'+strtrim(i,2)+'/DATA/DIRECT_GRISM/F160W_drz.fits  '+mastpath+'/par'+strtrim(i,2)+'/hlsp_wisp_hst_wfc3_par'+strtrim(i,2)+'-80mas_f160w_v6.2_drz.fits'
IF G141_OBS eq 'YES' then spawn, 'rsync -a $WISPDATA/aXe/Par'+strtrim(i,2)+'/DATA/DIRECT_GRISM/G141_drz.fits  '+mastpath+'par'+strtrim(i,2)+'/hlsp_wisp_hst_wfc3_par'+strtrim(i,2)+'-80mas_g141_v6.2_drz.fits'
   
;-------catalog files-------
IF J_OBS eq 'YES' then  spawn, 'rsync -a $WISPDATA/aXe/Par'+strtrim(i,2)+'/DATA/DIRECT_GRISM/fin_F110.cat  '+mastpath+'par'+strtrim(i,2)+'/hlsp_wisp_hst_wfc3_par'+strtrim(i,2)+'_f110w_v6.2_cat.txt'
IF H_OBS eq 'F140' then spawn, 'rsync -a $WISPDATA/aXe/Par'+strtrim(i,2)+'/DATA/DIRECT_GRISM/fin_F140.cat  '+mastpath+'/par'+strtrim(i,2)+'/hlsp_wisp_hst_wfc3_par'+strtrim(i,2)+'_f140w_v6.2_cat.txt'
IF H_OBS eq 'F160' then spawn, 'rsync -a $WISPDATA/aXe/Par'+strtrim(i,2)+'/DATA/DIRECT_GRISM/fin_F160.cat  '+mastpath+'par'+strtrim(i,2)+'/hlsp_wisp_hst_wfc3_par'+strtrim(i,2)+'_f160w_v6.2_cat.txt'


;-------1d spectra files-------
spawn, 'rsync -a $WISPDATA/aXe/Par'+strtrim(i,2)+'/Spectra/Par'+strtrim(i,2)+'*.dat '+mastpath+'par'+strtrim(i,2)+'/1dspectra'
;-------2dspectra files-------
spawn, 'rsync -a $WISPDATA/aXe/Par'+strtrim(i,2)+'/Stamps/Par'+strtrim(i,2)+'*.fits '+mastpath+'par'+strtrim(i,2)+'/2dstamp'



;------ rename 1dspectra files --------
IF G102_OBS eq 'YES' THEN BEGIN
 spawn, 'ls -1 '+mastpath+'par'+strtrim(i,2)+'/1dspectra/*G102*.dat', spec1
 for k=0, n_elements(spec1)-1 do begin
  name = strtrim(spec1[k], 2)
  fooa = strpos(name, '/', /reverse_search)
  foob = strpos(name, '_', /reverse_search)
  fooe = strpos(name, '.', /reverse_search)
  foo0 = strmid(name, 0, fooa+1)     ;path
  foo1 = strmid(name, foob+1,fooe-foob-1)   ;file name w/o grism info and '.dat'
  foo2 = strmid(name, fooa+1,fooe-fooa-1)   ;file name w/o '.dat'
  inputname = 'par'+strtrim(i,2)+'_BEAM_'+foo1
  result=chname(inputname,/oned,/g102,outputname)
  spawn, 'mv '+foo0+foo2+'.dat '+foo0+outputname+'.dat'
 endfor
ENDIF

IF G141_OBS eq 'YES' THEN BEGIN
 spawn, 'ls -1 '+mastpath+'par'+strtrim(i,2)+'/1dspectra/*G141*.dat', spec1
 for k=0, n_elements(spec1)-1 do begin
  name = strtrim(spec1[k], 2)
  fooa = strpos(name, '/', /reverse_search)
  foob = strpos(name, '_', /reverse_search)
  fooe = strpos(name, '.', /reverse_search)
  foo0 = strmid(name, 0, fooa+1)     ;path
  foo1 = strmid(name, foob+1,fooe-foob-1)   ;file name w/o grism info and '.dat'
  foo2 = strmid(name, fooa+1,fooe-fooa-1)   ;file name w/o '.dat'
  inputname = 'par'+strtrim(i,2)+'_BEAM_'+foo1
  result=chname(inputname,/oned,/g141,outputname)
  spawn, 'mv '+foo0+foo2+'.dat '+foo0+outputname+'.dat'
 endfor
ENDIF

spawn, 'ls -1 '+mastpath+'par'+strtrim(i,2)+'/1dspectra/Par'+strtrim(i,2)+'_BEAM*.dat', spec1
for k=0, n_elements(spec1)-1 do begin
   name = strtrim(spec1[k], 2)
   fooa = strpos(name, '/', /reverse_search)
   foob = strpos(name, '_', /reverse_search)
   fooe = strpos(name, '.', /reverse_search)
   foo0 = strmid(name, 0, fooa+1)     ;path
   foo1 = strmid(name, foob+1,fooe-foob-1)   ;file name w/o '.dat'
   inputname = 'par'+strtrim(i,2)+'_BEAM_'+foo1
   result=chname(inputname,/oned,/g102,/g141,outputname)
   spawn, 'mv '+foo0+inputname+'.dat '+foo0+outputname+'.dat'
endfor

;------ rename 2dstamp files --------
IF G102_OBS eq 'YES' THEN BEGIN
 spawn, 'ls -1 '+mastpath+'par'+strtrim(i,2)+'/2dstamp/*G102*.fits', spec1
 for k=0, n_elements(spec1)-1 do begin
  name = strtrim(spec1[k], 2)
  fooa = strpos(name, '/', /reverse_search)
  foob = strpos(name, '_', /reverse_search)
  fooe = strpos(name, '.', /reverse_search)
  foo0 = strmid(name, 0, fooa+1)     ;path
  foo1 = strmid(name, foob+1,fooe-foob-1) ;file name w/o grism info and '.fits'
  foo2 = strmid(name, fooa+1,fooe-fooa-1) ;file name w/o '.dat'
  inputname = 'par'+strtrim(i,2)+'_BEAM_'+foo1
  result=chname(inputname,/g102,outputname)
  spawn, 'mv '+foo0+foo2+'.fits '+foo0+outputname+'.fits'
 endfor
ENDIF

IF G141_OBS eq 'YES' THEN BEGIN
 spawn, 'ls -1 '+mastpath+'par'+strtrim(i,2)+'/2dstamp/*G141*.fits', spec1
 for k=0, n_elements(spec1)-1 do begin
  name = strtrim(spec1[k], 2)
  fooa = strpos(name, '/', /reverse_search)
  foob = strpos(name, '_', /reverse_search)
  fooe = strpos(name, '.', /reverse_search)
  foo0 = strmid(name, 0, fooa+1)     ;path
  foo1 = strmid(name, foob+1,fooe-foob-1)   ;file name w/o grism info and '.fits'
  foo2 = strmid(name, fooa+1,fooe-fooa-1)   ;file name w/o '.dat'
  inputname = 'par'+strtrim(i,2)+'_BEAM_'+foo1
  result=chname(inputname,/g141,outputname)
  spawn, 'mv '+foo0+foo2+'.fits '+foo0+outputname+'.fits'
 endfor
ENDIF


;------- below for drz file name udpate & generate readme files ----------
spawn, 'ls -1 $WISPDATA/aXe/Par'+strtrim(i,2)+'/DATA/DIRECT/*flt.fits', flt

ra_targ = dindgen(n_elements(flt))
dec_targ = dindgen(n_elements(flt))
dateobs = strarr(n_elements(flt))
timeobs = strarr(n_elements(flt))
expstart = dindgen(n_elements(flt))
expend = dindgen(n_elements(flt))
exptime = dindgen(n_elements(flt))

for idx=0, n_elements(flt)-1 do begin
 h=headfits(flt[idx])
 ;hprint,h
 ;filter=strcompress(sxpar(h,'FILTER'),/remove_all)
 ;if filter eq 'F110W' then print, flt[idx]
 ra_targ[idx]  = strcompress(sxpar(h,'RA_TARG'),/remove_all)
 dec_targ[idx]  = strcompress(sxpar(h,'DEC_TARG'),/remove_all)
 dateobs[idx] = strcompress(sxpar(h,'DATE-OBS'),/remove_all)
 timeobs[idx] =  strcompress(sxpar(h,'TIME-OBS'),/remove_all)
 expstart[idx] = strcompress(sxpar(h,'EXPSTART'),/remove_all)
 EXPEND[idx] = strcompress(sxpar(h,'EXPEND'),/remove_all)
 EXPTIME[idx] = sxpar(h,'EXPTIME')
endfor


openw, lun, mastpath+'par'+strtrim(i,2)+'/hlsp_wisp_hst_wfc3_par'+strtrim(i,2)+'_v6.2_catalog-header.txt',/get_lun
printf,lun,'#TELESCOP=         HST                / telescope used to acquire data'                 
printf,lun,'#INSTRUME=         WFC3/IR               / instrument used to acquire data'                              
printf,lun,'#RA_TARG =  '+string(min(ra_targ))+' / right ascension of target (deg) (J2000) '       
printf,lun,'#DEC_TARG=  '+string(min(dec_targ))+' / declination of target (deg) (J2000) '              
printf,lun,'#DATE-OBS=         '+min(dateobs)+'    / UT date of start of first exposure'             
printf,lun,'#TIME-OBS=         '+min(timeobs)+'      / UT start time of first exposure  '              

IF G102_OBS eq 'YES' and G141_OBS eq 'YES' then openw, lun2, mastpath+'par'+strtrim(i,2)+'/1dspectra/hlsp_wisp_hst_wfc3_par'+strtrim(i,2)+'_g102-g141_v6.2_1dspectra-header.txt',/get_lun
IF G102_OBS ne 'YES' and G141_OBS eq 'YES' then openw, lun2, mastpath+'par'+strtrim(i,2)+'/1dspectra/hlsp_wisp_hst_wfc3_par'+strtrim(i,2)+'_g141_v6.2_1dspectra-header.txt',/get_lun
IF G102_OBS eq 'YES' and G141_OBS ne 'YES' then openw, lun2, mastpath+'par'+strtrim(i,2)+'/1dspectra/hlsp_wisp_hst_wfc3_par'+strtrim(i,2)+'_g102_v6.2_1dspectra-header.txt',/get_lun
   
printf,lun2,'1D Spectra Header:'
printf,lun2,'                                                        '    
printf,lun2,'   // DATA DESCRIPTION KEYWORDS'
printf,lun2,'#TELESCOP =        HST'
printf,lun2,'#INSTRUME =        WFC3/IR'

IF UV_OBS eq 'YES' THEN printf,lun2,'#INSTRUME =        WFC3/UVIS'
printf,lun2,'#FILTER   =        MULTI'

IF G102_OBS eq 'YES' THEN printf,lun2,'#FILTER[0] =       G102'
IF G102_OBS ne 'YES' THEN printf,lun2,'#FILTER[0] =           '
IF G141_OBS eq 'YES' THEN printf,lun2,'#FILTER[1] =       G141'
IF G141_OBS ne 'YES' THEN printf,lun2,'#FILTER[1] =           '


printf,lun2,'#TARGNAME =        MULTI'
printf,lun2,'#RA_TARG  = '+string(min(ra_targ))+' / right ascension of target (deg) (J2000) '       
printf,lun2,'#DEC_TARG = '+string(min(dec_targ))+' / declination of target (deg) (J2000) '              
printf,lun2,'                                                        '    

printf,lun2,'       // DATE AND TIME KEYWORDS'
printf,lun2,'#DATE-OBS =           '+string(min(dateobs))+'    / UT date of start of first exposure'             
printf,lun2,'#TIME-OBS =           '+string(min(timeobs))+'      / UT start time of first exposure  '              
printf,lun2,'#EXPTIME  =            1 '
printf,lun2,'#EXPSTART =    '+string(min(expstart))+'    /	start time of observation, or first exposure if composite [MJD]'
printf,lun2,'#EXPEND   =    '+string(max(expend))+'    / end time of observation, or last exposure if composite [MJD]'
printf,lun2,'#EXPDEFN  =             SUM'
printf,lun2,'#EXPSUM   =    '+string(total(exptime))  +'  / total exposure time in seconds'
printf,lun2,'                                                        '    
printf,lun2,'                                                        '    

printf,lun2,'       // For Tabular Spectra: FITS BINARY/ ASCII TABLE EXTENSION KEYWORDS'
printf,lun2,'#XTENSION=   ASCIITABLE               '
printf,lun2,'#BITPIX  =   8                             '
printf,lun2,'#NAXIS   =   2                /Binary table                  '
printf,lun2,'#NAXIS1  =   1152000          /Number of bytes per row      '
printf,lun2,'#NAXIS2  =   1                /Number of rows              '
printf,lun2,'#PCOUNT  =   0                /Random parameter count     '
printf,lun2,'#GCOUNT  =   1                /Group count               '
printf,lun2,'#TFIELDS =   5                /Number of columns        '
printf,lun2,'#EXTNAME =   txt              /Extension name              '             
printf,lun2,'#EXTNO   =   1                /Extension number           '
printf,lun2,'#TFORM1  = 64000E             /Real*4 (floating point)           '
printf,lun2,'#TTYPE1  = WAVE               /Column 1: Wavelength             '
printf,lun2,'#TUNIT1  = Angstroms          /Units of column 1               '
printf,lun2,'#TFORM2  = 64000D             /Real*8 (double precision)      '
printf,lun2,'#TTYPE2  = FLUX               /Column 2: Flux Density        '
printf,lun2,'#TUNIT2  = erg/s/cm^2/A       /Units of column 2            '
printf,lun2,'#TFORM3  = 64000D             /Real*8 (double precision)   '
printf,lun2,'#TTYPE3  = ERROR              /Column 3: Photometric Error'
printf,lun2,'#TUNIT3  = erg/s/cm^2/A       /Units of column 3         '
printf,lun2,'#TFORM4  = 64000D             /Real*8 (double precision)'
printf,lun2,'#TTYPE4  = CONTAM             /Column 4: Contamination '
printf,lun2,'#TUNIT4  = erg/s/cm^2/A       /Units of column 4  '
printf,lun2,'#TFORM5  = 64000E             /Real*4 (floating point)'
printf,lun2,'#TTYPE5  = ZEROTH             /Column 5: Zeroth order, 0:no contamination, 1: contamination from zeroth order, 2: edge truncation'
printf,lun2,'#UNIT5   = unitless           /Units of column 5 '
printf,lun2,'#'
printf,lun2,'#COMMENT = Delivered to MAST from the WISP survey'
printf,lun2,'#END'



free_lun,lun
free_lun,lun2

eend:
endfor
end
