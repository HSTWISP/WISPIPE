;##############################################################
;# WISPIPE
;# Reduction Pipeline for the WISP program (Atek et al. 2010)
;# Hakim Atek 2009
;# Last edit by Sophia Dai 2014.06.17
;#
;# Purpose: modify the *_flt.fits file to *_flt_clean.fits by
;# rejecting low pixels along X & Y directions, bad pixels are
;# replaced with means of neighboring pixels. 
;# Calls: 
;# Input: 
;#       filename
;#       exten
;#       bad pixel mask: ~/WISPIPE/aXe/CONFIG/bp_mask_v5.fits
;# Output: cleaned image named dq
;
;
;###############################################################

function BADPIX,filename,exten,dq

image=MRDFITS(filename,exten,/silent)
dq_image=MRDFITS(dq,3,/silent)

;image=filename
;dq_image=dq

goto,bp

; These are pixel and sigma rejection used before the bad-pixel maps
;image=filename
;===============================================
;    LOW PIXELS REJECTION
;===============================================
mean_width=3
bw2=mean_width^2
mean_2d=( filter_image( image,SMO=3,ALL=all)*bw2 - image )/(bw2-1)


low_fact=1.2

mean_image=mean(image)
median_image=median(image)
minimum=min([mean_image,median_image])

ind_low=where(image LT (minimum/low_fact))
image(ind_low)=mean_2d(ind_low)


;===============================================
;    SIGMA FILTER 2D FRAMES
;===============================================

;image_clean=sigma_filter(image,5, N_sigma=5, /ALL,/MON,/ITERATE )
;image=image_clean


;===========================================================
;    SIGMA FILTER ALONG LINES (following dispersion axis)
;===========================================================

;========== LOOP ALONG X AXIS ============================

;************************* 
; width of box for mean and diff measurement, 
; and signal level used for stdev rejection
       width=3
       sig=6
;************************
diff=fltarr(n_elements(image[*,0]))
mean_val=fltarr(n_elements(image[*,0]))
dev=fltarr(n_elements(image[*,0]))
med=fltarr(n_elements(image[*,0]))


for i=0, n_elements(image[0,*])-1 do begin
     
    line=image[*,i]


        for j=width,n_elements(line)-width-1 do begin
        box=[line[j-width:j-1],line[j+1:j+width]]
        mean_val(j)=mean(box)       
        diff(j)=abs(line(j)-mean_val(j))
        dev(j)=stddev(box)
        med(j)=median(box) 

          ; low pixels rejection along x
          median_1d=median(line)
          IF (line(j) LT (median_1d/low_fact)) THEN BEGIN        
          line(j)=mean_val(j)
          ENDIF


        endfor
;sigma rejection along x
ind_sig=where(diff ge sig*dev)
if (ind_sig[0] NE -1) THEN BEGIN
line(ind_sig)=mean_val(ind_sig)
ENDIF


 image[*,i]=line

endfor

image_inter=image

;========== LOOP ALONG Y AXIS ============================



;*************************
       width=2
       sig=2.5
;************************

diff=fltarr(n_elements(image[0,*]))
mean_val=fltarr(n_elements(image[0,*]))
dev=fltarr(n_elements(image[0,*]))
med=fltarr(n_elements(image[0,*]))

for iter=0,1 DO BEGIN  ; begin iterations

;filtering image for y axis rejection
image=filter_image(image,SMOOTH =3, /ALL_PIXELS,/ITERATE)
;writefits,path+'test_filter.fits', image

for i=0, n_elements(image[*,0])-1 do begin
    col_inter=image_inter[i,*] 
    col=image[i,*]
  
        for j=width,n_elements(col)-width-1 do begin
        box=[col[j-width:j-1],col[j+1:j+width]]
        mean_val(j)=TOTAL(box)/N_elements(box)       
        diff(j)=abs(col(j)-mean_val(j))
        dev(j)=stddev(box)
        med(j)=median(box) 
        endfor

;sigma rejection along y
ind_sig=where(diff ge sig*dev)
if (ind_sig[0] NE -1) THEN BEGIN
col_inter(ind_sig)=med(ind_sig)
ENDIF

image_inter[i,*]=col_inter

endfor

image=image_inter

endfor ; end iterations

;MODFITS,path+list[k], image, EXTEN_NO = 1







bp:
; cleaning using the bad pixel mask
;********************************************

bp_mask=MRDFITS(getenv('WISPIPE')+'/aXe/CONFIG/bp_mask_v5.fits',0,/silent)
bp_mask_2=bp_mask


; for size 1014*1014
for i=3,1010 do begin
    for j=3,1010 do begin

       cross=mean([image[i-1,j],image[i+1,j],image[i,j+1],image[i,j-1],image[i,j]])
     ;  sky=median([reform(image[i-2,j-2:j+2]),reform(image[i+2,j-2:j+2]),reform(image[i-2:i+2,j-2]),reform(image[i-2:i+2,j+2])]) 
     ;  std=stddev([reform(image[i-2,j-2:j+2]),reform(image[i+2,j-2:j+2]),reform(image[i-2:i+2,j-2]),reform(image[i-2:i+2,j+2])])      
       sky=median([image[i-1,j-1],image[i-1,j+1],image[i+1,j-1],image[i+1,j+1]]) 
       std=stddev([image[i-1,j-1],image[i-1,j+1],image[i+1,j-1],image[i+1,j+1]])   
       side=median([image[i-3,j],image[i-2,j],image[i+2,j],image[i+3,j]]) 
       
     ;;  if ((bp_mask[i,j] eq 1 or bp_mask_2[i,j] eq 1) and ((cross- sky) gt (1.5*std)) and ((side-sky) lt (2.5*std)) ) then begin
;;          ; print,"stddev=",std  
;;           image[i,j]=sky
;;           image[i-1,j]=sky        
;;           image[i+1,j]=sky
;;           image[i,j-1]=sky
;;           image[i,j+1]=sky
;;       endif

      if ( (dq_image[i,j] gt 0 and dq_image[i,j] ne 512 ) and (abs(cross- sky) gt (1.8*std)) and (abs(side-sky) lt (1.3*std))) then begin
       image[i,j]=sky
       image[i-1,j]=sky        
       image[i+1,j]=sky
       image[i,j-1]=sky
       image[i,j+1]=sky
      endif

   endfor
endfor




for i=2,1011 do begin
     for j=2,1011 do begin
         med_box=median([reform(image[i-1,j-1:j+1]),reform(image[i+1,j-1:j+1]),reform(image[i-1:i+1,j-1]),reform(image[i-1:i+1,j+1])]) 
         std_box=stddev([reform(image[i-1,j-1:j+1]),reform(image[i+1,j-1:j+1]),reform(image[i-1:i+1,j-1]),reform(image[i-1:i+1,j+1])])
         side=median([image[i-2,j],image[i-1,j],image[i+1,j],image[i+2,j]]) 
         v1_side=mean([reform(image[i-1:i+1,j-1]),reform(image[i-1:i+1,j+1])]) 
         v2_side=mean([reform(image[i-1:i+1,j-2]),reform(image[i-1:i+1,j+2])]) 
 
         if ( (image[i,j]-med_box) gt (std_box*3.5) and ((side-med_box) lt (std_box*1.5)) and ( (v1_side-v2_side) lt (std_box*1.2) )) then begin        
         image[i,j]=med_box      
         endif   

      endfor
  endfor




;image=sigma_filter(image,2, N_sigma=7, /ALL,/MON,/ITERATE )

;endfor ; end loop on images
;==========================================================================================
fin:
;MODFITS,filename+'_clean'+'.fits', image, EXTEN_NO = 1

return, image

END



;===============================================================



