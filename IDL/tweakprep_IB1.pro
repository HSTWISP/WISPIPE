;##############################################################
;# WISPIPE
;# Reduction Pipeline for the WISP program
;# Generated by Sophia Dai 2014
;# Purpose: 
;#       generate catalogs with SEX on the flt images
;# Input:  
;#       F110_clean.list, F160_clean.list or F140_clean.list
;# Output:
;#   direct_clean.list, direct_crclean.list,
;#   direct_clean_catfile.list, direct_crclean_catfile.list 
;#   F110_crclean.list
;#   F160_crclean.list
;#   or F140_crclean.list
;#   tweakprep.py
;# 
;# Last edit, Ivano Baronchelli 2016
;#  -  Elaborations are made on grism and direct exposures in the
;#   GRISM and DIRECT folders and not DIRECT_GRISM anymore.
;#  -  in all the "imcalc" tasks, the rms image is computed using
;#   the following form:
;#   iraf.imcalc(input="input_wht.fits", output="output_rms.fits", equals="1.0/sqrt(im1)")
;#   This new solution is suggested by Marc Rafelsky and sobstitutes
;#   the old one:
;#   iraf.imcalc(input="input_wht.fits", output="output_rms.fits", equals="value/sqrt(im1)")
;###############################################################
pro tweakprep_IB1,field,path0

path = expand_path(path0)+'/aXe/'+field+'/'
;path = '/Volumes/Kudo/DATA/WISPS/aXe/Par288-full/'
;tweakprep,'Par288-full','/Volumes/Kudo/DATA/WISPS'


readcol,path+'DATA/DIRECT_GRISM/F110_clean.list',f110_list,format=('A')
readcol,path+'DATA/DIRECT_GRISM/F160_clean.list',f160_list,format=('A')
if  f160_list[0] eq 'none' then begin
   readcol,path+'DATA/DIRECT_GRISM/F140_clean.list',f140_list,format=('A')
endif

   spawn,'mkdir '+path+'DATA/DIRECT/DIRECT_orig/'
   spawn,'cp '+path+'DATA/DIRECT/*_flt_clean.fits '+path+'DATA/DIRECT/DIRECT_orig/'
   spawn,'cp '+path+'DATA/DIRECT_GRISM/F110_clean.list '+path+'DATA/DIRECT/'
   spawn,'cp '+path+'DATA/DIRECT_GRISM/F160_clean.list '+path+'DATA/DIRECT/'

 openw,1, path+'DATA/DIRECT/direct_clean.list'
 openw,11, path+'DATA/DIRECT/direct_crclean.list'
 openw,2, path+'DATA/DIRECT/direct_clean_catfile.list'
 openw,22, path+'DATA/DIRECT/direct_crclean_catfile.list'

 openw,3, path+'DATA/DIRECT/F110_crclean.list'
 openw,4, path+'DATA/DIRECT/F160_crclean.list'
 openw,5, path+'DATA/DIRECT/F140_crclean.list'

;   prepare the lists and generate
;   catalogs on the *flt files (no CR
;   removal yet)
   ;**************************************
for i = 0, n_elements(f110_list)-1 do begin
  ; Modified by I.B.  
; spawn,'cp '+path+'DATA/DIRECT_GRISM/'+f110_list[i]+' '+path+'SEX/'
    spawn,'cp '+path+'DATA/DIRECT/'+f110_list[i]+' '+path+'SEX/'
   printf,1,f110_list[i]
   printf,2,f110_list[i],'  ',f110_list[i]+'.coo'
   printf,11,strmid(f110_list[i],0,19),'_crclean.fits'
   printf,22,strmid(f110_list[i],0,19),'_crclean.fits   ',strmid(f110_list[i],0,19),'_crclean.fits.coo'
   printf,3,strmid(f110_list[i],0,9),'_flt_clean_crclean.fits'

    h1=headfits(path+'DATA/DIRECT/'+f110_list[i]) 
    exptime1=strcompress(sxpar(h1,'EXPTIME'),/remove_all)
    if exptime1 gt 1041 then det='1.9'
    if exptime1 le 1041 then det='2.3'
    spawn,'sex '+path+'DATA/DIRECT/'+f110_list[i]+'[0] -c '+path+'SEX/config.sex -catalog_name '+path+$
    'SEX/'+f110_list[i]+'.coo -mag_zeropoint 26.83  -WEIGHT_TYPE MAP_WEIGHT -weight_image '+path+'DATA/DIRECT/'+strmid(f110_list[i],0,19)+'.fits -parameters_name '+path+$
    'SEX/config.param -filter Y -filter_name '+path+'SEX/gauss_2.0_5x5.conv -detect_minarea 6 -detect_thresh '+det+$
    ' -ANALYSIS_THRESH 2.0 -DEBLEND_NTHRESH 64 -DEBLEND_MINCONT 0.005 -GAIN '+exptime1+' -STARNNW_NAME '+path+$
    'SEX/default.nnw -CHECKIMAGE_NAME '+path+'SEX/'+strmid(f110_list[i],0,19)+'_seg.fits'
   spawn,'cp '+path+'SEX/'+f110_list[i]+'.coo '+path+'DATA/DIRECT/'
endfor

if f160_list[0] ne 'none' then begin
   for i = 0, n_elements(f160_list)-1 do begin
; Modified by I.B.  
; spawn,'cp '+path+'DATA/DIRECT_GRISM/'+f160_list[i]+' '+path+'SEX/'
   spawn,'cp '+path+'DATA/DIRECT/'+f160_list[i]+' '+path+'SEX/'
   printf,1,f160_list[i]
   printf,2,f160_list[i],'  ',f160_list[i]+'.coo'
   printf,11,strmid(f160_list[i],0,19),'_crclean.fits'
   printf,22,strmid(f160_list[i],0,19),'_crclean.fits   ',strmid(f160_list[i],0,19),'_crclean.fits.coo'
   printf,4,strmid(f160_list[i],0,9),'_flt_clean_crclean.fits'
;   printf,41,strmid(f160_list[i],0,9),'_flt_clean_crclean_sci.fits'
;   printf,42,strmid(f160_list[i],0,9),'_flt_clean_crclean_wht.fits'
    h2=headfits(path+'DATA/DIRECT/'+f160_list[i]) 
    exptime2=strcompress(sxpar(h2,'EXPTIME'),/remove_all)
    det='2.3'
    spawn,'sex '+path+'DATA/DIRECT/'+strmid(f160_list[i],0,19)+'.fits[0] -c '+path+'SEX/config.sex -catalog_name '+path+$
      'SEX/'+f160_list[i]+'.coo  -mag_zeropoint 25.96  -WEIGHT_TYPE MAP_WEIGHT -weight_image '+path+'DATA/DIRECT/'+strmid(f160_list[i],0,19)+'.fits -parameters_name '+path+$
      'SEX/config.param -filter Y -filter_name '+path+'SEX/gauss_2.0_5x5.conv -detect_minarea 6 -detect_thresh '+det+$
      ' -ANALYSIS_THRESH 2.0 -DEBLEND_NTHRESH 64 -DEBLEND_MINCONT 0.005 -GAIN '+exptime2+' -STARNNW_NAME '+path+$
      'SEX/default.nnw -CHECKIMAGE_NAME '+path+'SEX/'+strmid(f160_list[i],0,19)+'_seg.fits'
    spawn,'cp '+path+'SEX/'+f160_list[i]+'.coo '+path+'DATA/DIRECT/'
  endfor
endif else begin
   for i = 0, n_elements(f140_list)-1 do begin
; spawn,'cp '+path+'DATA/DIRECT_GRISM/'+f140_list[i]+' '+path+'SEX/'
   spawn,'cp '+path+'DATA/DIRECT/'+f140_list[i]+' '+path+'SEX/'
   printf,1,f140_list[i]
   printf,2,f140_list[i],'  ',f140_list[i]+'.coo'
   printf,11,strmid(f140_list[i],0,19),'_crclean.fits'
   printf,22,strmid(f140_list[i],0,19),'_crclean.fits   ',strmid(f140_list[i],0,19),'_crclean.fits.coo'
   printf,5,strmid(f140_list[i],0,9),'_flt_clean_crclean.fits'
;   printf,51,strmid(f140_list[i],0,9),'_flt_clean_crclean_sci.fits'
;   printf,52,strmid(f140_list[i],0,9),'_flt_clean_crclean_wht.fits'
    h2=headfits(path+'DATA/DIRECT/'+strmid(f140_list[i],0,19)+'.fits') 
    exptime2=strcompress(sxpar(h2,'EXPTIME'),/remove_all)
    det='2.0'
   spawn,'sex '+path+'DATA/DIRECT/'+strmid(f140_list[i],0,19)+'.fits'+'[0] -c '+path+'SEX/config.sex -catalog_name '+path+$
         'SEX/'+f140_list[i]+'.coo -mag_zeropoint 26.46 -WEIGHT_TYPE MAP_WEIGHT -weight_image '+path+'DATA/DIRECT/'+strmid(f140_list[i],0,19)+'.fits'+$
         ' -parameters_name '+path+$
         'SEX/config.param -filter Y -filter_name '+path+'SEX/gauss_2.0_5x5.conv -detect_minarea 6 -detect_thresh '+det+$
         ' -ANALYSIS_THRESH 2 -DEBLEND_NTHRESH 64 -DEBLEND_MINCONT 0.005 -GAIN '+exptime2+' -STARNNW_NAME '+path+$
         'SEX/default.nnw -CHECKIMAGE_NAME  '+path+'SEX/'+strmid(f140_list[i],0,19)+'_seg.fits'
   spawn,'cp '+path+'SEX/'+f140_list[i]+'.coo '+path+'DATA/DIRECT/'
   endfor
endelse


;Now generate the python code to do the CR removal ONLY
openw,6,path+'DATA/DIRECT/'+'tweakprep.py'
   printf,6,'import os,string,time'
   printf,6,'import sys'
   printf,6,'import shutil'
   printf,6,'from pyraf import iraf'
   printf,6,'from iraf import stsdas, dither'
   printf,6,'from pyraf.irafpar import IrafParS'
   printf,6,'from stsci.tools import teal'
   printf,6,'import drizzlepac'
   printf,6,'from drizzlepac import tweakreg'
   printf,6,'from drizzlepac import astrodrizzle'
   printf,6,'from drizzlepac import tweakback'
   printf,6,'import glob'
   printf,6,'from stwcs import wcsutil'
   printf,6,'import stwcs.wcsutil.headerlet'
;   printf,6,'unlearn tweakreg force=yes'
;   printf,6,'unlearn imagefindpars force=yes'
;   printf,6,'unlearn astrodrizzle force=yes'
;   printf,6,'unlearn tweakback force=yes'

                                ; --------------- this step is to
                                ; align *flt.fits files in each
                                ; filter, separately, to the 1st
                                ; flt.fits image -------------

;first trial to check if the pos-shift makes sense, debug use only
;   printf,6,'tweakreg.TweakReg("@direct_clean.list",catfile="direct_clean_catfile.list", refimage="'+f110_list[0]+'",refcat="'+f110_list[0]+'.coo",updatehdr=False,updatewcs=False, xcol=2, ycol=3, fluxcol=12, fluxunits="mag", xyunits="pixels",  refxcol=7, refycol=8, refxyunits="degrees", rfluxcol=12, rfluxunits="mag", minobj=15, searchrad=2.0, sigma=3.0, nclip=3, shiftfile=True,outshifts="shift_pos_flt.txt")'

;apply the pos-shift in the header
   printf,6,'tweakreg.TweakReg("@direct_clean.list",catfile="direct_clean_catfile.list", refimage="'+f110_list[0]+'",refcat="'+f110_list[0]+'.coo",updatehdr=True, wcsname="shift0",updatewcs=False,xcol=2, ycol=3, fluxcol=12, fluxunits="mag", xyunits="pixels",  refxcol=7, refycol=8, refxyunits="degrees", rfluxcol=12, rfluxunits="mag", minobj=15, searchrad=1.0, sigma=4.0, nclip=3, shiftfile=True,outshifts="shift_pos_flt_final.txt",fitgeometry="shift")'

                                ; --------------- this step is for CR
                                ; removal to generate the
                                ; *crclean.fits*  --------------
   printf,6,'teal.unlearn("astrodrizzle")'
   num = n_elements(f110_list)
   if num gt 1 then begin
      printf,6,'astrodrizzle.AstroDrizzle("@F110_clean.list",output="F110W_orig",num_cores=5,final_wcs=False,final_wht_type="IVM",build=True,updatewcs=False,clean=True,preserve=False, driz_cr_corr=True, driz_combine=True)'
   endif else begin
      printf,6,'astrodrizzle.AstroDrizzle("@F110_clean.list",output="F110W_orig",num_cores=5,final_wcs=False,final_wht_type="IVM",build=True,updatewcs=False,clean=True,preserve=False,median=False,blot=False,driz_cr=False)'
            spawn,'cp '+path+'DATA/DIRECT/'+f110_list[0]+' '+path+'DATA/DIRECT/'+strmid(f110_list[0],0,19)+'_crclean.fits'
   endelse
   printf,6,'iraf.imcopy(input="F110W_orig_drz.fits[1]", output="F110W_orig_sci.fits")'
   printf,6,'iraf.imcopy(input="F110W_orig_drz.fits[2]", output="F110W_orig_wht.fits")'
   ;;; printf,6,'iraf.imcalc(input="F110W_orig_wht.fits", output="F110W_orig_rms.fits", equals="1.5/sqrt(im1)")'
   printf,6,'iraf.imcalc(input="F110W_orig_wht.fits", output="F110W_orig_rms.fits", equals="1.0/sqrt(im1)")' ; NEW solution (Suggeste by Marc)

   if f160_list[0] ne 'none' then begin
         num2 = n_elements(f160_list)
         if num2 gt 1 then begin
            printf,6,'astrodrizzle.AstroDrizzle("@F160_clean.list",output="F160W_orig",num_cores=5,final_wcs=False,final_wht_type="IVM",build=True,updatewcs=False,clean=True,preserve=False, driz_cr_corr=True, driz_combine=True)'
         endif else begin
            printf,6,'astrodrizzle.AstroDrizzle("@F160_clean.list",output="F160W_orig",num_cores=5,final_wcs=False,final_wht_type="IVM",build=True,updatewcs=False,clean=True,preserve=False,median=False,blot=False,driz_cr=False)'
            spawn,'cp '+path+'DATA/DIRECT/'+f160_list[0]+' '+path+'DATA/DIRECT/'+strmid(f160_list[0],0,19)+'_crclean.fits'
         endelse
         
   printf,6,'iraf.imcopy(input="F160W_orig_drz.fits[1]", output="F160W_orig_sci.fits")'
   printf,6,'iraf.imcopy(input="F160W_orig_drz.fits[2]", output="F160W_orig_wht.fits")'
   ;;; printf,6,'iraf.imcalc(input="F160W_orig_wht.fits", output="F160W_orig_rms.fits", equals="1.5/sqrt(im1)")'
   printf,6,'iraf.imcalc(input="F160W_orig_wht.fits", output="F160W_orig_rms.fits", equals="1.0/sqrt(im1)")' ; NEW solution (Suggeste by Marc)
endif else begin
   num3 = n_elements(f140_list)
      if num3 gt 1 then begin
         printf,6,'astrodrizzle.AstroDrizzle("@F140_clean.list",output="F140W_orig",num_cores=5,final_wcs=False,final_wht_type="IVM",build=True,updatewcs=False,clean=True,preserve=False, driz_cr_corr=True, driz_combine=True)'
      endif else begin
         printf,6,'astrodrizzle.AstroDrizzle("@F140_clean.list",output="F140W_orig",num_cores=5,final_wcs=False,final_wht_type="IVM",build=True,updatewcs=False,clean=True,preserve=False,median=False,blot=False,driz_cr=False)'
            spawn,'cp '+path+'DATA/DIRECT/'+f140_list[0]+' '+path+'DATA/DIRECT/'+strmid(f140_list[0],0,19)+'_crclean.fits'
      endelse      
   printf,6,'iraf.imcopy(input="F140W_orig_drz.fits[1]", output="F140W_orig_sci.fits")'
   printf,6,'iraf.imcopy(input="F140W_orig_drz.fits[2]", output="F140W_orig_wht.fits")'
   ;;; printf,6,'iraf.imcalc(input="F140W_orig_wht.fits", output="F140W_orig_rms.fits", equals="1.5/sqrt(im1)")'
   printf,6,'iraf.imcalc(input="F140W_orig_wht.fits", output="F140W_orig_rms.fits", equals="1.0/sqrt(im1)")' ; NEW solution (Suggeste by Marc)
endelse   


close,1,2,3,4,11,22,5,6  ;31,32,41,42,51,52
free_lun,1,2,3,4,11,22,5,6 ;;31,32,41,42,51,52

end
