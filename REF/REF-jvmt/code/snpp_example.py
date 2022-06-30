#
# Name:
#	snpp
# PURPOSE:
#	calculate the noise and generate a simulated spectrum with noise for any given template or spectrum.
#IUTPUTS:
#   targetmag      the surface brightness of the target you want to calculate the S/N (defult: 17 mag/arcsec^2)
#   galtpl         the filename of star-forming galaxy template you want to use.
#   filtera         the filter you chosed to estimate the S/N (defult: sdss_g)
#   readnoise      read noise, in e/pix. (defult: 5.0)
#   fovp           diameter of fiber (or spaxel) in arcsec (defult: 0.2 arcsec)
#   npixel_width   the width of the spectrum on the CCD (defult: 2.0)
#   obstime        in seconds, single integration time (defult: 300s)
#   repeatnum      repeat number (defult: 20.0)
#   skyr           r band sky brightness in Johnson r mag/arcsec^2 unit (defult: 21 mag/arcsec^2)
#   qinput         the throughput correct factor (defult: 1.0)
#   skyperpixel    a second way of estimating the Sky, if know the sky photon number per pixel
#
#IDL
# v5: 15 August 2018      writen by Lei Hao, rivised by Jun Yin
# v7: 10 Sep 2019  by Jun Yin
#
#python
# V1.04beta by Mengting Ju
##################################################################
##################################################################

from __future__ import print_function
from astropy.io import fits
import numpy as np
from snpp import *


#################################################################
##################################################################

def snpp_example(): 
    
    select=1 # 1 or 2
    resu=snpp_model(select)
    wavearr,galflux,targetmag,filtera=resu[0],resu[1],resu[2],resu[3]
    
    filename='mg_ell_17_30020_blue.fits'

    ss=snpp(wavearr=wavearr,galflux=galflux,
            filename=filename, dlambda=0.98,
            readnoise=5.5,fovp=3.2,npixel_width=5.0,
            obstime=300,repeatnum=20,skyr=21.0, qinput=0.05,
            snlimit=1,targetmaglimit=targetmag,filtera=filtera,
            aa=0.,bb=0.)
    
    '''
    ra=np.zeros(6*4*8).reshape(6,4,8)
    slim=[1,3,5,10] #4
    reqtime=[1,3,5,10,15,20,30,40] #8
    obt=[300,600,900,1200,1500,1800] #6
    galtpl='../obs/SFgal_tpl/SFgal_texp_FeH0_tau5_Ew10_AGN1.fits'
    filtera='../obs/filters/sdss_g0.par'
    for x in range(6):
        for i in range(4):
            for j in range(8):
    
                #x,i,j=3,1,5
                filename='mg_ang_17_30020.fits'

                if(os.path.exists(filename))==1:
                    os.remove(filename)

                ss=snpp(wavearr=wavearr,galflux=galflux,
                        filename=filename, dlambda=1.30,
                        readnoise=5.5,fovp=3.2,npixel_width=5.0,
                        obstime=obt[x],repeatnum=reqtime[j],skyr=21.0, qinput=0.05,
                        snlimit=slim[i],targetmaglimit=targetmag,filtera=filtera,
                        aa=0.,bb=0.)

                ra[x,i,j]=ss.printf()
               
    hun1=fits.PrimaryHDU()
    hun2=fits.ImageHDU(ra)
    hdulist = fits.HDUList([hun1,hun2])
    hdulist.writeto('snresult_5.5_red.fits')
    '''
   ################################################################

def snpp_model(s):
    
    if s==1:
    
        #select model and magnitude
        targetmag=17.
        galtpl='../obs/SFgal_tpl/SFgal_texp_FeH0_tau1_Ewd.fits'
        #galtpl='../obs/SFgal_tpl/SFgal_texp_FeH-2_tau10_Ew50.fits'
        filtera='../obs/filters/sdss_g0.par'

        result=input_mag_model(targetmag,galtpl,filtera,0.98,3490.71,5500)
        #result=input_mag_model(targetmag,galtpl,filtera,1.30,4579.30,7200)
        wavearr=result[0]   #A
        galflux=result[1]   #10^-12 erg/s/A/cm2   
        

    elif s==2:
        
        #select put in wave and flux    
        filee=fits.open('MockGal-M21Z0.01-W350n1000n.fits')
        fluxx=filee[1].data  #erg/s/A/cm2
        wavee=filee[2].data    #A

        result=input_wave_flux(wavee,fluxx)
        wavearr=result[0]  #A
        galflux=result[1]  #10^-12 erg/s/A/cm2
    
    return wavearr,galflux,targetmag,filtera


#########################################################

if __name__=='__main__':
    snpp_example()
    
