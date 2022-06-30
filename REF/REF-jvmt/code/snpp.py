#+
# Name:
#	snpp
# PURPOSE:
#	calculate the S/N per pixel for CSST and simulate a noisy spectrum for any given template.
# CALLING SEQUENCE:
#   snpp,limitmag, repeatnum=10,obstime=300,targetmag=18,/skyperpixel,$
#         galtpl=,wavearr=wavearr,mockgal=mockgal,galflux=galflux
#   plot, wavearr, galflux  ; the input galaxy template
#   plot, wavearr, mockgal  ; the output spectrum with noise
#     
# INPUTS:
# OPTIONAL IUTPUTS:
#   darkcurrent    dark current, in e/s/pix, (defult: 0.0017)
#   deltal         the delta lambda per pixel, in unit of nm (defult: 0.1755555 nm)
#   fovp           diameter of fiber (or spaxel) in arcsec (defult: 3.2 arcsec)
#   filtera         the filter you chosed to estimate the S/N (defult: bessell_V)
#   galtpl         the filename of star-forming galaxy template you want to use. 
#                  They are in the ../obs/SFgal_tpl/ folder (default: SFgal_texp_FeH0_tau5_Ew10.fits)
#   lambdac        the noise at this wavelength wanted (defult: 550 nm)
#   npixel_width   the width of the spectrum on the CCD (defult: 3.0)
#   obstime        in seconds, single integration time (defult: 300s)
#   outfile        the output file name (defult: '../results/noise.dat' )
#   qinput         the throughput correct factor (defult: 1.0)
#   readnoise      read noise, in e/pix. (defult: 4.0)
#   redshift       the redshift of the target spectrum. (defult: 0.0)
#   repeatnum      repeat number (defult: 1.0)
#   skyperpixel    a second way of estimating the Sky, if know the sky photon number per pixel
#   skyv           V band sky brightness in Johnson V mag/arcsec^2 unit (defult: 22.5 mag/arcsec^2)
#  slitwidth      suit to the slit case. the length assumed to be 0.15 arcsec
#   snlimit        S/N limit (defult: 1.0)
#   specsample     pixels per spectral resolution element (defult: 2)
#   targetmag      the surface brightness of the target you want to calculate the S/N (defult: 22 .5 mag/arcsec^2)
#   teld           diameter of the telescope, in cm unit. (defult: d=200 cm)
# OUTPUTS:
#   limitmag       the Vband surface brightness needed to achieve the S/N limit (defult: 1.0)
# OPTIONAL OUTPUTS:
#   limitemi       the medien of Flambda*dlambda*sampling value of Ha line
#   limitemif      the limit detection of Ha flux 
#   snmean         the median S/N of the whole input target spectrum (mag_v=targetmag)
#   wavearr        the wave array (nm)
#   galflux        the input galaxy flux  (1e-13 erg/s/cm2/nm)
#   mockgal        the mocked galaxy flux  with noise (1e-13 erg/s/cm2/nm)
#
# v5: 15 August 2018      writen by Lei Hao, rivised by Jun Yin
# v7: 10 Sep 2019  by Jun Yin
#     1) remove the function im_filtermag, so do not need the Kcorrect package anymore.
#     2) 
#python 
# v7: 22 Sep 2019 by Mengting Ju
#1.02alpha : 26 Sep 2019 
#1.03alpha: 29 Sep 2019
#1.03beta: 14 Oct 2019
#-
#####################################################################################
#####################################################################################

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
import matplotlib
import pandas as pd
from scipy import interpolate
from sympy import *
import os
from scipy.integrate import simps


####################################################################################
####################################################################################
def filteraa(filtera):
    #load the filters          
    filtersel=filtera #'../sdss_g0.par'            
    filterpath='./'
    filterfile=filterpath+filtersel   # ;fluxfilter: max=1, min=0, no particular unit
    print(filterfile)   

    ia=0
    with open(filterfile,'r') as fh:
        for line in fh:
            if line.startswith('#'):
                ia=ia+1
                continue

    band=pd.read_csv(filterfile,sep='\s+',header=None,skiprows=ia)
    wavefilter=np.array(band[0])
    fluxfilter=np.array(band[1])
    wavefilter=wavefilter  # A
    vmin=wavefilter[0]
    vmax=wavefilter[-1]

    # find the central wavelength, effective wavelength, and FWHM of the given filter
    filtermid=(vmax-vmin)*0.5  #A, central wavelength
    dwave=wavefilter[1:]-wavefilter[:-1]
    filtereff=np.nansum(dwave*wavefilter[1:]*fluxfilter[1:])/np.nansum(dwave*fluxfilter[1:]) #A, effective wavelength
    rmax=np.max(fluxfilter)
    nnn=np.where(fluxfilter > 0.5*rmax)[0]
    FWHMmin=wavefilter[nnn[0]]
    FWHMmax=wavefilter[nnn[-1]]
    filterwid=FWHMmax-FWHMmin  #A, FWHM

    return wavefilter,fluxfilter,vmin,vmax,filtereff,FWHMmin,FWHMmax
################################################################################
#select input parameters
def input_mag_model(targetmag,galtpl,filtera,dlambda,a,b):

    #filter
    resulte=filteraa(filtera)
    wavefilter=resulte[0]
    fluxfilter=resulte[1]
    vmin=resulte[2]
    vmax=resulte[3]
    filtereff=resulte[4]

    ####################################################################################
    # define wavelength array,
    #cover the range of 350nm to 1050nm, depend on the spectral resolution wanted. 

    delta_lambda=dlambda # has to be in unit of A
    print('delta_lambda:', delta_lambda)
    
    narray=int((b-a)/delta_lambda) 
    wavearr=a+delta_lambda*np.float64(np.arange(narray))
    
    # select out the array of V band filter
    nar=int((10000-3500)/delta_lambda) 
    wavearr1=3500+delta_lambda*np.float64(np.arange(nar))
    
    ii=np.logical_and(wavearr1 >= vmin, wavearr1 <= vmax)
    wavetmp2=wavearr1[ii]
    x=np.interp(wavetmp2,wavefilter,fluxfilter)
    integratef4=x*wavetmp2
    integconst=simps(integratef4,wavetmp2) # int(lambda*Rlambda*dlambda)

    lambdav=filtereff #A
    ###############################################################
    # define basic target brightness, parameters constantly change
    itarget=targetmag    # in Johnson V mag/arcsec^2 unit
    print('itarget:',itarget)

    itarget_jy=3631.0*10**(-itarget/2.5+3.0)  # target flux in V in mJy/arcsec^2 unit
    itarget_nm=itarget_jy*3.0/(lambdav/100.0)**2 #target flux in V in 10^(-12)erg/s/cm^2/A/arcsec^2 

    tplfile=galtpl
    print('tplfile:',tplfile)   

    sfgal=fits.open(tplfile)
    wavegal=sfgal[1].data['wavelength'] # A
    galflux2=sfgal[1].data['flux']
    galflux1=np.interp(wavearr1,wavegal,galflux2)

    #;normalize the galaxy spectrum to the V band magnitude specified.
    ii=np.logical_and(wavegal >= vmin, wavegal <= vmax)
    wavetmp=wavegal[ii]
    fluxtmp=galflux2[ii]
    x=np.interp(wavetmp,wavefilter,fluxfilter)
    vfluxtmp=x*wavetmp*fluxtmp #bandpass*lambda*F_gal_lambda
    galintegrate=simps(vfluxtmp,wavetmp)
    galnorm=itarget_nm*integconst/galintegrate
    galflux1=galnorm*galflux1   # the unit should now be in 10^(-12)erg/s/A/cm^2/arcsec^2 
    
    #galflux=np.interp(wavearr,wavearr1,galflux1)
    galflux=galflux1
    
    return wavearr, galflux, wavearr1


#######################################################################################################
# select put in wave and flu
def input_wave_flux(wave,flux,dlambda):

    wave=wave
    flux=flux
    delta_lambda=dlambda
    narray=int((wave[-1]-wave[0])/delta_lambda) 
    wavearr=wave[0]+delta_lambda*np.float64(np.arange(narray))

    galflux=np.interp(wavearr,wave,flux)      #erg/s/A/cm2

    return wavearr, galflux*10**12

################################################################################
#####################################################################################
class snpp(object):
    def __init__(self, wavearr=False,galflux=False,
                 filename=False,dlambda=False,
                 readnoise=5.5,fovp=3.2,npixel_width=5.0,
                 obstime=300,repeatnum=20.0,skyr=21.0, qinput=0.05,
                 snlimit=10.0,targetmaglimit=None,filtera=False,
                 aa=0.,bb=0.):

        #Do extensive checking of possible input errors
        #
         
        self.qinput=qinput 
        self.fovp=fovp
        self.obstime=obstime
        self.skyr=skyr
        self.repeatnum=repeatnum
        self.readnoise=readnoise
        self.npixel_width=npixel_width
        
        ###########################################################################
        #some basic unchanged parameters     
        d=240        # diameter of the telescope, in cm unit
        obscure=0.1  #effective central obscuration, no unit
        telarea=3.14159/4.0*d*d*(1.0-obscure)  #effective area of the telescope, cm^2
        darkc=0.017   #dark current, in e/s/pix
        planckh=6.626    # 10^{-27} erg*s
        cc=3.0   # speed of light, 10^{18} A/s
        #slitunit=0.074  # arcsec. the length of slit which conresponds to a pixel length on IFU CCD 
        
        rn=readnoise  #read noise, in e/pix
        print('readnoise:', rn)
                       
        npixw=npixel_width
        print('npixel_width:',npixw)
               
        obst=obstime  # in seconds, single integration time
        print('obstime:',obst)
        
        repn=repeatnum   # repeating time
        print('repeatnum:',repn)

        iskyr0=skyr  # in Johnson V mag/arcsec^2 unit
        print('skyr:', iskyr0)
        
        sampling=1.0
        #print('specsample:',sampling)
        
        delta_lambda=dlambda
        #####################################################################
        
        
        
        # some less basic parameters, may change, but not often
        
        #throughput=0.05
        #lambdaq=np.array(throughput[0])*10 # A
        qtot=0.1 #; throughput of the whole system,
 

        #;assuming the total throughput cannot reach the theory value, 0.3 is the upper limit. 
        #qtot[qtot>=0.3]=0.3 
        
        qinput=qinput
        print('qinput:', qinput)
        
        q=qtot*qinput #*qe ;qtot of CSST already includes the CCD efficiency 
        

        fov2=((fovp)**2)*3.14159/4.0
        print('fov2:', fov2)               
        ##############################################################################
        
        # SKY

        
        resulte=filteraa('../obs/filters/sdss_g0.par')
        wavefilter=resulte[0]
        fluxfilter=resulte[1]
        vmin=resulte[2]
        vmax=resulte[3]
        filtereff=resulte[4]
        FWHMmin = resulte[5]
        FWHMmax = resulte[6]

        # select out the array of r band filter
        nar=int((10000-3500)/delta_lambda) 
        wavearr1=3500+delta_lambda*np.float64(np.arange(nar))
        
        ii=np.logical_and(wavearr1 >= vmin, wavearr1 <= vmax)
        
        wavetmp2=wavearr1[ii]
        x=np.interp(wavetmp2,wavefilter,fluxfilter)
        integratef4=x*wavetmp2    
        integconst=simps(integratef4,wavetmp2) # int(lambda*Rlambda*dlambda)

        #####################################################################

        #define r band sky brightness

        lambdar=filtereff   #in A

        #sky brightness corresponding to this sky magnitude
        iskyr0_jy=3631.0*10**(-iskyr0/2.5+3.0)  # sky flux in V in mJy/arcsec^2 unit
        iskyr0_nm=iskyr0_jy*3.0/(lambdar/100.0)**2 #sky flux in V in 10^(-12)erg/s/cm^2/A /arcsec^2 

        #readin the ground sky spectrum 
        skybg_50=pd.read_csv('../obs/skybg_50_10.dat',sep='\s+',header=None,skiprows=14)
        wavesky=np.array(skybg_50[0])*10 #in A
        fluxsky1=np.array(skybg_50[1])/10 #phot/s/A/arcsec^2/m^2
        fluxsky2=fluxsky1/wavesky*1.98 #change the sky flux unit to 10^(-12)erg/s/cm^2/A/arcsec^2


        #This fluxsky is in unit of phot/s/A/arcsec^2/m^2, to convert it to F_lambda/arcsec^2, 
        #need to do fluxsky(phot/s/A/arcsec^2/m^2)*h(6.625*10^{-27}erg.s)*nu(1/s)*10{-4}(m^2/cm^2)
        #=fluxsky*c(3.0*10^{18}A/s)/lambda(A)*6.6*10{-31} erg/s/cm^2/A/arcsec^2
        #=fluxsky/lambda*1.98*10^{-12}erg/s/cm^2/A/arcsec^2 

        #find out the normalization of the sky,
        ii=np.logical_and(wavesky >= vmin, wavesky <= vmax)
        wavetmp=wavesky[ii]
        fluxtmp=fluxsky1[ii]

        x=np.interp(wavetmp,wavefilter,fluxfilter)
        vfluxtmp=x*fluxtmp*1.98  #bandpass*lambda*F_sky_lambda(fluxsky2)=bandpass*fluxsky*1.98
        skyintegrate=simps( vfluxtmp ,wavetmp)
        skynorm=iskyr0_nm*integconst/skyintegrate 
        fluxsky3=np.interp(wavearr1,wavesky,fluxsky2)
        fluxsky4=fluxsky3*skynorm   
        # get the sky spectrum in wavearr grid, the unit should now be the same as fluxvega: 10^(-12) erg/s/A/cm^2/arcsec^2 
        
        #fluxskypp=np.interp(wavearr,wavearr1,fluxsky4)                
        fluxskypp=fluxsky4
        ##########################################################################
        
        #define observation information, parameters constantly change
        
        narray=len(wavearr1)
        expf2=np.zeros(narray)
        mockgal=np.zeros(narray)
        tmp=np.zeros(narray)
        lista=np.zeros(narray*10).reshape(narray,10)
        
        for i in range(narray):
            lambda0=wavearr1[i]
            qlambda=q
            hv=planckh*cc/lambda0 #;10^{-9}erg
            delta_hz=cc*delta_lambda/lambda0/lambda0 #;10^18 1/s
            delta_shz=delta_hz*sampling
            
            #now that many fluxes are in 10^(-12)erg/s/A/cm^2, to convert it to Jy, need to multiple: 
            #lambda0^2/c(in A)=lambda0^2(A)/(3.*10^(18))*10^(-12)erg/s/Hz/cm^2
            #=lambda^2(A)*3.33*10^(-31)erg/s/Hz/cm^2=lambda^2(A)*3.33*10^(-8)Jy
            #=lambda^2(A)*0.0333uJy

            #find out sky value at lambda0    
            #calculate n_sky/pixel
            isky=fluxskypp[i]*lambda0**2*0.0333*fov2   #in uJy/spaxel unit
            iskyall=isky*telarea/1000.0   #in 10-26 erg/s/Hz /spaxel
            fsky=qlambda*iskyall*delta_hz   #10^{-8} erg/s /spaxel
            fsky=fsky*sampling
            nsky=fsky/hv*10.0   #in unit of #e/s /spaxel
                
            #calculate n_source/pixel
            isource=galflux[i]*lambda0**2*0.0333*fov2   #in uJy/spaxel unit
            isall=isource*telarea/1000.0   #in 10-26 erg/s/Hz /spaxel
            fs=qlambda*isall*delta_hz   #10^{-8} erg/s /spaxel
            fs=fs*sampling
            ns=fs/hv*10.0   #in unit of #e/s /spaxel

            darkn=(darkc*repn*obst*npixw*sampling)
            rnn2=rn**2*(repn*npixw*sampling)
            sourcenn=(ns*repn*obst)
            skynn=(nsky*repn*obst)
            tmp[i]=skynn

            nn1=np.sqrt(rnn2+darkn+skynn+sourcenn)  #total noise
            sn1=repn*ns*obst/nn1  #S/N
            nn=np.sqrt(rnn2+darkn+skynn)  #system noise
            
            
            mockgal[i]=galflux[i]+galflux[i]/sn1*np.random.randn(1,1)[0][0]  #in 10^{-12} erg/s/A/cm^2
            
            lista[i,:]=[lambda0, sn1, galflux[i]*10**5, nn1,\
                        np.sqrt(sourcenn), nn, np.sqrt(rnn2),np.sqrt(darkn), \
                        np.sqrt(skynn), mockgal[i]*10**5]
            
            #set the detection limit 
            detlimit=snlimit
            
            nntmp=detlimit**2+np.sqrt(detlimit**4+4.0*detlimit**2*nn**2)
            nntmp=nntmp/2.0            
            #calculate detection limit in uJy and mag
            fnn=nntmp
            f1=fnn/obst/repn    #in e/s
            f2=f1*hv       #in 10^{-9} erg/s
            f3=f2/delta_lambda #in 10^{-9} erg/s/A
            f1=f3/telarea #in 10^{-9} erg/s/A/cm^2
            f2=f1/qlambda    #in 10^{-9} erg/s/A/cm^2
            expf2[i]=f2/fov2*100000.   # in 10^{-14} erg/s/A/cm^2/arcsec^2
    ############################################################################################
        self.bb=lista
       
        
        lambda1=wavearr
        snarray=np.interp(wavearr,wavearr1,lista[:,1])
        galflux=np.interp(wavearr,wavearr1,lista[:,2]) 
        tot_noise=np.interp(wavearr,wavearr1,lista[:,3]) 
        sc_noise=np.interp(wavearr,wavearr1,lista[:,4])            
        sys_noise=np.interp(wavearr,wavearr1,lista[:,5])
        readnoise=np.interp(wavearr,wavearr1,lista[:,6])
        dark_noise=np.interp(wavearr,wavearr1,lista[:,7])
        sky_noise=np.interp(wavearr,wavearr1,lista[:,8])
        mockgal=np.interp(wavearr,wavearr1,lista[:,9])  
        
        
        # write file
        namedat=np.array(['lambda','S/N','tar_flux','tot_noise','sc_noise', \
                          'sys_noise', 'readnoise','dark_noise', 'sky_noise', 'mockgal'])
        unit=np.array(['A', ' ','1e-17 erg/s/cm^2/A',\
                       '#e','#e','#e','#e','#e','#e', '1e-17 erg/s/cm^2/A'])
        
        hdr=fits.Header()
        for i in range(len(namedat)):
            hdr[str(i)]=unit[i]
        hun1=fits.PrimaryHDU(header=hdr)
        hun2=fits.ImageHDU(wavearr,name='lambda')
        hun3=fits.ImageHDU(snarray,name='S/N')
        hun4=fits.ImageHDU(galflux,name='tar_flux')
        hun5=fits.ImageHDU(tot_noise,name='tot_noise')
        hun6=fits.ImageHDU(sc_noise,name='sc_noise')
        hun7=fits.ImageHDU(sys_noise,name='sys_noise')
        hun8=fits.ImageHDU(readnoise,name='readnoise')
        hun9=fits.ImageHDU(dark_noise,name='dark_noise')
        hun10=fits.ImageHDU(sky_noise,name='sky_noise')
        hun11=fits.ImageHDU(mockgal,name='mockgal')
        hdulist = fits.HDUList([hun1,hun2,hun3,hun4,hun5,hun6,hun7,hun8,hun9,hun10,hun11])

        print('output filt:',filename)
        hdulist.writeto('./'+filename)
        
        #####################################################################
        resulte=filteraa(filtera)
        wavefilter=resulte[0]
        fluxfilter=resulte[1]
        vmin=resulte[2]
        vmax=resulte[3]
        filtereff=resulte[4]
        FWHMmin = resulte[5]
        FWHMmax = resulte[6]
        
        ii=np.logical_and(wavearr1 >= FWHMmin , wavearr1 <= FWHMmax)
        wavetmp=wavearr1[ii]
        snm=sorted(lista[:,1][ii])[len(wavetmp)//2]
        #the median SN of FWHM range to acchieve the sn limit
        im=np.where(lista[:,1][ii] == snm)[0]
        fact=expf2[ii][im]*0.01/((lista[:,2][ii][im])*10**-5)
        limitmag=-2.5*np.log10(fact)+targetmaglimit
        self.aa=limitmag
        print('limitmag:',limitmag)
    ############################################################################
    
    def printf(self):
        return self.aa

    #############################################################################
    def fits(self):
        return self.bb
    ###############################################################################
    print('The END!')
        
#######################################################################################################

