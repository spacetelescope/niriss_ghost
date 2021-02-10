import glob
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

import sys
DIR_CODE = './code/'
sys.path.append(DIR_CODE)
from utils import get_gap,get_ghost,str2bool

from astropy.io import fits,ascii
import astropy.coordinates as coord
import astropy.units as u
from astropy.wcs import WCS
from astroquery.vizier import Vizier


def tweak_imaege2(id_gst, infile, segfile, outfile=None, DQ_SET=1):
    '''
    This script tweaks image2 DQ arrays, by making positions of detected ghost DQ=1 ("Do not use").

    Return:
    =======
    Save fits files with _gsf.fits
    '''
    if not outfile == None:
        outfile = infile.replace('.fits','_gst.fits')

    os.system('cp %s %s'%(infile, outfile))
    
    seg = fits.open(segfile)[0].data
    with fits.open(outfile, mode='update') as hdul:
        for id in id_gst:
            con = np.where(seg == id)
            hdul['DQ'].data[con] = DQ_SET
        
        hdul.flush()  # changes are written back to original.fits


if __name__ == "__main__":
    '''
    '''
    parser = argparse.ArgumentParser(description='Run the NIRISS ghost detection script.')
    parser.add_argument('--f_verbose',default=False,help='Print comments.', type=str2bool)
    parser.add_argument('--rlim',default=10,help='Ghost search radius (pixel)', type=float)
    parser.add_argument('--frac_ghost',default=0.01,help='Flux fraction for ghosts', type=float)
    parser.add_argument('--f_tweak_imaege2',default=False,help='Tweak Image2 products.', type=str2bool)
    parser.add_argument('--o',default='./',help='Output directory.', type=str)
    #parser.add_argument('--filt',default='f150w',help='Flux fraction for ghosts', type=float)
    args = parser.parse_args()
    f_mirage = True

    f_verbose = args.f_verbose
    rlim = args.rlim
    frac_ghost = args.frac_ghost
    f_tweak_imaege2 = args.f_tweak_imaege2
    DIR_OUT = args.o    
    if not os.path.exists(DIR_OUT):
        os.mkdir(DIR_OUT)

    '''
    f_verbose = True
    rlim = 10 # pixel
    frac_ghost = 0.01
    f_tweak_imaege2 = True
    #filt = 'f150w'
    '''

    # Ghost detection must be done in distortion corrected frame
    '''
    DIR_DATA = './reduced/'
    file_root = 'jw01085001001_01101'
    infiles = ['%s/%s_00001_nis_cal.fits'%(DIR_DATA, file_root), \
            '%s/%s_00006_nis_cal.fits'%(DIR_DATA, file_root), \
            '%s/%s_00007_nis_cal.fits'%(DIR_DATA, file_root), \
            '%s/%s_00008_nis_cal.fits'%(DIR_DATA, file_root)]
    '''
    DIR_DATA = './'
    file_root = 'l3_nis_test_f150w'
    infiles = ['%s/%s_i2d.fits'%(DIR_DATA, file_root)]
    # Catalog;
    files_cat = [infile.replace('_i2d.fits','_cat_man.ecsv') for infile in infiles]

    for ff,infile in enumerate(infiles[:1]):
        '''
        '''
        
        # Read header;
        hd = fits.open(infile)[0].header
        hd1 = fits.open(infile)[1].header
        filt = hd['FILTER']
        pupil = hd['PUPIL']
        XOFFSET = hd['XOFFSET']
        YOFFSET = hd['YOFFSET']
        CDELT1 = hd1['CDELT1']
        
        if f_mirage:
            # image size difference;
            # May be specific only to mirage data.        
            xshift = (2048-hd1['NAXIS1'])/2.
            yshift = (2048-hd1['NAXIS2'])/2.
        else:
            xshift = 0
            yshift = 0

        if not filt == 'CLEAR':
            print('Grism data are currently not supported. Skipped.')
            continue

        file_cat = files_cat[ff]
        fd_cat = ascii.read(file_cat)

        # Empty Array for ghost
        flag_gst = np.zeros(len(fd_cat['id']), 'int')
        prob_gst = np.zeros(len(fd_cat['id']), 'float')
        id_gst = np.zeros(len(fd_cat['id']), 'int')

        fw = open('%s/ghs_%s.reg'%(DIR_OUT,file_root), 'w')
        fw.write('# xgst ygst\n')

        if True:
            # No1; Root method
            # Look into detected sources in Image2 catalog;
            gst = get_ghost(fd_cat['xcentroid'][:].value, fd_cat['ycentroid'][:].value, \
                            flux=fd_cat['source_sum'][:], filt=pupil,\
                            xshift=xshift, yshift=yshift)

            for ii in range(len(fd_cat['id'])):
                # This is because photutils is 0-based, while ds9 is not.
                shift = 1.0 
                if fd_cat['source_sum'][ii]>0:
                    #fw.write('%.3f %.3f\n'%(fd_cat['xcentroid'][ii].value, fd_cat['ycentroid'][ii].value))
                    fw.write('%.3f %.3f\n'%(gst[0][ii]+shift, gst[1][ii]+shift))

                # Check the position:
                rtmp = np.sqrt( (fd_cat['xcentroid'].value-gst[0][ii])**2 + (fd_cat['ycentroid'].value-gst[1][ii])**2 )
                iiy = np.argmin(rtmp)

                # Check if source at the predicted positions is brighter than object of [ii].
                if rtmp[iiy] < rlim and (fd_cat['source_sum'][ii]-fd_cat['source_sum'][iiy]) > 0:

                    id_gst[ii] = fd_cat['id'][iiy]
                    flag_gst[ii] = 1
                    
                    #if f_verbose:
                    #    print('Found an object at r=%.2f pixel'%(rtmp[iiy]))
                    
                    prob_pos = np.exp(-0.5 * rtmp[iiy]**2)
                    residual = np.abs(fd_cat['source_sum'][ii]*frac_ghost - fd_cat['source_sum'][iiy])
                    expectation = fd_cat['source_sum'][ii]*frac_ghost
                    prob_flux = np.exp(-0.5 * residual**2/expectation)
                    
                    prob_gst[ii] = prob_pos * prob_flux

                #else:
                #if f_verbose and fd_cat['id'][ii] == 490:
                #print(gst)
                #print('No counterpart is for %d found: r=%.2f>%.2f pixel'%(fd_cat['id'][ii], rtmp[iiy], rlim))

        if True:
            # No2; Catalog method
            fd = fits.open(infile)[1].data
            hd = fits.open(infile)[1].header        
            RA = hd['CRVAL1']
            DEC = hd['CRVAL2']

            print('Target coordinates are:',RA,DEC)
            
            # Astro query;
            Vizier.ROW_LIMIT = 1000
            # To cover pickoff mirror;
            wid_extra = 00 # This is from e-mail exchange with Kevin; to be confirmed after flight operation.
            width_detector = 2048 + wid_extra
            pixel = 0.06 # arcsec / pixel

            Vizier.ROW_LIMIT = 1000
            result = Vizier.query_region(coord.SkyCoord(ra=RA, dec=DEC,
                                                        unit=(u.deg, u.deg),
                                                        frame='icrs'),
                                        width='%.1fs'%(width_detector*pixel), 
                                        catalog=["I/345/gaia2"], 
                                        column_filters={'Gmag': '<20'})
                                        #catalog=["GSC"]) # GSC has some issues..

            # Get ghost positions;
            # Input magnitude;
            GSmag = result[0]['Gmag']
            magzp = 25.0
            flux_GS = 10**(-(GSmag-magzp)/(2.5))

            RA_key = 'RA_ICRS'
            DE_key = 'DE_ICRS'
            #RA_key = 'RAICRS'
            #DE_key = 'DEICRS'

            if False:
                # Save;
                result[0].write('2mass.dat', format='ascii')
            
            # Get pixel position of 2mass sources;
            #w = WCS(hd,naxis=2)
            #x, y = w.wcs_world2pix(result[0][RA_key],result[0][DE_key],0)
            # This needs to be distortion-uncorrected format;
            from jwst import datamodels
            im = datamodels.open(infile)
            x, y = im.meta.wcs.backward_transform(result[0][RA_key],result[0][DE_key])
            id_pub = np.arange(0,len(x),1)

            if True:
                # Add manually added sources to the query result;
                ra_ext = np.asarray([82.74528081, 82.74549915])
                de_ext = np.asarray([-68.79706061, -68.79157048])
                x_ext, y_ext = im.meta.wcs.backward_transform(ra_ext, de_ext)
                plt.scatter(x_ext, y_ext, marker='o', s=10, edgecolor='r', facecolor='none', label='')
                flux_GS_ext = [100,100]
                #gst_ext = get_ghost(x_ext, y_ext, flux=flux_GS_ext, filt=pupil, xshift=xshift, yshift=yshift)
                x = np.append(x,x_ext)
                y = np.append(y,y_ext)
                flux_GS = np.append(flux_GS, flux_GS_ext)

            gst_cat = get_ghost(x, y, flux=flux_GS, filt=pupil, xshift=xshift, yshift=yshift)
            #plt.scatter(gst_cat[0], gst_cat[1], marker='o', s=30, edgecolor='lightgreen', color='', label='GSC sources')
            
            for ii in range(len(flag_gst[:])):
                if flag_gst[ii] != 1: # If not ghost flagged above.
                    # This is because photutils is 0-based, while ds9 is not.
                    shift = 1.0 
                    #if fd_cat['source_sum'][ii]>100:
                    #    fw.write('%.3f %.3f\n'%(gst[0]+shift, gst[1]+shift))

                    # Check the position:
                    rtmp = np.sqrt( (fd_cat['xcentroid'].value[ii] - gst_cat[0])**2 + (fd_cat['ycentroid'].value[ii] - gst_cat[1])**2 )
                    iiy = np.argmin(rtmp)

                    if rtmp[iiy] < rlim: # and (fd_cat['source_sum'][ii] - fd_cat['source_sum'][iiy])>0 :
                        id_gst[ii] = fd_cat['id'][ii]
                        flag_gst[ii] = 1 
                        prob_pos = np.exp(-0.5 * rtmp[iiy]**2)
                        residual = np.abs(fd_cat['source_sum'][ii]*frac_ghost - flux_GS[iiy])
                        expectation = fd_cat['source_sum'][ii]*frac_ghost
                        prob_flux = np.exp(-0.5 * residual**2/expectation)
                        
                        prob_gst[ii] = prob_pos * prob_flux


        # Close file
        fw.close()

        # Save result;
        if True:
            suff = infile.split('/')[-1].replace('.fits','')

            fd_sci = fits.open(infile)[1].data
            plt.imshow(fd_sci, vmin=0, vmax=1)

            # All sources;
            plt.scatter(fd_cat['xcentroid'].value, fd_cat['ycentroid'].value, marker='o', s=30, edgecolor='cyan', color='', label='i2d sources')
            # Public;
            if False:
                plt.scatter(x, y, marker='o', s=30, edgecolor='lightgreen', color='', label='GSC sources')
            
            fw_cat = open('%s/ghost_detected_cat_%s.txt'%(DIR_OUT,suff),'w')
            fw_cat.write('# id raghs decghs xghs yghs\n')

            fw = open('%s/ghost_detected_%s.reg'%(DIR_OUT,suff),'w')
            fw.write('# xghs yghs\n')

            f_label = True
            for ii in range(len(flag_gst)):
                if flag_gst[ii] == 1:
                    iix = np.where(fd_cat['id']==id_gst[ii])
                    yghs = fd_cat['ycentroid'].value[iix]
                    xghs = fd_cat['xcentroid'].value[iix]

                    if f_label:
                        label = 'IDed ghost'
                        f_label = False
                    else:
                        label = ''
                    plt.scatter(xghs, yghs, marker='s', s=30, edgecolor='r', color='', label=label)
                    shift = 1.0 # This is because photutils is 0-based, while ds9 is not.
                    fw.write('%.3f %.3f\n'%(xghs+shift, yghs+shift))
                    fw_cat.write('%d %.7f %.7f %.7f %.7f\n'\
                                %(fd_cat['id'][iix], \
                                fd_cat['sky_centroid'][iix].ra.value, fd_cat['sky_centroid'][iix].dec.value, xghs, yghs))

            # Plot GAP:
            gaps = get_gap(pupil)
            plt.scatter(gaps[0], gaps[1], marker='x', s=30, color='orange', label='GAP')#, color='')
            
            fw_cat.close()
            fw.close()
            plt.legend(loc=0)
            plt.savefig('%s/results_%s.png'%(DIR_OUT,suff), dpi=300)
            plt.close()
        

        if f_tweak_imaege2:
            con = (flag_gst==1)
            print(id_gst[con])
            segfile = infile.replace('.fits', '_seg.fits')
            outfile = infile.replace('.fits', '_gst.fits')
            tweak_imaege2(id_gst[con], infile, segfile, outfile=outfile)


