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


if __name__ == "__main__":
    '''
    Purpose:
    ========
    This script is for detecting ghost in Image3 step products.
    There are two algorithms used to flag possible ghost images.

    1. Root method: Look into detected sources in Image3 catalog, predict positions of ghosts based on GAP, 
    and then confirm if the source is consistent with the prediction.;

    2. Catalog method: Retrieve bright sources from all sky catalogs (here GAIA, through astroquery),
    then apply the same method as for root method to predict and confirm any sources in image3 catalog as ghosts.

    Return:
    =======
    A subset of the input source catalog, with a ghost flag column, "is_this_ghost", added.
    For ghosts without original sources in the input catalog but in a catalog from astroquery, idsrc is set to >idpub.

    '''

    parser = argparse.ArgumentParser(description='Run the NIRISS ghost detection script.')
    parser.add_argument('input_image', metavar='input_image', type=str, nargs=1, help='Input image to be assessed.')
    parser.add_argument('input_catalog', metavar='input_catalog', type=str, nargs=1, help='The source catalog derived with input_image.')
    parser.add_argument('--f_verbose',default=False,help='Print comments.', type=str2bool)
    parser.add_argument('--rlim',default=10,help='Search radius for ghost around the predicted position (in pixel).', type=float)
    parser.add_argument('--frac_ghost',default=0.01,help='Flux fraction for ghosts.', type=float)
    parser.add_argument('--o',default='./',help='Output directory.', type=str)
    parser.add_argument('--f_tweak_imaege2',default=False,help='Tweak Image2 products.', type=str2bool)
    parser.add_argument('--f_mirage',default=True,help='Is dataset created by Mirage?', type=str2bool)
    args = parser.parse_args()
    
    # ghosts with idsrc greater than the following number are identified through the Catalog method.
    idpub = 100000

    # Flags for analysis method;
    f_rootmethod = True
    f_gsc = True

    f_verbose = args.f_verbose
    rlim = args.rlim
    frac_ghost = args.frac_ghost
    f_tweak_imaege2 = args.f_tweak_imaege2
    DIR_OUT = args.o    
    if not os.path.exists(DIR_OUT):
        os.mkdir(DIR_OUT)

    # Ghost detection must be done in distortion corrected frame
    DIR_DATA = './'

    # Image;
    infiles = args.input_image
    # Catalog;
    #files_cat = [infile.replace('_i2d.fits','_cat_man.ecsv') for infile in infiles]
    files_cat = args.input_catalog
    
    for ff,infile in enumerate(infiles):
        file_root = infile.split('/')[-1].replace('.fits','')

        # Read header;
        hd = fits.open(infile)[0].header
        hd1 = fits.open(infile)[1].header
        filt = hd['FILTER']
        pupil = hd['PUPIL']
        XOFFSET = hd['XOFFSET']
        YOFFSET = hd['YOFFSET']
        CDELT1 = hd1['CDELT1']
        
        if args.f_mirage:
            # Add shift for image size difference.
            # This is due to the fact that ghosts were added in the seed image dimention, whereas analysis is done in i2d image.
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
        id_src = np.zeros(len(fd_cat['id']), 'int')

        ##################
        # No1; Root method
        ##################
        if f_rootmethod:
            gst = get_ghost(fd_cat['xcentroid'][:].value, fd_cat['ycentroid'][:].value, \
                            flux=fd_cat['source_sum'][:], filt=pupil,\
                            xshift=xshift, yshift=yshift)

            for ii in range(len(fd_cat['id'])):

                # Check the position:
                rtmp = np.sqrt( (fd_cat['xcentroid'].value-gst[0][ii])**2 + (fd_cat['ycentroid'].value-gst[1][ii])**2 )
                iiy = np.argmin(rtmp)

                # Check if source at the predicted positions is brighter than object of [ii].
                if rtmp[iiy] < rlim and (fd_cat['source_sum'][iiy]-fd_cat['source_sum'][ii]) > 0:

                    id_src[ii] = fd_cat['id'][iiy]
                    flag_gst[ii] = 1
                    
                    #if f_verbose:
                    #    print('Found an object at r=%.2f pixel'%(rtmp[iiy]))
                    
                    prob_pos = np.exp(-0.5 * rtmp[iiy]**2)
                    residual = np.abs(fd_cat['source_sum'][ii]*frac_ghost - fd_cat['source_sum'][iiy])
                    expectation = fd_cat['source_sum'][ii]*frac_ghost
                    prob_flux = np.exp(-0.5 * residual**2/expectation)
                    
                    prob_gst[ii] = prob_pos * prob_flux

        #####################
        # No2; Catalog method
        #####################
        if f_gsc:
            from jwst import datamodels
            fd = fits.open(infile)[1].data
            hd = fits.open(infile)[1].header        
            RA = hd['CRVAL1']
            DEC = hd['CRVAL2']

            print('Target coordinates are:',RA,DEC)
            
            # Astro query;
            Vizier.ROW_LIMIT = 1000

            # To cover pickoff mirror;
            wid_extra = 200 # This is from e-mail exchange with Kevin; to be confirmed after flight operation.
            width_detector = 2048 + wid_extra
            Vizier.ROW_LIMIT = 1000

            result = Vizier.query_region(coord.SkyCoord(ra=RA, dec=DEC,
                                                        unit=(u.deg, u.deg),
                                                        frame='icrs'),
                                        width='%.1fs'%(width_detector*CDELT1*3600),
                                        catalog=["I/345/gaia2"], 
                                        column_filters={'Gmag': '<20'})
                                        #catalog=["GSC"]) # GSC has some issues..

            # Get ghost positions;
            # Input magnitude;
            mvega2ab = -0.08
            GSmag = result[0]['Gmag'] + mvega2ab
            magzp = 25.0
            flux_GS = 10**(-(GSmag-magzp)/(2.5)) # Fnu with magzp

            RA_key = 'RA_ICRS'
            DE_key = 'DE_ICRS'
            #RA_key = 'RAICRS'
            #DE_key = 'DEICRS'
            
            # Get pixel position of GAIA sources;
            im = datamodels.open(infile)
            x, y = im.meta.wcs.backward_transform(result[0][RA_key],result[0][DE_key])
            id_pub = np.arange(0,len(x),1)

            ra_src_pub = np.zeros(len(flag_gst[:]),'float')
            dec_src_pub = np.zeros(len(flag_gst[:]),'float')
            
            # Get ghost xy;
            gst_cat = get_ghost(x, y, flux=flux_GS, filt=pupil, xshift=xshift, yshift=yshift)

            # Cross match;
            for ii in range(len(flag_gst[:])):
                if flag_gst[ii] != 1: # If the source has not been flagged above yet.

                    # Check the position:
                    rtmp = np.sqrt( (fd_cat['xcentroid'].value[ii] - gst_cat[0])**2 + (fd_cat['ycentroid'].value[ii] - gst_cat[1])**2 )
                    iiy = np.argmin(rtmp)

                    if rtmp[iiy] < rlim and (flux_GS[iiy] - fd_cat['source_sum'][ii])>0 :
                        id_src[ii] = idpub + ii
                        flag_gst[ii] = 1
                        prob_pos = np.exp(-0.5 * rtmp[iiy]**2)
                        residual = np.abs(fd_cat['source_sum'][ii]*frac_ghost - flux_GS[iiy])
                        expectation = fd_cat['source_sum'][ii]*frac_ghost
                        prob_flux = np.exp(-0.5 * residual**2/expectation)
                        prob_gst[ii] = prob_pos * prob_flux

                        ra_src_pub[ii] = result[0][RA_key][iiy]
                        dec_src_pub[ii] = result[0][DE_key][iiy]


        # Save result;
        if True:
            fig = plt.figure(figsize=(6,4))
            fig.subplots_adjust(top=0.93, bottom=0.1, left=0.05, right=0.9, hspace=0.05, wspace=0.2)
            ax = plt.subplot(111)

            fd_sci = fits.open(infile)[1].data
            ax.imshow(fd_sci, vmin=0, vmax=1)

            # All sources;
            ax.scatter(fd_cat['xcentroid'].value, fd_cat['ycentroid'].value, marker='o', s=30, edgecolor='cyan', color='', label='i2d sources')

            # Public;
            if False:
                ax.scatter(x, y, marker='o', s=30, edgecolor='lightgreen', color='', label='GSC sources')
            
            fw_cat = open('%s/ghost_detected_cat_%s.txt'%(DIR_OUT, file_root),'w')
            fw_cat.write('# id ra dec x y is_this_ghost id_src ra_src dec_src\n')

            f_label = True
            for ii in range(len(flag_gst)):
                yghs = fd_cat['ycentroid'].value[ii]
                xghs = fd_cat['xcentroid'].value[ii]
                if flag_gst[ii] == 1:
                    iix = np.where(fd_cat['id']==id_src[ii])
                    if len(iix[0])>0:
                        ysrc = fd_cat['ycentroid'].value[iix]
                        xsrc = fd_cat['xcentroid'].value[iix]
                        rasrc = fd_cat['sky_centroid'].ra.value[iix]
                        decsrc = fd_cat['sky_centroid'].dec.value[iix]

                        if f_label:
                            label = 'IDed ghost'
                            f_label = False
                        else:
                            label = ''

                        ax.scatter(xghs, yghs, marker='s', s=30, edgecolor='r', color='', label=label)
                        shift = 1.0 # This is because photutils is 0-based, while ds9 is not.
                        fw_cat.write('%d %.7f %.7f %.7f %.7f %s %d %.7f %.7f\n'\
                                    %(fd_cat['id'][ii], \
                                    fd_cat['sky_centroid'][ii].ra.value, fd_cat['sky_centroid'][ii].dec.value, \
                                    xghs, yghs, 'True', fd_cat['id'][iix], rasrc, decsrc))

                    else: # Ghost from outside FoV;
                        fw_cat.write('%d %.7f %.7f %.7f %.7f %s %d %.7f %.7f\n'\
                                    %(fd_cat['id'][ii], \
                                    fd_cat['sky_centroid'][ii].ra.value, fd_cat['sky_centroid'][ii].dec.value, \
                                    xghs, yghs, 'True', id_src[ii], ra_src_pub[ii], dec_src_pub[ii]))

                else:
                    fw_cat.write('%d %.7f %.7f %.7f %.7f %s %d %.7f %.7f\n'\
                                %(fd_cat['id'][ii], \
                                fd_cat['sky_centroid'][ii].ra.value, fd_cat['sky_centroid'][ii].dec.value, \
                                xghs, yghs, 'False', -99, np.nan, np.nan))


            # Plot GAP:
            gaps = get_gap(pupil)
            ax.scatter(gaps[0], gaps[1], marker='x', s=30, color='orange', label='GAP')
            fw_cat.close()
            ax.legend(bbox_to_anchor=(1., 1.05))
            plt.savefig('%s/results_%s.png'%(DIR_OUT,file_root), dpi=300)
            plt.close()


        '''
        from utils import tweak_imaege2
        # This is for future implementation;
        if f_tweak_imaege2:
            con = (flag_gst==1)
            segfile = infile.replace('.fits', '_seg.fits')
            outfile = infile.replace('.fits', '_gst.fits')
            tweak_imaege2(id_gst[con], infile, segfile, outfile=outfile)
        '''