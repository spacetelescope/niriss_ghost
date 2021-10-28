import glob
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys

from astropy.io import fits,ascii
import astropy.coordinates as coord
import astropy.units as u
from astropy.wcs import WCS
from astroquery.vizier import Vizier
from jwst import datamodels

# This repository;
from niriss_ghost.utils import get_gap,get_ghost,str2bool,tweak_dq,check_keyword


def run(infiles, files_cat, f_verbose=True, rlim=10, frac_ghost=0.01, f_tweak_dq=True, DIR_OUT='./output/',
    f_mirage=True, segmap=None, idarx=100000, keyword_id='label', keyword_flux='source_sum', 
    keyword_xcent='xcentroid', keyword_ycent='ycentroid', keyword_coord='sky_centroid', 
    f_save_result=True, out_cat=None, file_gap=None, DQ_KEY='DQ'):
    '''
    Parameters
    ----------
    infiles : list
        List of input fits image files.
    files_cat : list
        List of input catalog files. The number of the elements must be same as those in infiles.
    idarx : int
        ghosts with idsrc greater than this number are those identified through the Catalog method.
    keyword_id : str
        keyword for the object id column in files_cat. (can be `label`,`id`,`number`, depending on the version of photutils or sourcextractor)
    file_gap : str
        file name for gap summary file. If none, this will be niriss_ghost_gap_summary.txt
    '''
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\nRunning ghost detection script\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    

    # Flags for analysis method;
    f_rootmethod = True
    f_gsc = True

    if not os.path.exists(DIR_OUT):
        os.mkdir(DIR_OUT)
    
    for ff,infile in enumerate(infiles):
        file_root = infile.split('/')[-1].replace('.fits','')

        # Read header;
        hd = fits.open(infile)[0].header
        hd1 = fits.open(infile)[1].header
        filt = hd['FILTER']
        pupil = hd['PUPIL']
        XOFFSET = hd['XOFFSET']
        YOFFSET = hd['YOFFSET']
        try:
            # _cal.fits
            CDELT1 = np.abs(hd1['CD1_1'])
        except:
            CDELT1 = hd1['CDELT1']
            print('CAUTION : Your input seems to be i2d products.')
            print('Use results with causion. Pixel scale is currently set to a default value.\n')

        try:
            # Magnitude zeropoint:
            PIXAR_SR = hd1['PIXAR_SR'] # str / pixel
            PHOTMJSR = hd1['PHOTMJSR'] # Mjy / str
            fluxzp = PIXAR_SR * PHOTMJSR * 1e6
            magzp = -2.5 * np.log10(fluxzp) + 8.9
        except:
            magzp = 25.0
            print('Magzp cannot be calculated from header. Set to %.1f'%magzp)

        if f_mirage:
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

        # Check keyword:
        keys_used = [keyword_id,keyword_flux,keyword_xcent,keyword_ycent,keyword_coord]
        keys_str = ['keyword_id', 'keyword_flux', 'keyword_xcent', 'keyword_ycent', 'keyword_coord']
        key_result = check_keyword(file_cat,keys_used,keys_str=keys_str)
        if not key_result:
            print('Exiting.')
            sys.exit()

        # Empty Array for ghost
        flag_gst = np.zeros(len(fd_cat[keyword_id]), 'int')
        prob_gst = np.zeros(len(fd_cat[keyword_id]), 'float')
        id_src = np.zeros(len(fd_cat[keyword_id]), 'int')

        # Check flux column:
        flux_tmp = fd_cat[keyword_flux]
        flux_cat = fd_cat[keyword_flux]
        xcent = fd_cat[keyword_xcent]
        ycent = fd_cat[keyword_ycent]
        try:# Remove unit, if it has any.
            flux_cat = fd_cat[keyword_flux].value
            xcent = fd_cat[keyword_xcent].value
            ycent = fd_cat[keyword_ycent].value
        except:
            pass

        ##################
        # No1; Root method
        ##################
        if f_rootmethod:
            gst = get_ghost(xcent, ycent, \
                            flux=flux_cat[:], filt=pupil,\
                            xshift=xshift, yshift=yshift)

            for ii in range(len(fd_cat[keyword_id])):

                # Check the position:
                rtmp = np.sqrt( (xcent-gst[0][ii])**2 + (ycent-gst[1][ii])**2 )
                iiy = np.argmin(rtmp)

                # Check if source at the predicted positions is brighter than object of [ii].
                if rtmp[iiy] < rlim and (flux_cat[iiy]-flux_cat[ii]) > 0:

                    id_src[ii] = fd_cat[keyword_id][iiy]
                    flag_gst[ii] = 1
                    
                    #if f_verbose:
                    #    print('Found an object at r=%.2f pixel'%(rtmp[iiy]))
                    
                    prob_pos = np.exp(-0.5 * rtmp[iiy]**2)
                    residual = np.abs(flux_cat[ii]*frac_ghost - flux_cat[iiy])
                    expectation = flux_cat[ii]*frac_ghost
                    prob_flux = np.exp(-0.5 * residual**2/expectation)
                    
                    prob_gst[ii] = prob_pos * prob_flux

        #####################
        # No2; Catalog method
        #####################
        if f_gsc:
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

            # Retrieve source catalog in a box.
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
            flux_GS = 10**(-(GSmag-magzp)/(2.5)) # Fnu with the same magzp as input image.

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
                    rtmp = np.sqrt( (xcent[ii] - gst_cat[0])**2 + (ycent[ii] - gst_cat[1])**2 )
                    iiy = np.argmin(rtmp)

                    if rtmp[iiy] < rlim and (flux_GS[iiy] - flux_cat[ii]) > 0:
                        id_src[ii] = idarx + ii
                        flag_gst[ii] = 1
                        prob_pos = np.exp(-0.5 * rtmp[iiy]**2)
                        residual = np.abs(flux_cat[ii]*frac_ghost - flux_GS[iiy])
                        expectation = flux_cat[ii]*frac_ghost
                        prob_flux = np.exp(-0.5 * residual**2/expectation)
                        prob_gst[ii] = prob_pos * prob_flux

                        ra_src_pub[ii] = result[0][RA_key][iiy]
                        dec_src_pub[ii] = result[0][DE_key][iiy]


        # Save result;
        if f_save_result:
            fig = plt.figure(figsize=(4,4))
            fig.subplots_adjust(top=0.98, bottom=0.1, left=0.1, right=0.98, hspace=0.05, wspace=0.2)
            ax = plt.subplot(111)

            fd_sci = fits.open(infile)[1].data
            ax.imshow(fd_sci, vmin=0, vmax=1, origin='lower')

            # All sources;
            ax.scatter(xcent, ycent, marker='o', s=30, edgecolor='cyan', color='none', label='Catalog sources')

            # Source in retrieved catalog;
            if False:
                ax.scatter(x, y, marker='o', s=30, edgecolor='lightgreen', color='none', label='GSC sources')

            # Write to an ascii; 
            if out_cat == None:
                out_cat = '%s/ghost_detected_cat_%s.txt'%(DIR_OUT, file_root)
            fw_cat = open(out_cat,'w')
            fw_cat.write('# id ra dec x y is_this_ghost id_src ra_src dec_src\n')

            f_label = True
            for ii in range(len(flag_gst)):
                yghs = ycent[ii]
                xghs = xcent[ii]
                if flag_gst[ii] == 1:
                    iix = np.where(fd_cat[keyword_id]==id_src[ii])
                    if len(iix[0])>0:
                        ysrc = ycent[iix]
                        xsrc = xcent[iix]
                        rasrc = fd_cat[keyword_coord].ra.value[iix]
                        decsrc = fd_cat[keyword_coord].dec.value[iix]

                        if f_label:
                            label = 'IDed ghost'
                            f_label = False
                        else:
                            label = ''

                        ax.scatter(xghs, yghs, marker='s', s=30, edgecolor='r', color='none', label=label)
                        shift = 1.0 # This is because photutils is 0-based, while ds9 is not.
                        fw_cat.write('%d %.7f %.7f %.7f %.7f %s %d %.7f %.7f\n'\
                                    %(fd_cat[keyword_id][ii], \
                                    fd_cat[keyword_coord][ii].ra.value, fd_cat[keyword_coord][ii].dec.value, \
                                    xghs, yghs, 'True', fd_cat[keyword_id][iix], rasrc, decsrc))

                    else: # Ghost from outside FoV;
                        fw_cat.write('%d %.7f %.7f %.7f %.7f %s %d %.7f %.7f\n'\
                                    %(fd_cat[keyword_id][ii], \
                                    fd_cat[keyword_coord][ii].ra.value, fd_cat[keyword_coord][ii].dec.value, \
                                    xghs, yghs, 'True', id_src[ii], ra_src_pub[ii], dec_src_pub[ii]))

                else:
                    fw_cat.write('%d %.7f %.7f %.7f %.7f %s %f %f %f\n'\
                                %(fd_cat[keyword_id][ii], \
                                fd_cat[keyword_coord][ii].ra.value, fd_cat[keyword_coord][ii].dec.value, \
                                xghs, yghs, 'False', np.nan, np.nan, np.nan))

            fw_cat.close()
            print('Catalog saved to : %s'%(out_cat))

            # Plot GAP:
            outplot = '%s/results_%s.png'%(DIR_OUT,file_root)
            gaps = get_gap(pupil, file_gap=file_gap)
            ax.scatter(gaps[0], gaps[1], marker='x', s=30, color='orange', label='GAP (%s)'%pupil)
            ax.legend()#bbox_to_anchor=(1., 1.05))
            plt.savefig(outplot, dpi=300)
            plt.close()
            print('Plot saved to : %s\n'%outplot)

        # Tweak DQ array;
        if f_tweak_dq:
            print('Tweaking DQ array...')
            con = (flag_gst==1)
            if segmap != None:
                segfile = segmap
            else:
                print('No segmtation map specified by --segmap. Guessing from the input catalog name...\n')
                segfile = infile.replace('.fits', '_seg.fits')
            
            outfile = DIR_OUT + infile.split('/')[-1].replace('.fits', '_gst.fits')
            if not os.path.exists(segfile):
                print('Segmentation file (%s) is missing. No DQ tweaking.\n'%segfile)
            else:
                results_qd = tweak_dq(fd_cat[keyword_id][con], infile, segfile, outfile=outfile, DQ_SET=1, DQ_KEY=DQ_KEY)
                if results_qd:
                    print('New image with updated DQ saved to: %s\n'%outfile)
        
        print('Successfully done! (%d/%d)\n'%(ff+1,len(infiles)))


if __name__ == "__main__":
    '''
    Purpose
    -------
    This script is for detecting ghost in Image2 step products.
    There are two algorithms used to flag possible ghost images.

    1. Root method: Look into detected sources in Image2 catalog, predict positions of ghosts based on GAP, 
    and then confirm if the source is consistent with the prediction.;

    2. Catalog method: Retrieve bright sources from all sky catalogs (here GAIA, through astroquery),
    then apply the same method as for root method to predict and confirm any sources in image2 catalog as ghosts.


    Arguments
    ---------
    f_mirage : bool
        If input images are real data, turn this off. If images are from Mirage, turn this on. 
        This is due to the fact that ghosts were added in the seed image dimention, whereas analysis is done in i2d image.
    f_tweak_dq : bool
        Ghost detection in image2 products. Currently not supported.


    Notes
    -----
    Ghost detection must be done on the distortion corrected frame, as the GAP coordinates were calculated so.


    Returns
    -------
    A subset of the input source catalog, with a ghost flag column, "is_this_ghost", added.
    For ghosts without original sources in the input catalog but in a catalog from astroquery, idsrc is set to > idarx.

    '''

    parser = argparse.ArgumentParser(description='Run the NIRISS ghost detection script.')
    parser.add_argument('input_image', metavar='input_image', type=str, nargs=1, help='Input image to be assessed.')
    parser.add_argument('input_catalog', metavar='input_catalog', type=str, nargs=1, help='The source catalog derived with input_image.')
    parser.add_argument('--f_verbose',default=False,help='Print comments.', type=str2bool)
    parser.add_argument('--rlim',default=10,help='Search radius for ghost around the predicted position (in pixel).', type=float)
    parser.add_argument('--frac_ghost',default=0.01,help='Flux fraction for ghosts.', type=float)
    parser.add_argument('--o',default='./',help='Output directory.', type=str)
    parser.add_argument('--f_tweak_dq',default=True,help='Tweak DQ array of the input image.', type=str2bool)
    parser.add_argument('--f_mirage',default=True,help='Is input image created by Mirage?', type=str2bool)
    parser.add_argument('--keyword_flux',default='source_sum',help='Keyword for a flux column in input_catalog', type=str)
    parser.add_argument('--segmap',default=None,help='Segmentation map associated with input_catalog', type=str)
    parser.add_argument('--keyword_id',default='id',help='keyword for the object id column in files_cat. (can be `label`,`id`,`number`)', type=str)
    parser.add_argument('--keyword_xcent',default='xcentroid',help='keyword for object x-position.', type=str)
    parser.add_argument('--keyword_ycent',default='ycentroid',help='keyword for object y-position.', type=str)
    parser.add_argument('--keyword_coord',default='sky_centroid',help='keyword for object skyposition (ie radec).', type=str)
    parser.add_argument('--file_gap',default=None,help='file name for gap summary.', type=str)
    args = parser.parse_args()

    f_verbose = args.f_verbose
    rlim = args.rlim
    frac_ghost = args.frac_ghost
    f_tweak_dq = args.f_tweak_dq
    DIR_OUT = args.o    

    input_images = args.input_image[0].split(',')
    input_catalogs = args.input_catalog[0].split(',')
    run(input_images, input_catalogs, f_verbose=f_verbose, rlim=rlim, frac_ghost=frac_ghost, \
        f_tweak_dq=f_tweak_dq, DIR_OUT=DIR_OUT, f_mirage=args.f_mirage, keyword_flux=args.keyword_flux, \
        segmap=args.segmap, keyword_id=args.keyword_id, keyword_xcent=args.keyword_xcent, keyword_ycent=args.keyword_ycent, \
        keyword_coord=args.keyword_coord, file_gap=args.file_gap)