import numpy as np
import argparse
import os
import astropy.wcs as wcs
from astropy.io import fits,ascii


def run_photutils(image_cal, file_out=None, nsigma=1.5, npixels=20):
    '''
    image_cal : str
        file name for _cal.fits file.
    file_out : str
        name for output catalog.
    '''
    from photutils import Background2D, MedianBackground, detect_sources, deblend_sources, source_properties
    from astropy.stats import gaussian_fwhm_to_sigma
    from astropy.convolution import Gaussian2DKernel
    if file_out == None:
        file_out = image_cal.replace('_cal.fits', '_cat_man.ecsv')
        if file_out == image_cal:
            print('Input file is not `_cal.fits` file. No photutils.')
            return None

    hdu = fits.open(image_cal)
    data = hdu[1].data
    imwcs = wcs.WCS(hdu[1].header, hdu)
    err = hdu[2].data

    # Measure background and set detection threshold
    bkg_estimator = MedianBackground()
    bkg = Background2D(data, (50, 50), filter_size=(3, 3), bkg_estimator=bkg_estimator)
    #threshold = bkg.background + (15. * bkg.background_rms)
    threshold = bkg.background + (10. * bkg.background_rms)

    # Before detection, smooth image with Gaussian FWHM = n pixels
    sigma = nsigma * gaussian_fwhm_to_sigma  
    kernel = Gaussian2DKernel(sigma)#, x_size=5, y_size=5)
    kernel.normalize()

    # Detect and deblend
    segm_detect = detect_sources(data, threshold, npixels=npixels)#, filter_kernel=kernel)

    # Save segmentation map of detected objects
    segm_hdu = fits.PrimaryHDU(segm_detect.data.astype(np.uint32), header=imwcs.to_header())
    segm_hdu.writeto(image_cal.replace('_cal.fits','_cal_seg.fits'), overwrite=True)
    
    # Save cat;
    cat = source_properties(data-bkg.background, segm_detect, wcs=imwcs, background=bkg.background, error=err)
    tbl = cat.to_table()#columns=columns)

    tbl.write(file_out, format='ascii.ecsv',overwrite=True)
    return file_out


def check_keyword(file_cat, keywords, keys_str=None):
    '''
    file_cat : str
        Ascii catalog to be used.
    keywords : list
        list of strings for keywords to be checked.
    '''
    flag = True
    fd_cat = ascii.read(file_cat)
    for kk,key in enumerate(keywords):
        try:
            value_tmp = fd_cat[key]
        except:
            print('\n!!! Warning !!!\n`%s` column is not found in the input catalog.'%key)
            flag = False
            if not keys_str == None:
                print('Specify the column name by adding --%s argument.'%keys_str[kk])
                print('e.g.,\n python detect_ghost_image2.py image catalog --%s column-that-exists\n'%keys_str[kk])
    return flag


def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def get_gap_wfss(pupil, filt, file_gap=None):
    '''
    Purpose
    -------
    To get GAP coordinates corresponding to input filter.
    
    Parameters
    ----------
    pupil : str
        element in the pupil wheel.
    filt : str
        element in the filter wheel.

   
    Returns
    -------
    x, y, and fraction of ghost flux over the source flux.
    
    Note
    ----
    GAP is based on CV3 data.
    We currently do not know fractional flux, tab_gap['frac_50'], i.e. there may be positional dependence too.
    '''
    if file_gap == None:
        this_path = os.path.realpath(__file__)
        file_gap = '%s/niriss_ghost_gap_summary.txt'%(this_path.replace('utils.py',''))
        print('Using gap summary file: %s'%file_gap)

    tab_gap = ascii.read(file_gap)
    iix0 = np.where(tab_gap['pupil']==pupil.upper())
    iix1 = np.where(tab_gap['filt'][iix0]==filt.upper())

    if len(iix0[0])>0 and len(iix1[0])>0:
        xgap, ygap = tab_gap['gapx_50'][iix0][iix1], tab_gap['gapy_50'][iix0][iix1]
        frac = tab_gap['frac_50'][iix0][iix1]
        #print(xgap, ygap)
    else:
        print('%s-%s combination is not found in the table.'%(pupil,filt))
        return np.nan

    return xgap, ygap, frac


def get_ghost_wfss(x, y, flux=None, pupil='F200W',filt='CLEAR', shift=0, xshift=0, yshift=0, gap_tmp=[None,None,None], file_gap=None):
    '''
    Purpose
    -------
    A function that gives expected ghost positions
    given positions of a source.
    
    Parameters
    ----------
    x,y : arrays
        coordinates for sources (arrays).
    flux : array
        fluxes for sources (array)
    gap_tmp : 
        Manual input for coordinates of GAP and flux fraction. (x,y,f_flux)
    xshift : float
        Shift may be used, because photutils is 0-based, while ds9 is not.
    yshift : float
        Shift may be used, because photutils is 0-based, while ds9 is not.

    Returns
    -------
    xgs, ygs :
        coordinates for ghosts
    flux_gs :
        fluxes for ghosts    
    x,y,flux :
        input source coordinates and fluxes.
    '''
    if file_gap == None:
        this_path = os.path.realpath(__file__)
        file_gap = '%s/niriss_ghost_gap_summary.txt'%(this_path.replace('utils.py',''))
        #print('Using gap summary file: %s'%file_gap)
    
    if gap_tmp[0] == None:
        xgap, ygap, frac = get_gap_wfss(pupil, filt, file_gap=file_gap)
    else:
        xgap, ygap, frac = gap_tmp
    
    xgap += shift
    ygap += shift
    
    xgs = 2 * xgap - (x+xshift)
    ygs = 2 * ygap - (y+yshift)
    
    xgs -= xshift
    ygs -= yshift

    xgs.name = 'x_ghost'
    ygs.name = 'y_ghost'
    
    try:
        flux_gs = flux * frac/100.
    except:
        flux_gs = np.zeros(len(xgap),'float')
        
    return xgs,ygs,flux_gs,x,y,flux


def get_corners(x0, y0, filt, C, ORD='0', delx=40, dely=40, f_ds9=False, oversample_factor=1, f_verbose=False, nele=100):
    '''
    Purpose is to get extraction box, based on the potision of the source in the direct image.
    
    Parameters
    ----------
    x0, y0 : float
        position of source in direct image
    f_ds9 : bool
        is input coordinate retrieved directly from DS9?
    filt : str
        Target Filter.
        
    Returns
    -------
    x0s,x1s,y0s,y1s : float
        coordinates of the trace at wmin (0) and wmax (1)
    ''' 
    # This should match filts;
    lambdas = [0.90, 1.15, 1.405, 1.498, 1.587, 1.984]
    filts = ['F090W','F115W','F140M','F150W','F158M','F200W']
    try:
        iiw = np.where(np.asarray(filts) == filt.upper())[0][0]
    except:
        print('%s not found in '%(filt.upper()),filts)
        print('Returning None.')
        return None

    minx0 = x0 - delx/2.
    miny0 = y0 - dely/2.

    maxx0 = x0 + delx/2.
    maxy0 = y0 + dely/2.

    order = ORD
    #s = C.SENS[order]

    # Figuring out a few things about size of order, dispersion and wavelengths to use
    wmin = C.WRANGE[order][0]
    wmax = C.WRANGE[order][1]

    t0 = C.INVDISPL(order,x0,y0,wmin)
    t1 = C.INVDISPL(order,x0,y0,wmax)    
    
    ts = np.linspace(t0,t1,nele)
    x0s = x0 + C.DISPX(order,x0,y0,t0)
    x1s = x0 + C.DISPX(order,x0,y0,t1)
    y0s = y0 + C.DISPY(order,x0,y0,t0)
    y1s = y0 + C.DISPY(order,x0,y0,t1)
    xss = x0 + C.DISPX(order,x0,y0,ts)
    yss = y0 + C.DISPY(order,x0,y0,ts)
    
    #dx0 = C.DISPX(order,x0,y0,t0) - C.DISPX(order,x0,y0,t1)
    #dx1 = C.DISPY(order,x0,y0,t0) - C.DISPY(order,x0,y0,t1)

    # Get central point of the trace;
    dS = C.INVDISPL(order,x0,y0,lambdas[iiw])
    dXs = C.DISPX(order,x0,y0,dS)
    dYs = C.DISPY(order,x0,y0,dS)
    if f_verbose:
        print(dXs,dYs)

    x0s_diff = x0 + dXs 
    y0s_diff = y0 + dYs 
    #return x0s,x1s,y0s,y1s 
    return xss,yss,x0s_diff,y0s_diff


def rotate(origin, point, angle):
    '''
    For CV3 DMS conversion.
    Rotate a point counterclockwise by a given angle around a given origin.

    The angle should be given in radians.
    '''
    ox, oy = origin
    px, py = point

    qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
    qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
    return qx, qy


def flip(origin, point, direction='y'):
    '''
    For CV3 DMS conversion.
    '''
    ox, oy = origin
    px, py = point


    if direction == 'y':
        qy = oy + (oy - py)
        qx = px
    else:
        qx = ox + (ox - px)
        qy = py

    return qx, qy
    

def get_gap(filt, pupil='CLEAR', file_gap=None):
    '''
    Purpose
    -------
    To get GAP coordinates corresponding to input filter.
    
    Input
    -----
    filt : str
        filter of interest.
   
    Returns
    -------
    x, y, and fraction of ghost flux over the source flux.
    
    Notes
    -----
    GAP is based on CV3 data.
    We currently do not know fractional flux, tab_gap['frac_50'], i.e. there may be positional dependence too.
    input pupil and filter keywords are reversed. Sorry for confusion...

    '''
    if file_gap == None:
        this_path = os.path.realpath(__file__)
        file_gap = '%s/niriss_ghost_gap_summary.txt'%(this_path.replace('utils.py',''))
        print('Using gap summary file: %s'%file_gap)

    tab_gap = ascii.read(file_gap)
    iix = np.where( (tab_gap['pupil']==filt.upper()) & (tab_gap['filt']==pupil.upper()))

    if len(iix[0])>0:
        xgap, ygap = tab_gap['gapx_50'][iix], tab_gap['gapy_50'][iix]
        frac = tab_gap['frac_50'][iix]
    else:
        print('%s is not found in the table.'%filt)
        return np.nan

    return xgap, ygap, frac


def get_ghost(x, y, flux=None, filt='F200W', shift=0, xshift=0, yshift=0, gap_tmp=[None,None,None], file_gap=None):
    '''
    Purpose
    -------
    A function that gives expected ghost positions
    given positions of a source.
    
    Parameters
    ----------
    x,y : arrays
        coordinates for sources (arrays).
    flux : array
        fluxes for sources (array)
    gap_tmp : 
        Manual input for coordinates of GAP and flux fraction. (x,y,f_flux)
    xshift : float
        Shift may be used, because photutils is 0-based, while ds9 is not.
    yshift : float
        Shift may be used, because photutils is 0-based, while ds9 is not.

    Returns
    -------
    xgs, ygs :
        coordinates for ghosts
    flux_gs :
        fluxes for ghosts    
    x,y,flux :
        input source coordinates and fluxes.
    '''
    if file_gap == None:
        this_path = os.path.realpath(__file__)
        file_gap = '%s/niriss_ghost_gap_summary.txt'%(this_path.replace('utils.py',''))
        #print('Using gap summary file: %s'%file_gap)
    
    if gap_tmp[0] == None:
        xgap, ygap, frac = get_gap(filt, file_gap=file_gap)
    else:
        xgap, ygap, frac = gap_tmp
    
    xgap += shift
    ygap += shift
    
    xgs = 2 * xgap - (x+xshift)
    ygs = 2 * ygap - (y+yshift)
    
    xgs -= xshift
    ygs -= yshift

    xgs.name = 'x_ghost'
    ygs.name = 'y_ghost'

    try:
        flux_gs = flux * frac/100.
    except:
        flux_gs = np.zeros(len(xgap),'float')
        
    return xgs,ygs,flux_gs,x,y,flux


def tweak_dq(id_gst, infile, segfile, outfile=None, DQ_SET=1, DQ_KEY='DQ'):
    '''
    Purpose
    -------
    To make a copy of input fits file and tweak its DQ array.
    '''
    if not outfile == None:
        outfile = infile.replace('.fits','_gst.fits')

    try:
        dq_array = fits.open(infile)[DQ_KEY]
    except:
        print('DQ array, `%s`, not found. No DQ tweaking.'%DQ_KEY)
        return False

    os.system('cp %s %s'%(infile, outfile))
    
    seg = fits.open(segfile)[0].data
    with fits.open(outfile, mode='update') as hdul:
        for id in id_gst:
            con = np.where(seg == id)
            # 1073741824
            #hdul[DQ_KEY].data[con] += DQ_SET
            hdul[DQ_KEY].data[con] = DQ_SET
        
        hdul.flush()  # changes are written back to original.fits

    return True