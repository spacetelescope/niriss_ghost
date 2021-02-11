import numpy as np
import argparse
from astropy.io import fits,ascii

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def get_gap(filt, file_gap='/ifs/jwst/wit/niriss/tmorishita/GhostDetection/ghost_analysis/gap_summary.txt'):
    '''
    Purpose:
    ========
    To get GAP coordinates corresponding to input filter.
    
    Input:
    ======
    filt : filter of interest.
   
    Return:
    =======
    x, y, and fraction of ghost flux over the source flux.
    
    Note:
    =====
    GAP is based on CV3 data.
    We currently do not know fractional flux, tab_gap['frac_50'], i.e. there may be positional dependence too.
    '''

    tab_gap = ascii.read(file_gap)
    iix = np.where(tab_gap['filt']==filt.upper())

    if len(iix[0])>0:
        xgap, ygap = tab_gap['gapx_50'][iix], tab_gap['gapy_50'][iix]
        frac = tab_gap['frac_50'][iix]
    else:
        print('%s is not found in the table.'%filt)
        return np.nan

    return xgap, ygap, frac



def get_ghost(x, y, flux=None, filt='F200W', shift=0, xshift=0, yshift=0, file_gap='/ifs/jwst/wit/niriss/tmorishita/GhostDetection/ghost_analysis/gap_summary.txt'):
    '''
    Purpose:
    ========
    A function that gives expected ghost positions
    given positions of a source.
    
    Input:
    ======
    x,y : coordinates for sources (arrays).
    flux : fluxes for sources (array)

    *shift may be used, because photutils is 0-based, while ds9 is not.

    xshift :
    yshift :


    Return:
    =======
    xgs, ygs : coordinates for ghosts
    flux_gs : fluxes for ghosts
    
    x,y,flux : input source coordinates and fluxes.
    '''

    xgap, ygap, frac = get_gap(filt, file_gap=file_gap)
    
    xgap += shift
    ygap += shift
    
    xgs = 2 * xgap - (x+xshift)
    ygs = 2 * ygap - (y+yshift)
    
    xgs -= xshift
    ygs -= yshift
    
    try:
        flux_gs = flux * frac/100.
    except:
        flux_gs = np.zeros(len(xgap),'float')
        
    return xgs,ygs,flux_gs,x,y,flux


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

