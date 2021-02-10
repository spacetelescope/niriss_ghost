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




def run_mirage(pfile, DIR_INP='./files/', f_plot=True):
    '''
    Purpose:
    ========
    Create NIRISS scene with ghost.

    Input:
    ======
    pfile : yaml file. No configuration for ghost is needed.
    f_plot : Bool, for plot.
    '''
    #from utils import get_gap,get_ghost

    ghost_image = '%s/ghost_cen.fits'%DIR_INP

    # 1.1 Run;
    print('Running catalog_seed_image for {}'.format(pfile))
    imseeds = []
    cat = catalog_seed_image.Catalog_seed(offline=False)
    cat.paramfile = pfile

    cat.make_seed()
    imseeds.append(cat.seed_file)

    # 1.2. Get output source catalog;

    # Taken From Mirage;
    source_type = 'pointsources'
    psfile = cat.params['Output']['file'][0:-5] + '_{}.list'.format(source_type)
    
    # This should be a list of real sources;
    tab_source = ascii.read(psfile, comment = '#', header_start=2) #, data_start=0)


    # _s represents source, as opposed to _gs for ghost;
    mag_s = tab_source['magnitude']

    # Let's focus on bright sources
    con = (mag_s<15)

    id_s = tab_source['Index'][con]
    xpix_s = tab_source['pixel_x'][con]
    ypix_s = tab_source['pixel_y'][con]
    ra_s = tab_source['RA_degrees'][con]
    dec_s = tab_source['DEC_degrees'][con]
    count_s = tab_source['counts/sec'][con]


    # Use median for GAP coordinate and flux ratio;
    filter = cat.params['Readout']['pupil'] #'F150W'

    # Coordinates and flux for ghosts:
    xgs,ygs,flux_gs,xin,yin,fluxin = get_ghost(xpix_s, ypix_s, \
                                            flux=count_s, filt=filter,\
                                            file_gap=file_gap)

    # Assign new ids for Ghost:
    # Here, id+10000.
    t = Table()
    t['idgs']= id_s + 100000
    t['xgs'] = xgs
    t['ygs'] = ygs
    t['counts/sec'] = flux_gs


    if f_plot:
        plt.style.use('dark_background')

        fig = plt.figure(figsize=(8,8))
        fig.subplots_adjust(top=0.94, bottom=0.12, left=0.14, right=0.98, hspace=0.12, wspace=0.2)

        # Source;
        plt.scatter(xpix_s,ypix_s,marker='o',c='w',alpha=0.5,s=count_s/1000)

        xdet = 2048
        ydet = 2048

        plt.xlim(0,xdet)
        plt.ylim(0,ydet)
        plt.title('Source')

        plt.savefig('ps.png', dpi=150)
