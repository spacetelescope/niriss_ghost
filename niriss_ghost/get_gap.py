import glob
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys
import emcee
from scipy.optimize import minimize

from astropy.io import fits,ascii
import astropy.coordinates as coord
import astropy.units as u
from astropy.wcs import WCS
from astroquery.vizier import Vizier
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.convolution import Gaussian2DKernel
from photutils import Background2D, MedianBackground, detect_sources, deblend_sources, source_properties

from niriss_ghost.utils import get_gap,get_ghost,str2bool


def get_rms_img(file_image, segfile_out=None, f_detect=False):
    '''
    '''
    data = fits.open(file_image)['SCI'].data

    # Measure background and set detection threshold
    bkg_estimator = MedianBackground()
    bkg = Background2D(data, (50, 50), filter_size=(3, 3), bkg_estimator=bkg_estimator)
    bkg_rms = bkg.background_rms
    bkg_rms_med = bkg.background_rms_median

    if f_detect:
        threshold = bkg.background + (10. * bkg.background_rms)
        # Before detection, smooth image with Gaussian FWHM = n pixels
        sigma = 1.5 * gaussian_fwhm_to_sigma  

        kernel = Gaussian2DKernel(sigma)#, x_size=5, y_size=5)
        kernel.normalize()

        # Detect and deblend
        segm_detect = detect_sources(data, threshold, npixels=20)#, filter_kernel=kernel)

        # Save segmentation map of detected objects
        segm_hdu = fits.PrimaryHDU(segm_detect.data.astype(np.uint32), header=imwcs.to_header())

        if segfile_out == None:
            segfile_out = file_image.split('/')[-1].replace('.fits','_seg.fits')
        segm_hdu.writeto(segfile_out, overwrite=True)

    return bkg_rms_med

def log_prior(theta):
    '''
    Currently, no prior is set.
    '''
    x, y, f = theta
    #if f>-2.2 and f<-1.8: # Prior, or range?
    if True: # Prior, or range?
        return 0.0
    return -np.inf


def get_lnlike(gap_tmp, fd_cat, pupil, xshift, yshift, rlim, check_flux=False):
    '''
    Purpose
    -------
    Get log likelihood to evaluate current set of params, gap_tmp.
    
    Parameters
    ----------
    rlim : float
        Ghost will be identified if source is within this radius, in pixel.
    '''

    # Check if params are within range;
    lp = log_prior(gap_tmp)
    if not np.isfinite(lp):
        return -np.inf

    # Empty Array for ghost
    flag_gst = np.zeros(len(fd_cat['id']), 'int')
    prob_gst = np.zeros(len(fd_cat['id']), 'float')
    id_src = np.zeros(len(fd_cat['id']), 'int')
    
    # Get ghost position with the root method;
    gst = get_ghost(fd_cat['xcentroid'][:].value, fd_cat['ycentroid'][:].value, \
                    flux=fd_cat['source_sum'][:], filt=pupil,\
                    xshift=xshift, yshift=yshift, gap_tmp=gap_tmp)

    # Convert to linear scale;
    frac_ghost = 10**gap_tmp[2]

    # Check each obj;
    for ii in range(len(fd_cat['id'])):
        # Check the position:
        rtmp = np.sqrt( (fd_cat['xcentroid'].value-gst[0][ii])**2 + (fd_cat['ycentroid'].value-gst[1][ii])**2 )
        iiy = np.argmin(rtmp)

        # Check if source at the predicted positions is brighter than object of [ii].
        if rtmp[iiy] < rlim and (fd_cat['source_sum'][iiy]-fd_cat['source_sum'][ii]) > 0:

            id_src[ii] = fd_cat['id'][iiy]
            flag_gst[ii] = 1
            
            # Likelihood at the position;
            prob_pos = np.exp(-0.5 * rtmp[iiy]**2)
            prob_gst[ii] = prob_pos

            # Likelihood at the flux ratio.
            if check_flux:
                residual = np.abs(fd_cat['source_sum'][iiy]*frac_ghost - fd_cat['source_sum'][ii])
                expectation = fd_cat['source_sum'][iiy]*frac_ghost
                prob_flux = np.exp(-0.5 * residual**2 / expectation)
                prob_gst[ii] *= prob_flux

    #lnlike = np.exp(-np.abs(len(fd_cat['id'])-flag_gst.sum()))
    lnlike = np.log(np.sum(prob_gst))

    return lp + lnlike


def calc_gap(file_image, file_catalog, xshift=0, yshift=0, rlim=10, f_mirage=True, check_flux=False, nmc=3000, nwalkers=20,
    xini=1162.600, yini=937.900, logfini=-2):
    '''
    Parameters
    ----------
    check_flux : bool
        If include flux into log likeligood calculation. This depends on the catalog flux quality.
    xini, yini, logfini : float
        Initial guess for [x,y] position of GAP (in pixel) and flux ratio (in log).
    '''

    import emcee
    import corner
    fd_cat = ascii.read(file_catalog)

    # Read image header;
    hd = fits.open(file_image)[0].header
    hd1 = fits.open(file_image)[1].header
    filt = hd['FILTER']
    pupil = hd['PUPIL']
    XOFFSET = hd['XOFFSET']
    YOFFSET = hd['YOFFSET']
    try:
        # _cal.fits
        CDELT1 = np.abs(hd1['CD1_1'])
    except:
        CDELT1 = hd1['CDELT1']
        print('CAUTION : Your input seems to be IMAGE3 products.')

    if f_mirage:
        xshift = (2048-hd1['NAXIS1'])/2.
        yshift = (2048-hd1['NAXIS2'])/2.

    # Initial params;
    gap_true = [xini, yini, logfini] # x,y,log(f_flux)
    gap_tmp = gap_true
    lnlike = get_lnlike(gap_tmp, fd_cat, pupil, xshift, yshift, rlim)

    # Minimization;
    nll = lambda *args: -get_lnlike(*args)

    initial = np.array(gap_tmp)
    # Random fluctuation in xy
    initial[:2] +=  0.1 * np.random.randn(len(gap_tmp[:2]))
    # Random fluctuation in logf
    initial[2] +=  1e-2 * np.random.randn(1)

    soln = minimize(nll, initial, args=(fd_cat, pupil, xshift, yshift, rlim))
    m_ml, b_ml, log_f_ml = soln.x
    print("parameters at the minimization are;",soln.x)

    pos = soln.x + 1e-1 * np.random.randn(nwalkers, len(gap_tmp))
    ndim = pos.shape[1]

    sampler = emcee.EnsembleSampler(nwalkers, ndim, get_lnlike, args=(fd_cat, pupil, xshift, yshift, rlim), kwargs={'check_flux': check_flux})
    sampler.run_mcmc(pos, nmc, progress=True)

    # Make chain flat;
    flat_samples = sampler.get_chain(discard=int(nmc/2), thin=1, flat=True)

    # Corner;
    labels = ["$x_\mathrm{gap}$", "$y_\mathrm{gap}$", "$\log f$"]
    fig = corner.corner(
        flat_samples, labels=labels, truths=gap_true
    )
    plt.savefig('corner.png', dpi=100)

    # Get params
    for i in range(ndim):
        mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
        q = np.diff(mcmc)
        txt = "\mathrm{{{3}}} = {0:.3f}_{{-{1:.3f}}}^{{{2:.3f}}}"
        txt = txt.format(mcmc[1], q[0], q[1], labels[i])
        print(txt)


if __name__ == "__main__":
    '''
    Purpose
    -------
    Estimate GAP from image3 source catalog.

    Returns
    -------

    '''
    parser = argparse.ArgumentParser(description='Run the NIRISS gap analyzer.')
    parser.add_argument('file_image', metavar='file_image', type=str, nargs=1, help='Input image to be assessed.')
    parser.add_argument('file_catalog', metavar='file_catalog', type=str, nargs=1, help='The source catalog derived with input_image.')
    parser.add_argument('--nmc',default=3000,help='Number of iterations.', type=int)
    parser.add_argument('--nwalkers',default=20,help='Number of walkers.', type=int)
    parser.add_argument('--check_flux',default=True,help='Is flux ratio included in posterior calculation?', type=str2bool)
    args = parser.parse_args()
    
    for ii, image in enumerate(args.file_image):
        calc_gap(args.file_image[ii], args.file_catalog[ii], check_flux=args.check_flux, nmc=args.nmc, nwalkers=args.nwalkers)
