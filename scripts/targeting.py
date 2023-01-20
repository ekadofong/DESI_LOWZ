import desispec
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd 
from astropy.table import Table, vstack
from astropy.io import fits
import scipy.constants as spc
from astropy.io import fits

def build_targets (filename="LOWZ_TARGETS_FUJI_HEALPIX.CSV"):
    target_cat = pd.read_csv(filename, index_col=0).set_index('TARGETID')

    # \\ build catalogs
    specfit_fuji = '/global/cfs/cdirs/desi/public/edr/vac/fastspecfit/fuji/v1.0/catalogs/fastspec-fuji.fits'
    photo_fuji = '/global/cfs/cdirs/desi/public/edr/vac/fastspecfit/fuji/v1.0/catalogs/fastphot-fuji.fits'
    
    with fits.open(specfit_fuji) as hdu:
        specfit_fuji1 = Table(hdu[1].data)
        specfit_fuji2 = Table(hdu[2].data).to_pandas()
        
    with fits.open(photo_fuji) as hdu:
        photo_fuji1 = Table(hdu[1].data)

    del specfit_fuji1['CONTINUUM_COEFF']
    del photo_fuji1['CONTINUUM_COEFF']

    specfit_fuji1 = specfit_fuji1.to_pandas()
    photo_fuji1 = photo_fuji1.to_pandas()    
    
    # \\ remove duplicates
    sp = specfit_fuji1.set_index('TARGETID')
    ph = photo_fuji1.set_index("TARGETID")

    rankings = np.zeros(sp.shape[0], dtype=int)
    rankings[sp['PROGRAM']=='dark'] = 4
    rankings[sp['PROGRAM']=='other'] = 3
    rankings[sp['PROGRAM']=='backup'] = 2
    rankings[sp['PROGRAM']=='bright'] = 1

    sp['ranking'] = rankings
    ph['ranking'] = rankings
    sp = sp.sort_values('ranking', ascending=False)
    ph = ph.sort_values('ranking', ascending=False)

    keep = ~sp.index.duplicated(keep='first')
    assert np.isclose(sp.index,ph.index).all()
    sp = sp.loc[keep]
    ph = ph.loc[keep]
    assert not sp.index.duplicated().any()
    assert not ph.index.duplicated().any()    
    
    lowz_sp = sp.reindex(target_cat.index)
    lowz_ph = ph.reindex(target_cat.index)
    return lowz_sp, lowz_ph