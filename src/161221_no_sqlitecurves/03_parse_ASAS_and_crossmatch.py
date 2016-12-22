'''
Having generated a list of HAT IDs and their correspond 2MASS IDs, we want to
find objects in the field that are known EBs.
Here we do this by crossmatching with the ASAS3 Catalog of Variable Stars.
'''

from astropy.io import ascii
col_str = 'ID,PER,HJD0,VMAX,VAMP,TYPE,GCVS_ID,GCVS_TYPE,IR12,IR25,IR60,'+\
    'IR100,J,H,K,V_IR12,V_J,V_H,V_K,J_H,H_K'
col_names = col_str.split(',')

catalog_path = '../data/catalogs/'
tab = ascii.read(catalog_path+'ACVS.1.1',names=col_names)

# works, but I'm being silly b/c ASAS3 is in the south.
