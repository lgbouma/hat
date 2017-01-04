from astrobase import hatlc

LC_path = '../data/LCs/'
field = 'G199'+'/'

this_sqlite = 'HAT-199-0087236-V0-DR0-hatlc.sqlite.gz'

lcd, msg = hatlc.read_and_filter_sqlitecurve(LC_path+field+this_sqlite)
