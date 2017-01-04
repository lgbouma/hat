'''
The rules for data which TOPCAT understands data are as follows:
* Each row must have the same number of comma-separated fields.
* Whitespace (space or tab) adjacent to a comma is ignored.
* Adjacent commas, or a comma at the start or end of a line (whitespace apart) 
    indicates a null field.
* Lines are terminated by any sequence of carriage-return or newline characters 
    ('\r' or '\n') (a corollary of this is that blank lines are ignored).
* Cells may be enclosed in double quotes; quoted values may contain linebreaks 
    (or any other character); a double quote character within a quoted value is 
    represented by two adjacent double quotes.
* The first line may be a header line containing column names rather than a row 
    of data. Exactly the same syntactic rules are followed for such a row as for 
    data rows.
* You cannot use a "#" character (or anything else) to introduce "comment" lines.
'''

def split_HAT_format_to_TOPCAT_readable(lc_paths):
    for lc in lc_paths:
        print lc
        file_contents = []
        for line in open(lc, 'r'):
            file_contents.append(line)
        comments = [l for l in file_contents if ('#' in l) or ('\n' == l)]
        assert file_contents[len(comments)-1] == '\n'

        lc_data = file_contents[len(comments):]

        col_names = [c[11:14] for c in comments if '# col' in c]
        header_str = ''
        for ix, c in enumerate(col_names):
            if ix != len(col_names)-1:
                header_str += c+', '
            if ix == len(col_names)-1:
                header_str += c+'\n'

        comment_fname = lc.split('/')[1].split('.')[0] + '.head'
        data_fname = lc.split('/')[1]
        with open('headers_to_data/'+comment_fname, 'w') as f:
            f.writelines(comments)
            f.close()
        lc_data.insert(0, header_str)
        with open('data/'+data_fname, 'w') as f:
            f.writelines(lc_data)
            f.close()

import numpy as np, pandas as pd, matplotlib.pyplot as plt
import os

if __name__ == '__main__':
    lc_dir = 'light_curves/'
    lc_names = np.sort(os.listdir(lc_dir))
    lc_paths = [lc_dir+name for name in lc_names]

    split_HAT_format_to_TOPCAT_readable(lc_paths)
