# microscript: gzip a file

import gzip
import shutil

with open( "/path/to/file", "rb" ) as f_in:
    with open( "/path/to/file.gz", "wb" ) as f_out:
        shutil.copyfileobj( f_in, f_out )
