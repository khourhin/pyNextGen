import sqlite3
import glob
import os
from mylog import get_logger

logger = get_logger(__file__, __name__)


pups_fs = sorted(glob.glob('/home/ekornobis/analysis/batsche/methyleu/methylome/test/*.pup'))

logger.info(pups_fs)

# Initiate the db
conn = sqlite3.connect('demo_data/example.db')
c = conn.cursor()

for pup in pups_fs:
    # Feed the db
    table_id = (os.path.basename(pup))
    table_id = os.path.splitext(table_id)[0]
    logger.info(table_id)
    c.execute('CREATE TABLE {0} (chr text, position int, base text, nreads int)'.format(table_id))


    
    with open(pup) as f:
        for line in f:
            logger.info(line)
            c.execute('INSERT INTO {0} VALUES ("{1}", {2}, "{3}", {4})'.format(table_id,
                                                                               line.split()[0],
                                                                               line.split()[1],
                                                                               line.split()[2],
                                                                               line.split()[3]))
    conn.commit()
    break

# Query the db
c.execute('SELECT * FROM {0} WHERE position > 60908 AND position < 61000'.format('SRR2028035_chr11_filt_H_sub'))
print(c.fetchall())

conn.close()
