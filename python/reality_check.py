import subprocess
import sqlite3
from optparse import OptionParser
from pandas import DataFrame
import matplotlib.pyplot as plt
import seaborn

import decrypt_lib as decrypt
import os
import re

# options
parser = OptionParser()
parser.add_option("-d", type="str", dest="database", help="path to database")
parser.add_option("-l", type="int", dest="sequence_size", help="sequence_size")
parser.add_option("-s", type="float", dest="scale_tree", help="scale tree branch length")
parser.add_option("-b", type="str", dest="bpp", help="path to bpp executable")
parser.add_option("-c", type="str", dest="bpp_ctl", help="path to bpp config file")

(options, args) = parser.parse_args()
print(options.sequence_size)
conn = sqlite3.connect(options.database)

cursor = conn.execute("SELECT id, genealogies from results")
d = {'lon': [], 'lat': [], 'p1': [], 'p2': []}

for row in cursor:
    print('ID: ', row[0])
    trees = row[1].split("\n\n")
    imap = row[2]
    d['lon'].append(row[3])
    d['lat'].append(row[4])
    sequences = decrypt.simulate_sequences(trees, options.sequence_size, options.scale_tree)
    bpp_output = decrypt.run(options.bpp, options.bpp_ctl, sequences, imap)
    p1, p2 = decrypt.extract_probabilities(bpp_output)
    d['p1'].append(p1)
    d['p2'].append(p2)

df = DataFrame(data=d)
df.to_csv('data.txt')
