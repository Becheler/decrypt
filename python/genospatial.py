import subprocess
import sqlite3
from optparse import OptionParser
from pandas import DataFrame
import matplotlib.pyplot as plt
import seaborn

import genospatial_lib as gensp
import os
import re

# options
parser = OptionParser()
parser.add_option("-c", action="store_true", dest="recompile", help="recompile the demographic expansion cpp file")
parser.add_option("-s", action="store_true", dest="resimulate", help="resimulate demographic history and sampling schemes")
parser.add_option("-d", action="store_true", dest="rebuild_database", help="reinitialize the results database")
parser.add_option("-e", action="store_true", dest="evolve", help="evolve sequences along genealogies")
parser.add_option("-t", action="store_true", dest="show_tree", help="show tree and sequences for each locus")
parser.add_option("-b", action="store_true", dest="delimit", help="rerun bpp species delimitation analysis")

(options, args) = parser.parse_args()

# database
conn = sqlite3.connect('data/results.db')

# rebuild
if(options.rebuild_database == True):
    gensp.drop_results_table(conn)
    gensp.create_results_table(conn)
    print("Results database reinitialized")

# compile quetzal
if(options.recompile == True):
    try:
        subprocess.call([
        'g++', 'main.cpp', '-Wall', '-Ofast', '-std=c++17', '-I/usr/include/gdal',
        '-lboost_program_options', '-L/usr/lib/gdal',
        '-lgdal', '-lsqlite3', '-I/home/becheler/dev'
        ])
        print("Quetzal simulation compiled")
    except:
        sys.exit("Failed to compile main.cpp")

# run quetzal
if(options.resimulate == True):
    tmp = subprocess.call(
    ['./a.out',
    '--landscape', 'australia_precipitation_6032.tif',
    '--n_sim_demo', '1',
    '--n_sim_gen', '300',
    '--n_sample', '30',
    '--n_loci', '5'
    ])
    print("Quetzal simulation done")

cursor = conn.execute("SELECT id, genealogies, imap from results")

p1 = []
p2 = []
d = {'lon': [], 'lat': [], 'dist': [], 'area':[], 'p1': [], 'p2': []}

cursor = conn.execute("SELECT id, genealogies, imap from results")
for row in cursor:
    ID = row[0]
    newicks = row[1].split("\n\n")
    imap = row[2]
    print("ID = ", ID)
    if(options.evolve == True):
        sequence_size = 150
        scale_tree = 0.00001
        phylip = gensp.evolve(newicks, sequence_size, scale_tree, options.show_tree)
        with open(phylip) as f:
            s = f.read()
            conn.execute("""UPDATE results SET sequences = ? WHERE id = ?""", (s, ID))
            conn.commit()
        os.remove(phylip)
        print("Sequences inserted")
    if(options.delimit == True):
        # Run BPP
        cursor2 = conn.execute("""SELECT sequences from results WHERE id = ?""", (ID,))
        tup = cursor2.fetchone()
        sequences = tup[0]
        output = gensp.run_bpp(sequences, imap)
        conn.execute("""UPDATE results SET bpp = ? WHERE id  = ?""", (output, ID))
        conn.commit()
        print("BPP output inserted")
    # Retrieve probabilities and number of species
    cursor3 = conn.execute("""SELECT bpp, lon, lat, distance FROM results WHERE id = ?""", (ID,))
    tup = cursor3.fetchone()
    content = tup[0]
    d['lon'].append(tup[1])
    d['lat'].append(tup[2])
    d['dist'].append(tup[3])
    p = re.compile('P\[(\d)\] = (\d+\.\d+)')
    iterator = p.finditer(content)
    for match in iterator:
        if(match.group(1) == "1"):
            d['p1'].append(float(match.group(2)))
        elif(match.group(1) == "2"):
            d['p2'].append(float(match.group(2)))
        else:
            raise ValueError("Number of species should be 1 or 2")

d['area'] = [ 'Origin' if item <= 137 else 'York' for item in d['lon'] ]
print(d)
df = DataFrame(data=d)
print(df)

#create unique list of names
UniqueNames = df.area.unique()
#create a data frame dictionary to store your data frames
DataFrameDict = {elem : DataFrame for elem in UniqueNames}

for key in DataFrameDict.keys():
    DataFrameDict[key] = df[:][df.area == key]

df1 = DataFrameDict['Origin']
plt.scatter(df1['dist'], df1['p2'])
plt.show()

df2 = DataFrameDict['York']
plt.scatter(df2['dist'], df2['p2'])
plt.show()

df.to_csv('data.txt')
