import pyvolve
import sys, os, random, string, re
import ete3
import subprocess

class Sequence(object):
    """The Sequence object has a string *header* and
    various representations."""

    def __init__(self, header, seq):
        self.header = re.findall('^>(\S+)', header)[0]
        self.seq = seq

    def __len__(self):
        return len(self.seq)

    @property
    def phylip(self):
        return self.header + " " + self.seq.replace('.','-') + "\n"

    @property
    def fasta(self):
        return ">" + self.header + "\n" + self.seq + "\n"


def fasta_parse(path):
    """Reads the file at *path* and yields
       Sequence objects in a lazy fashion"""
    header = ''
    seq = ''
    with open(path) as f:
        for line in f:
            line = line.strip('\n')
            if line.startswith('>'):
                if header: yield Sequence(header, seq)
                header = line
                seq = ''
                continue
            seq += line
    yield Sequence(header, seq)

def fasta_to_phyl(fasta_file, phylip_file):
    # Check that the path is valid #
    if not os.path.exists(fasta_file): raise Exception("No file at %s." % fasta_file)
    # Use our two functions #
    seqs = fasta_parse(fasta_file)
    # Write the output to temporary file #
    temp_path = phylip_file + '.' + ''.join(random.choice(string.ascii_letters) for i in range(10))
    # Count the sequences #
    count = 0
    with open(temp_path, 'w') as f:
        for seq in seqs:
            f.write("^" + seq.phylip)
            count += 1
    # Add number of entries and length at the top #
    with open(temp_path, 'r') as old, open(phylip_file, 'w') as new:
        new.write(" " + str(count) + " " + str(len(seq)) + "\n")
        new.writelines(old)
    # Clean up #
    os.remove(temp_path)
    return

def create_results_table(conn):
    cursor = conn.cursor()
    cursor.execute("""
    CREATE TABLE results(
    id INTEGER PRIMARY KEY AUTOINCREMENT UNIQUE,
    sim_demo INTEGER,
    sim_gen INTEGER,
    n1 INTEGER,
    n2 INTEGER,
    n_loci INTEGER,
    genealogies TEXT,
    sequences TEXT,
    imap TEXT,
    lat REAL,
    lon REAL,
    distance REAL,
    bpp TEXT
    )
    """)
    conn.commit()
    return

def drop_results_table(conn):
    cursor = conn.cursor()
    cursor.execute("""DROP TABLE results""")
    conn.commit()
    return

def evolve(newicks, sequence_size, scale_tree, show_tree):
    temp = "temporary_sequences.fasta"
    phy_files = []

    my_model = pyvolve.Model("nucleotide")
    partition = pyvolve.Partition(models = my_model, size = sequence_size)

    for i in range(0, len(newicks)):

        newick = newicks[i]
        tree = pyvolve.read_tree(tree = newick, scale_tree = scale_tree)
        my_evolver = pyvolve.Evolver(tree = tree, partitions = partition)

        fasta_seqfile = "temp" + str(i) + ".fasta"
        phylip_seqfile = "temp" + str(i) + ".phyl"
        phy_files.append(phylip_seqfile)

        my_evolver(seqfile=fasta_seqfile, seqfmt = "fasta", ratefile = None, infofile = None)
        fasta_to_phyl(fasta_seqfile, phylip_seqfile)

        if(show_tree == True):
            t = ete3.PhyloTree(newick)
            t.link_to_alignment(alignment=open(fasta_seqfile, "r").read(), alg_format="fasta")
            t.show()

        os.remove(fasta_seqfile)

    phyl_output = "temp_seq.phyl"

    with open(phyl_output, 'w') as outfile:
        for fname in phy_files:
            with open(fname) as infile:
                outfile.write(infile.read())
                outfile.write("\n")
            os.remove(fname)

    return phyl_output

def run_bpp(sequences, imap):

    imap_f  = open("temp_imap.txt", "w")
    imap_f.write(imap)
    imap_f.close()

    seq_f = open("temp_seq.phyl", "w")
    seq_f.write(sequences)
    seq_f.close()

    tmp = subprocess.call(['./bpp','--cfile', 'bpp.ctl'])

    try:
        output = ""
        with open("out.txt") as infile:
            output = infile.read()
            os.remove("out.txt")
            os.remove("mcmc.txt")
            return output
    except IOError:
        os.remove("temp_imap.txt")
        os.remove("temp_seq.phyl")
        sys.exit("BPP outfiles not accessible")

    return
