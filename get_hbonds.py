import sys
import string
import math

import Bio.PDB
from matplotlib import pyplot



def get_acceptors(model):

    acceptors = []

    for chain in model:
        for residue in chain:
            if residue.has_id("C") and residue.has_id("O"):
                acceptors.append(residue)

    return acceptors



def get_donors(model):

    donors = []

    for chain in model:
        for residue in chain:
            if (residue.has_id("H") or residue.has_id("HN")) and residue.has_id("N"):
                donors.append(residue)

    return donors



def get_model(filename):

    parser = Bio.PDB.PDBParser()
    structure = parser.get_structure("pdb", filename)

    for model in structure:

        # Return first model
        return model



def calc_energy(donor, acceptor):

    o = acceptor['O'].get_vector()
    c = acceptor['C'].get_vector()
    try:
        h = donor['H'].get_vector()
    except:
        h = donor['HN'].get_vector()
    n = donor['N'].get_vector()

    r_on = (o - n).norm()
    r_ch = (c - h).norm()
    r_oh = (o - h).norm()
    r_cn = (c - n).norm()

    e = 0.084 * ( 1.0/r_on +  1.0/r_ch - 1.0/r_oh - 1.0/r_cn) * 332.0

    return e



def has_hbond(donor, acceptor):

    if donor.get_id()[1] == acceptor.get_id()[1] or \
            donor.get_id()[1] == acceptor.get_id()[1] + 1:
        return False

    e = calc_energy(donor, acceptor)

    if e < -0.5:
        return True
    else:
        return False



def parse_data(filename):

    data = []

    f = open(filename, "r")
    lines = f.readlines()

    for line in lines:
        if line[0] == "#" or len(line) < 4:
            continue

        tokens = string.split(line)

        start = int(tokens[0])
        end = int(tokens[1])
        uptake = float(tokens[2])/100.0

        residues = range(start, end + 1)
        length = len(residues)

        data_point = dict()

        data_point['residues'] = residues
        data_point['uptake'] = uptake
        data_point['length'] = length

        data.append(data_point)

    f.close()

    return data



def calc_hbonds(model):

    donors = get_donors(model)
    acceptors = get_acceptors(model)

    hbonds = []

    for donor in donors:
        for acceptor in acceptors:
            if has_hbond(donor, acceptor):

                hbonds.append(donor.get_id()[1])
                # print donor.get_id()[1], acceptor.get_id()[1]

    return hbonds



if __name__ == "__main__":

    pdb_filename = sys.argv[1]
    data_filename = sys.argv[2]

    model = get_model(pdb_filename)

    hbonds = calc_hbonds(model)
    data = parse_data(data_filename)

    print hbonds

    x = []
    y = []

    for d in data:

        n_bonds = 0.0

        for residue in d['residues']:
            if residue in hbonds:
                n_bonds += 1

        if d['length'] == 2:
            print d['uptake'], n_bonds/d['length'], d['length']

            x.append(d['uptake'])
            y.append(n_bonds/d['length'])


    pyplot.plot(x, y, 'ko')
    pyplot.xlabel("Deuterium uptake")
    pyplot.ylabel("Avg. hbonds")
    pyplot.grid(True)
    pyplot.savefig("hGH.png")

