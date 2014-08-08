"""
Copyright (c) 2014, Anders S. Christensen,
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

__author__ = "Anders S. Christensen"
__copyright__ = "Copyright 2013, Anders S. Christensen"
__credits__ = ["Anders S. Christensen (2014) ProHXMS, https://github.com/andersx/prohxms/"]
__license__ = "BSD 2-clause"
__version__ = "1.0.0"
__maintainer__ = "Anders S. Christensen"
__email__ = "andersx@nano.ku.dk"
__status__ = "Development"


import sys
import string
import math
import os.path
import cPickle
import Bio.PDB
from scipy.special import binom, gamma

# Energy in kcal/mol before an electrostatic 
# interaction is considered a H-bond.
HB_BOND_CUT_OFF = -0.5

DEBUG = False

def get_donors(model):

    donors = []

    for chain in model:
        for residue in chain:
            if (residue.has_id("H") or residue.has_id("HN")) and residue.has_id("N"):
                donors.append(residue)

    return donors



def get_bb_acceptors(model):

    acceptors = []

    for chain in model:
        for residue in chain:
            if residue.has_id("C") and residue.has_id("O"):acceptors.append(residue)


    return acceptors



def get_sc_acceptors(model):

    donors = []

    for chain in model:
        for residue in chain:

            if residue.get_resname() == "ASP":
                if residue.has_id("CG") and residue.has_id("OD1") \
                    and ( residue.has_id("CG") and residue.has_id("OD2")):
                    donors.append(residue)

            elif residue.get_resname() == "GLU":
                if residue.has_id("CD") and residue.has_id("OE1") \
                    and ( residue.has_id("CD") and residue.has_id("OE2")):
                    donors.append(residue)

            elif residue.get_resname() == "ASN":
                if residue.has_id("CG") and residue.has_id("OD1"):
                    donors.append(residue)

            elif residue.get_resname() == "GLN":
                if residue.has_id("CD") and residue.has_id("OE1"):
                    donors.append(residue)

            elif residue.get_resname() == "SER":
                if residue.has_id("CB") and residue.has_id("OG"):
                    donors.append(residue)

            elif residue.get_resname() == "THR":
                if residue.has_id("CB") and residue.has_id("OG1"):
                    donors.append(residue)

            elif residue.get_resname() == "TYR":
                if residue.has_id("CZ") and residue.has_id("OH"):
                    donors.append(residue)
    return donors



def get_model(filename):

    # parser = Bio.PDB.PDBParser(QUIET=False)
    parser = Bio.PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("pdb", filename)

    for model in structure:

        # Return first model
        return model



def calc_energy(n, h, o, c):

    r_on = (o - n).norm()
    r_ch = (c - h).norm()
    r_oh = (o - h).norm()
    r_cn = (c - n).norm()

    e = 0.084 * ( 1.0/r_on +  1.0/r_ch - 1.0/r_oh - 1.0/r_cn) * 332.0

    return e



def has_bb_hbond(donor, acceptor):

    if donor.get_id()[1] == acceptor.get_id()[1] or \
            donor.get_id()[1] == acceptor.get_id()[1] + 1:
        return False

    o = acceptor['O'].get_vector()
    c = acceptor['C'].get_vector()
    try:
        h = donor['H'].get_vector()
    except:
        h = donor['HN'].get_vector()
    n = donor['N'].get_vector()


    e = calc_energy(n, h, o, c)

    if e < HB_BOND_CUT_OFF:
        if (DEBUG): print donor.get_id()[1], " -> ", acceptor.get_id()[1], "%4.2f" % (e)
        return True
    else:
        return False



def has_sc_hbond(donor, acceptor):

    try:
        h = donor['H'].get_vector()
    except:
        h = donor['HN'].get_vector()
    n = donor['N'].get_vector()


    e1 = float("inf")
    e2 = float("inf")


    if acceptor.get_resname() == "ASP":
        if acceptor.has_id("CG") and acceptor.has_id("OD1"):
            o = acceptor['OD1'].get_vector()
            c = acceptor['CG'].get_vector()

            e1 = calc_energy(n, h, o, c)

        if acceptor.has_id("CG") and acceptor.has_id("OD2"):
            o = acceptor['OD2'].get_vector()
            c = acceptor['CG'].get_vector()

            e2 = calc_energy(n, h, o, c)

    elif acceptor.get_resname() == "GLU":
        if acceptor.has_id("CD") and acceptor.has_id("OE1"):
            o = acceptor['OE1'].get_vector()
            c = acceptor['CD'].get_vector()

            e1 = calc_energy(n, h, o, c)

        if acceptor.has_id("CD") and acceptor.has_id("OE2"):
            o = acceptor['OE2'].get_vector()
            c = acceptor['CD'].get_vector()

            e2 = calc_energy(n, h, o, c)

    elif acceptor.get_resname() == "ASN":
        if acceptor.has_id("CG") and acceptor.has_id("OD1"):
            o = acceptor['OD1'].get_vector()
            c = acceptor['CG'].get_vector()

            e1 = calc_energy(n, h, o, c)

    elif acceptor.get_resname() == "GLN":
        if acceptor.has_id("CD") and acceptor.has_id("OE1"):
            o = acceptor['OE1'].get_vector()
            c = acceptor['CD'].get_vector()

            e1 = calc_energy(n, h, o, c)

    elif acceptor.get_resname() == "SER":
        if acceptor.has_id("CB") and acceptor.has_id("OG"):
            o = acceptor['OG'].get_vector()
            c = acceptor['CB'].get_vector()

            e1 = calc_energy(n, h, o, c)

    elif acceptor.get_resname() == "THR":
        if acceptor.has_id("CB") and acceptor.has_id("OG1"):
            o = acceptor['OG1'].get_vector()
            c = acceptor['CB'].get_vector()

            e1 = calc_energy(n, h, o, c)

    elif acceptor.get_resname() == "TYR":
        if acceptor.has_id("CZ") and acceptor.has_id("OH"):
            o = acceptor['OH'].get_vector()
            c = acceptor['CZ'].get_vector()

            e1 = calc_energy(n, h, o, c)


    if (e1 < HB_BOND_CUT_OFF) or (e2 < HB_BOND_CUT_OFF):
        if (DEBUG): print donor.get_id()[1], " -> ", acceptor.get_id()[1], min(e1, e2)
        return True
    else:
        return False


def calc_hbonds(model):

    donors = get_donors(model)
    bb_acceptors = get_bb_acceptors(model)
    sc_acceptors = get_sc_acceptors(model)

    hbonds = []

    for donor in donors:
        for acceptor in bb_acceptors:
            if has_bb_hbond(donor, acceptor):
                hbonds.append(donor.get_id()[1])

        for acceptor in sc_acceptors:
            if has_sc_hbond(donor, acceptor):
                hbonds.append(donor.get_id()[1])

    return hbonds


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


def load_pickle(filename):
    f = open(filename,"rb")
    p = cPickle.load(f)
    f.close()
    return(p)


def factorial(x):
    """
        Returns the factorial of x
    """
    try:
        return gamma(x + 1)
    except OverflowError:
        print "Overflow, x =",x
        exit(1)


def B(x, y):
    """
        Returns the factorial of x
    """
    return factorial(x - 1) * factorial(y - 1) / factorial(x + y - 1)


def beta_binom(k, n, alpha, beta):
    """
        Returns the factorial of x
    """
    return binom(n, k) * B(k + alpha, n - k + beta) / B(alpha, beta)

