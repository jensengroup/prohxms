#!/usr/bin/python
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

import sys
import os

from numpy import log

import prohxms_api as prohxms

from parameters import model1 as parameters


# Ideal gas constant in kcal/mol
R_GAS = 0.001987204118 

# 300 kelvin in units of kelvin.
T300K = 300.0

# R * T at 300K (i.e. kB * T in units of kcal/mol)
RT = R_GAS * T300K



def get_log_lik(k, n, uptake):

    alpha = None
    beta = None

    for parameter in parameters:
        if uptake < parameter[0]:

            (max_uptake, alpha, beta) = parameter

    if alpha is None or beta is None:
        print "An error occurred determining alpha and beta"
        exit(1)

    p = prohxms.beta_binom(k, n, alpha, beta)

    log_lik = -1.0 * log(p) 

    return log_lik



if __name__ == "__main__":

    pdb_filename = sys.argv[1]
    data_filename = sys.argv[2]

    model = prohxms.get_model(pdb_filename)
    hbonds = prohxms.calc_hbonds(model)

    data = prohxms.parse_data(data_filename)

    log_lik_sum = 0.0

    for d in data:

        #print d
        n_bonds = 0.0

        for residue in d['residues']:
            if residue in hbonds:
                n_bonds += 1

        # print d, d['uptake'], n_bonds/d['length'], d['length']
        k = n_bonds
        n = d['length']
        uptake = d['uptake']

        log_lik = get_log_lik(k, n, uptake)
        log_lik_sum += log_lik



        #print d, d['uptake'], n_bonds/d['length'], log_lik

    print " E = %4.3f kcal/mol" % (log_lik_sum * RT)

