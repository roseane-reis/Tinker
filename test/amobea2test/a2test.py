#!/usr/bin/env python

import os
from glob import glob
import subprocess
import sys
import numpy as np

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

CSI = "\x1B["
PASS = CSI + "32;1m" + "PASSED" + CSI + "0m"
FAIL = CSI + "31;1m" + "FAILED" + CSI + "0m"
POOP = u'\U0001F4A9'
BEER = u'\U0001F37A'

def make_green(s):
    return CSI + "32;1m" + s + CSI + "0m"

def make_red(s):
    return CSI + "31;1m" + s + CSI + "0m"

# remove all previous test output files
for filename in glob('*.out'):
    os.remove(filename)

print "Testing Tinker..."
my_bin = str("../../source/")

print "Testing binaries located at %s" % ''.join(my_bin)

# Find this script's directory and go there, to allow execution from any location
scriptpath = os.path.dirname(os.path.realpath(__file__))
os.chdir(scriptpath)
print "Running from %s" % scriptpath

def run_test(testname,my_exec,keyfile,xyzfile,options):
    outfile = open("%s.out"%testname, 'w')
    do_me = [my_bin + my_exec,'-k', keyfile, xyzfile] + options
#    print subprocess.list2cmdline(do_me)
    process = subprocess.Popen(do_me, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    stdout,stderr = process.communicate()
    # Echo stdout
    outfile.write(stdout)
    # Echo stderr
    outfile.write(stderr + "\n")
    outfile.close()

def get_energy(component,type,testname):
    if type == "reference":
        file = testname + ".out"
    elif type == "test":
        file = testname + ".ref"
    else:
        print "must input reference or test"
    for line in open(file, 'r'):
        if component in line:
            line_s = line.split("               ")
            energy = float(line_s[1])
    return energy

def get_gradient(testname):
    file = testname + ".out"
    for line in open(file, 'r'):
        if "Anlyt      RMS Gradient over All Atoms" in line:
            line_s = line.split("        ")
            ana = float(line_s[1])
        if "Numer      RMS Gradient over All Atoms" in line:
            line_s = line.split("        ")
            num = float(line_s[1])
    return ana,num

def get_virial(testname):
    file = testname + ".out"
    counter = 0
    for line in open(file, 'r'):
        if "Internal Virial Tensor :" in line:
            line_s = line.split()
            xx = float(line_s[4])
            counter = 1
        elif counter == 1:
            line_s = line.split()
            yy = float(line_s[1])
            counter = 2
        elif counter == 2:
            line_s = line.split()
            zz = float(line_s[2])
            ana = [xx,yy,zz]
            counter = 0
        if "Numerical Virial Diagonal :          " in line:
            line_s = line.split()
            xx0 = float(line_s[4])
            yy0 = float(line_s[5])
            zz0 = float(line_s[6])
            num = [xx0,yy0,zz0]
    return ana,num


def test_energy(component,method):
    name = method + "_energy"
    run_test(name,"analyze.x",method + ".key","dimer01.xyz",["e"])
    ref_ene = get_energy(component,"reference",name) 
    test_ene = get_energy(component,"test",name)
    if np.allclose(ref_ene,test_ene):
        print method + ' ' + component + ' Energy         PASSED'
    else:
        print method + ' ' + component + ' Energy         FAILED'

def test_gradient(method):
    name = method + "_gradient"
    options = ["y","y",str(0.00001)]
    run_test(name,"testgrad.x",method + ".key","dimer01.xyz",options)
    analytic,numeric = get_gradient(name)
    if np.allclose(analytic,numeric):
        print method + ' ' + ' Gradient         PASSED'
    else:
        print method + ' ' + ' Gradient         FAILED'

def test_virial(method):
    name = method + "_virial"
    options = ["v"]
    run_test(name,"analyze.x",method + ".key","dimer01.xyz",options)
    analytic, numeric = get_virial(name)
    if np.allclose(analytic,numeric):
        print method + ' ' + ' Virial         PASSED'
    else:
        print method + ' ' + ' Virial         FAILED'


test_energy("Multipole","pairwise")
test_energy("Pauli Repulsion","pairwise")
test_energy("Dispersion       ","pairwise")
test_energy("Polarization     ","pairwise")
test_energy("Extra ","pairwise")

test_energy("Multipole","neighborlist-ewald")
test_energy("Pauli Repulsion","neighborlist-ewald")
test_energy("Dispersion       ","neighborlist-ewald")
test_energy("Polarization     ","neighborlist-ewald")
test_energy("Extra ","neighborlist-ewald")

test_gradient("pairwise")

test_gradient("neighborlist-ewald")

test_virial("neighborlist-ewald")
