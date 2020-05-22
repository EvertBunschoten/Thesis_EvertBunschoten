import numpy as np
import matplotlib.pyplot as plt
import pdb
import sys
import os
from Meangen2Parablade import Meangen2Parablade
from Parablade2UMG2 import WriteUMG
import time

BLADE_HOME = '/home/evert/Documents/TU_Delft_administratie/Thesis/parablade'
# BLADE_HOME = os.environ["BLADE_HOME"]

sys.path.append(BLADE_HOME+'/bin/')
sys.path.append(BLADE_HOME+'/src/')
sys.path.append(BLADE_HOME+'/common/')

#---------------------------------------------------------------------------------------------#
# Importing ParaBlade classes and functions
#---------------------------------------------------------------------------------------------#
from common import *
from config import *

DIR = os.getcwd() + '/'
# try:
#     INFile = DIR + sys.argv[-1]
# except:
#     INFile = DIR + 'M2P.cfg'  # Default File name

try:
    IN = ReadUserInput(DIR + 'M2P.cfg')
except:
    raise Exception('\n\n\n''Something went wrong when reading the configuration file,exiting the program...'
                    '\n\nTo call MakeBlade.py from terminal type:'
                    '\n\tMakeBlade.py <configuration file name>')
t_start = time.time()

M = Meangen2Parablade(IN)

n_stage = int(IN["N_stage"][0])
n_rows = 2*n_stage

# for i in range(n_rows):
#     os.system("mkdir Bladerow_"+str(i+1))
#     os.system("mv "+os.getcwd() + "/Bladerow_"+str(i+1)+".cfg " + os.getcwd() + "/Bladerow_"+str(i+1)+"/")
for i in range(n_stage):
    if os.path.isdir("Stage_"+str(i+1)):
        if os.path.isdir("Stage_"+str(i+1) + "/Bladerow_1"):
            os.system("mv Bladerow_"+str(2*i + 1) + ".cfg" + " Stage_"+str(i+1)+"/Bladerow_1/Bladerow.cfg")
        else:
            os.system("mkdir "+"Stage_"+str(i+1)+"/Bladerow_1")
            os.system("mv Bladerow_" + str(2 * i + 1) + ".cfg" + " Stage_" + str(i + 1) + "/Bladerow_1/Bladerow.cfg")
        if os.path.isdir("Stage_"+str(i+1) + "/Bladerow_2"):
            os.system("mv Bladerow_" + str(2 * i + 2) + ".cfg" + " Stage_" + str(i + 1) + "/Bladerow_2/Bladerow.cfg")
        else:
            os.system("mkdir " + "Stage_" + str(i + 1) + "/Bladerow_2")
            os.system("mv Bladerow_" + str(2 * i + 2) + ".cfg" + " Stage_" + str(i + 1) + "/Bladerow_2/Bladerow.cfg")
    else:
        os.system("mkdir Stage_"+str(i+1))
        os.system("mkdir " + "Stage_" + str(i + 1) + "/Bladerow_1")
        os.system("mv Bladerow_" + str(2 * i + 1) + ".cfg" + " Stage_" + str(i + 1) + "/Bladerow_1/Bladerow.cfg")
        os.system("mkdir " + "Stage_" + str(i + 1) + "/Bladerow_2")
        os.system("mv Bladerow_" + str(2 * i + 2) + ".cfg" + " Stage_" + str(i + 1) + "/Bladerow_2/Bladerow.cfg")

row = 1

for i in range(n_stage):
    for j in [1, 2]:
        os.chdir(os.getcwd()+"/Stage_"+str(i+1)+"/Bladerow_"+str(j)+"/")
        if IN['PLOT_BLADE'] == 'YES':
            print("plotting blade")
            os.system("PlotBlade.py Bladerow.cfg")
        os.system("MakeBlade.py Bladerow.cfg")

        if IN['N_dim'][0] == 2:
            WriteUMG(row, i+1, M)
            print("Starting UMG2 meshing...", end='                 ')
            os.system("HYMESH.sh > mesher.log")
            print("Done!")
        row += 1
        os.chdir(DIR)

print("Total geometry and mesh generation took "+str(format(time.time() - t_start, ".2f")) + " seconds")

