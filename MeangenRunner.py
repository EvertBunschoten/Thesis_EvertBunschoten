from Meangen2Parablade import Meangen2Parablade
from Parablade2UMG2 import WriteUMG
import os
import sys
import math
import time

R = 1.0
phi = 0.5
psi = 0.8
r_m = 0.35

start_time = time.time()
BLADE_HOME = '/home/evert/Documents/TU_Delft_administratie/Thesis/parablade'
sys.path.append(BLADE_HOME+'/src/')
sys.path.append(BLADE_HOME+'/common/')
sys.path.append(BLADE_HOME+'/bin/')

M = Meangen2Parablade(R, phi, psi, r_m)
N_b = M.N_b
pitch_1 = 2*math.pi*r_m/N_b[0, 0]
pitch_2 = 2*math.pi*r_m/N_b[0, 1]

os.system("mkdir Bladerow1")
os.system("mv Bladerow_1.cfg ./Bladerow1/")
os.chdir("./Bladerow1")
os.system("MakeBlade.py Bladerow_1.cfg")
WriteUMG(1, M)
os.system("HYMESH.sh")
os.chdir("./..")

os.system("mkdir Bladerow2")
os.system("mv Bladerow_2.cfg ./Bladerow2/")
os.chdir("./Bladerow2")
os.system("MakeBlade.py Bladerow_2.cfg")
WriteUMG(2, M)
os.system("HYMESH.sh")

T = time.time() - start_time

print("Total meshing took "+ str(T) + "seconds")


