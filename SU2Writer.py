import numpy as np
import os
import sys

class writeBFMinput:
    n_sec = None
    n_points = None
    n_rows = 0

    def __init__(self, Meangen):
        self.M = Meangen
        self.dir = os.getcwd()
        self.BFMfile = open(self.dir + "/BFM_stage_input", "w+")
        self.getBFMinputs()
        self.BFMfile.write("%i\t%i\t%i\n" % (self.n_rows, self.n_sec, self.n_points))
        self.writeInputFile()
        self.BFMfile.close()

    def getBFMinputs(self):
        for i in range(self.M.n_stage):
            for j in [1, 2]:
                with open(self.dir + "/Stage_"+str(i+1)+"/Bladerow_"+str(j)+"/output/mesh_files/BFM_input", "r") as file:
                    first_line = file.readline().split('\t', 2)
                    first_line[-1] = first_line[-1].strip()
                    self.n_sec = int(first_line[0])
                    self.n_points = int(first_line[1])
                    self.n_rows += 1
                file.close()

    def writeInputFile(self):

        for i in range(self.M.n_stage):
            for j in [1, 2]:
                with open(self.dir + "/Stage_"+str(i+1)+"/Bladerow_"+str(j)+"/output/mesh_files/BFM_input", "r") as file:
                    lines = file.readlines()[1:]
                    if self.M.machineType == 'C':
                        rotFac = [1, 0]
                        blades = [self.M.N_b_R[i], self.M.N_b_S[i]]
                    else:
                        rotFac = [0, 1]
                        blades = [self.M.N_b_S[i], self.M.N_b_R[i]]
                    for line in lines:
                        self.BFMfile.write(line.strip()+"\t"+str(rotFac[j-1])+"\t"+str(int(blades[j-1]))+"\n")
                    #self.BFMfile.writelines(lines)
                file.close()






