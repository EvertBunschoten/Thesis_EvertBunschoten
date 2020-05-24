import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os

class WriteUMG:
    def __init__(self, rowNumber, stage, Meangen, bodyForce, blade):

        self.n_blade = 1
        self.rowNumber = rowNumber
        self.stage = stage
        self.Meangen = Meangen
        self.name = "test"

        thisDir = os.getcwd()
        self.makeGeom()

        self.n_bl = 5
        self.thc_bl = 1.4 * 10**-5
        self.h_max_inflow = 0.0025
        self.h_min_inflow = 0.00025
        self.radCrv_inflow = 5

        self.h_max_perio = 0.0025
        self.h_min_perio = 0.00025
        self.radCrv_perio = 5

        self.h_max_blade = 0.0003
        self.h_min_blade = 0.00003
        self.radCrv_blade = 5

        self.h_max_outflow = 0.0025
        self.h_min_outflow = 0.00025
        self.radCrv_outflow = 5

        os.system("cp ../../createmesh.template ./")
        if bodyForce:
            print("Starting BFM mesh computation")
            self.makeBFMMesh()
        else:
            print("No BFM mesh creation requested")
        os.chdir(thisDir)
        if blade:
            print("Starting blade mesh computation")
            self.makeBladeMesh()
        else:
            print("No blade mesh creation requested")
    def makeGeom(self):

        x_le = self.Meangen.X_LE
        x_te = self.Meangen.X_TE
        coordDir = os.getcwd() + "/output/coordinates/surface_coordinates.csv"
        p, x, y, u, v = np.loadtxt(coordDir, unpack=True, delimiter=',\t', skiprows=1)
        i_min = list(x).index(min(x))
        i_max = list(x).index(max(x))
        if i_min < i_max:
            x_down = x[i_min:i_max+1]
            y_down = y[i_min:i_max+1]
            x_up = np.concatenate((x[i_max:], x[:i_min+1]), axis=0)
            y_up = np.concatenate((y[i_max:], y[:i_min + 1]), axis=0)
        else:
            x_down = np.concatenate((x[i_min:], x[:i_max+1]), axis=0)
            y_down = np.concatenate((y[i_min:], y[:i_max+1]), axis=0)
            x_up = x[i_max:i_min+1]
            y_up = y[i_max:i_min+1]
        Y_up = interp1d(x_up, y_up, kind='linear')
        y_up = Y_up(x_down)
        x_av = x_down
        y_av = 0.5*(y_up + y_down)

        x_0 = x_av[0]
        y_0 = y_av[0]
        x_1 = x_av[-1]
        y_1 = y_av[-1]

        if self.Meangen.machineType == 'C':
            if self.rowNumber % 2 != 0:
                pitch = (2 * np.pi * self.Meangen.r_m[self.stage - 1])/self.Meangen.N_b_R[self.stage - 1]
            else:
                pitch = (2 * np.pi * self.Meangen.r_m[self.stage - 1]) / self.Meangen.N_b_S[self.stage - 1]
        else:
            if self.rowNumber % 2 != 0:
                pitch = (2 * np.pi * self.Meangen.r_m[self.stage - 1])/self.Meangen.N_b_S[self.stage - 1]
            else:
                pitch = (2 * np.pi * self.Meangen.r_m[self.stage - 1]) / self.Meangen.N_b_R[self.stage - 1]
        self.Pitch = pitch
        if self.stage == 1 and self.rowNumber == 1:
            x_2 = x_0 - (x_1 - x_0)
            x_6 = 0.5 * (x_te[1, self.rowNumber - 1] + self.Meangen.X_LE[1, self.rowNumber])
        elif self.stage == self.Meangen.n_stage and self.rowNumber == 2:
            x_2 = 0.5 * (x_le[1, self.rowNumber - 1] + self.Meangen.X_TE[1, self.rowNumber - 2])
            x_6 = x_1 + (x_1 - x_0)
        else:
            x_2 = 0.5 * (x_le[1, self.rowNumber - 1] + self.Meangen.X_TE[1, self.rowNumber - 2])
            x_6 = 0.5 * (x_te[1, self.rowNumber - 1] + self.Meangen.X_LE[1, self.rowNumber])

        y_2 = y_0 - 0.5 * pitch

        x_3 = x_2
        y_3 = y_0 + (0.5 + self.n_blade - 1) * pitch

        x_4 = x_0
        y_4 = y_3

        x_5 = x_1
        y_5 = y_1 + (0.5 + self.n_blade - 1) * pitch

        y_6 = y_5
        x_7 = x_6
        y_7 = y_1 - 0.5 * pitch

        x_8 = x_1
        y_8 = y_7

        x_9 = x_0
        y_9 = y_0 - 0.5 * pitch

        self.names = ["BLADE", "BLADE", "INFLOW", "PER_INF_DOWN", "PER_CHANNEL_DOWN", "PER_GAP_DOWN",
                      "PER_INF_UP", "PER_CHANNEL_UP", "PER_GAP_UP", "OUTFLOW"]
        self.types = [8, 8, 1, 4, 4, 4, 5, 5, 5, 3]
        self.order = [1, 2, 4, 5, 6, 10, 9, 8, 7, 3]
        self.periodic = [0, 0, 0, 7, 8, 9, 4, 5, 6, 0]
        self.X_curve = [x_down, x_down, [x_2, x_3], [x_2, x_9], x_av, [x_8, x_7], [x_3, x_4], x_av, [x_5, x_6], [x_7, x_6]]
        self.Y_curve = [y_down, y_up, [y_2, y_3], [y_2, y_9], y_av - 0.5 * pitch, [y_8, y_7], [y_3, y_4],
                        y_av + (0.5 + self.n_blade - 1) * pitch, [y_5, y_6], [y_7, y_6]]


    def makeBFMMesh(self):
        if os.path.exists("BFMMesh"):
            print("Directory already exists, rewriting files")
        else:
            os.system("mkdir BFMMesh")
        meshDir = os.getcwd() + "/BFMMesh"
        os.chdir(meshDir)
        if not os.path.exists("Db"):
            os.system("mkdir Db")
        fileDir = os.getcwd() + "/Db/"
        os.chdir(fileDir)

        geomFile = open("geometry."+self.name, "w+")
        geomFile.write("Number of surfaces\n%i\n" % (len(self.names) - 2))
        for i in range(2, len(self.names)):
            geomFile.write("%s\n' S '\n" % (self.names[i]))
            geomFile.write("dim\tnp\n2\t%i\nx\ty\n" % (len(self.X_curve[i])))
            for j in range(len(self.X_curve[i])):
                geomFile.write("%+.5e\t%+.5e\n" % (self.X_curve[i][j], self.Y_curve[i][j]))
        geomFile.close()

        topoFile = open("topology."+self.name, "w+")
        topoFile.write("curve type\tperiodic curve\tModifiable curve\n")
        for i in range(2, len(self.types)):
            j = i - 1
            topoFile.write("%i\t%i\t%i\n" % (self.types[i], max([0, self.periodic[i]-2]), 0))
        topoFile.write("Number of ZONE\n")
        topoFile.write("1\nZONE 1\n")
        for i in range(2, len(self.order)-1):
            topoFile.write(" %i\n" % (self.order[i] - 2))
        topoFile.write("%i\n" % (-(self.order[-1]-2)))
        topoFile.close()

        spacingFile = open("spacingcontrol."+self.name, "w+")
        spacingFile.write("thk_bl\tn\tBC\tGEOM\tCV\n")
        spacingFile.write("%+.5e\t%i\taxl\t0\n\n" % (self.thc_bl, 5))
        spacingFile.write("PITCH\txc\tyc\n")
        spacingFile.write("%+.5e\t1.0\t1.0\n\n" % (self.Pitch))
        spacingFile.write("1\tINFLOW\th_min\th_max\tNode per RadCRv\n")
        spacingFile.write("%+.5e\t%+.5e\t%i\n" % (self.h_min_inflow, self.h_max_inflow, self.radCrv_inflow))
        spacingFile.write("8\tBLADE\th_min\th_max\tNode per RadCRv\n")
        spacingFile.write("%+.5e\t%+.5e\t%i\n" % (self.h_min_blade, self.h_max_blade, self.radCrv_blade))
        spacingFile.write("3\tOUTFLOW\th_min\th_max\tNode per RadCRv\n")
        spacingFile.write("%+.5e\t%+.5e\t%i\n" % (self.h_min_outflow, self.h_max_outflow, self.radCrv_outflow))
        spacingFile.write("4\tPERIO\th_min\th_max\tNode per RadCRv\n")
        spacingFile.write("%+.5e\t%+.5e\t%i\n" % (self.h_min_perio, self.h_max_perio, self.radCrv_perio))
        spacingFile.write("5\tPERIO\th_min\th_max\tNode per RadCRv\n")
        spacingFile.write("%+.5e\t%+.5e\t%i\n\n" % (self.h_min_perio, self.h_max_perio, self.radCrv_perio))
        spacingFile.write("NZONES\n1\n")
        spacingFile.close()

        optionsFile = open("options", "w+")
        optionsFile.write("fmt\tname\n")
        optionsFile.write("'grd'\t'%s'\n" % (self.name))
        optionsFile.write("optimization\n1\nmax element deformation\n1.0\nlayer of background grid\n3\n")
        optionsFile.write("Periodic geometry\n.true\t1.e-6\nScaling for SU2 file\n1.0\nnumber of boundary layers\n")
        optionsFile.write("0\t%+.5e\n" % self.thc_bl)
        optionsFile.write("Graph for hybrid mesh construction\n.true\nKind of radial basis function(1-11)\n11\n")
        optionsFile.write("Support radius for compact basis functions\n0.05\n")
        optionsFile.close()

        os.chdir(meshDir)
        os.system("HYMESH.sh")

        os.system("cp ../createmesh.template ./createmesh.cfg")
        os.system("sed -i 's/PITCH/%+.5e/g' createmesh.cfg" % self.Pitch)
        os.system("SU2_PERIO < createmesh.cfg")



    def makeBladeMesh(self):
        if os.path.exists("BladeMesh"):
            print("Directory already exists, rewriting files")
        else:
            os.system("mkdir BladeMesh")
        meshDir = os.getcwd() + "/BladeMesh"
        os.chdir(meshDir)
        if not os.path.exists("Db"):
            os.system("mkdir Db")
        fileDir = os.getcwd() + "/Db/"
        os.chdir(fileDir)

        geomFile = open("geometry." + self.name, "w+")
        geomFile.write("Number of surfaces\n%i\n" % (len(self.names)))
        for i in range(len(self.names)):
            geomFile.write("%s\n' S '\n" % (self.names[i]))
            geomFile.write("dim\tnp\n2\t%i\nx\ty\n" % (len(self.X_curve[i])))
            for j in range(len(self.X_curve[i])):
                geomFile.write("%+.5e\t%+.5e\n" % (self.X_curve[i][j], self.Y_curve[i][j]))
        geomFile.close()

        topoFile = open("topology." + self.name, "w+")
        topoFile.write("curve type\tperiodic curve\tModifiable curve\n")
        for i in range(len(self.types)):
            j = i - 1
            topoFile.write("%i\t%i\t%i\n" % (self.types[i], max([0, self.periodic[i]]), 0))
        topoFile.write("Number of ZONE\n")
        topoFile.write("1\nZONE 1\n")
        for i in range(len(self.order) - 1):
            topoFile.write(" %i\n" % (self.order[i]))
        topoFile.write("%i\n" % (-(self.order[-1])))
        topoFile.close()

        spacingFile = open("spacingcontrol." + self.name, "w+")
        spacingFile.write("thk_bl\tn\tBC\tGEOM\tCV\n")
        spacingFile.write("%+.5e\t%i\taxl\t0\n\n" % (self.thc_bl, self.n_bl))
        spacingFile.write("PITCH\txc\tyc\n")
        spacingFile.write("%+.5e\t1.0\t1.0\n\n" % (self.Pitch))
        spacingFile.write("1\tINFLOW\th_min\th_max\tNode per RadCRv\n")
        spacingFile.write("%+.5e\t%+.5e\t%i\n" % (self.h_min_inflow, self.h_max_inflow, self.radCrv_inflow))
        spacingFile.write("8\tBLADE\th_min\th_max\tNode per RadCRv\n")
        spacingFile.write("%+.5e\t%+.5e\t%i\n" % (self.h_min_blade, self.h_max_blade, self.radCrv_blade))
        spacingFile.write("3\tOUTFLOW\th_min\th_max\tNode per RadCRv\n")
        spacingFile.write("%+.5e\t%+.5e\t%i\n" % (self.h_min_outflow, self.h_max_outflow, self.radCrv_outflow))
        spacingFile.write("4\tPERIO\th_min\th_max\tNode per RadCRv\n")
        spacingFile.write("%+.5e\t%+.5e\t%i\n" % (self.h_min_perio, self.h_max_perio, self.radCrv_perio))
        spacingFile.write("5\tPERIO\th_min\th_max\tNode per RadCRv\n")
        spacingFile.write("%+.5e\t%+.5e\t%i\n\n" % (self.h_min_perio, self.h_max_perio, self.radCrv_perio))
        spacingFile.write("NZONES\n1\n")
        spacingFile.close()

        optionsFile = open("options", "w+")
        optionsFile.write("fmt\tname\n")
        optionsFile.write("'grd'\t'%s'\n" % (self.name))
        optionsFile.write("optimization\n1\nmax element deformation\n1.0\nlayer of background grid\n3\n")
        optionsFile.write("Periodic geometry\n.true\t1.e-6\nScaling for SU2 file\n1.0\nnumber of boundary layers\n")
        optionsFile.write("%i\t%+.5e\n" % (self.n_bl, self.thc_bl))
        optionsFile.write("Graph for hybrid mesh construction\n.true\nKind of radial basis function(1-11)\n11\n")
        optionsFile.write("Support radius for compact basis functions\n0.05\n")
        optionsFile.close()

        os.chdir(meshDir)
        os.system("HYMESH.sh")
        os.system("cp ../createmesh.template ./createmesh.cfg")
        os.system("sed -i 's/PITCH/%+.5e/g' createmesh.cfg" % self.Pitch)
        os.system("SU2_PERIO < createmesh.cfg")

