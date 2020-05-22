import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os

class WriteUMG:
    def __init__(self, rowNumber, stage, Meangen):
        if os.path.exists("Db"):
            print("Directory already exists")
        else:
            os.system("mkdir Db")
        self.rowNumber = rowNumber
        self.stageNumber = stage
        self.Meangen = Meangen
       # self.pitch = 2*np.pi*Meangen.r_m/Meangen.N_b
        self.bladeName = "test"
        self.BL_count = 10
        self.t_BL = 0.1
        self.n_x = 1.0
        self.n_y = 1.0
        self.N_blade = 1

        self.CURVES = []
        self.NP = []
        self.X = []
        self.Y = []
        self.TYPE = []
        self.PERIO = []

        self.cdir = os.getcwd()
        self.WriteGeom()
        self.WriteTopo()
        self.WriteOptions()
        self.WriteSpacing()

    def WriteSpacing(self):
        f = open(self.cdir + "/Db/spacingcontrol."+self.bladeName, "w+")
        f.write("thk_bl\tn BC\tGEOM\tCV\n")
        f.write("4e-5\t5\taxl\t0\n\n")
        f.write("PITCH\txc\tyc\n")
        pitch = 2*np.pi*self.Meangen.r_m[self.stageNumber - 1]/self.Meangen.N_b[self.rowNumber-1]
        f.write(str(format(pitch, ".5e"))+"\t"+str(self.n_x) + "\t" + str(self.n_y) + "\n\n")
        f.write("1 INFLOW\t h_min\th_max\tNode per RadCRv\n")
        f.write(str(0.0030)+"\t"+str(0.0030)+"\t5.\n")
        f.write("3 OUTFLOW\t h_min\th_max\tNode per RadCRv\n")
        f.write(str(0.0015)+"\t"+str(0.0015)+"\t5.\n")
        f.write("5 PERIO\t h_min\th_max\tNode per RadCRv\n")
        f.write(str(0.00125)+"\t"+str(0.00125)+"\t5.\n")
        f.write("4 PERIO\t h_min\th_max\tNode per RadCRv\n")
        f.write(str(0.00125) + "\t" + str(0.00125) + "\t5.\n")
        f.write("8 BLADE\t h_min\th_max\tNode per RadCRv\n")
        f.write(str(0.00025) + "\t"+str(0.0005)+"\t15.\n\n")
        f.write("NZONES\n")
        f.write("0\n")
        f.write("RADIUS\tXC\tYC\th\n")
        f.write("0.005\t"+str(0.0)+"\t"+str(0.0)+"\t"+str(0.000025)+"\n")

        f.close()
    def WriteOptions(self):
        f = open(self.cdir + "/Db/options", "w+")
        f.write("fmt \t name\n")
        f.write("'grd'\t '"+self.bladeName+"'\n")
        f.write("optimization\n")
        f.write("1.\n")
        f.write("max element deformation\n")
        f.write("1.\n")
        f.write("layer of the background grid\n")
        f.write("3\n")
        f.write("Periodic geometry\t tol perio\n")
        f.write(".true\t1.e-6\n")
        f.write("Scaling for SU2 file\n")
        f.write("1.\n")
        f.write("number of boundary layers\tBL thickness\n")
        f.write(str(self.BL_count) + "\t" + str(self.t_BL) + "\n")
        f.write("Graph for hybrid mesh construction\n")
        f.write(".true\n")
        f.write("Kind of radial basis function(1-11)\n")
        f.write("11\n")
        f.write("Support radius for compact basis functions\n")
        f.write("0.05\n")
        f.close()

    def WriteTopo(self):
        f = open(self.cdir + "/Db/topology.test", "w+")
        f.write("\t\t  curve type\tperiodic curve\tModifiable curve\n")
        for i in range(self.N_blade):
            f.write("\t\t\t\t  " + str(self.Types[0]) + "\t\t\t\t\t" + str(self.Perio[0]) + "\t\t\t\t\t0\n")
            f.write("\t\t\t\t  " + str(self.Types[0]) + "\t\t\t\t\t" + str(self.Perio[0]) + "\t\t\t\t\t0\n")
        for i in range(1, len(self.Names)):
            f.write("\t\t\t\t  " + str(self.Types[i]) + "\t\t\t\t\t" + str(self.Perio[i]) + "\t\t\t\t\t0\n")
        f.write("Number of ZONE\n")
        f.write("1\n")
        f.write("ZONE 1\n")
        for i in range(1, 2 * self.N_blade + 8):
            f.write(" " + str(i) + "\n")
        f.write(str(-2 * self.N_blade - 8) + "\n")
        f.close()
    def WriteGeom(self):

        N_blade = self.N_blade
        pitch = 2*np.pi*self.Meangen.r_m[self.stageNumber - 1]/self.Meangen.N_b[self.rowNumber-1]

        rowGap = 0.25

        p, x, y, u, v = np.loadtxt(self.cdir + "/output/coordinates/surface_coordinates.csv", unpack=True, delimiter=',',
                                   skiprows=1)

        self.Names = ["BLADE", "PER_INF_DOWN", "INFLOW", "PER_INF_UP", "PER_PAS_UP", "PER_GAP_UP", "OUTFLOW",
                 "PER_GAP_DOWN", "PER_PAS_DOWN"]
        self.Types = [8, 4, 1, 5, 5, 5, 3, 4, 4]
        self.Perio = [0, 7, 0, 4, 5, 6, 0, 8, 9]

        side_switch = int(len(u)/2)

        c_ax = max(x) - min(x)
        i_min = list(x).index(min(x))
        i_max = list(x).index(max(x))

        if i_min < i_max:
            x_s1 = x[i_min:i_max + 1]
            y_s1 = y[i_min:i_max + 1]

            x_s2 = np.concatenate((x[i_max:], x[:i_min+1]), axis=0)
            y_s2 = np.concatenate((y[i_max:], y[:i_min+1]), axis=0)
        else:
            x_s1 = np.append(x[i_min:], x[:i_max + 1])
            y_s1 = np.append(y[i_min:], y[:i_max + 1])
            x_s2 = x[i_max:i_min + 1]
            y_s2 = y[i_max:i_min + 1]

        x_0 = min(x)
        y_0 = y[i_min]
        x_1 = max(x)
        y_1 = y[i_max]


        x_2 = x_0
        y_2 = y_0 - 0.5*pitch

        if self.Meangen.machineType == 'C':
            chords = np.array([self.Meangen.chord_R, self.Meangen.chord_S])
        else:
            chords = np.array([self.Meangen.chord_S, self.Meangen.chord_R])

        if self.rowNumber == 1 and self.stageNumber == 1:
            x_3 = x_2 - chords[0, self.stageNumber - 1]
            x_7 = 0.5*(self.Meangen.X_TE[1, self.rowNumber - 1] + self.Meangen.X_LE[1, self.rowNumber])

        elif self.rowNumber == self.Meangen.n_stage * 2:
            x_3 = 0.5*(self.Meangen.X_LE[1, self.rowNumber - 1] + self.Meangen.X_TE[1, self.rowNumber - 2])
            x_7 = x_1 + chords[1, self.stageNumber - 1]
        else:
            x_3 = 0.5*(self.Meangen.X_LE[1, self.rowNumber - 1] + self.Meangen.X_TE[1, self.rowNumber - 2])
            x_7 = 0.5*(self.Meangen.X_TE[1, self.rowNumber - 1] + self.Meangen.X_LE[1, self.rowNumber])

        y_3 = y_2

        x_4 = x_3
        y_4 = y_0 + pitch*(0.5 + N_blade - 1)

        x_5 = x_0
        y_5 = y_4

        # x_6 = x_1 + 0.5 * rowGap * c_ax
        x_6 = x_1
        y_6 = y_1 + pitch*(0.5 + N_blade - 1)

        y_7 = y_6

        x_8 = x_7
        y_8 = y_1 - 0.5*pitch

        x_9 = x_1
        y_9 = y_8

        y_up = interp1d(x_s2, y_s2)
        y_av = 0.5 * (y_s1 + y_up(x_s1))

        # plt.figure(0)
        # plt.title("Stage "+str(self.stageNumber)+", row "+str(self.rowNumber))
        # plt.plot([x_2, x_3, x_4, x_5], [y_2, y_3, y_4, y_5], 'k-*')
        # plt.plot(x_s1, y_av + (0.5 + self.N_blade - 1)*pitch, 'k')
        # plt.plot(x_s1, y_av - pitch * 0.5, 'k')
        # plt.plot([x_6, x_7, x_8, x_9], [y_6, y_7, y_8, y_9], 'k-*')
        #
        # for i in range(N_blade):
        #     plt.plot(x_s1, y_s1 + pitch * i, 'b')
        #     plt.plot(x_s2, y_s2 + pitch * i, 'r')
        # plt.axis('equal')
        # plt.show()

        X = np.array([[x_2, x_3], [x_3, x_4], [x_4, x_5], x_s1, [x_6, x_7], [x_7, x_8], [x_8, x_9], x_s1[::-1]])
        Y = np.array([[y_2, y_3], [y_3, y_4], [y_4, y_5], y_av + (0.5 + self.N_blade - 1)*pitch, [y_6, y_7], [y_7, y_8], [y_8, y_9],
                      y_av[::-1] - 0.5*pitch])

        f = open(self.cdir + "/Db/geometry."+self.bladeName, "w+")

        f.write("Number of surfaces\n")
        f.write("\t\t\t\t " + str(int(8 + 2 * N_blade)) + "\n")
        for i in range(N_blade):
            f.write(self.Names[0] + "\n")
            f.write("\t\t\t\t' S '\n")
            f.write("\t\t\t\tdim\t\t\t\tnp\n")
            f.write("\t\t\t\t  2\t\t\t\t" + str(len(x_s1)) + "\n")
            f.write("\t\t\t\t  x\t\t\t\t  y\n")
            for j in range(len(x_s1)):
                f.write(str(format(x_s1[j], ".5e")) + "\t" + str(format(y_s1[j] + pitch * i, ".5e")) + "\n")
            f.write(self.Names[0] + "\n")
            f.write("\t\t\t\t' S '\n")
            f.write("\t\t\t\tdim\t\t\t\tnp\n")
            f.write("\t\t\t\t  2\t\t\t\t" + str(len(x_s2)) + "\n")
            f.write("\t\t\t\t  x\t\t\t\t  y\n")
            for j in range(len(x_s2)):
                f.write(str(format(x_s2[j], ".5e")) + "\t" + str(format(y_s2[j] + pitch * i, ".5e")) + "\n")

        for i in range(1, len(self.Names)):
            f.write(self.Names[i] + "\n")
            f.write("\t\t\t\t' S '\n")
            f.write("\t\t\t\tdim\t\t\t\tnp\n")
            f.write("\t\t\t\t  2\t\t\t\t" + str(len(X[i - 1])) + "\n")
            f.write("\t\t\t\t  x\t\t\t\t  y\n")
            for j in range(len(X[i - 1])):
                f.write(str(format(X[i - 1][j], ".5e")) + "\t" + str(format(Y[i - 1][j], ".5e")) + "\n")
        f.close()



