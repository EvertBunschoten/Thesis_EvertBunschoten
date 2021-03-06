#!/home/evert/anaconda3/bin/python3
import numpy as np
import os
import time
from StagenReader import StagenReader
from scipy.interpolate import interp1d
import math



class Meangen2Parablade:

    Compressor = True   # Machine type(if False, it's an axial turbine)
    n_stage =       1
    Dimension =     3
    N_b =           [30, 30]
    R =             0.5
    phi =           0.5
    psi =           0.4
    r_m =           0.35
    chord =         0.04
    rowGap =        0.25
    stageGap =      0.5
    QO_LE =         90
    QO_TE =         90
    tc_max =        0.1
    x_tmax =        0.5


    def __init__(self, IN):
        # Getting the input parameters
        self.IN = IN

        self.n_stage = int(IN["N_stage"][0])    # Stage count
        self.Dimension = int(IN["N_dim"][0])    # Parablade output dimension
        self.Omega = IN["Omega"][0]             # Rotation speed
        self.mdot = IN["mass_flow"][0]          # Mass flow rate
        self.R_gas = IN["R_gas"][0]
        self.gamma = IN["gamma"][0]
        self.P_tin = IN["P_t_in"][0]
        self.T_tin = IN["T_t_in"][0]
        self.machineType = IN["TYPE"]           # Machine type
        self.N_b_R = IN["N_blade_R"]            # Blade count in rotor row
        self.N_b_S = IN["N_blade_S"]            # Blade count in stator row
        self.R = IN["R"]                        # Degree of reaction for each respective stage
        self.phi = IN["phi"]                    # Flow coefficient for each respective stage
        self.psi = IN["psi"]                    # Work coefficient for each respective stage
        self.r_m = IN["r_m"]                    # Mean stage radius for each respective stage
        self.chord_R = IN["chord_R"]            # Rotor chords
        self.chord_S = IN["chord_S"]            # Stator chords
        self.rowGap = IN["rowGap"]              # Row gap, normalized wrt first bladerow chord
        self.stageGap = IN["stageGap"]          # Stage gap, normalized wrt first bladerow chord
        self.Twist = IN["twist"]                # Meangen twist value
        self.QO_LE_R = IN["QO_LE_R"]            # Leading edge QO angle for rotor rows
        self.QO_TE_R = IN["QO_TE_R"]            # Trailing edge QO angle for rotor rows
        self.QO_LE_S = IN["QO_LE_S"]            # Leading edge QO angle for stator rows
        self.QO_TE_S = IN["QO_TE_S"]            # Trailing edge QO angle for stator rows
        self.eta_guess = IN["eta_guess"]        # Efficiency guess for each respective stage
        self.Dev_R_LE = IN["dev_R_LE"]          # Rotor incidence angles
        self.Dev_R_TE = IN["dev_R_TE"]          # Rotor deviation angles
        self.Dev_S_LE = IN["dev_S_LE"]          # Stator incidence angles
        self.Dev_S_TE = IN["dev_S_TE"]          # Stator deviation angles
        self.t_LE_R = IN["t_le_R"]
        self.t_TE_R = IN["t_te_R"]
        self.t_LE_S = IN["t_le_S"]
        self.t_TE_S = IN["t_te_S"]
        self.D_1_R = IN["d_1_R"]
        self.D_2_R = IN["d_2_R"]
        self.D_1_S = IN["d_1_S"]
        self.D_2_S = IN["d_2_S"]
        self.T_R = np.zeros([self.n_stage, 6])
        self.T_S = np.zeros([self.n_stage, 6])
        for i in range(self.n_stage):
            for j in range(6):
                self.T_R[i, j] = IN["T_"+str(j+1)+"_R"][i]
                self.T_S[i, j] = IN["T_" + str(j + 1) + "_S"][i]
        # Defining the machine type
        self.machineDefinition()

        # Writing Meangen input file
        self.meangenWriter()

        # Getting directory to Meangen executable and initializing Meangen
        print("Starting Meangen...", end='                 ')
        start_time = time.time()
        # meangen_dir = "/home/evert/Documents/TU_Delft_administratie/Thesis/Pythonscripts"
        #os.system(meangen_dir + "/a.out < "+meangen_dir+"/input > Output")
        HOME = os.environ["M2BFM"]
        os.system("Meangen.exe < "+HOME+"templates/input > Output")
        print("Done!")
        print("Meangen took "+str(time.time() - start_time) + " seconds")
        # Reading stagen.dat file for blade row geometry
        S = StagenReader()

        # Storing blade row inlet metal angles
        self.theta_in = S.theta_in
        # Storing blade row outlet metal angles
        self.theta_out = S.theta_out

        # Storing leading edge and trailing edge coordinates
        self.X_LE = S.X_LE
        self.X_TE = S.X_TE
        self.Z_LE = S.Z_LE
        self.Z_TE = S.Z_TE

        # Writing Parablade input file
        print("Writing Parablade input file...", end='          ')
        start_time = time.time()
        self.ParabladeWriter()
        print("Done!")
        print("Parablade file writing took "+ str(time.time() - start_time) + " seconds")

        # Storing Meangen output files
        self.storeFiles()

    def storeFiles(self):
        Dir = os.getcwd()
        os.system("mkdir "+Dir+"/MeangenOutput")
        os.system("mv meangen.in "+Dir+"/MeangenOutput/")
        os.system("mv meangen.out " + Dir + "/MeangenOutput/")
        os.system("mv meandesign.out " + Dir + "/MeangenOutput/")
        os.system("mv stagen.dat " + Dir + "/MeangenOutput/")
        os.system("mv Output " + Dir + "/MeangenOutput/")
    def meangenWriter(self):
        # Defining meangen input file
        save_path = os.getcwd()
        filename = "meangen"
        completename = os.path.join(save_path, filename + ".in")
        self.f = open(completename, 'wt')  # Opening input file
        self.f.write(self.machineType + "\n")  # Defining machine type
        self.f.write(self.flowPath + "\n")  # Defining flow path
        self.f.write(str(self.R_gas) + "  " + str(self.gamma) + "\n")  # Defining working fluid properties
        self.f.write(str(self.P_tin) + "  " + str(self.T_tin) + "\n")  # Defining inlet conditions
        self.f.write(str(self.n_stage) + "\n")  # Defining number of stages
        self.f.write(str(self.designPoint) + "\n")  # Defining blade design point
        self.f.write(str(self.Omega) + "\n")  # Input of rotation speed
        self.f.write(str(self.mdot) + "\n")  # Input of target mass flow rate
        # Looping over each stage, inputting the duty coefficients and radius of each stage
        for i in range(self.n_stage):
            radius = self.r_m[i]  # Calculating mean radius
            if i > 0:
                self.f.write("N\nN\n")  # Enabling new input for next stage
            self.f.write("A\n")  # Choice for flow angle calculation
            # Inputting duty coefficients
            self.f.write(str(self.R[i]) + "  " + str(self.phi[i]) + "  " + str(self.psi[i]) + "\n")
            self.f.write(str(self.radType) + "\n")  # Inputting radius input type
            self.f.write(str(radius) + "\n")  # Inputting mean radius
            if self.machineType == 'C':
                self.f.write(str(self.chord_R[i]) + "  " + str(self.chord_S[i]) + "\n")  # Inputting chords
            else:
                self.f.write(str(self.chord_S[i]) + "  " + str(self.chord_R[i]) + "\n")  # Inputting chords
            self.f.write(
                str(self.rowGap[i]) + "  " + str(self.stageGap[i]) + "\n")  # Inputting row and stage gaps
            self.f.write(
                str(0.00) + "  " + str(0.00) + "\n")  # Inputting blockage factors
            self.f.write(str(self.eta_guess[i]) + "\n")  # Inputting guess values of stage isentropic efficiency
            if self.machineType == 'C':
                self.f.write(str(self.Dev_R_TE[i]) + "  " + str(self.Dev_S_TE[i]) + "\n")  # Inputting deviation angles
                self.f.write(str(self.Dev_R_LE[i]) + "  " + str(self.Dev_S_LE[i]) + "\n")  # Inputting incidence angles
            else:
                self.f.write(str(self.Dev_S_TE[i]) + "  " + str(self.Dev_R_TE[i]) + "\n")  # Inputting deviation angles
                self.f.write(str(self.Dev_S_LE[i]) + "  " + str(self.Dev_R_LE[i]) + "\n")  # Inputting incidence angles
            self.f.write(str(self.Twist[i]) + "\n")  # Inputting blade twist value
            self.f.write("N\n")
            if self.machineType == 'C':
                self.f.write(str(self.QO_LE_R[i]) + "  " + str(self.QO_TE_R[i]) + "\n")  # Inputting rotor QO angles
                self.f.write(str(self.QO_LE_S[i]) + "  " + str(self.QO_TE_S[i]) + "\n")  # Inputting stator QO angles
            else:
                self.f.write(str(self.QO_LE_S[i]) + "  " + str(self.QO_TE_S[i]) + "\n")  # Inputting rotor QO angles
                self.f.write(str(self.QO_LE_R[i]) + "  " + str(self.QO_TE_R[i]) + "\n")  # Inputting stator QO angles
        self.f.write("N\nY\n")
        for i in range(self.n_stage):
            self.f.write("Y\nY\n")
        self.f.close()

    def Variables(self, R, phi, psi, r_m, chords, rowGap, Twist, QO_LE, QO_TE, tc_m, x_tm):
        self.R = R * np.ones([self.n_stage, 1])
        self.phi = phi * np.ones([self.n_stage, 1])
        self.psi = psi * np.ones([self.n_stage, 1])
        self.r_m = r_m * np.ones([self.n_stage, 1])
        self.chords = chords * np.ones([self.n_stage, 2])
        self.rowGap = rowGap * np.ones([self.n_stage, 1])
        self.Twist = Twist*np.ones([self.n_stage, 1])
        self.QO_LE = QO_LE*np.ones([self.n_stage, 2])
        self.QO_TE = QO_TE*np.ones([self.n_stage, 2])
        self.tc_m = tc_m * np.ones([self.n_stage, 3])
        self.x_tm = x_tm * np.ones([self.n_stage, 3])
        self.eta_guess = 1.0*np.ones([self.n_stage, 1])



    def machineDefinition(self):

        self.flowPath = 'AXI'
        self.designPoint = 'M'
        self.inType = 'A'
        self.radType = 'A'


    def ParabladeWriter(self):
        n_rows = len(self.theta_in[0, :])
        if self.Dimension == 2:
            n_start = 1
            n_end = 2
            n_sec = 1
            CASCADE_TYPE = "LINEAR"
            sec_count = 30
            point_count = 30
        else:
            n_start = 0
            n_end = 3
            n_sec = 3
            CASCADE_TYPE = "ANNULAR"
            sec_count = 30
            point_count = 30

        HOME = os.environ["M2BFM"]

        template_dir = HOME + "templates/"

        N_b = np.zeros(n_rows)
        D1 = np.zeros(n_rows)
        D2 = np.zeros(n_rows)
        R_LE = np.zeros(n_rows)
        R_TE = np.zeros(n_rows)
        T = np.zeros([6, n_rows])

        if self.machineType == 'C':

            for i in range(self.n_stage):
                N_b[2*i] += self.N_b_R[i]
                N_b[2*i + 1] += self.N_b_S[i]
                D1[2*i] += self.D_1_R[i]
                D1[2*i + 1] += self.D_1_S[i]
                D2[2 * i] += self.D_2_R[i]
                D2[2 * i + 1] += self.D_2_S[i]
                R_LE[2 * i] += self.t_LE_R[i]
                R_LE[2 * i + 1] += self.t_LE_S[i]
                R_TE[2 * i] += self.t_TE_R[i]
                R_TE[2 * i + 1] += self.t_TE_S[i]
                T[:, 2 * i] += self.T_R[i, :]
                T[:, 2 * i + 1] += self.T_S[i, :]
            self.N_b = N_b
        else:
            for i in range(self.n_stage):
                N_b[2 * i] += self.N_b_S[i]
                N_b[2 * i + 1] += self.N_b_R[i]
                D1[2 * i] += self.D_1_S[i]
                D1[2 * i + 1] += self.D_1_R[i]
                D2[2 * i] += self.D_2_S[i]
                D2[2 * i + 1] += self.D_2_R[i]
                R_LE[2 * i] += self.t_LE_S[i]
                R_LE[2 * i + 1] += self.t_LE_R[i]
                R_TE[2 * i] += self.t_TE_S[i]
                R_TE[2 * i + 1] += self.t_TE_R[i]
                T[:, 2 * i] += self.T_S[i, :]
                T[:, 2 * i + 1] += self.T_R[i, :]
            self.N_b = N_b
        for i in range(n_rows):
            os.system("cp " + template_dir + "/template_turbine.cfg ./Bladerow_" + str(i+1) + ".cfg")

            os.system("sed -i 's/CAS_type/"+CASCADE_TYPE+"/g' Bladerow_"+str(i+1)+ ".cfg")
            os.system("sed -i 's/N_sec/" + str(sec_count) + "/g' Bladerow_" + str(i + 1) + ".cfg")
            os.system("sed -i 's/N_point/" + str(point_count) + "/g' Bladerow_" + str(i + 1) + ".cfg")
            os.system("sed -i 's/N_dim/"+str(int(self.Dimension))+"/g' Bladerow_"+str(i+1) + ".cfg")

            os.system("sed -i 's/N_blade/" + str(int(self.N_b[i])) + "/g' Bladerow_"+str(i+1)+ ".cfg")
            tan_stagger = np.transpose(np.array(0.5 * (np.tan(self.theta_in[n_start:n_end, i]*np.pi/180.0) +
                                                       np.tan(self.theta_out[n_start:n_end, i]*np.pi/180.0))))
            stagger = []
            for m in range(len(tan_stagger)):
                stagger.append(math.atan(tan_stagger[m])*180.0/np.pi)
            os.system("sed -i 's/STAGGER/"+", ".join([str(s) for s in stagger])+"/g' Bladerow_"+str(i+1)+".cfg")

            X_hub = [0.75*self.X_LE[0, i] + 0.25*self.X_TE[0, i], 0.75*self.X_TE[0, i] + 0.25*self.X_LE[0, i]]
            os.system("sed -i 's/X_HUB/" + ", ".join([str(s) for s in X_hub]) + "/g' Bladerow_" + str(i + 1) + ".cfg")
            Z_hub = [0.75*self.Z_LE[0, i] + 0.25*self.Z_TE[0, i], 0.75*self.Z_TE[0, i] + 0.25*self.Z_LE[0, i]]
            os.system("sed -i 's/Z_HUB/" + ", ".join([str(s) for s in Z_hub]) + "/g' Bladerow_" + str(i + 1) + ".cfg")
            X_shroud = [0.75 * self.X_LE[-1, i] + 0.25 * self.X_TE[-1, i], 0.75 * self.X_TE[-1, i] + 0.25 *
                        self.X_LE[-1, i]]
            os.system("sed -i 's/X_SHROUD/" + ", ".join([str(s) for s in X_shroud]) +
                      "/g' Bladerow_" + str(i + 1) + ".cfg")
            Z_shroud = [0.75 * self.Z_LE[-1, i] + 0.25 * self.Z_TE[-1, i], 0.75 * self.Z_TE[-1, i] + 0.25 *
                        self.Z_LE[-1, i]]
            os.system("sed -i 's/Z_SHROUD/" + ", ".join([str(s) for s in Z_shroud]) +
                      "/g' Bladerow_" + str(i + 1) + ".cfg")
            X_le = np.transpose(self.X_LE[n_start:n_end, i])
            os.system("sed -i 's/X_LE/" + ", ".join([str(s) for s in X_le]) + "/g' Bladerow_" + str(i + 1) + ".cfg")
            Z_le = np.transpose(self.Z_LE[:, i])
            Z_te = np.transpose(self.Z_TE[:, i])
            if (i+1) % 2 == 0:
                os.system("sed -i 's/Z_LE/" + ", ".join([str(s) for s in Z_le]) + "/g' Bladerow_" + str(i + 1) + ".cfg")
                os.system("sed -i 's/Z_TE/" + ", ".join([str(s) for s in Z_te]) + "/g' Bladerow_" + str(i + 1) + ".cfg")
            else:
                z_le = []
                z_te = []
                for j in range(len(Z_le)):
                    z_le.append(Z_le[0] + 0.99 * (Z_le[j] - Z_le[0]))
                    z_te.append(Z_te[0] + 0.99 * (Z_te[j] - Z_te[0]))
                os.system("sed -i 's/Z_LE/" + ", ".join([str(s) for s in z_le]) + "/g' Bladerow_" + str(i + 1) + ".cfg")
                os.system("sed -i 's/Z_TE/" + ", ".join([str(s) for s in z_te]) + "/g' Bladerow_" + str(i + 1) + ".cfg")
            X_te = np.transpose(self.X_TE[n_start:n_end, i])
            os.system("sed -i 's/X_TE/" + ", ".join([str(s) for s in X_te]) + "/g' Bladerow_" + str(i + 1) + ".cfg")


            for n in range(6):
                os.system("sed -i 's/T"+str(n+1)+"/"+", ".join([str(T[n, i]) for j in range(n_sec)]) +"/g' Bladerow_" + str(i + 1) + ".cfg")
            os.system("sed -i 's/D1/"+", ".join([str(D1[i]) for j in range(n_sec)]) + "/g' Bladerow_" + str(i + 1) + ".cfg")
            os.system("sed -i 's/D2/" + ", ".join([str(D2[i]) for j in range(n_sec)]) + "/g' Bladerow_" + str(i + 1) + ".cfg")
            os.system("sed -i 's/R_LE/" + ", ".join([str(R_LE[i]) for j in range(n_sec)]) + "/g' Bladerow_" + str(i + 1) + ".cfg")
            os.system("sed -i 's/R_TE/" + ", ".join([str(R_TE[i]) for j in range(n_sec)]) + "/g' Bladerow_" + str(i + 1) + ".cfg")
            os.system("sed -i 's/THETA_IN/" + ", ".join([str(d) for d in self.theta_in[n_start:n_end, i]]) +
                      "/g' Bladerow_" + str(i + 1) + ".cfg")
            os.system("sed -i 's/THETA_OUT/" + ", ".join([str(d) for d in self.theta_out[n_start:n_end, i]]) +
                      "/g' Bladerow_" + str(i + 1) + ".cfg")
