
import os
import numpy as np

class meangenWriter:
    
    def __init__(self, Variables, Parameters):
        Pt_in = 2.507   # Inlet stagnation pressure[bar]
        Tt_in = 340     # Inlet stagnation temperature[K]
        R = 287         # Gas constant[J/kg/K]
        gamma = 1.4     # Specific heat ratio
        Omega = 13700   # Rotation speed[rpm]
        mdot = 88.27    # Mass flow rate[kg/s]

        # Possible design variables
        self.n_stage = n_stage    # Number of stages[integer]
        self.c_ax = phi*Omega*(np.pi/30)*r_m     # Axial velocity[m/s]
        self.BLOCK_LE = 0.002*np.ones([self.n_stage, 1])    # Leading edge blockage factor
        self.BLOCK_TE = 0.002*np.ones([self.n_stage, 1])    # Trailing edge blockage factor
        self.phi = 0.5*np.ones([self.n_stage, 1])   # Flow coefficient for each stage
        self.psi = 0.28*np.ones([self.n_stage, 1])  # Work coefficient for each stage
        self.R = 0.8*np.ones([self.n_stage, 1])     # Stage reaction
        self.EFF_GUESS = 0.67*np.ones([self.n_stage, 1])    # Guess of isentropic efficiency for each stage
        self.DEV_R = 2*np.ones([self.n_stage, 1])   # Rotor deviation angle[deg]
        self.DEV_S = 2*np.ones([self.n_stage, 1])   # Stator deviation angle[deg]
        self.INC_R = -1.5*np.ones([self.n_stage, 1])    # Rotor incidence angle[deg]
        self.INC_S = -1.5*np.ones([self.n_stage, 1])    # Stator incidence angle[deg]
        self.TWIST = 0.66*np.ones([self.n_stage, 1])    # Blade twist for controlled vortex design
        self.QO_LE_R = 90*np.ones([self.n_stage, 1])    # Rotor leading edge QO angle[deg]
        self.QO_TE_R = 90*np.ones([self.n_stage, 1])    # Rotor trailing edge QO angle[deg]
        self.QO_LE_S = 90*np.ones([self.n_stage, 1])    # Stator leading edge QO angle[deg]
        self.QO_TE_S = 90*np.ones([self.n_stage, 1])    # Stator trailing edge QO angle[deg]
        self.tc_MAX1 = 0.075*np.ones([self.n_stage, 1])     # Maximum thickness-to-chord ratio of first blade row
        self.xc_MAX1 = 0.5*np.ones([self.n_stage, 1])       # Chordwise maximum thickness location of first blade row
        self.tc_MAX2 = 0.075*np.ones([self.n_stage, 1])     # Maximum thickness-to-chord ratio of second blade row
        self.xc_MAX2 = 0.5*np.ones([self.n_stage, 1])       # Chordwise maximum thickness location of second blade row

        self.machineDefinition()    # Defining the machine type
        self.bladeDesign()          # Defining blade design

        # Defining meangen input file
        save_path = os.getcwd()
        filename = "meangen"
        completename = os.path.join(save_path, filename+".IN")

        self.f = open(completename, 'wt')       # Opening input file
        self.f.write(self.machineType+"\n")     # Defining machine type
        self.f.write(self.flowPath+"\n")        # Defining flow path
        self.f.write(str(R) + "  " + str(gamma) + "\n")     # Defining working fluid properties
        self.f.write(str(Pt_in) + "  " + str(Tt_in) + "\n") # Defining inlet conditions
        self.f.write(str(self.n_stage) + "\n")  # Defining number of stages
        self.f.write(str(self.designPoint) + "\n")  # Defining blade design point
        self.f.write(str(Omega) + "\n")     # Input of rotation speed
        self.f.write(str(mdot) + "\n")      # Input of target mass flow rate
        # Looping over each stage, inputting the duty coefficients and radius of each stage
        for i in range(0, self.n_stage):
            radius = self.c_ax/(self.phi[i, 0]*Omega*np.pi/30)  # Calculating mean radius
            if i > 0:
                self.f.write("N\nN\n")  # Enabling new input for next stage
            self.f.write("A\n")     # Choice for flow angle calculation
            # Inputting duty coefficients
            self.f.write(str(self.R[i, 0]) + "  " + str(self.phi[i, 0]) + "  " + str(self.psi[i, 0]) + "\n")
            self.f.write(str(self.radType) + "\n")  # Inputting radius input type
            self.f.write(str(radius) + "\n")        # Inputting mean radius
            self.f.write(str(self.rotorChord[i, 0]) + "  " + str(self.statorChord[i, 0]) + "\n")    # Inputting chords
            self.f.write(str(self.rowGap[i, 0]) + "  " + str(self.stageGap[i, 0]) + "\n")   # Inputting row and stage gaps
            self.f.write(str(self.BLOCK_LE[i, 0]) + "  " + str(self.BLOCK_TE[i, 0]) + "\n") # Inputting blockage factors
            self.f.write(str(self.EFF_GUESS[i, 0]) + "\n")      # Inputting guess values of stage isentropic efficiency
            self.f.write(str(self.DEV_R[i, 0]) + "  " + str(self.DEV_S[i, 0]) + "\n")   # Inputting deviation angles
            self.f.write(str(self.INC_R[i, 0]) + "  " + str(self.INC_S[i, 0]) + "\n")   # Inputting incidence angles
            self.f.write(str(self.TWIST[i, 0]) + "\n")      # Inputting blade twist value
            self.f.write(str(self.bladeRotation) + "\n")
            self.f.write(str(self.QO_LE_R[i, 0]) + "  " + str(self.QO_TE_R[i, 0]) + "\n")   # Inputting rotor QO angles
            self.f.write(str(self.QO_LE_S[i, 0]) + "  " + str(self.QO_TE_S[i, 0]) + "\n")   # Inputting stator QO angles
        self.f.write("N\nY\n")
        for i in range(0, self.n_stage):
            self.f.write("Y\n")
            # Inputting first blade row thickness-to-chord ratio and location
            self.f.write(str(self.tc_MAX1[i, 0]) + "  " + str(self.xc_MAX1[i, 0]) + "\n")
            self.f.write(str(self.tc_MAX1[i, 0]) + "  " + str(self.xc_MAX1[i, 0]) + "\n")
            self.f.write(str(self.tc_MAX1[i, 0]) + "  " + str(self.xc_MAX1[i, 0]) + "\n")
            self.f.write("Y\n")
            # Inputting second blade row thickness-to-chord ratio and location
            self.f.write(str(self.tc_MAX2[i, 0]) + "  " + str(self.xc_MAX2[i, 0]) + "\n")
            self.f.write(str(self.tc_MAX2[i, 0]) + "  " + str(self.xc_MAX2[i, 0]) + "\n")
            self.f.write(str(self.tc_MAX2[i, 0]) + "  " + str(self.xc_MAX2[i, 0]) + "\n")
        self.f.close()


    def machineDefinition(self):
        self.machineType = 'C'
        self.flowPath = 'AXI'
        self.designPoint = 'M'
        self.inType = 'A'
        self.radType = 'A'
        self.stageGap = 0.25*np.ones([self.n_stage, 1])
        self.rowGap = 0.5*np.ones([self.n_stage, 1])

    def bladeDesign(self):
        self.rotorChord = 0.05*np.ones([self.n_stage, 1])
        self.statorChord = 0.05 * np.ones([self.n_stage, 1])
        self.bladeRotation = 'N'

