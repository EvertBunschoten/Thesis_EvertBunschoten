import gmsh
import numpy as np
import matplotlib.pyplot as plt

class Gmesh3D:
    wedge = 10
    pointID = None
    xCoords = None
    yCoords = None
    zCoords = None

    def __init__(self, Meangen):
        self.M = Meangen

        self.model = gmsh.model
        self.factory = self.model.geo
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 1)

        self.model.add("3DBFM")

    def makePoints(self):
        X_LE = self.M.X_LE
        R_LE = self.M.Z_LE
        X_TE = self.M.X_TE
        R_TE = self.M.Z_TE
        n_rows = len(X_LE[0, :])
        n_stage = self.M.n_stage
        X_p = []
        Y_p = []
        Z_p = []
        points = []
        i_point = 1
        X_p.append(X_LE[0, 0] - 2 * self.M.chord_R[0])
        Y_p.append(R_LE[0, 0] * np.sin(0.5 * self.wedge * np.pi / 180))
        Z_p.append(R_LE[0, 0] * np.cos(0.5 * self.wedge * np.pi / 180))
        points.append(i_point)
        i_point += 1
        X_p.append(X_LE[-1, 0] - 2 * self.M.chord_R[0])
        Y_p.append(R_LE[-1, 0] * np.sin(0.5 * self.wedge * np.pi / 180))
        Z_p.append(R_LE[-1, 0] * np.cos(0.5 * self.wedge * np.pi / 180))
        points.append(i_point)
        i_point += 1
        for j in range(n_rows):
                X_p.append(X_LE[0, j])
                Y_p.append(R_LE[0, j] * np.sin(0.5 * self.wedge * np.pi/180))
                Z_p.append(R_LE[0, j] * np.cos(0.5 * self.wedge * np.pi / 180))
                points.append(i_point)
                i_point += 1

                X_p.append(X_LE[-1, j])
                Y_p.append(R_LE[-1, j] * np.sin(0.5 * self.wedge * np.pi / 180))
                Z_p.append(R_LE[-1, j] * np.cos(0.5 * self.wedge * np.pi / 180))
                points.append(i_point)
                i_point += 1

                X_p.append(X_TE[0, j])
                Y_p.append(R_TE[0, j] * np.sin(0.5 * self.wedge * np.pi / 180))
                Z_p.append(R_TE[0, j] * np.cos(0.5 * self.wedge * np.pi / 180))
                points.append(i_point)
                i_point += 1

                X_p.append(X_TE[-1, j])
                Y_p.append(R_TE[-1, j] * np.sin(0.5 * self.wedge * np.pi / 180))
                Z_p.append(R_TE[-1, j] * np.cos(0.5 * self.wedge * np.pi / 180))
                points.append(i_point)
                i_point += 1

        X_p.append(X_p[-2] + 2 * self.M.chord_S[-1])
        Y_p.append(Y_p[-2])
        Z_p.append(Z_p[-2])
        points.append(i_point)
        i_point += 1
        X_p.append(X_p[-1] + 2 * self.M.chord_S[-1])
        Y_p.append(Y_p[-1])
        Z_p.append(Z_p[-1])
        points.append(i_point)
        i_point += 1

        self.pointID = points
        self.xCoords = X_p
        self.yCoords = Y_p
        self.zCoords = Z_p

        plt.plot(X_p[:-1:2], Z_p[:-1:2], 'k')
        plt.plot(X_p[1::2], Z_p[1::2], 'k')
        plt.show()