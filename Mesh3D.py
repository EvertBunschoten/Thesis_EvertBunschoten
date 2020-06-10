import gmsh
import numpy as np
import os

class Gmesh3D:
    wedge = 10
    n_sec = 8
    pointID = None
    xCoords = None
    yCoords = None
    zCoords = None
    lC = None

    lineID = None
    lineStart = None
    lineEnd = None

    Rev = None

    def __init__(self, Meangen):
        self.M = Meangen

        self.model = gmsh.model
        self.factory = self.model.geo
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 1)

        self.model.add("3DBFM")
        self.makePoints()
        self.makeLines()
        self.makePlane()
        self.revolve()
        self.nameBoundaries()

        gmsh.model.geo.synchronize()
        gmsh.model.mesh.generate(3)
        gmsh.write("3DBFM.su2")
        # gmsh.fltk.run()

    def nameBoundaries(self):
        periodic_1 = self.model.addPhysicalGroup(2, [1])
        self.model.setPhysicalName(2, periodic_1, "periodic_1")
        periodic_2 = self.model.addPhysicalGroup(2, [self.Rev[0][1]])
        self.model.setPhysicalName(2, periodic_2, "periodic_2")
        self.model.addPhysicalGroup(3, [1], 1)
        self.model.setPhysicalName(3, 1, "FlowField")

        i = 2
        inlet = self.model.addPhysicalGroup(2, [self.Rev[i][1]])
        self.model.setPhysicalName(2, inlet, "inlet")
        i += 1
        shroud_list = []
        for j in range(i, i + 3 * self.M.n_stage + 2):
            shroud_list.append(self.Rev[j][1])
            i += 1
        shroud = self.model.addPhysicalGroup(2, shroud_list)
        self.model.setPhysicalName(2, shroud, "shroud")
        outlet = self.model.addPhysicalGroup(2, [self.Rev[i][1]])
        self.model.setPhysicalName(2, outlet, "outlet")
        i += 1
        hub_list = []
        for j in range(i, i + 3 * self.M.n_stage + 2):
            hub_list.append(self.Rev[j][1])
            i += 1
        hub = self.model.addPhysicalGroup(2, hub_list)
        self.model.setPhysicalName(2, hub, "hub")

    def revolve(self):
        self.Rev = self.factory.revolve([(2, 1)], 0, 0, 0, 1, 0, 0, self.wedge*np.pi/180, [self.n_sec])

    def makePlane(self):
        for i in range(len(self.pointID)):
            self.factory.addPoint(self.xCoords[i], self.yCoords[i], self.zCoords[i], self.lC[i], self.pointID[i])
        for i in range(len(self.lineID)):
            self.factory.addLine(self.lineStart[i], self.lineEnd[i], self.lineID[i])
        self.factory.addCurveLoop(self.lineID, 1)
        self.factory.addPlaneSurface([1], 1)

    def makeLines(self):
        i_line = 1
        start_line = []
        end_line = []
        lines = []
        line_names = []

        start_line.append(self.pointID[0])
        end_line.append(self.pointID[1])
        lines.append(i_line)
        i_line += 1
        for i in range(1, len(self.pointID)-1, 2):
            start_line.append(self.pointID[i])
            end_line.append(self.pointID[i+2])
            lines.append(i_line)
            i_line += 1
        start_line.append(self.pointID[-1])
        end_line.append(self.pointID[-2])
        lines.append(i_line)
        i_line += 1

        for i in range(len(self.pointID)-2, 0, -2):
            start_line.append(self.pointID[i])
            end_line.append(self.pointID[i-2])
            lines.append(i_line)
            i_line += 1

        self.lineID = lines
        self.lineStart = start_line
        self.lineEnd = end_line

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
        L_c = []

        points = []
        i_point = 1
        X_p.append(X_LE[0, 0] - 2 * self.M.chord_R[0])
        Y_p.append(R_LE[0, 0] * np.sin(0.5 * self.wedge * np.pi / 180))
        Z_p.append(R_LE[0, 0] * np.cos(0.5 * self.wedge * np.pi / 180))
        points.append(i_point)
        L_c.append(2 * self.M.chord_R[0] / 30)
        i_point += 1
        X_p.append(X_LE[-1, 0] - 2 * self.M.chord_R[0])
        Y_p.append(R_LE[-1, 0] * np.sin(0.5 * self.wedge * np.pi / 180))
        Z_p.append(R_LE[-1, 0] * np.cos(0.5 * self.wedge * np.pi / 180))
        points.append(i_point)
        L_c.append(2 * self.M.chord_R[0] / 30)
        i_point += 1
        for j in range(n_rows):
                X_p.append(X_LE[0, j])
                Y_p.append(R_LE[0, j] * np.sin(0.5 * self.wedge * np.pi/180))
                Z_p.append(R_LE[0, j] * np.cos(0.5 * self.wedge * np.pi / 180))
                points.append(i_point)
                L_c.append((X_TE[0, j] - X_LE[0, j]) / 30)
                i_point += 1

                X_p.append(X_LE[-1, j])
                Y_p.append(R_LE[-1, j] * np.sin(0.5 * self.wedge * np.pi / 180))
                Z_p.append(R_LE[-1, j] * np.cos(0.5 * self.wedge * np.pi / 180))
                points.append(i_point)
                L_c.append((X_TE[-1, j] - X_LE[-1, j]) / 30)
                i_point += 1

                X_p.append(X_TE[0, j])
                Y_p.append(R_TE[0, j] * np.sin(0.5 * self.wedge * np.pi / 180))
                Z_p.append(R_TE[0, j] * np.cos(0.5 * self.wedge * np.pi / 180))
                points.append(i_point)
                L_c.append((X_TE[0, j] - X_LE[0, j]) / 30)
                i_point += 1

                X_p.append(X_TE[-1, j])
                Y_p.append(R_TE[-1, j] * np.sin(0.5 * self.wedge * np.pi / 180))
                Z_p.append(R_TE[-1, j] * np.cos(0.5 * self.wedge * np.pi / 180))
                points.append(i_point)
                L_c.append((X_TE[-1, j] - X_LE[-1, j]) / 30)
                i_point += 1

        X_p.append(X_p[-2] + 2 * self.M.chord_S[-1])
        Y_p.append(Y_p[-2])
        Z_p.append(Z_p[-2])
        points.append(i_point)
        L_c.append(2 * self.M.chord_S[-1] / 30)
        i_point += 1
        X_p.append(X_p[-2] + 2 * self.M.chord_S[-1])
        Y_p.append(Y_p[-2])
        Z_p.append(Z_p[-2])
        points.append(i_point)
        L_c.append(2 * self.M.chord_S[-1] / 30)
        i_point += 1

        self.pointID = points
        self.xCoords = X_p
        self.yCoords = Y_p
        self.zCoords = Z_p
        self.lC = L_c
        # plt.plot(X_p[:-1:2], Z_p[:-1:2], 'k')
        # plt.plot(X_p[1::2], Z_p[1::2], 'k')
        # plt.show()