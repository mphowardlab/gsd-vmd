"""File courtesy of Joshua A. Anderson to test different precisions in GSD files."""

import gsd.hoomd

s = gsd.hoomd.Frame()
s.particles.N = 4
s.particles.types = ["A", "B"]
s.particles.typeid = [0, 0, 1, 1]
s.particles.position = [[0, 0, 0], [1, 1, 1], [-1, -1, -1], [1, -1, -1]]
s.configuration.box = [3, 3, 3, 0, 0, 0]
traj = gsd.hoomd.open(name="single.gsd", mode="w", precision="single")
traj.append(s)

traj = gsd.hoomd.open(name="double.gsd", mode="w", precision="double")
traj.append(s)

traj = gsd.hoomd.open(name="mixed.gsd", mode="w", precision="single")
traj.append(s)

s = gsd.hoomd.Frame()
s.particles.N = 4
s.particles.position = [[0, 0, 0], [2, 2, 2], [-2, -2, -2], [2, -2, -2]]
s.configuration.box = [6, 6, 6, 0, 0, 0]
traj._precision = "double"
traj.append(s)
