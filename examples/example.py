import matplotlib.pyplot as plt

from tinycwrap import CModule


cg = CModule("cgeom.c")

origin=cg.G2DPoint()
p_circ = cg.geom2d_circle_get_points(origin,1)

# make axes with equal aspect ratio
fig, ax = plt.subplots()
ax.set_aspect('equal')
ax.plot(p_circ['x'], p_circ['y'])


p_rect = cg.geom2d_rect_get_points(-1,-1,1,1)
ax.plot(p_rect['x'], p_rect['y'])
