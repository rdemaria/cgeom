from tinycwrap import CModule
import numpy as np



cm = CModule("base.c","base.h")

assert cm.geom2d_norm(3,4)==5.0
assert np.allclose(cm.geom2d_points_distance(1,2,4,6),5.0)

p=cm.G2DPoint(x=3,y=4)
print(p._data)