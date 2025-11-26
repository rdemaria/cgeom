import cgeom as cg
from cgeom import clib

cg.Path2D.from_rectellipse(1,1,1.2,1.2)

cg.Path2D.from_circle(1.2).plot()
cg.Path2D.from_rectangle(1,1).plot()
cg.Path2D.from_rectellipse(1,1,1.2,1.2).plot()

cg.Path2D.from_ellipse(1.2,0.8).plot()
p=cg.Path2D.from_circle(0.5).plot()

cg.geom2d_path_get_length(p)
steps=cg.geom2d_path_get_steps(p,0.1)
pts=cg.geom2d_path_get_points_at_steps(p, steps)
