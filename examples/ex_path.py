from tinycwrap import CModule
import numpy as np



cm = CModule("base.c", "path.c")

p=cm.geom2d_rectangle_to_path(1,1)
cm.geom2d_path_get_length(p)
steps=cm.geom2d_path_get_steps(p,0.1)
pts=cm.geom2d_path_get_points_at_steps(p, steps)