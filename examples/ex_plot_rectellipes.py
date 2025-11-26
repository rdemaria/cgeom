import cgeom as cg
from cgeom import clib

hw,hh,rr=1,1,1.2

cg.Path2D.from_circle(rr).plot()
cg.Path2D.from_rectangle(hw,hh).plot()
cg.Path2D.from_rectellipse(hw,hh,rr,rr).plot()
