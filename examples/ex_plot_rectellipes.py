import cgeom as cg
#from cgeom import clib
import matplotlib.pyplot as plt


args=[
    (0.9,0.7,1.0,1.0,"4 corners"),
    (1.2,0.7,1.0,1.0,"flat top/bottom"),
    (0.7,1.2,1.0,1.0,"flat sides"),
    (1.2,1.2,1.0,1.0, "rect outside"),
    (0.5,0.5,1.0,1.0, "rect inside"),
    (1.0,1.0,1.0,1.0, "tangent" ),
    (0.9,0.9,1.2,1.0,"4 corners ellipse"),
    (0.9,0.9,1.0,1.2,"4 corners ellipse"),
    (0.9,1.2,0.9,1.2,"4 tangent ellipse" ),
]

fig,axes=plt.subplots(3,3, num="rectellipse", figsize=(9,9))

for ax,arg in zip(axes.flatten(),args):
    hw,hh,rx,ry,lbl=arg
    cg.Path2D.from_ellipse(rx, ry).plot(ax=ax)
    cg.Path2D.from_rectangle(hw,hh).plot(ax=ax)
    p=cg.Path2D.from_rectellipse(hw,hh,rx,ry).plot(ax=ax)
    ax.set_aspect('equal')
    ax.text(0.5, 0.95, lbl, ha='center', va='top', transform=ax.transAxes)
    ax.text(0.5, 0.15, f"{p.len_segments} segments", ha='center', va='top', transform=ax.transAxes)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(-2,2)
    ax.set_ylim(-2,2)


#plt.tight_layout()


#cg.Path2D.from_rectellipse(1.2,0.9,rr,rr).plot()