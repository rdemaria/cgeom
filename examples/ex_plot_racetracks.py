
import cgeom as cg
#from cgeom import clib
import matplotlib.pyplot as plt


args=[
    (1.0,1.0,0.5,"4 corners"),
    (1.0,0.5,0.5,"flat top/bottom"),
    (0.5,1.0,0.5,"flat sides"),
    (1.0,1.0,0.0, "just rect"),
    (0.5,0.5,0.5, "just circle"),
    (0.7,1.0,0.5, "4 corners" ),
]


fig,axes=plt.subplots(2,3, num="rectellipse", figsize=(9,6))

for ax,arg in zip(axes.flatten(),args):
    hw,hh,rr,lbl=arg
    cg.Path2D.from_circle(rr).plot(ax=ax)
    cg.Path2D.from_rectangle(hw,hh).plot(ax=ax)
    p=cg.Path2D.from_racetrack(hw,hh,rr,rr).plot(ax=ax)
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