import cgeom as cg
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv


def comulative_distance(p_aa,p_bb):
    dx=p_aa['x']-p_bb['x']
    dy=p_aa['y']-p_bb['y']
    dd=np.sum(dx*dx+dy*dy)
    return dd

def distance_points(p_aa,p_bb):
    return np.sqrt((p_aa['x']-p_bb['x'])**2+(p_aa['y']-p_bb['y'])**2)

#def best_shift_points(p_aa,p_bb):
#    best_dd=comulative_distance(p_aa,p_bb)
#    best_shift=0
#    for shift in range(len(p_bb)):
#        p_bb_shifted={}
#        p_bb_shifted['x']=np.roll(p_bb['x'],shift)
#        p_bb_shifted['y']=np.roll(p_bb['y'],shift)
#        dd=comulative_distance(p_aa,p_bb_shifted)
#        if dd<best_dd:
#            best_dd=dd
#            best_shift=shift
#    return best_shift

def best_shift_points(p_aa,p_bb):
    best_shift=0
    best_dd=distance_points(p_aa[0],p_bb[0])
    for shift in range(len(p_bb)):
        dd=distance_points(p_aa[0],p_bb[shift])
        if dd<best_dd:
            best_dd=dd
            best_shift=shift
    return -best_shift




aa=cg.Path2D.from_rectellipse(1,1,1.2,1.2)
bb=cg.Path2D.from_circle(1.5)


nsteps=155
c_aa=aa.get_corner_steps()[1:-1]/aa.length
c_bb=bb.get_corner_steps()[1:-1]/bb.length

ss=np.arange(nsteps)/nsteps
#ss=np.r_[ss,c_aa,c_bb]
#ss=np.sort(ss)
p_aa=aa.get_points_at_steps(ss*aa.length)
p_bb=bb.get_points_at_steps(ss*bb.length)


best_shift=best_shift_points(p_aa,p_bb)

#best_shift=-1
#best_shift+=1
p_bb_aligned={}
p_bb_aligned['x']=np.roll(p_bb['x'],best_shift)
p_bb_aligned['y']=np.roll(p_bb['y'],best_shift)

plt.clf()
aa.plot()
c_aa_pp=aa.get_points_at_steps(aa.get_corner_steps()[:-1])
plt.plot(c_aa_pp['x'],c_aa_pp['y'],'ro')
bb.plot()
x1=p_aa['x']
x2=p_bb_aligned['x']
y1=p_aa['y']
y2=p_bb_aligned['y']
out=[]
for xx1,xx2,yy1,yy2 in zip(x1,x2,y1,y2):
    plt.plot([xx1,xx2],[yy1,yy2],'k--',alpha=0.3)
    out.append( ((xx1+xx2)/2,(yy1+yy2)/2) )
out.append(out[0])
out=np.array(out)
plt.plot(out[:,0],out[:,1],'g-')
plt.axis('equal')
plt.show()

# render in 3
def plot_surface(x1,y1,x2,y2):
    n=len(x1)
    points = np.zeros((n*2, 3))
    for i in range(n):
        points[i, 0] = x1[i]
        points[i, 1] = y1[i]
        points[i, 2] = 0.0
        points[n + i, 0] = x2[i]
        points[n + i, 1] = y2[i]
        points[n + i, 2] = 0.0

    # Create the faces (quads)
    faces = []
    for i in range(n):
        next_i = (i + 1) % n
        faces.append([4, i, next_i, n + next_i, n + i])  # Quad face

    # Convert faces to a flat array
    faces_flat = np.hstack(faces)

    # Create the PolyData object
    mesh = pv.PolyData(points, faces_flat)

    # Plot the surface
    plotter = pv.Plotter()
    plotter.add_mesh(mesh, color='lightblue', show_edges=True, opacity=0.7)
    plotter.show()

plot_surface(p_aa['x'],p_aa['y'],p_bb_aligned['x'],p_bb_aligned['y'])


