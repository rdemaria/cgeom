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




aa=cg.Path2D.from_rectellipse(0.8,0.8,0.9,0.9)
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

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
x1=np.r_[x1, x1[0]]
y1=np.r_[y1, y1[0]]
x2=np.r_[x2, x2[0]]
y2=np.r_[y2, y2[0]]
ax.fill_between(x1, y1,0, x2,y2, z2=1, alpha=0.5)
ax.plot(out[:,0],out[:,1], zs=0.5, zdir='z', label='interpolated', color='g')

def intersect_line_plane(line,plane):
    x1,y1,z1=line[0]
    x2,y2,z2=line[1]
    a,b,c,d=plane[0]
    t=(d - a*x1 - b*y1 - c*z1)/(a*(x2 - x1) + b*(y2 - y1) + c*(z2 - z1))
    x=x1 + t*(x2 - x1)
    y=y1 + t*(y2 - y1)
    z=z1 + t*(z2 - z1)
    return (x,y,z)

out=[]
for xx1,xx2,yy1,yy2 in zip(x1,x2,y1,y2):
    p1=(xx1,yy1,0)
    p2=(xx2,yy2,1)
    line=(p1,p2)
    plane=((0,0.3,1,0.5),)
    p_int=intersect_line_plane(line,plane)
    out.append(p_int)

out=np.array(out)
ax.plot(out[:,0],out[:,1],out[:,2], zdir='z', label='interpolated', color='r')