# TODO





# Naming conventiosn

methods: 
- geom2d_
- geom3d_

2d profile functions:
- generate points for plotting
- generate points for transitions

profile list:
- circle
- ellipse
- rectangle
- octagon
- polygon
- rectellipse
- racetrack
- Bcurve (straight and arcs)
- SVG Path (from string)

2d primitives
- 2d_intersect_ray_polygon
- 2d_distance_point_to_polygon

3d primitives
- 3d_transform_points
- 3d_intersect_segment_plane





Conventions

Use double everywhere for geometry.
Use size_t for sizes and counts.

AoS (Array of Structures) for collections of points:
2D: points[2*i + 0] = x_i, points[2*i + 1] = y_i
3D: points[3*i + 0] = x_i, points[3*i + 1] = y_i, points[3*i + 2] = z_i

Row-major for matrices:
M[row*cols + col]

g2d_
g3d_

p or p2 point
pts, pts points


#define GEOM2D_POINT(pts, i) ((pts) + 2*(i))




