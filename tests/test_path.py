from cgeom import clib
import numpy as np


def test_geom2d_line_segment_from_start_length():
    seg = clib.geom2d_line_segment_from_start_length(1.0, 2.0, 0.5, 0.5, 2.0)
    assert seg.type == 0
    assert seg.data[0] == 1.0
    assert seg.data[1] == 2.0
    assert seg.data[2] == 1.0 + 2.0 * 0.5
    assert seg.data[3] == 2.0 + 2.0 * 0.5

def test_geom2d_line_segment_get_length():
    seg = clib.geom2d_line_segment_from_start_length(0.0, 0.0, 1.0, 0.0, 5.0)
    length = clib.geom2d_segment_get_length(seg)
    assert length == 5.0

def test_geom2d_line_segment_get_points_at_steps():
    seg = clib.geom2d_line_segment_from_start_end(0.0, 0.0, 3.0, 4.0)
    steps = np.array([0.0, 2.5, 5])
    points = clib.geom2d_line_segment_get_points_at_steps(seg, steps)
    assert len(points['x']) == 3
    assert points[0][0] == 0.0 and points[0][1] == 0.0
    assert points[1][0] == 1.5 and points[1][1] == 2.0
    assert points[2][0] == 3.0 and points[2][1] == 4.0

def test_geom2d_arc_segment_get_length():
    seg = clib.geom2d_arc_segment_from_center_radius_angles(0.0, 0.0, 1.0, 0.0, np.pi/2)
    length = clib.geom2d_segment_get_length(seg)
    assert length == np.pi/2

def test_geom2d_arc_segment_get_points_at_steps():
    seg = clib.geom2d_arc_segment_from_center_radius_angles(0.0, 0.0, 1.0, 0.0, np.pi/2)
    steps = np.array([0.0, np.pi/4, np.pi/2])
    points = clib.geom2d_arc_segment_get_points_at_steps(seg, steps)
    assert len(points['x']) == 3
    assert np.isclose(points[0][0], 1.0) and np.isclose(points[0][1], 0.0)
    assert np.isclose(points[1][0], np.sqrt(2)/2) and np.isclose(points[1][1], np.sqrt(2)/2)
    assert np.isclose(points[2][0], 0.0) and np.isclose(points[2][1], 1.0)

def test_geom2d_segments_from_rectangle():
    halfwidth = 2.0
    halfheight = 1.0
    segments = clib.geom2d_segments_from_rectangle(halfwidth, halfheight)
    assert len(segments) == 4
    # Check first segment
    seg0 = segments[0]
    assert seg0['type'] == 0  # Line segment
    assert seg0['data'][0] == -2.0 and seg0['data'][1] == -1.0
    assert seg0['data'][2] == 2.0 and seg0['data'][3] == -1.0

def test_geom2d_segments_from_circle():
    radius = 1.0
    segments = clib.geom2d_segments_from_circle(radius)
    assert len(segments) > 0
    # Check that all segments are arc segments
    for seg in segments:
        assert seg['type'] == 1  # Arc segment
