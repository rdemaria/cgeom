from typing import Callable, Dict, Union, Type, Tuple
from itertools import pairwise

import numpy as np

from xtrack.beam_elements import apertures
from xtrack.progress_indicator import progress

from . import clib
from .path import Path2D


LimitTypes = Union[
    apertures.LimitRect,
    apertures.LimitEllipse,
    apertures.LimitRectEllipse,
    apertures.LimitRacetrack,
    apertures.LimitPolygon,
]


def transform_matrix(dx, dy, ds, theta, phi, psi):
    """Generate a 3D transformation matrix.

    Parameters
    ----------
    dx, dy, ds : float
        Shifts in x, y, and s directions
    theta : float
        Rotation around the y-axis (positive s to x) in radians
    phi
        Rotation around the x-axis (positive s to y) in radians
    psi
        Rotation around the s-axis (positive y to x) in radians
    """
    s_phi, c_phi = np.sin(phi), np.cos(phi)
    s_theta, c_theta = np.sin(theta), np.cos(theta)
    s_psi, c_psi = np.sin(psi), np.cos(psi)
    matrix = np.array(
        [
            [-s_phi * s_psi * s_theta + c_psi * c_theta,
                -c_psi * s_phi * s_theta - c_theta * s_psi, c_phi * s_theta, dx],
            [c_phi * s_psi, c_phi * c_psi, s_phi, dy],
            [-c_theta * s_phi * s_psi - c_psi * s_theta,
                -c_psi * c_theta * s_phi + s_psi * s_theta, c_phi * c_theta, ds],
            [0, 0, 0, 1],
        ]
    )
    return matrix


class Aperture:
    def __init__(self, line):
        self.line = line
        self._twiss = line.twiss()

        s_positions, twiss_data, apertures, types = self.build_aperture_data()
        self.s_positions = s_positions
        self.twiss_data = twiss_data
        self.apertures = apertures
        self.types = types

    def build_aperture_data(self):
        line_table = self.line.get_table()
        s_positions, twiss_data, paths, types = [], [], [], []
        names = self.line.element_names

        for idx, name in progress(enumerate(names), desc="Building aperture data", total=len(names)):
            element = self.line.element_dict[name]
            if not isinstance(element, LimitTypes):
                continue

            aperture = self.get_aperture_data_from_limit_element(name)
            if element.transformations_active:
                self.apply_transform(aperture, element)
            twiss = self.get_twiss_data_at_element(idx)

            s_positions.append(line_table.s[idx])
            twiss_data.append(twiss)
            paths.append(aperture)
            types.append(type(element))

        return s_positions, twiss_data, paths, types

    def get_aperture_data_from_limit_element(self, name):
        """
        Get 2d points from limit element

        Parameters
        ----------
        name : str
            Name of the limit element
        Returns: G2DBeamApertureData
            2D beam aperture data
        """
        element: LimitTypes = self.line.element_dict[name]
        path, centre_x, centre_y = self.CONVERTER_FUNCTIONS[type(element)](element)
        points = clib.geom2d_path_get_points(path)
        points['x'] += centre_x
        points['y'] += centre_y
        beam_aperture = clib.G2DBeamApertureData(
            points=points,
            n_points=len(points),
            tol_r=0,
            tol_x=0,
            tol_y=0,
        )
        return beam_aperture

    def get_twiss_data_at_element(self, name):
        """
        Get twiss data at element 
        """
        row = self._twiss.rows[name]
        c_twiss_data = clib.G2DTwissData(
            x=row['x'][0],
            y=row['y'][0],
            betx=row['betx'][0],
            bety=row['bety'][0],
            dx=row['dx'][0],
            dy=row['dy'][0],
            delta=row['delta'][0],
            gamma=row['gamma0'],
        )
        return c_twiss_data

    @staticmethod
    def get_path2d_from_limit_rect(element: apertures.LimitRect) -> Tuple[Path2D, float, float]:
        half_width = (element.max_x - element.min_x) / 2
        half_height = (element.max_y - element.min_y) / 2
        x = (element.min_x + element.max_x) / 2
        y = (element.min_y + element.max_y) / 2
        return Path2D.from_rectangle(half_width, half_height), x, y

    @staticmethod
    def get_path2d_from_limit_ellipse(element: apertures.LimitEllipse) -> Tuple[Path2D, float, float]:
        rx = element.a
        ry = element.b
        return Path2D.from_ellipse(rx, ry), 0, 0

    @staticmethod
    def get_path2d_from_limit_rect_ellipse(element: apertures.LimitRectEllipse) -> Tuple[Path2D, float, float]:
        max_x = element.max_x
        max_y = element.max_y
        rx = element.a
        ry = element.b
        return Path2D.from_rectellipse(max_x, max_y, rx, ry), 0, 0

    @staticmethod
    def get_path2d_from_limit_racetrack(element: apertures.LimitRacetrack) -> Tuple[Path2D, float, float]:
        half_width = (element.max_x - element.min_x) / 2
        half_height = (element.max_y - element.min_y) / 2
        x = (element.min_x + element.max_x) / 2
        y = (element.min_y + element.max_y) / 2
        rx = element.a
        ry = element.b
        return Path2D.from_racetrack(half_width, half_height, rx, ry), x, y

    @staticmethod
    def get_path2d_from_limit_polygon(element: apertures.LimitPolygon) -> Tuple[Path2D, float, float]:
        xs = element.x_vertices + [element.x_vertices[0]]
        ys = element.y_vertices + [element.y_vertices[0]]
        segments = []
        for (x0, y0), (x1, y1) in pairwise(zip(xs, ys)):
            segment = clib.geom2d_line_segment_from_start_end(x0, y0, x1, y1)
            segments.append(segment)

        return Path2D(segments, len_segments=len(segments))

    def apply_transform(self, aperture, element):
        points = aperture.points
        matrix = transform_matrix(
            dx=element.shift_x,
            dy=element.shift_y,
            ds=element.shift_s,
            theta=element.rot_y_rad,
            phi=element.rot_x_rad,
            psi=element.rot_s_rad_no_frame,
        )
        point_vecs = np.vstack([
            points['x'],
            points['y'],
            np.zeros(len(points)),
            np.ones(len(points)),
        ])
        transformed_points = matrix @ point_vecs

        # Project back on 2D plane
        points['x'] = transformed_points[0]
        points['y'] = transformed_points[1]

    @property
    def extents(self):
        for aperture in self.apertures:
            points = aperture.points
            yield max(points['x']), min(points['y']), min(points['x']), max(points['y'])

    @property
    def extents_on_axes(self):
        for aperture in self.apertures:
            points = aperture.points
            xs, ys = self._find_axis_intersections(points['x'], points['y'])
            yield xs[-1], ys[0], xs[0], ys[-1]

    @staticmethod
    def _find_axis_intersections(x, y):
        """Finds intersections of a closed polygon with x = 0 and y = 0.

        Returns:
            x_axis_points: sorted array of x-values where polygon crosses y = 0
            y_axis_points: sorted array of y-values where polygon crosses x = 0
        """
        # X Axis
        mask_xaxis = (y[:-1] * y[1:] <= 0) & (y[:-1] != y[1:])

        t_xaxis = y[:-1][mask_xaxis] / (
                    y[:-1][mask_xaxis] - y[1:][mask_xaxis])
        x_xaxis = x[:-1][mask_xaxis] + t_xaxis * (
                    x[1:][mask_xaxis] - x[:-1][mask_xaxis])

        # Y Axis
        mask_yaxis = (x[:-1] * x[1:] <= 0) & (x[:-1] != x[1:])

        t_yaxis = x[:-1][mask_yaxis] / (
                    x[:-1][mask_yaxis] - x[1:][mask_yaxis])
        y_yaxis = y[:-1][mask_yaxis] + t_yaxis * (
                    y[1:][mask_yaxis] - y[:-1][mask_yaxis])

        return np.sort(x_xaxis), np.sort(y_yaxis)

    CONVERTER_FUNCTIONS: Dict[Type[LimitTypes], Callable[[LimitTypes], Path2D]] = {
        apertures.LimitRect: get_path2d_from_limit_rect,
        apertures.LimitEllipse: get_path2d_from_limit_ellipse,
        apertures.LimitRectEllipse: get_path2d_from_limit_rect_ellipse,
        apertures.LimitRacetrack: get_path2d_from_limit_racetrack,
        apertures.LimitPolygon: get_path2d_from_limit_polygon,
    }
