import matplotlib.pyplot as plt

from .cwrapper import clib


class Path2D(clib.G2DPath):
    @classmethod
    def from_circle(cls, radius):
        """Create a circular Path2D centered at the origin.

        Parameters:
            radius: Radius of the circle.

        Returns:
            Path2D: The created circular path.
        """
        segments = clib.geom2d_segments_from_circle(radius)
        return cls(segments=segments, len_segments=len(segments))

    @classmethod
    def from_ellipse(cls, rx, ry):
        """Create an elliptical Path2D centered at the origin.

        Parameters:
            rx: Radius along x-axis of the ellipse.
            ry: Radius along y-axis of the ellipse.

        Returns:
            Path2D: The created elliptical path.
        """
        segments = clib.geom2d_segments_from_ellipse(rx, ry)
        return cls(segments=segments, len_segments=len(segments))

    @classmethod
    def from_rectellipse(cls, halfwidth, halfheight, rx, ry):
        """Create a rectangular-ellipse Path2D centered at the origin.

        Parameters:
            halfwidth: Half the width of the rectangle.
            halfheight: Half the height of the rectangle.
            rx: Radius along x-axis of the ellipse.
            ry: Radius along y-axis of the ellipse.
        """
        segments, _ = clib.geom2d_segments_from_rectellipse(
            halfwidth, halfheight, rx, ry
        )
        return cls(segments=segments, len_segments=len(segments))

    @classmethod
    def from_rectangle(cls, halfwidth, halfheight):
        """Create a rectangular Path2D centered at the origin.

        Parameters:
            halfwidth: Half the width of the rectangle.
            halfheight: Half the height of the rectangle.

        Returns:
            Path2D: The created rectangular path.
        """
        segments = clib.geom2d_segments_from_rectangle(halfwidth, halfheight)
        return cls(segments=segments, len_segments=len(segments))

    @classmethod
    def from_racetrack(cls, halfhside, halfvside, rx, ry):
        """Create a racetrack Path2D centered at the origin.

        Parameters:
            halfhside: Half the length of the straight horizontal sides.
            halfvside: Half the length of the straight vertical sides.
            rx: Radius along x-axis of the ellipse segments.
            ry: Radius along y-axis of the ellipse segments.

        Returns:
            Path2D: The created racetrack path.
        """
        segments, _ = clib.geom2d_segments_from_racetrack(
            halfhside, halfvside, rx, ry
        )
        return cls(segments=segments, len_segments=len(segments))
    
    @classmethod
    def from_octagon(cls, halfwidth, halfheight, halfdgap):
        """Create an octagonal Path2D centered at the origin.

        Parameters:
            halfwidth: Half the width of the octagon.
            halfheight: Half the height of the octagon.
            halfdgap: Half the gap along 45-degree edges.

        Returns:
            Path2D: The created octagonal path.
        """
        segments, ret = clib.geom2d_segments_from_octagon(halfwidth, halfheight, halfdgap)
        if ret == 0:
            raise ValueError("Invalid parameters for octagon")
        return cls(segments=segments, len_segments=len(segments))

    def get_corner_steps(self):
        """Get the corner points of the path.

        Returns:
            numpy.ndarray: Array of corner points.
        """
        return clib.geom2d_path_get_corner_steps(self)

    def get_points(self, ds_min):
        """Get points along the path with a minimum step size.

        Parameters:
            ds_min: Minimum step size.

        Returns:
            numpy.ndarray: Array of points along the path.
        """
        steps = clib.geom2d_path_get_steps(self, ds_min)
        return clib.geom2d_path_get_points_at_steps(self, steps)

    def get_points_at_steps(self, steps):
        """Get points along the path at specified steps.

        Parameters:
            steps: Array of step positions along the path.

        Returns:
            numpy.ndarray: Array of points at the specified steps.
        """
        return clib.geom2d_path_get_points_at_steps(self, steps)

    def get_steps(self, ds_min):
        """Get steps along the path with a minimum step size.

        Parameters:
            ds_min: Minimum step size.

        Returns:
            numpy.ndarray: Array of step positions along the path.
        """
        return clib.geom2d_path_get_steps(self, ds_min)

    @property
    def length(self):
        """Get the total length of the path."""
        return clib.geom2d_path_get_length(self)

    def plot(self, ax=None, ds_min=0.1, **kwargs):
        """Plot the path using Matplotlib.

        Parameters:
            ax: Matplotlib Axes
                If None, uses the current axes.
            ds_min: float | None
                Minimum step size. If None, uses adaptive point distribution,
                instead of uniform step size.
            **kwargs: dict
                Additional keyword arguments passed to the plot function.
        """

        if ax is None:
            ax = plt.gca()
            ax.set_aspect("equal")

        if ds_min is not None:
            steps = clib.geom2d_path_get_steps(self, ds_min)
            points = clib.geom2d_path_get_points_at_steps(self, steps)
        else:
            points = clib.geom2d_path_get_points(self)

        ax.plot(points["x"], points["y"], **kwargs)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        return self


clib._struct_dtypes["G2DPath"] = Path2D
