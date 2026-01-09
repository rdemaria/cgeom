import xobjects as xo


class Circle(xo.Struct):
    radius= xo.Float32


class Rectangle(xo.Struct):
    half_width= xo.Float32
    half_height= xo.Float32


class Ellipse(xo.Struct):
    half_major= xo.Float32
    half_minor= xo.Float32


class Octagon(xo.Struct):
    half_width= xo.Float32
    half_height= xo.Float32
    half_diagonal= xo.Float32


class Polygon(xo.Struct):
    vertices= xo.Float32[:, 2]

class SVGShape(xo.Struct):
    svg_data = xo.String()


class ProfileUnion(xo.UnionRef):
    _ref_types = [Circle, Rectangle, Ellipse, Octagon, Polygon, SVGShape]


class ProfileList(xo.Struct):
    profiles = ProfileUnion[:]


class Profile:
    def __init__(self, profiles):
        self.names = {}  # dict of profile names to indices
        self.profiles = profiles  # list of profile objects


class ProfilePosition(xo.Struct):
    profile_index = xo.Int64
    s_position = xo.Float32  # position along the type axiswhere this profile sits
    shift_x = xo.Float32  # horizontal shift of the profile center from the type axis
    shift_y = xo.Float32  # vertical shift of the profile center from the type axis
    rot_x = xo.Float32  # rotation of the profile around the type axis in radians
    rot_y = xo.Float32  # rotation of the profile around the vertical axis in radians
    rot_z = xo.Float32  # rotation of the profile around the horizontal axis in radians


class ApertureType(xo.Struct):
    curvature = (
        xo.Float32
    )  # curvature of the type axis assumed to be in the horizontal plane
    profiles = ProfilePosition[:]


class ApertureTypeList(xo.Struct):
    aperture_types = ApertureType[:]


class ApertureTypes:
    def __init__(self, aperture_types):
        self.names = {}  # dict of aperture type names to indices
        self.aperture_types = aperture_types  # list of pipes in the lab frame


class TypePosition(xo.Struct):
    type_index = xo.Int32
    ref_position = xo.String  # identify a point in survey
    idx_position = xo.Int32  # index of the point in the survey
    transformation = xo.Float32[
        4, 4
    ]  # 3D rigid transformation matrix from the ref point to the center of the aperture type


class TypePositionList(xo.Struct):
    type_positions = TypePosition[:]


class ApertureModel:
    def __init__(self, line_name, type_positions):
        self.names = {}  # dict of aperture model names to indices
        self.line_name = line_name
        self.type_positions = type_positions  # positioning of types in line frame


class Aperture:
    def __init__(self, env, profiles, aperture_types, aperture_model, cross_sections):
        self.env = env
        self.profiles = profiles  # list of profile objects
        self.aperture_types = aperture_types  # list of pipes in the lab frame
        self.aperture_model = aperture_model  # positioning of types in line frame
        self.cross_sections = cross_sections

