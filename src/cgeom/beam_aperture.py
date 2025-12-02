from .cwrapper import clib



class Aperture:

    def __init__(self, line):
        self.line = line

    def get_2dpoints_from_limit_element(self,name):
        """
        Get 2d points from limit element

        Parameters
        ----------
        name : str
            Name of the limit element
        Returns: G2DBeamApertureData
            2D beam aperture data
        """


    def get_twiss_data_at_element(self, name):
        """
        Get twiss data at element 
        """