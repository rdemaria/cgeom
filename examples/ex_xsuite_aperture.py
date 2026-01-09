import xtrack as xt
import cgeom

lhc=xt.Line.from_json('../test_data/line_and_particle/line_lhc_no_errors.json')

aperture=cgeom.Aperture.from_line_with_limit(lhc)
aperture.get_aperture_margin_mm(line="b1",element="vmdqb.a1r1.a.b1")# return a vector at each s position
aperture.get_aperture_sigma(line="b1",element="vmdqb.a1r1.a.b1")# return a vector at each s position
aperture.get_aperture_sigma_hv(line="b1",element="vmdqb.a1r1.a.b1")# return a vector at each s position


aperture=cgeom.Aperture.from_line_with_aperture(lhc)
aperture.get_aperture_sigma(line="b1",element="mbrc.4r1", resolution=0.1) # return a vector at each s position
aperture.get_aperture_sigma(line="b1",element="mb.a18r1.b1", resolution=0.1) # return a vector at each s position

# for the MBRC lofting is needed
# look into ex_interpolate.py