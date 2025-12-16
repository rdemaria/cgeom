import numpy as np

import xtrack as xt
import matplotlib.pyplot as plt

from cgeom import clib, Path2D
from cgeom.beam_aperture import Aperture, LimitTypes

url = 'https://raw.githubusercontent.com/xsuite/xcoll/refs/heads/main/examples/machines/lhc_run3_b1.json'
b1_with_aper = xt.load(url)


# Clean up the line by removing 10m apertures
def clean_up_line(line):
    def _delete_if(name):
        element = line.element_dict[name]
        if not isinstance(element, xt.LimitRectEllipse):
            return False

        min_qty = np.min([element.max_x, element.max_y, element.a, element.b])
        if min_qty >= 9:
            return True

        return False

    clean_names = [name for name in line.element_names if not _delete_if(name)]
    line.element_names = clean_names


clean_up_line(b1_with_aper)
aperture_model = Aperture(b1_with_aper, tol_x=0.001, tol_y=0.002, tol_r=0.001)
aperture_model.set_beam_data(
    emitt_x_norm=2.5e-6,
    emitt_y_norm=2.5e-6,
    delta_rms=0,
    tol_co=0,
    tol_disp=0,
    tol_disp_ref_dx=0,
    tol_disp_ref_beta=0,
    tol_energy=0,
    tol_betabeating=0,
    halo_x=5,
    halo_y=5,
    halo_r=5.8,
    halo_primary=0,
)

# Compute envelopes
nemitt_x = 2.5e-6
nemitt_y = 2.5e-6
sigmas = 11
tw = b1_with_aper.twiss4d()
one_sigma_envelope_x = np.sqrt(tw.betx * nemitt_x / tw.gamma0)
one_sigma_envelope_y = np.sqrt(tw.bety * nemitt_y / tw.gamma0)

table = b1_with_aper.get_table()
table['one_sigma_envelope_x'] = one_sigma_envelope_x
table['one_sigma_envelope_y'] = one_sigma_envelope_y
table['x'] = tw.x
table['y'] = tw.y

# name_pattern = 'mbrc.4r1.b1.*'
name_pattern = 'mqy.4r1.b1..*'
section_to_plot = table.rows[name_pattern]

ax = plt.gca()
ax.set_aspect('equal')

for name in section_to_plot.name:
    elem = b1_with_aper.element_dict[name]

    this_row = section_to_plot.rows[name]
    sig_x = this_row.one_sigma_envelope_x
    sig_y = this_row.one_sigma_envelope_y

    hr = aperture_model.beam_data.halo_r
    hx = aperture_model.beam_data.halo_x
    hy = aperture_model.beam_data.halo_y
    tmp = np.sqrt(2 * (hr - hx) * (hr - hy))
    sh = hr - hy + tmp
    sv = hr - hx + tmp
    sr = hx + hy - hr - tmp

    bh, bv, brx, bry = sh * sig_x, sh * sig_y, sr * sig_x, sr * sig_y

    # Plot the beam halo
    beam_re = Path2D.from_racetrack(
        halfhside=bh + brx,
        halfvside=bv + bry,
        rx=brx,
        ry=bry,
    )
    beam_points = beam_re.get_points(ds_min=None)
    beam_points['x'] += section_to_plot.rows[name].x
    beam_points['y'] += section_to_plot.rows[name].y
    ax.plot(beam_points['x'], beam_points['y'], color='g', label=fr'beam with halo = (x = {hx}$\sigma$, y = {hy}$\sigma$, r = {hr}$\sigma$))')

    if not isinstance(elem, LimitTypes):
        continue

    aperture = aperture_model.get_aperture(name)
    twiss = aperture_model.get_twiss(name)
    beam_data = aperture_model.beam_data

    # Plot the aperture
    points = aperture.points
    ax.plot(points["x"], points["y"], color='k', label='aperture')

    # Plot the beam envelope with tolerances
    envelope_points = clib.geom2d_get_beam_envelope(beam_data, twiss, aperture, len_points=30)
    envelope_points['x'] += section_to_plot.rows[name].x
    envelope_points['y'] += section_to_plot.rows[name].y
    ax.plot(envelope_points["x"], envelope_points["y"], color='r', label=f'envelope for aperture tol = (x = {aperture.tol_x} m, y = {aperture.tol_y} m, r = {aperture.tol_r} m)')

handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys())

plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.title(fr'Envelopes for {name_pattern}: $s \in [{section_to_plot.s[0]}, {section_to_plot.s[-1]}]$')
plt.show()
