import numpy as np

import xtrack as xt
import matplotlib.pyplot as plt

from cgeom import clib
from cgeom.beam_aperture import Aperture

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
aperture_model = Aperture(b1_with_aper)

# Compute envelopes
nemitt_x = 3e-6
nemitt_y = 3e-6
sigmas = 11
tw = b1_with_aper.twiss4d()
one_sigma_envelope_x = np.sqrt(tw.betx * nemitt_x / tw.gamma0)
one_sigma_envelope_y = np.sqrt(tw.bety * nemitt_y / tw.gamma0)

# Obtain extents
extents = np.array(list(aperture_model.extents))
extents_on_axes = np.array(list(aperture_model.extents_on_axes))

# Plot
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8, 6))

ax1.plot(aperture_model.s_positions, extents[:, 2], label=r'$x_{min}$')
ax1.plot(aperture_model.s_positions, extents[:, 0], label=r'$x_{max}$')
ax1.plot(aperture_model.s_positions, extents_on_axes[:, 2], linestyle='--', label=r'$x_{0,min}$')
ax1.plot(aperture_model.s_positions, extents_on_axes[:, 0], linestyle='--', label=r'$x_{0,max}$')
ax1.fill_between(tw.s, -one_sigma_envelope_x + tw.x, one_sigma_envelope_x + tw.x, alpha=0.8, label=r'$1\sigma$ envelope')
ax1.fill_between(tw.s, -sigmas * one_sigma_envelope_x + tw.x, sigmas * one_sigma_envelope_x + tw.x, alpha=0.2, label=rf'${sigmas}\sigma$ envelope')
ax1.set_ylabel("x [m]")
ax1.legend()

ax2.plot(aperture_model.s_positions, extents[:, 1], label=r'$y_{min}$')
ax2.plot(aperture_model.s_positions, extents[:, 3], label=r'$y_{max}$')
ax2.plot(aperture_model.s_positions, extents_on_axes[:, 1], linestyle='--', label=r'$y_{0,min}$')
ax2.plot(aperture_model.s_positions, extents_on_axes[:, 3], linestyle='--', label=r'$y_{0,max}$')
ax2.fill_between(tw.s, -one_sigma_envelope_y + tw.y, one_sigma_envelope_y + tw.y, alpha=0.8, label=r'$1\sigma$ envelope')
ax2.fill_between(tw.s, -sigmas * one_sigma_envelope_y + tw.y, sigmas * one_sigma_envelope_y + tw.y, alpha=0.2, label=rf'${sigmas}\sigma$ envelope')
ax2.set_ylabel("y [m]")
ax2.set_xlabel("s [m]")
ax2.legend()

plt.setp(ax1.get_xticklabels(), visible=False)
plt.tight_layout()
plt.show()

# Show some interesting shapes
debug_elements = [
    'vaehn.4r1.b.b1_aper',
    'vctnc.4r1.b.b1_aper',
    'vctcw.6l2.b.b1_aper',
    'msib.c6l2.b1_mken_aper',
    'tdisa.a4l2.a.b1_aper',
    'btvss.6l2.b1_mken_aper',
]
for name in debug_elements:
    elem = b1_with_aper.element_dict[name]
    print(f'==> {name} is {elem}')
    path, cx, cy = Aperture.CONVERTER_FUNCTIONS[type(elem)](elem)
    print(f'    * len_points = {clib.geom2d_path_get_len_points(path)}')
    print(f'    * segments = {path.segments}')
    print()

    path.plot(label='adaptive point distribution', marker='.')
    path.plot(ds_min=None, label='uniform point distribution', linestyle='--', marker='x')
    plt.title(f'{name} at s = {b1_with_aper.get_table().rows[name].s[0]}')
    plt.legend()
    plt.show()
