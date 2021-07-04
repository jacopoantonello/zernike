#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This file is part of zernike.
#
# zernike is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# zernike is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with zernike.  If not, see <http://www.gnu.org/licenses/>.

import argparse
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from zernike import RZern
matplotlib.rcParams['mathtext.fontset'] = 'cm'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Plot the pyramid of Zernike polynomials',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--max-radial-order',
                        default=6,
                        type=int,
                        metavar='NMAX')
    args = parser.parse_args()

    plt.close('all')

    fs = 7
    inset_width = 0.9
    inset_height = 0.9

    nradial = args.max_radial_order
    cart = RZern(nradial)
    L, K = 100, 100
    ddx = np.linspace(-1.0, 1.0, K)
    ddy = np.linspace(-1.0, 1.0, L)
    xv, yv = np.meshgrid(ddx, ddy)
    cart.make_cart_grid(xv, yv)
    c = np.zeros(cart.nk)

    fig, ax = plt.subplots()
    ax.set(xlim=(-nradial - 1, nradial + 1), ylim=(nradial + 0.5, 0))
    ax.set_xlabel('$m$')
    ax.set_ylabel('$n$')

    # Styling: no frame, no ticks, tex rendered tick labels.
    ax.set_frame_on(False)
    ax.tick_params(left=False, bottom=False)
    ax.xaxis.set_major_formatter('${x:n}$')
    ax.yaxis.set_major_formatter('${x:n}$')

    for k in range(1, cart.nk + 1):
        # Plot the phase for the `k`-th Zernike polynomial and add a label
        # next to it.
        n, m = cart.noll2nm(k)
        left = m - inset_width / 2
        bottom = n - inset_height / 2
        # Create inset in data coordinates using ax.transData as transform
        axins = inset_axes(ax, width="100%", height="100%",
                           bbox_to_anchor=(left, bottom, inset_width,
                                           inset_height),
                           bbox_transform=ax.transData, borderpad=0)

        c *= 0.0
        c[k - 1] = 1.0
        Phi = cart.eval_grid(c, matrix=True)
        axins.imshow(Phi, origin='lower', extent=[0, 1, 0, 1])
        axins.axis('off')

        axins.text(0.9, 0.9, '$\\#' + str(k) + '$', fontsize=fs)

    plt.show()
