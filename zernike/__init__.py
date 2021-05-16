#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Complex- and real-valued Zernike polynomials.

"""

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

from __future__ import division, print_function  # Python 2

try:
    from abc import ABC, abstractmethod
except ImportError:  # Probably Python 2
    from abc import ABCMeta, abstractmethod
    ABC = ABCMeta('ABC', (object, ), {'__slots__': ()})
from math import factorial

import h5py
import numpy as np
from numpy.linalg import lstsq, matrix_rank, norm

from zernike import version

__author__ = 'J Antonello'
__copyright__ = 'Copyright 2016-2021, J. Antonello'
__license__ = 'GPLv3+'
__email__ = 'jacopo@antonello.org'
__status__ = 'Production'
__version__ = version.__version__
__date__ = version.__date__
__commit__ = version.__commit__
__doc__ = """
Python code for Zernike polynomials.

date:    {}
version: {}
commit:  {}
""".format(
    __date__,
    __version__,
    __commit__,
)
__docformat__ = 'restructuredtext'

HDF5_options = {
    'chunks': True,
    'shuffle': True,
    'fletcher32': True,
    'compression': 'gzip',
    'compression_opts': 9
}


class Zern(ABC):
    """Shared code for `RZern` and `CZern`.

    This is an abstract class, use `RZern` and `CZern` instead. Only
    `NORM_NOLL` is implemented. The polynomials are ordered and normalised
    according to [N1976]_ and [M1994]_, see also Appendix A in [A2015]_.

    References
    ----------
    ..  [N1976] R. Noll, "Zernike polynomials and atmospheric turbulence,"
        J. Opt. Soc. Am.  66, 207-211 (1976). `doi
        <http://dx.doi.org/10.1364/JOSA.66.000207>`__.
    ..  [M1994] V. N. Mahajan, "Zernike circle polynomials and optical
        aberrations of systems with circular pupils," Appl. Opt. 33, 8121–8124
        (1994). `doi <http://dx.doi.org/10.1364/AO.33.008121>`__.
    ..  [A2015] Jacopo Antonello and Michel Verhaegen, "Modal-based phase
        retrieval for adaptive optics," J. Opt. Soc. Am. A 32, 1160-1170
        (2015). `url <http://dx.doi.org/10.1364/JOSAA.32.001160>`__.

    """

    # unimplemented, ENZ papers of Janssen, Braat, etc.
    NORM_ZERNIKE = 0
    # Noll/Mahajan's normalisation, unit variance over the unit disk
    NORM_NOLL = 1

    # FIXME
    def _print_rhotab(self):
        for i in range(self.nk):
            for j in range(self.n + 1):
                print('{:< 3.3f} '.format(self.rhotab[i, j]), end='')
            print('')

    def _print_nmtab(self):
        for i in range(self.nk):
            print('{:< d} {:< d}'.format(self.ntab[i], self.mtab[i]))

    # FIXME
    def _make_rhotab_row(self, c, n, m):
        # col major, row i, col j
        self.coefnorm[c] = self.ck(n, m)
        for s in range((n - m) // 2 + 1):
            self.rhotab[c, self.n - (n - 2 * s)] = (
                ((-1)**s) * factorial(n - s) /
                (factorial(s) * factorial((n + m) // 2 - s) *
                 factorial((n - m) // 2 - s)))

    def __init__(self, n, normalise=NORM_NOLL):
        """Initialise Zernike polynomials up to radial order `n`.

        This is an abstract class, use `RZern` and `CZern` instead. Only
        `NORM_NOLL` is implemented.

        """
        self.shape = None
        self.numpy_dtype = 'undefined'

        nk = (n + 1) * (n + 2) // 2
        self.n = n
        self.nk = nk
        self.normalise = normalise
        assert (self.normalise == self.NORM_NOLL)

        # coefficients of R_n^m(\rho), see [N1976]_
        rhotab = np.zeros((nk, n + 1), order='F')
        # coefficients of \int R_n^m(\rho) \rho \,d\rho
        rhoitab = np.zeros((nk, n + 3), order='F')
        ntab = np.zeros(nk, dtype=int)
        mtab = np.zeros(nk, dtype=int)
        coefnorm = np.zeros(nk)

        self.rhotab = rhotab
        self.rhoitab = rhoitab
        self.ntab = ntab
        self.mtab = mtab
        self.coefnorm = coefnorm

        self.rhotab[0, n] = 1.0
        self.coefnorm[0] = 1.0

        for ni in range(1, n + 1):
            for mi in range(-ni, ni + 1, 2):
                k = self.nm2noll(ni, mi) - 1
                self._make_rhotab_row(k, ni, abs(mi))
                ntab[k], mtab[k] = ni, mi

        # make rhoitab
        for ci in range(nk):
            for ni in range(n + 1):
                self.rhoitab[ci, ni] = self.rhotab[ci, ni] / (n + 2 - ni)

    @abstractmethod
    def ck(self, n, m):
        r"""Normalisation coefficient for the `k`-th Zernike polynomial.

        For real-valued Zernike polynomials,

        .. math::

            c_n^m =
            \begin{cases}
                \sqrt{n + 1} & m = 0\\
                \sqrt{2(n + 1)} & m \neq 0
            \end{cases}

        For complex-valued Zernike polynomials, :math:`c_n^m = \sqrt{n + 1}`.

        """
        pass

    def Rnm(self, k, rho):
        r"""Compute the `k`-th radial polynomial :math:`R_n^m(\rho)`.

        The radial polynomial is defined in Eq. (2) of [N1976]_ and Eq. (2) of
        [M1994]_.

        References
        ----------
        ..  [N1976] R. Noll, "Zernike polynomials and atmospheric turbulence,"
            J. Opt. Soc. Am.  66, 207-211 (1976). `doi
            <http://dx.doi.org/10.1364/JOSA.66.000207>`__.
        ..  [M1994] V. N. Mahajan, "Zernike circle polynomials and optical
            aberrations of systems with circular pupils," Appl. Opt.  33,
            8121–8124 (1994). `doi
            <http://dx.doi.org/10.1364/AO.33.008121>`__.

        """
        return np.polyval(self.rhotab[k, :], rho)

    def I_Rnmrho(self, k, rho):
        r"""Compute :math:`\int R_n^m(\rho)\rho`."""
        return sum([(rho**(self.n + 2 - i)) * self.rhoitab[k, i]
                    for i in range(self.n + 3)])

    @abstractmethod
    def angular(self, k, theta):
        r"""Compute the angular function for the `k`-th Zernike polynomial.

        For real-valued polynomials, the angular function is

        .. math::

            \Theta_n^m(\theta) =
            \begin{cases}
                \cos(m\theta) & m \ge 0\\
                -\sin(m\theta) & m < 0
            \end{cases}

        For complex-valued Zernike polynomials, :math:`\Theta_n^m(\theta) =
        \exp(i m\theta)`.

        References
        ----------
        ..  [A2015] Jacopo Antonello and Michel Verhaegen, "Modal-based phase
            retrieval for adaptive optics," J. Opt. Soc. Am. A 32, 1160-1170
            (2015). `url <http://dx.doi.org/10.1364/JOSAA.32.001160>`__.

        """
        pass

    def radial(self, k, rho):
        r"""Compute the radial function for the `k`-th Zernike polynomial.

        The radial function is :math:`c_n^m R_n^{|m|}(\rho)`, see Appendix A in
        [A2015]_.

        References
        ----------
        ..  [A2015] Jacopo Antonello and Michel Verhaegen, "Modal-based phase
            retrieval for adaptive optics," J. Opt. Soc. Am. A 32, 1160-1170
            (2015). `url <http://dx.doi.org/10.1364/JOSAA.32.001160>`__.

        """
        return self.coefnorm[k] * self.Rnm(k, rho)

    def Zk(self, k, rho, theta):
        r"""Compute the `k`-th Zernike polynomial.

        For real-valued Zernike polynomials, :math:`\mathcal{Z}_n^m(\rho,
        \theta) = c_n^m R_n^{|m|}(\rho) \Theta_n^m(\theta)`. For complex-valued
        Zernike polynomials, :math:`\mathcal{N}_n^m(\rho, \theta) = c_n^m
        R_n^{|m|}(\rho)\exp(i m\theta)`. See Appendix A in [A2015]_.

        References
        ----------
        ..  [A2015] Jacopo Antonello and Michel Verhaegen, "Modal-based phase
            retrieval for adaptive optics," J. Opt. Soc. Am. A 32, 1160-1170
            (2015). `url <http://dx.doi.org/10.1364/JOSAA.32.001160>`__.

        """
        return self.radial(k, rho) * self.angular(k, theta)

    def eval_a(self, a, rho, theta):
        r"""Compute the sum of `self.nk` Zernike polynomials at a point.

        This evaluates :math:`\sum_k a_k \mathcal{Z}_k(\rho, \theta)` or
        :math:`\sum_k a_k \mathcal{N}_k(\rho, \theta)` at `rho` and `theta`,
        where :math:`a_k` are the elements of `a`.

        """
        return sum([a[j] * self.Zk(j, rho, theta) for j in range(self.nk)])

    @staticmethod
    def nm2noll(n, m):
        """Convert indices `(n, m)` to the Noll's index `k`.

        Note that Noll's index `k` starts from one and Python indexing is
        zero-based.

        """
        k = n * (n + 1) // 2 + abs(m)
        if (m <= 0 and n % 4 in (0, 1)) or (m >= 0 and n % 4 in (2, 3)):
            k += 1
        return k

    def noll2nm(self, k):
        """Convert Noll's index `k` to the indices `(n, m)`.

        Note that Noll's index `k` starts from one and Python indexing is
        zero-based.

        """
        n = self.ntab[k - 1]
        m = self.mtab[k - 1]
        return n, m

    def vect(self, Phi):
        r"Reshape `Phi` into a vector"
        return Phi.ravel(order='F')

    def matrix(self, Phi):
        r"Reshape `Phi` into a matrix"
        if self.shape is None:
            raise ValueError('Use make_cart_grid() to define the shape first')
        elif self.shape[0] * self.shape[1] != Phi.size:
            raise ValueError('Phi.shape should be {}'.format(self.shape))
        return Phi.reshape(self.shape, order='F')

    def make_cart_grid(self, xx, yy, unit_circle=True):
        r"""Make a cartesian grid to evaluate the Zernike polynomials.

        Parameters
        ----------
        - `xx`: `numpy` array generated with `numpy.meshgrid()`.
        - `yy`: `numpy` array generated with `numpy.meshgrid()`.
        - `unit_circle`: set `np.nan` for points where :math:`\rho > 1`.

        Examples
        --------

        .. code:: python

            import numpy as np
            from zernike import RZern

            cart = RZern(6)
            dd = np.linspace(-1.0, 1.0, 200)
            xv, yv = np.meshgrid(dd, dd)
            cart.make_cart_grid(xv, yv)

        Notes
        -----
        `ZZ` is stored with `order='F'`.

        """
        self.ZZ = np.zeros((xx.size, self.nk),
                           order='F',
                           dtype=self.numpy_dtype)
        self.shape = xx.shape
        rho = np.sqrt(np.square(xx) + np.square(yy))
        theta = np.arctan2(yy, xx)
        for k in range(self.nk):
            prod = self.radial(k, rho) * self.angular(k, theta)
            if unit_circle:
                prod[rho > 1.0] = np.nan
            self.ZZ[:, k] = self.vect(prod)

    def eval_grid(self, a, matrix=False):
        """Evaluate the Zernike polynomials using the coefficients in `a`.

        Parameters
        ----------
        - `a`: Zernike coefficients.

        Returns
        -------
        -   `Phi`: a `numpy` vector in column major order. Use the `vect()` or
            `matrix()` methods to flatten or unflatten the matrix.
        -   `matrix`: return a matrix instead of a vector

        Examples
        --------

        .. code:: python

            import numpy as np
            import matplotlib.pyplot as plt
            from zernike import RZern

            cart = RZern(6)
            L, K = 200, 250
            ddx = np.linspace(-1.0, 1.0, K)
            ddy = np.linspace(-1.0, 1.0, L)
            xv, yv = np.meshgrid(ddx, ddy)
            cart.make_cart_grid(xv, yv)

            c = np.zeros(cart.nk)
            plt.figure(1)
            for i in range(1, 10):
                plt.subplot(3, 3, i)
                c *= 0.0
                c[i] = 1.0
                Phi = cart.eval_grid(c, matrix=True)
                plt.imshow(Phi, origin='lower', extent=(-1, 1, -1, 1))
                plt.axis('off')

            plt.show()

        """
        if a.size != self.nk:
            raise ValueError('a.size = {} but self.nk = {}'.format(
                a.size, self.nk))
        Phi = np.dot(self.ZZ, a)
        if matrix:
            return self.matrix(Phi)
        else:
            return Phi

    def fit_cart_grid(self, Phi, rcond=None):
        """Fit a cartesian grid using least-squares.

        Parameters
        ----------
        - `Phi`: cartesian grid, e.g., generated with make_cart_grid().
        - `rcond`: rcond supplied to `lstsq`

        Returns
        -------
        -   `a`, `numpy` vector of Zernike coefficients
        -   `res`, see `lstsq`
        -   `rnk`, see `lstsq`
        -   `sv`, see `lstsq`

        Examples
        --------

        .. code:: python

            import numpy as np
            import matplotlib.pyplot as plt
            from zernike import RZern

            cart = RZern(6)
            L, K = 200, 250
            ddx = np.linspace(-1.0, 1.0, K)
            ddy = np.linspace(-1.0, 1.0, L)
            xv, yv = np.meshgrid(ddx, ddy)
            cart.make_cart_grid(xv, yv)

            c0 = np.random.normal(size=cart.nk)
            Phi = cart.eval_grid(c0, matrix=True)
            c1 = cart.fit_cart_grid(Phi)[0]
            plt.figure(1)
            plt.subplot(1, 2, 1)
            plt.imshow(Phi, origin='lower', extent=(-1, 1, -1, 1))
            plt.axis('off')
            plt.subplot(1, 2, 2)
            plt.plot(range(1, cart.nk + 1), c0, marker='.')
            plt.plot(range(1, cart.nk + 1), c1, marker='.')

            plt.show()

        """
        vPhi = self.vect(Phi)
        zfm = np.logical_and(np.isfinite(self.ZZ[:, 0]), np.isfinite(vPhi))
        zfA = self.ZZ[zfm, :]
        Phi1 = vPhi[zfm]

        a, res, rnk, sv = lstsq(np.dot(zfA.T, zfA),
                                np.dot(zfA.T, Phi1),
                                rcond=rcond)

        return a, res, rnk, sv

    def make_pol_grid(self, rho_j, theta_i):
        r"""Make a polar grid to evaluate the Zernike polynomials.

        Parameters
        ----------
        - `rho_j`: `numpy` vector
        - `theta_i`: `numpy` vector

        Examples
        --------

        .. code:: python

            import numpy as np
            from zernike import RZern

            pol = RZern(6)
            pol.make_pol_grid(np.linspace(0.0, 1.0), np.linspace(0.0, 2*np.pi))

        """
        L, K = theta_i.size, rho_j.size
        self.ZZ = np.zeros((K * L, self.nk), order='F', dtype=self.numpy_dtype)
        self.shape = (L, K)
        for k in range(self.nk):
            rad = self.radial(k, rho_j).reshape((1, K), order='F')
            ang = self.angular(k, theta_i).reshape((L, 1), order='F')
            prod = rad * ang
            self.ZZ[:, k] = prod.ravel(order='F')

    def save(self,
             filename,
             prepend=None,
             params=HDF5_options,
             libver='latest'):
        """Save object into an HDF5 file."""
        f = h5py.File(filename, 'w', libver=libver)
        self.save_h5py(f, prepend=prepend, params=params)
        f.close()

    def save_h5py(self, f, prepend=None, params=HDF5_options):
        """Dump object contents into an opened HDF5 file object."""
        prefix = self.__class__.__name__ + '/'

        if prepend is not None:
            prefix = prepend + prefix

        params['data'] = self.coefnorm
        f.create_dataset(prefix + 'coefnorm', **params)

        params['data'] = self.ntab
        f.create_dataset(prefix + 'ntab', **params)

        params['data'] = self.mtab
        f.create_dataset(prefix + 'mtab', **params)

        f.create_dataset(prefix + 'n', data=np.array([self.n], dtype=int))

        f.create_dataset(prefix + 'nk', data=np.array([self.nk], dtype=int))

        f.create_dataset(prefix + 'normalise',
                         data=np.array([self.normalise], dtype=int))

        params['data'] = self.rhoitab
        f.create_dataset(prefix + 'rhoitab', **params)

        params['data'] = self.rhotab
        f.create_dataset(prefix + 'rhotab', **params)

        f.create_dataset(prefix + 'numpy_dtype',
                         data=np.array(self.numpy_dtype.encode('utf-8'),
                                       dtype=h5py.string_dtype(
                                           'utf-8', len(self.numpy_dtype))))

        try:
            params['data'] = self.ZZ
            f.create_dataset(prefix + 'ZZ', **params)
            f.create_dataset(prefix + 'shape', data=self.shape)
        except AttributeError:
            pass

    def make_rotation(self, alpha):
        r"""Make an orthogonal matrix to rotate the pupil.

        Parameters
        ----------
        - `alpha`: `float` rotation angle in degrees

        """
        alpha *= np.pi / 180
        nml = list(zip(self.ntab.tolist(), self.mtab.tolist()))
        R = np.zeros((self.nk, self.nk))
        for i, nm in enumerate(nml):
            n, m = nm[0], nm[1]
            if m == 0:
                R[i, i] = 1.0
            elif m > 0:
                R[i, i] = np.cos(m * alpha)
                R[i, nml.index((n, -m))] = np.sin(m * alpha)
            else:
                R[i, nml.index((n, -m))] = -np.sin(abs(m) * alpha)
                R[i, i] = np.cos(abs(m) * alpha)

        # checks
        assert (matrix_rank(R) == R.shape[0])
        assert (norm((np.dot(R, R.T) - np.eye(self.nk)).ravel()) < 1e-11)
        assert (norm((np.dot(R.T, R) - np.eye(self.nk)).ravel()) < 1e-11)

        return R

    def make_yflip(self):
        r"Make an orthogonal matrix to flip the pupil along y."
        nml = list(zip(self.ntab.tolist(), self.mtab.tolist()))
        R = np.zeros((self.nk, self.nk))
        for i, nm in enumerate(nml):
            m = nm[1]
            if m < 0:
                R[i, i] = -1.0
            else:
                R[i, i] = 1.0

        # checks
        assert (matrix_rank(R) == R.shape[0])
        assert (norm((np.dot(R, R.T) - np.eye(self.nk)).ravel()) < 1e-11)
        assert (norm((np.dot(R.T, R) - np.eye(self.nk)).ravel()) < 1e-11)

        return R

    def make_xflip(self):
        r"Make an orthogonal matrix to flip the pupil along x."
        nml = list(zip(self.ntab.tolist(), self.mtab.tolist()))
        R = np.zeros((self.nk, self.nk))
        for i, nm in enumerate(nml):
            m = nm[1]
            if abs(m) % 2 == 0 and m < 0:
                R[i, i] = -1.0
            elif abs(m) % 2 == 1 and m > 0:
                R[i, i] = -1.0
            else:
                R[i, i] = 1.0
        # checks
        assert (matrix_rank(R) == R.shape[0])
        assert (norm((np.dot(R, R.T) - np.eye(self.nk)).ravel()) < 1e-11)
        assert (norm((np.dot(R.T, R) - np.eye(self.nk)).ravel()) < 1e-11)

        return R

    def make_permutation(self, noll):
        P = np.eye(self.nk)
        return P[:, np.asarray(noll) - 1]

    @classmethod
    def load(cls, filename, prepend=None):
        """Load object from an HDF5 file."""
        f = h5py.File(filename, 'r')
        z = cls.load_h5py(f, prepend=prepend)
        f.close()

        return z

    @classmethod
    def load_h5py(cls, f, prepend=None):
        """Load object contents from an opened HDF5 file object."""
        z = cls(1)

        prefix = cls.__name__ + '/'

        if prepend is not None:
            prefix = prepend + prefix

        z.coefnorm = f[prefix + 'coefnorm'][()]
        z.ntab = f[prefix + 'ntab'][()]
        z.mtab = f[prefix + 'mtab'][()]
        z.n = int(f[prefix + 'n'][0])
        z.nk = int(f[prefix + 'nk'][0])
        z.normalise = int(f[prefix + 'normalise'][0])
        z.rhoitab = f[prefix + 'rhoitab'][()]
        z.rhotab = f[prefix + 'rhotab'][()]
        z.numpy_dtype = f[prefix + 'numpy_dtype'][()]
        if isinstance(z.numpy_dtype, bytes):
            z.numpy_dtype = z.numpy_dtype.decode('utf-8')
        try:
            z.ZZ = f[prefix + 'ZZ'][()]
            z.shape = f[prefix + 'shape'][()]
        except KeyError:
            pass

        return z


class CZern(Zern):
    r"""Complex-valued Zernike polynomials.

    .. math::

        \mathcal{N}_k(\rho, \theta) = \mathcal{N}_n^m(\rho, \theta)

        \mathcal{N}_n^m(\rho, \theta) = c_n^m R_n^{|m|}(\rho)\exp(i m\theta)

        c_n^m = \sqrt{n + 1}

        \int |mathcal{N}_k(rho, theta)|^2 \rho d\rho d\theta = \pi, \;
        \text{for} \; k > 1

    See Eq. (A5) in [A2015]_.

    References
    ----------
    ..  [N1976] R. Noll, "Zernike polynomials and atmospheric turbulence,"
        J. Opt. Soc. Am.  66, 207-211 (1976). `doi
        <http://dx.doi.org/10.1364/JOSA.66.000207>`__.
    ..  [M1994] V. N. Mahajan, "Zernike circle polynomials and optical
        aberrations of systems with circular pupils," Appl. Opt. 33, 8121–8124
        (1994). `doi <http://dx.doi.org/10.1364/AO.33.008121>`__.
    ..  [A2015] Jacopo Antonello and Michel Verhaegen, "Modal-based phase
        retrieval for adaptive optics," J. Opt. Soc. Am. A 32, 1160-1170
        (2015). `url <http://dx.doi.org/10.1364/JOSAA.32.001160>`__.

    """
    def __init__(self, n, normalise=Zern.NORM_NOLL):
        super(CZern, self).__init__(n, normalise)
        self.numpy_dtype = 'complex'

    def ck(self, n, m):
        return np.sqrt(n + 1.0)

    def angular(self, j, theta):
        m = self.mtab[j]
        return np.exp(1j * m * theta)


class RZern(Zern):
    r"""Real-valued Zernike polynomials.

    .. math::

        \mathcal{Z}_k(\rho, \theta) = \mathcal{Z}_n^m(\rho, \theta)

        \mathcal{Z}_n^m(\rho, \theta) = c_n^m R_n^{|m|}(\rho)
        \Theta_n^m(\theta)

        c_n^m =
        \begin{cases}
            \sqrt{n + 1} & m = 0\\
            \sqrt{2(n + 1)} & m \neq 0
        \end{cases}

        \Theta_n^m(\theta) =
        \begin{cases}
            \cos(m\theta) & m \ge 0\\
            -\sin(m\theta) & m < 0
        \end{cases}

        \int |\mathcal{Z}_k(rho, theta)|^2 \rho d\rho d\theta = \pi, \;
        \text{for} \; k > 1

    See Eq. (A1) in [A2015]_.

    References
    ----------
    ..  [N1976] R. Noll, "Zernike polynomials and atmospheric turbulence,"
        J. Opt. Soc. Am.  66, 207-211 (1976). `doi
        <http://dx.doi.org/10.1364/JOSA.66.000207>`__.
    ..  [M1994] V. N. Mahajan, "Zernike circle polynomials and optical
        aberrations of systems with circular pupils," Appl. Opt. 33, 8121–8124
        (1994). `doi <http://dx.doi.org/10.1364/AO.33.008121>`__.
    ..  [A2015] Jacopo Antonello and Michel Verhaegen, "Modal-based phase
        retrieval for adaptive optics," J. Opt. Soc. Am. A 32, 1160-1170
        (2015). `url <http://dx.doi.org/10.1364/JOSAA.32.001160>`__.

    """
    def __init__(self, n, normalise=Zern.NORM_NOLL):
        super(RZern, self).__init__(n, normalise)
        self.numpy_dtype = 'float'

    def ck(self, n, m):
        if self.normalise == self.NORM_NOLL:
            if m == 0:
                return np.sqrt(n + 1.0)
            else:
                return np.sqrt(2.0 * (n + 1.0))
        else:
            return 1.0

    def angular(self, j, theta):
        m = self.mtab[j]
        if m >= 0:
            return np.cos(m * theta)
        else:
            return np.sin(-m * theta)


class FitZern:
    r"""Compute the approximate inner products using Riemann sums in polar
    coordinates.

    The grid in :math:`\mathbb{R}^2` is defined as follows:

    .. math::

        \{\theta_i = 2 \pi i/L, \; i = 0, \ldots, L - 1\}

        \{\rho_j = \cos((K - j - 1/2)\pi/(2K)), \; j = 0, \ldots, K - 1\}.

    See Appendix B in [A2015]_.

    References
    ----------
    ..  [A2015] Jacopo Antonello and Michel Verhaegen, "Modal-based phase
        retrieval for adaptive optics," J. Opt. Soc. Am. A 32, 1160-1170
        (2015). `url <http://dx.doi.org/10.1364/JOSAA.32.001160>`__.

    """
    def __init__(self, z, L, K):
        r"""Initialise a grid with `L` points for theta and `K` for rho.

        Parameters
        ----------
        - `z`: `RZern` or `CZern` object
        - `L`: number of points for sampling :math:`\theta`
        - `K`: number of points for sampling :math:`\rho`

        """
        self.z = z
        self.L = L
        self.K = K

        theta_i = np.array([2 * np.pi * i / L for i in range(L)])
        theta_a = np.array([np.pi * (2 * i - 1) / L for i in range(L)])
        theta_b = np.array([np.pi * (2 * i + 1) / L for i in range(L)])
        self.theta_i = theta_i
        self.theta_a = theta_a
        self.theta_b = theta_b

        rho_j = np.array(
            [np.cos(np.pi * (K - (j + 0.5)) / (2 * K)) for j in range(K)])
        rho_a = np.array(
            [np.cos(np.pi * (K - (j + 0.0)) / (2 * K)) for j in range(K)])
        rho_b = np.array(
            [np.cos(np.pi * (K - (j + 1.0)) / (2 * K)) for j in range(K)])
        self.rho_j = rho_j
        self.rho_a = rho_a
        self.rho_b = rho_b

        # L x (m -1) = L x (n - 1), col major, excluding m = 0
        I_cosm = np.array([
            (1 / m) * (np.sin(m * theta_b[thi]) - np.sin(m * theta_a[thi]))
            for m in range(1, z.n + 1) for thi in range(L)
        ])
        I_sinm = np.array([
            (1 / m) * (-np.cos(m * theta_b[thi]) + np.cos(m * theta_a[thi]))
            for m in range(1, z.n + 1) for thi in range(L)
        ])
        self.I_cosm = I_cosm
        self.I_sinm = I_sinm

        # K x nk, repeated entries for m < 0
        I_Rnmrho = np.array([
            z.I_Rnmrho(j, rho_b[rhi]) - z.I_Rnmrho(j, rho_a[rhi])
            for j in range(z.nk) for rhi in range(K)
        ])
        self.I_Rnmrho = I_Rnmrho

        A = np.zeros((K * L, z.nk), order='F', dtype=z.numpy_dtype)
        for k in range(z.nk):
            m = z.mtab[k]
            offK = self.K * k
            if m == 0:
                av = np.zeros((K * L, ), order='F', dtype=z.numpy_dtype)
                for j in range(K):
                    tmp1 = L * j
                    tmp2 = offK + j
                    for i in range(L):
                        av[tmp1 + i] = (2.0 / L) * I_Rnmrho[tmp2]
            else:
                av = np.zeros((K * L, ), order='F', dtype=z.numpy_dtype)
                if type(z) is RZern:
                    if m > 0:
                        I_cs = I_cosm
                        offL = L * (m - 1)
                    else:
                        I_cs = I_sinm
                        offL = L * (-m - 1)
                    for j in range(K):
                        tmp1 = L * j
                        tmp2 = offK + j
                        for i in range(L):
                            av[tmp1 + i] = ((1.0 / np.pi) * I_Rnmrho[tmp2] *
                                            I_cs[offL + i])
                elif type(z) is CZern:
                    offL = 0
                    sgn = 0.0
                    if m > 0:
                        offL = L * (m - 1)
                        sgn = -1.0
                    else:
                        offL = L * (-m - 1)
                        sgn = 1.0
                    for j in range(K):
                        tmp1 = L * j
                        tmp2 = offK + j
                        for i in range(L):
                            av[tmp1 + i] = ((1.0 / np.pi) * I_Rnmrho[tmp2] *
                                            (I_cosm[offL + i] +
                                             1j * sgn * I_sinm[offL + i]))
            A[:, k] = self.z.coefnorm[k] * av
        self.A = A.T

    def fit_ak(self, k, Phi):
        r"""Compute the inner product for the `k`-th coefficient.

        See Eq. (B2) and Eq. (B4) in Appendix B of [A2015]_.

        Parameters
        ----------
        -   `k`: `k`-th Zernike coefficient.
        -   `Phi`: `L` x `K` matrix stored in column major order into
            a list object.

        Returns
        ----------
        - `ak`: `k`-th Zernike coefficient

        References
        ----------
        ..  [A2015] Jacopo Antonello and Michel Verhaegen, "Modal-based phase
            retrieval for adaptive optics," J. Opt. Soc. Am. A 32, 1160-1170
            (2015). `url <http://dx.doi.org/10.1364/JOSAA.32.001160>`__.

        """
        m = self.z.mtab[k]
        offK = self.K * k
        offL = 0
        if type(self.z) is RZern:
            ak = tmp1 = tmp2 = 0.0
        else:
            ak = 0.0 + 1j * 0.0
            tmp1 = tmp2 = 0.0
        I_cs = None
        if m == 0:
            for j in range(self.K):
                tmp1 = self.L * j
                tmp2 = offK + j
                for i in range(self.L):
                    ak += Phi[tmp1 + i] * self.I_Rnmrho[tmp2]
            ak *= (2.0 / self.L)
        else:
            if type(self.z) is RZern:
                if m > 0:
                    I_cs = self.I_cosm
                    offL = self.L * (m - 1)
                else:
                    I_cs = self.I_sinm
                    offL = self.L * (-m - 1)
                for j in range(self.K):
                    tmp1 = self.L * j
                    tmp2 = offK + j
                    # FIXME maybe implement kron(B', A)*vec(X)
                    for i in range(self.L):
                        ak += ((1.0 / np.pi) * Phi[tmp1 + i] *
                               self.I_Rnmrho[tmp2] * I_cs[offL + i])
            elif type(self.z) is CZern:
                offL = 0
                sgn = 0.0
                if m > 0:
                    offL = self.L * (m - 1)
                    sgn = -1.0
                else:
                    offL = self.L * (-m - 1)
                    sgn = 1.0
                for j in range(self.K):
                    tmp1 = self.L * j
                    tmp2 = offK + j
                    for i in range(self.L):
                        ak += ((1.0 / np.pi) * Phi[tmp1 + i] *
                               self.I_Rnmrho[tmp2] *
                               (self.I_cosm[offL + i] +
                                1j * sgn * self.I_sinm[offL + i]))
            else:
                raise NotImplementedError()
        return self.z.coefnorm[k] * ak

    def _fit_slow(self, Phi):
        r"""Compute all the inner products. `Phi` is a phase grid in polar
        coordinates.

        Parameters
        ----------
        -   `Phi`: `L` x `K` matrix stored in column major order into
            a list object.

        Returns
        ----------
        `a`: list of Zernike coefficients

        """
        assert (len(Phi) == self.K * self.L)

        a = list()
        for k in range(self.z.nk):
            a.append(self.fit_ak(k, Phi))
        return a

    def fit(self, Phi):
        """Compute all the inner products. `Phi` is a phase grid in polar
        coordinates.

        Parameters
        ----------
        -   `Phi`: `L` x `K` matrix stored in column major order into
            a `numpy` array.

        Returns
        ----------
        `a`: `numpy` vector of Zernike coefficients

        Examples
        --------

        .. code:: python

            import numpy as np
            from zernike import RZern, FitZern

            pol = RZern(6)
            L, K = 200, 250
            ip = FitZern(pol, L, K)

            pol.make_pol_grid(ip.rho_j, ip.theta_i)
            c_true = np.random.normal(size=pol.nk)
            Phi = pol.eval_grid(c_true)
            c_hat = ip.fit(Phi)
            R = np.zeros((pol.nk, 3))
            R[:, 0] = c_true
            R[:, 1] = c_hat
            R[:, 2] = np.abs(c_true - c_hat)
            print(R)
            np.linalg.norm(c_true - c_hat)/np.linalg.norm(c_true)

        """
        assert (Phi.size == self.K * self.L and Phi.size == Phi.shape[0])
        return np.dot(self.A, Phi)

    def save(self,
             filename,
             prepend=None,
             params=HDF5_options,
             libver='latest'):
        """Save object into an HDF5 file."""
        f = h5py.File(filename, 'w', libver=libver)
        self.save_h5py(f, prepend=prepend, params=params)
        f.close()

    def save_h5py(self, f, prepend=None, params=HDF5_options):
        """Dump object contents into an opened HDF5 file object."""
        self.z.save_h5py(f, prepend=prepend)

        prefix = self.__class__.__name__ + '/'

        if prepend is not None:
            prefix = prepend + prefix

        params['data'] = self.A
        f.create_dataset(prefix + 'A', **params)

        params['data'] = self.I_cosm
        f.create_dataset(prefix + 'I_cosm', **params)

        params['data'] = self.I_sinm
        f.create_dataset(prefix + 'I_sinm', **params)

        params['data'] = self.I_Rnmrho
        f.create_dataset(prefix + 'I_Rnmrho', **params)

        f.create_dataset(prefix + 'K', data=np.array([self.K], dtype=int))

        f.create_dataset(prefix + 'L', data=np.array([self.L], dtype=int))

        params['data'] = self.rho_a
        f.create_dataset(prefix + 'rho_a', **params)

        params['data'] = self.rho_b
        f.create_dataset(prefix + 'rho_b', **params)

        params['data'] = self.rho_j
        f.create_dataset(prefix + 'rho_j', **params)

        params['data'] = self.theta_a
        f.create_dataset(prefix + 'theta_a', **params)

        params['data'] = self.theta_b
        f.create_dataset(prefix + 'theta_b', **params)

        params['data'] = self.theta_i
        f.create_dataset(prefix + 'theta_i', **params)

    @classmethod
    def load(cls, filename, prepend=None):
        """Load object from an HDF5 file."""
        f = h5py.File(filename, 'r')
        z = cls.load_h5py(f, prepend=prepend)
        f.close()

        return z

    @classmethod
    def load_h5py(cls, f, prepend=None):
        """Load object contents from an opened HDF5 file object."""
        myfit = cls(RZern(1), 1, 1)

        if prepend is None:
            ck = f
        else:
            ck = f[prepend]

        if 'RZern' in ck:
            myfit.z = RZern.load_h5py(f, prepend=prepend)
        elif 'CZern' in ck:
            myfit.z = CZern.load_h5py(f, prepend=prepend)
        else:
            raise NotImplementedError()

        prefix = cls.__name__ + '/'

        if prepend is not None:
            prefix = prepend + prefix

        myfit.A = f[prefix + 'A'][()]
        myfit.I_cosm = f[prefix + 'I_cosm'][()]
        myfit.I_sinm = f[prefix + 'I_sinm'][()]
        myfit.I_Rnmrho = f[prefix + 'I_Rnmrho'][()]
        myfit.K = int(f[prefix + 'K'][0])
        myfit.L = int(f[prefix + 'L'][0])
        myfit.rho_a = f[prefix + 'rho_a'][()]
        myfit.rho_b = f[prefix + 'rho_b'][()]
        myfit.rho_j = f[prefix + 'rho_j'][()]
        myfit.theta_a = f[prefix + 'theta_a'][()]
        myfit.theta_b = f[prefix + 'theta_b'][()]
        myfit.theta_i = f[prefix + 'theta_i'][()]

        return myfit
