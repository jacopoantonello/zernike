#!/usr/bin/env python
# -*- coding: utf-8 -*-

# zernike - Zernike polynomials implementation for Python
# Copyright 2016-2021 J. Antonello <jacopo@antonello.org>
#
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

import logging
import os
import sys
import unittest
from tempfile import NamedTemporaryFile
from time import time

import numpy as np
from numpy.linalg import norm
from numpy.random import normal
from zernike import CZern, FitZern, RZern


class TestZern(unittest.TestCase):
    def setUp(self):
        self.pupil = RZern(4)
        self.max_enorm = 1e-9
        self.max_ip_err = 5e-2

    def test_ntab(self):
        self.assertTrue(self.pupil.ntab.size == 15)
        expect = np.array([0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4])
        self.assertTrue(norm(self.pupil.ntab - expect) < self.max_enorm)

    def test_mtab(self):
        self.assertTrue(self.pupil.mtab.size == 15)
        expect = np.array([0, 1, -1, 0, -2, 2, -1, 1, -3, 3, 0, 2, -2, 4, -4])
        self.assertTrue(norm(self.pupil.mtab - expect) < self.max_enorm)

    def test_nm2noll(self):
        noll_indices = ((0, 0), (1, 1), (1, -1), (2, 0), (2, -2), (2, 2),
                        (3, -1), (3, 1), (3, -3), (3, 3), (4, 0), (4, 2),
                        (4, -2), (4, 4), (4, -4))
        for k, (n, m) in enumerate(noll_indices):
            self.assertEqual(RZern.nm2noll(n, m), k + 1)

    def test_rhotab(self):
        expect = [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 13.416407864999, 12.649110640674,
            12.649110640674, 3.1622776601684, 3.1622776601684, 0, 0, 0, 0, 0,
            0, 8.4852813742386, 8.4852813742386, 2.8284271247462,
            2.8284271247462, 0, 0, 0, 0, 0, 0, 0, 0, 3.4641016151378,
            2.4494897427832, 2.4494897427832, 0, 0, 0, 0, -13.416407864999,
            -9.4868329805051, -9.4868329805051, 0, 0, 0, 2, 2, 0, 0, 0,
            -5.6568542494924, -5.6568542494924, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
            -1.7320508075689, 0, 0, 0, 0, 0, 0, 2.2360679774998, 0, 0, 0, 0
        ]
        for c, v in enumerate(expect):
            i, j = c % 15, c // 15
            err = abs(self.pupil.coefnorm[i] * self.pupil.rhotab[i, j] - v)
            self.assertTrue(err < self.max_enorm)

    def test_Rnm(self):
        inrho = [
            0.19343375193909, 0.75442458148825, 0.34626071666477,
            0.41862541430271, 0.15571983035489, 0.81900058392684,
            0.62492352784628, 0.73856042161065, 0.80511242070695,
            0.067222657563179, 0.95079031650557, 0.49757701102216,
            0.75514590470052, 0.74240507067694, 0.83112957571011,
            0.1565018467502, 0.45730872055141, 0.61810049609438,
            0.93218335504705, 0.83508823555935, 0.89542351634235,
            0.58251852466803, 0.58274680636179, 0.85492593896995,
            0.034865700047931, 0.8854200954848, 0.40773081135598,
            0.036382030447832, 0.74614794273772, 0.15482877075111
        ]
        out3 = [
            -1.6024358463019, 0.2395649672106, -1.3167172040235,
            -1.1249765690963, -1.6480509660176, 0.59153676922049,
            -0.37921722803669, 0.157517884017, 0.5134006785331,
            -1.716396928352, 1.3995047634699, -0.87439854650643,
            0.24333698669463, 0.177241760149, 0.66086873705921,
            -1.6472051624093, -1.0075988516153, -0.40859694027142,
            1.2781350494377, 0.68371791511127, 1.0454079255648,
            -0.55658471812649, -0.55566323680867, 0.79985538570113,
            -1.7278397866178, 0.98369658989472, -1.1561632626913,
            -1.7274655420545, 0.19654187580572, -1.6490095429101
        ]
        out4 = [
            0.091651618055082, 1.3941428842409, 0.29368520752156,
            0.42926631070763, 0.0593968575795, 1.6430245322286,
            0.95659779790408, 1.3361268353382, 1.5877739726481,
            0.011068964146113, 2.214344179944, 0.60645172969723,
            1.3968101047969, 1.3500737219024, 1.6920496368403,
            0.059994931046482, 0.51226489069861, 0.93582320415359,
            2.1285228321212, 1.7082064455855, 1.9639599046646,
            0.83118004289954, 0.83183162858811, 1.7903280385894,
            0.0029776414702216, 1.9203234007362, 0.4072139881838,
            0.0032422723387354, 1.3637209645609, 0.058719041358534
        ]
        out10 = [
            1.7528544047959, -1.0538684824995, 0.82035094947331,
            0.29691883664534, 1.9186268081509, -0.72681628778909,
            -0.95725432162119, -1.0902934593531, -0.82334230297547,
            2.1757147314369, 1.0717625678433, -0.2632155912034,
            -1.0518319661711, -1.082911351763, -0.62973545308711,
            1.915510691432, 0.017056106891689, -0.93137658244662,
            0.7084189870086, -0.59538643046011, 0.10384045987001,
            -0.7716820983231, -0.77282798892389, -0.40275090649042,
            2.2197785892442, -0.036158358505776, 0.37645712628727,
            2.2183328268091, -1.0748457792076, 1.9221603389389
        ]
        out13 = [
            0.0044271987866637, 1.0243852641201, 0.0454582689999,
            0.09711858840943, 0.0018594122493614, 1.4227770316249,
            0.48228916269284, 0.9409014176698, 1.3286974659506,
            6.4574746661877e-05, 2.5842766270514, 0.19383902995264,
            1.0283086425533, 0.96064675024684, 1.5089503417078,
            0.0018970460208737, 0.13830501642735, 0.4615687191036,
            2.3878408401262, 1.5379048343826, 2.0328904888789,
            0.36411532970706, 0.3646864342091, 1.6893279835042,
            4.6729760834851e-06, 1.9435579666012, 0.08739651710603,
            5.540484342947e-06, 0.98016633844227, 0.0018172164647309
        ]
        outi = [3, 4, 10, 13]
        outl = [out3, out4, out10, out13]
        for c, i in enumerate(outi):
            for j in range(len(inrho)):
                a = self.pupil.coefnorm[i] * self.pupil.Rnm(i, inrho[j])
                b = outl[c][j]
                err = abs(a - b)
                self.assertTrue(err < self.max_enorm)

    def test_Zj(self):
        rho = [
            0, 0.11111111111111, 0.22222222222222, 0.33333333333333,
            0.44444444444444, 0.55555555555556, 0.66666666666667,
            0.77777777777778, 0.88888888888889, 1
        ]
        theta = [
            0, 0.69813170079773, 1.3962634015955, 2.0943951023932,
            2.7925268031909, 3.4906585039887, 4.1887902047864, 4.8869219055841,
            5.5850536063819, 6.2831853071796
        ]
        c = [
            0.54031374569774, 0.99075098653197, 0.98937811670199,
            -0.68884118291816, -0.85683760660055, 0.048390017659746,
            -0.66485340286758, 1.4527419609898, 1.3798457521444,
            0.095136418969694, -0.42712706398395, 0.51080780984162,
            -0.65627362438046, -0.12501677594865, -0.53046670030726
        ]
        z = [
            0.77833652275983, 0.085651919598506, -0.53679726827283,
            -0.98349975448282, -1.1477173085559, -0.92148475591224,
            -0.19560997786786, 1.1403260883655, 3.1979694496802,
            6.090193057073, 0.77833652275983, 0.74923429103286,
            0.92966992705808, 1.3229751043261, 1.8849792604988, 2.524009597409,
            3.1008910810612, 3.4289464416307, 3.2739961734642, 2.3543585350795,
            0.77833652275983, 1.3880071769557, 2.1381716926469,
            2.8957139599519, 3.4767967033161, 3.6468614815119, 3.1206286876391,
            1.5620975491243, -1.4154538722785, -6.2494686804882,
            0.77833652275983, 1.6882064144334, 2.5555787288029,
            3.2751193339259, 3.7304210304525, 3.7940035516249, 3.3273135632778,
            2.1807246638382, 0.19353738432559, -2.8060208116485,
            0.77833652275983, 1.5982868830036, 2.2744826664559,
            2.7560978140881, 3.0147604484415, 3.0445528736273, 2.8620115753266,
            2.5061272207906, 2.0383446588402, 1.5425629198664,
            0.77833652275983, 1.2214350933728, 1.5980879771486,
            1.7569084557018, 1.5178403949563, 0.67215824514518,
            -1.0175329591891, -3.8172975991948, -8.0218694717109,
            -13.954651789267, 0.77833652275983, 0.68350110425162,
            0.77861386147368, 0.97581235178894, 1.1341939908368,
            1.0598160525329, 0.50569566906943, -0.82819016908496,
            -3.2949046131853, -7.3005509562103, 0.77833652275983,
            0.12208636455387, -0.30877246905582, -0.48094227889849,
            -0.39896357657152, -0.10521508444048, 0.32008626436092,
            0.75888532593079, 1.0552887455991, 1.0155649579277,
            0.77833652275983, -0.16765298701795, -1.0634255875761,
            -1.8055846686938, -2.2742293038136, -2.3329542500415,
            -1.8288499481469, -0.59250252256298, 1.5620062186141,
            4.8370987836243, 0.77833652275983, 0.085651919598506,
            -0.53679726827283, -0.98349975448283, -1.1477173085559,
            -0.92148475591224, -0.19560997786786, 1.1403260883655,
            3.1979694496802, 6.090193057073
        ]
        count = 0
        for th in theta:
            for rh in rho:
                zv = 0.0
                for j, ci in enumerate(c):
                    zv += ci * self.pupil.Zk(j, rh, th)
                err = abs(zv - z[count])
                self.assertTrue(err <= self.max_enorm)
                count += 1

    def test_eval_a(self):
        rho = [
            0, 0.11111111111111, 0.22222222222222, 0.33333333333333,
            0.44444444444444, 0.55555555555556, 0.66666666666667,
            0.77777777777778, 0.88888888888889, 1
        ]
        theta = [
            0, 0.69813170079773, 1.3962634015955, 2.0943951023932,
            2.7925268031909, 3.4906585039887, 4.1887902047864, 4.8869219055841,
            5.5850536063819, 6.2831853071796
        ]
        c = [
            0.54031374569774, 0.99075098653197, 0.98937811670199,
            -0.68884118291816, -0.85683760660055, 0.048390017659746,
            -0.66485340286758, 1.4527419609898, 1.3798457521444,
            0.095136418969694, -0.42712706398395, 0.51080780984162,
            -0.65627362438046, -0.12501677594865, -0.53046670030726
        ]
        phi = [
            0.77833652275983, 0.085651919598506, -0.53679726827283,
            -0.98349975448282, -1.1477173085559, -0.92148475591224,
            -0.19560997786786, 1.1403260883655, 3.1979694496802,
            6.090193057073, 0.77833652275983, 0.74923429103286,
            0.92966992705808, 1.3229751043261, 1.8849792604988, 2.524009597409,
            3.1008910810612, 3.4289464416307, 3.2739961734642, 2.3543585350795,
            0.77833652275983, 1.3880071769557, 2.1381716926469,
            2.8957139599519, 3.4767967033161, 3.6468614815119, 3.1206286876391,
            1.5620975491243, -1.4154538722785, -6.2494686804882,
            0.77833652275983, 1.6882064144334, 2.5555787288029,
            3.2751193339259, 3.7304210304525, 3.7940035516249, 3.3273135632778,
            2.1807246638382, 0.19353738432559, -2.8060208116485,
            0.77833652275983, 1.5982868830036, 2.2744826664559,
            2.7560978140881, 3.0147604484415, 3.0445528736273, 2.8620115753266,
            2.5061272207906, 2.0383446588402, 1.5425629198664,
            0.77833652275983, 1.2214350933728, 1.5980879771486,
            1.7569084557018, 1.5178403949563, 0.67215824514518,
            -1.0175329591891, -3.8172975991948, -8.0218694717109,
            -13.954651789267, 0.77833652275983, 0.68350110425162,
            0.77861386147368, 0.97581235178894, 1.1341939908368,
            1.0598160525329, 0.50569566906943, -0.82819016908496,
            -3.2949046131853, -7.3005509562103, 0.77833652275983,
            0.12208636455387, -0.30877246905582, -0.48094227889849,
            -0.39896357657152, -0.10521508444048, 0.32008626436092,
            0.75888532593079, 1.0552887455991, 1.0155649579277,
            0.77833652275983, -0.16765298701795, -1.0634255875761,
            -1.8055846686938, -2.2742293038136, -2.3329542500415,
            -1.8288499481469, -0.59250252256298, 1.5620062186141,
            4.8370987836243, 0.77833652275983, 0.085651919598506,
            -0.53679726827283, -0.98349975448283, -1.1477173085559,
            -0.92148475591224, -0.19560997786786, 1.1403260883655,
            3.1979694496802, 6.090193057073
        ]
        count = 0
        for th in theta:
            for rh in rho:
                phi2 = self.pupil.eval_a(c, rh, th)
                err = abs(phi2 - phi[count])
                self.assertTrue(err <= self.max_enorm)
                count += 1

    def test_rhoitab(self):
        z = self.pupil
        for k in range(z.nk):
            for i in range(z.n + 1):
                a = (z.rhoitab[k, i]) * (z.n + 2 - i)
                b = z.rhotab[k, i]
                self.assertTrue(abs(a - b) < self.max_enorm)

    def test_eval_grid(self):
        log = logging.getLogger('TestZern.test_eval_grid')
        z = self.pupil
        rho_j = np.array([
            0, 0.11111111111111, 0.22222222222222, 0.33333333333333,
            0.44444444444444, 0.55555555555556, 0.66666666666667,
            0.77777777777778, 0.88888888888889, 1
        ])
        theta_i = np.array([
            0, 0.69813170079773, 1.3962634015955, 2.0943951023932,
            2.7925268031909, 3.4906585039887, 4.1887902047864, 4.8869219055841,
            5.5850536063819, 6.2831853071796
        ])
        a = normal(size=z.nk)

        t1 = time()
        Phi1 = [z.eval_a(a, rh, th) for rh in rho_j for th in theta_i]
        t2 = time()
        log.debug('list eval {:.6f}'.format(t2 - t1))

        z.make_pol_grid(rho_j, theta_i)
        t1 = time()
        Phi2 = z.eval_grid(np.array(a))
        t2 = time()
        log.debug('numpy eval {:.6f}'.format(t2 - t1))

        Phi1 = np.array(Phi1).ravel(order='F')
        Phi2 = np.array(Phi2).ravel(order='F')
        log.debug('norm(Phi1 - Phi2) {}'.format(norm(Phi1 - Phi2)))
        self.assertTrue(norm(Phi1 - Phi2) < self.max_enorm)

    def test_save_load(self):
        rho_j = np.array([
            0, 0.11111111111111, 0.22222222222222, 0.33333333333333,
            0.44444444444444, 0.55555555555556, 0.66666666666667,
            0.77777777777778, 0.88888888888889, 1
        ])
        theta_i = np.array([
            0, 0.69813170079773, 1.3962634015955, 2.0943951023932,
            2.7925268031909, 3.4906585039887, 4.1887902047864, 4.8869219055841,
            5.5850536063819, 6.2831853071796
        ])

        def do_test(cls, complex_a):
            rz1 = cls(4)

            if complex_a:
                a = normal(size=rz1.nk) + 1j * normal(size=rz1.nk)
            else:
                a = normal(size=rz1.nk)

            rz1.make_pol_grid(rho_j, theta_i)
            PhiA = rz1.eval_grid(a)

            # create tmp path
            tmpfile = NamedTemporaryFile()
            tmppath = tmpfile.name
            tmpfile.close()

            rz1.save(tmppath)

            rz2 = cls.load(tmppath)
            PhiB = rz2.eval_grid(a)

            self.assertTrue(norm(PhiA - PhiB) < self.max_enorm)
            self.assertTrue(type(rz1) == type(rz2))
            self.assertTrue(norm(rz1.coefnorm - rz2.coefnorm) == 0)
            self.assertTrue(norm(rz1.ntab - rz2.ntab) == 0)
            self.assertTrue(norm(rz1.mtab - rz2.mtab) == 0)
            self.assertTrue(rz1.n == rz2.n)
            self.assertTrue(rz1.nk == rz2.nk)
            self.assertTrue(rz1.normalise == rz2.normalise)
            self.assertTrue(norm(rz1.rhoitab - rz2.rhoitab) == 0)
            self.assertTrue(norm(rz1.rhotab - rz2.rhotab) == 0)
            self.assertTrue(rz1.numpy_dtype == rz2.numpy_dtype)
            self.assertTrue(norm(rz1.ZZ - rz2.ZZ) == 0)

            os.unlink(tmppath)
            del rz1, rz2, a, PhiA, PhiB

        do_test(RZern, False)
        do_test(CZern, True)

    def test_normalisations_real(self):
        log = logging.getLogger('TestZern.test_normalisations_real')
        n_alpha = 6
        L, K = 400, 357

        # polar grid
        pol = RZern(n_alpha)
        fitAlpha = FitZern(pol, L, K)
        t1 = time()
        pol.make_pol_grid(fitAlpha.rho_j, fitAlpha.theta_i)
        t2 = time()
        log.debug('make pol grid {:.6f}'.format(t2 - t1))

        # cartesian grid
        cart = RZern(n_alpha)
        dd = np.linspace(-1.0, 1.0, max(L, K))
        xx, yy = np.meshgrid(dd, dd)
        t1 = time()
        cart.make_cart_grid(xx, yy)
        t2 = time()
        log.debug('make cart grid {:.6f}'.format(t2 - t1))

        smap = np.isfinite(cart.eval_grid(np.zeros(cart.nk)))
        scale = (1.0 / np.sum(smap))
        log.debug('')
        log.debug('{} modes, {} x {} grid'.format(n_alpha, L, K))
        for i in range(pol.nk):
            a = np.zeros(pol.nk)
            a[i] = 1.0
            Phi_a = cart.eval_grid(a)
            for j in range(pol.nk):
                b = np.zeros(pol.nk)
                b[j] = 1.0
                Phi_b = cart.eval_grid(b)
                ip = scale * np.sum(Phi_a[smap] * Phi_b[smap])
                if i == j:
                    eip = 1.0
                else:
                    eip = 0.0
                iperr = abs(ip - eip)
                log.debug('<{:02},{:02}> = {:+e} {:+e}'.format(
                    i + 1, j + 1, ip, iperr))
                self.assertTrue(iperr < self.max_ip_err)

    def test_normalisations_complex(self):
        log = logging.getLogger('TestZern.test_normalisations_complex')
        n_beta = 6
        L, K = 400, 393

        # polar grid
        pol = CZern(n_beta)
        fitBeta = FitZern(pol, L, K)
        t1 = time()
        pol.make_pol_grid(fitBeta.rho_j, fitBeta.theta_i)
        t2 = time()
        log.debug('make pol grid {:.6f}'.format(t2 - t1))

        # cartesian grid
        cart = CZern(n_beta)
        dd = np.linspace(-1.0, 1.0, max(L, K))
        xx, yy = np.meshgrid(dd, dd)
        t1 = time()
        cart.make_cart_grid(xx, yy)
        t2 = time()
        log.debug('make cart grid {:.6f}'.format(t2 - t1))

        smap = np.isfinite(cart.eval_grid(np.zeros(cart.nk)))
        scale = (1.0 / np.sum(smap))
        log.debug('')
        log.debug('{} modes, {} x {} grid'.format(n_beta, L, K))
        for i in range(pol.nk):
            a = np.zeros(pol.nk)
            a[i] = 1.0
            Phi_a = cart.eval_grid(a)
            for j in range(pol.nk):
                b = np.zeros(pol.nk)
                b[j] = 1.0
                Phi_b = cart.eval_grid(b)
                ip = scale * np.sum(Phi_a[smap] * (Phi_b[smap].conj()))
                if i == j:
                    eip = 1.0
                else:
                    eip = 0.0
                iperr = abs(ip - eip)
                log.debug('<{:02},{:02}> = {:+e} {:+e}'.format(
                    i + 1, j + 1, ip, iperr))
                self.assertTrue(iperr < self.max_ip_err)


class TestFitZern(unittest.TestCase):
    L = 200
    K = 253
    max_fit_norm = 5e-2
    max_enorm = 1e-9

    def test_fit_real(self):
        log = logging.getLogger('TestFitZern.test_fit_real')
        z = RZern(4)
        F = FitZern(z, self.L, self.K)
        theta_i = F.theta_i
        rho_j = F.rho_j

        c = normal(size=z.nk)
        time1 = time()
        Phi = [z.eval_a(c, rh, th) for rh in rho_j for th in theta_i]
        time2 = time()
        log.debug('eval Phi {:.4f}'.format(time2 - time1))

        time1 = time()
        ce = F._fit_slow(Phi)
        time2 = time()
        log.debug('elapsed time {:.4f}'.format(time2 - time1))

        err1 = norm(np.array(c) - np.array(ce))
        max1 = max([abs(c[i] - ce[i]) for i in range(z.nk)])

        log.debug('err1 {:e} max1 {:e}'.format(err1, max1))
        self.assertTrue(err1 < self.max_fit_norm)

    def test_fit_complex(self):
        log = logging.getLogger('TestFitZern.test_fit_complex')
        z = CZern(4)
        F = FitZern(z, self.L, self.K)
        theta_i = F.theta_i
        rho_j = F.rho_j

        c = normal(size=z.nk) + 1j * normal(size=z.nk)
        time1 = time()
        Phi = [z.eval_a(c, rh, th) for rh in rho_j for th in theta_i]
        time2 = time()
        log.debug('eval Phi {:.4f}'.format(time2 - time1))

        time1 = time()
        ce = F._fit_slow(Phi)
        time2 = time()
        log.debug('elapsed time {:.4f}'.format(time2 - time1))

        err1 = np.sqrt(sum([abs(c[i] - ce[i])**2 for i in range(z.nk)]))
        max1 = max([abs(c[i] - ce[i]) for i in range(z.nk)])

        log.debug('err1 {:e} max1 {:e} max {:e}'.format(
            err1, max1, self.max_fit_norm))
        self.assertTrue(err1 < self.max_fit_norm)

    def test_fit_real_numpy(self):
        log = logging.getLogger('TestFitZern.test_fit_real_numpy')
        z = RZern(4)
        F = FitZern(z, self.L, self.K)
        theta_i = F.theta_i
        rho_j = F.rho_j

        c = normal(size=z.nk)
        Phi = [z.eval_a(c, rh, th) for rh in rho_j for th in theta_i]

        time1 = time()
        ce = F._fit_slow(Phi)
        time2 = time()
        log.debug('elapsed FIT_LIST {:.6f}'.format(time2 - time1))

        time1 = time()
        ce2 = F.fit(np.array(Phi, order='F'))
        time2 = time()
        log.debug('elapsed FIT_NUMPY {:.6f}'.format(time2 - time1))

        enorm = norm(ce2 - np.array(ce, order='F'))
        log.debug('enorm {:e}'.format(enorm))
        self.assertTrue(enorm < self.max_enorm)

    def test_fit_complex_numpy(self):
        log = logging.getLogger('TestFitZern.test_fit_complex_numpy')
        z = CZern(4)
        F = FitZern(z, self.L, self.K)
        theta_i = F.theta_i
        rho_j = F.rho_j

        c = normal(size=z.nk) + 1j * normal(size=z.nk)
        Phi = [z.eval_a(c, rh, th) for rh in rho_j for th in theta_i]

        time1 = time()
        ce = F._fit_slow(Phi)
        time2 = time()
        log.debug('elapsed FIT_LIST {:.6f}'.format(time2 - time1))

        PhiN = np.array(Phi, order='F')
        time1 = time()
        ce2 = F.fit(PhiN)
        time2 = time()
        log.debug('elapsed FIT_NUMPY {:.6f}'.format(time2 - time1))

        enorm = norm(ce2 - np.array(ce, order='F'))
        log.debug('enorm {:e}'.format(enorm))
        self.assertTrue(enorm < self.max_enorm)

    def test_fit_save_load(self):
        def do_test(cls, complex_a):
            z = cls(4)
            F = FitZern(z, self.L, self.K)
            theta_i = F.theta_i
            rho_j = F.rho_j
            z.make_pol_grid(rho_j, theta_i)

            if complex_a:
                c = normal(size=z.nk) + 1j * normal(size=z.nk)
            else:
                c = normal(size=z.nk)

            PhiN = np.array(
                [z.eval_a(c, rh, th) for rh in rho_j for th in theta_i],
                order='F')

            ce1 = F.fit(PhiN)

            # create tmp path
            tmpfile = NamedTemporaryFile()
            tmppath = tmpfile.name
            tmpfile.close()

            F.save(tmppath)

            F2 = FitZern.load(tmppath)

            PhiN1 = np.array(
                [F.z.eval_a(c, rh, th) for rh in rho_j for th in theta_i],
                order='F')
            PhiN2 = np.array(
                [F2.z.eval_a(c, rh, th) for rh in rho_j for th in theta_i],
                order='F')

            self.assertTrue(isinstance(F, FitZern))
            self.assertTrue(isinstance(F2, FitZern))
            if complex_a:
                self.assertTrue(isinstance(F.z, CZern))
                self.assertTrue(isinstance(F2.z, CZern))
            else:
                self.assertTrue(isinstance(F.z, RZern))
                self.assertTrue(isinstance(F2.z, RZern))
            self.assertTrue(norm(PhiN - PhiN1) == 0)
            self.assertTrue(norm(PhiN - PhiN2) == 0)

            self.assertTrue(norm(F.z.coefnorm - F2.z.coefnorm) == 0)
            self.assertTrue(norm(F.z.ntab - F2.z.ntab) == 0)
            self.assertTrue(norm(F.z.mtab - F2.z.mtab) == 0)
            self.assertTrue(F.z.n == F2.z.n)
            self.assertTrue(F.z.nk == F2.z.nk)
            self.assertTrue(F.z.normalise == F2.z.normalise)
            self.assertTrue(norm(F.z.rhoitab - F2.z.rhoitab) == 0)
            self.assertTrue(norm(F.z.rhotab - F2.z.rhotab) == 0)
            self.assertTrue(F.z.numpy_dtype == F2.z.numpy_dtype)
            self.assertTrue(norm(F.z.ZZ - F2.z.ZZ) == 0)

            self.assertTrue(norm(F.A - F2.A) == 0)
            self.assertTrue(norm(F.I_cosm - F2.I_cosm) == 0)
            self.assertTrue(norm(F.I_sinm - F2.I_sinm) == 0)
            self.assertTrue(norm(F.I_Rnmrho - F2.I_Rnmrho) == 0)
            self.assertTrue(F.K == F2.K)
            self.assertTrue(F.L == F2.L)
            self.assertTrue(norm(F.rho_a - F2.rho_a) == 0)
            self.assertTrue(norm(F.rho_b - F2.rho_b) == 0)
            self.assertTrue(norm(F.rho_j - F2.rho_j) == 0)
            self.assertTrue(norm(F.theta_a - F2.theta_a) == 0)
            self.assertTrue(norm(F.theta_b - F2.theta_b) == 0)
            self.assertTrue(norm(F.theta_i - F2.theta_i) == 0)

            ce2 = F2.fit(PhiN)
            self.assertTrue(norm(ce2 - ce1) < self.max_enorm)

            os.unlink(tmppath)
            del z, F, F2, c, ce1, ce2, PhiN, PhiN1, PhiN2

        do_test(RZern, False)
        do_test(CZern, True)


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    unittest.main()
