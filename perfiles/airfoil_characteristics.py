# -*- coding: iso-8859-1 -*-
from scipy import interpolate, polyfit
from scipy.optimize import fmin, fsolve
from scipy.misc import derivative
import logging
import numpy as np

logger = logging.getLogger(__name__)
ch = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s -%(funcName)s: %(message)s ')
ch.setFormatter(formatter)
ch.setLevel(logging.WARNING)
logger.addHandler(ch)

class Airfoil:

    class DataNotAvailableError(Exception):
        pass

    def __init__(self, file_name, reynold=4e6):
        logger.debug("Initializing airfoil %r with Reynolds %r", file_name, reynold)
        self.file_name = file_name
        self.reynolds = {'Re3':3e6, 'Re6':6e6, 'Re9':9e6, 'std':0}
        self.characteristics = ['CL_max', 'beta_max',
         'CL_betamax', 'alpha_betamax', 'dbeta_dalpha', 'b_max',
         'CD_min', 'CL_CD_max', 'cuspide']
        self.AIRFOIL_DATA = dict()
        for re in self.reynolds.keys():
            self.AIRFOIL_DATA[re] = {'AoA_Cl': [], 'AoA_Cm': [], 'Cl_Cd': []}
        self.lectura_perfiles(file_name)
        self.build_airfoil_functions()
        self.reynold_key = self.get_reynold_key(reynold)
        self.reynold_number = reynold
        self.get_airfoil_characteristics()

    def lectura_perfiles(self, file_name):
        aux_file = open(file_name)
        polar_data = aux_file.readlines()
        aux_file.close()
        i=0
        block = 1
        while block < 3:
            b = 0
            if polar_data[i].split()[0] == 'Data':
                i += 1
                data_set = float(polar_data[i])
                i += 1
                while b == 0:
                    if polar_data[i].split()[-1] == 'exist.':
                        i += 1
                    elif polar_data[i].split()[0] == 'Data':
                        b = 1
                    elif polar_data[i].split()[0] == 'Found':
                        block += 1
                        i += 1
                        b = 1
                    else:
                        self.build_data(block, data_set, polar_data[i])
                        i += 1
            else:
                i += 1

    def build_airfoil_functions(self):
        self._build_cl_aoa()
        self._build_cd_cl()
        self._build_cm_aoa()

    def _build_cl_aoa(self):
        aoa_l = [[aoa for aoa, cl in self.AIRFOIL_DATA["Re{}".format(re)]["AoA_Cl"]] for re in [3, 6, 9]]
        cl_l = [[cl for aoa, cl in self.AIRFOIL_DATA["Re{}".format(re)]["AoA_Cl"]] for re in [3, 6, 9]]
        re_l = [[re for i in range(len(aoa))] for aoa, re in zip(aoa_l, [3e6, 6e6, 9e6])]
        aoa_points = aoa_l[0] + aoa_l[1] + aoa_l[2]
        cl_points = cl_l[0] + cl_l[1] + cl_l[2]
        re_points = re_l[0] + re_l[1] + re_l[2]
        bisplrep = interpolate.bisplrep(aoa_points, re_points, cl_points, s=0.5, kx=4, ky=2)
        self.cl_aoa_func = bisplrep

    def _build_cd_cl(self):
        cl_l = [[cl for cl, cd in self.AIRFOIL_DATA["Re{}".format(re)]["Cl_Cd"]] for re in [3, 6, 9]]
        cd_l = [[cd for cl, cd in self.AIRFOIL_DATA["Re{}".format(re)]["Cl_Cd"]] for re in [3, 6, 9]]
        re_l = [[re for i in range(len(cl))] for cl, re in zip(cl_l, [3e6, 6e6, 9e6])]
        cl_points = cl_l[0] + cl_l[1] + cl_l[2]
        cd_points = cd_l[0] + cd_l[1] + cd_l[2]
        re_points = re_l[0] + re_l[1] + re_l[2]
        bisplrep = interpolate.bisplrep(cl_points, re_points, cd_points, s=0.5, kx=4, ky=2)
        self.cd_cl_func = bisplrep

    def _build_cm_aoa(self):
        aoa_l = [[aoa for aoa, cl in self.AIRFOIL_DATA["Re{}".format(re)]["AoA_Cm"]] for re in [3, 6, 9]]
        cm_l = [[cl for aoa, cl in self.AIRFOIL_DATA["Re{}".format(re)]["AoA_Cm"]] for re in [3, 6, 9]]
        re_l = [[re for i in range(len(aoa))] for aoa, re in zip(aoa_l, [3e6, 6e6, 9e6])]
        aoa_points = aoa_l[0] + aoa_l[1] + aoa_l[2]
        cm_points = cm_l[0] + cm_l[1] + cm_l[2]
        re_points = re_l[0] + re_l[1] + re_l[2]
        bisplrep = interpolate.bisplrep(aoa_points, re_points, cm_points, s=0.5, kx=4, ky=2)
        self.cm_aoa_func = bisplrep

    def cl_aoa(self, alpha, reynold=None):
        if reynold == None:
            reynold = self.reynold_number
        return interpolate.bisplev(alpha, reynold, self.cl_aoa_func)

    def cd_cl(self, cl, reynold=None):
        if reynold == None:
            reynold = self.reynold_number
        return interpolate.bisplev(cl, reynold, self.cd_cl_func)

    def cm_aoa(self, alpha, reynold=None):
        if reynold == None:
            reynold = self.reynold_number
        return interpolate.bisplev(alpha, reynold, self.cm_aoa_func)

    def cl_aoa_data(self, alpha):
        logger.debug("cl_sos called with aoa = %r deg", alpha)
        aoa_l = [aoa for aoa, cl in self.AIRFOIL_DATA[self.reynold_key]['AoA_Cl']]
        cl_l = [cl for aoa, cl in self.AIRFOIL_DATA[self.reynold_key]['AoA_Cl']]
        if len(cl_l) == 0 or len(aoa_l) == 0:
            raise self.DataNotAvailableError
        aoa_l.reverse()
        cl_l.reverse()
        try:
            spline = interpolate.splrep(aoa_l, cl_l, s=0)
            cl_out = interpolate.splev(alpha, spline, der=0)
        except ValueError:
            logger.warning("Fail to interpolate cl vs alpha curve with spline for Re=%r",self.reynold_number)
            fit = interpolate.interp1d(aoa_l, cl_l, fill_value="extrapolate")
            cl_out = fit(alpha)
        return cl_out

    def cd_cl_data(self, cl_input):
        cl_l = [cl for cl, cd in self.AIRFOIL_DATA[self.reynold_key]['Cl_Cd']]
        cd_l = [cd for cl, cd in self.AIRFOIL_DATA[self.reynold_key]['Cl_Cd']]
        if len(cl_l) == 0 or len(cd_l) == 0:
            raise self.DataNotAvailableError
        cl_l.reverse()
        cd_l.reverse()
        try:
            spline = interpolate.splrep(cl_l, cd_l, s=0)
            cd_out = interpolate.splev(cl_input, spline, der=0)
        except ValueError:
            logger.warning("Fail to interpolate cd vs cl curve with spline for Re=%r", self.reynold_number)
            fit = interpolate.interp1d(cl_l, cd_l, fill_value="extrapolate")
            cd_out = fit(cl_input)
        return cd_out

    def cd_aoa(self, alpha):
        return self.cd_cl(self.cl_aoa(alpha))

    def cm_aoa_data(self, alpha):
        aoa_l = [aoa for aoa, cm in self.AIRFOIL_DATA[self.reynold_key]['AoA_Cm']]
        cm_l = [cm for aoa, cm in self.AIRFOIL_DATA[self.reynold_key]['AoA_Cm']]
        if len(aoa_l) == 0 or len(cm_l) == 0:
            raise self.DataNotAvailableError
        aoa_l.reverse()
        cm_l.reverse()
        spline = interpolate.splrep(aoa_l, cm_l, s=0)
        return interpolate.splev(alpha, spline, der=0)

    def beta(self, cl):
        return cl / self.cd_cl(cl)

    def b(self, cl):
        return cl ** 1.5 / self.cd_cl(cl)

    def _cl_max_func(self):
        alpha_0 = 0.0
        f = lambda alpha: -self.cl_aoa(alpha, self.reynold_number)
        return fmin(f, x0=alpha_0, full_output=True, disp=False)

    @property
    def a2d(self):
        alphas = np.arange(start=-5, stop=10, step=0.5)
        cls = np.array([self.cl_aoa(alpha) for alpha in alphas])
        rta = polyfit(alphas, cls, 1)
        return rta[0]

    @property
    def alpha0(self):
        rta = fsolve(self.cl_aoa, x0=0.0)
        return rta[0]

    @property
    def cl_max(self):
        return -self._cl_max_func()[1]

    @property
    def alpha_cl_max(self):
        return self._cl_max_func()[0]

    def _beta_max_func(self):
        cl0 = 1.0
        f = lambda cl: -1.0 * self.beta(cl)
        return fmin(f, x0=cl0, full_output=True, disp=False)

    @property
    def beta_max(self):
        return -self._beta_max_func()[1]

    @property
    def cl_beta_max(self):
        return float(self._beta_max_func()[0])

    @property
    def alpha_beta_max(self):
        f = lambda alpha: self.cl_aoa(alpha) - self.cl_beta_max
        rta = fsolve(f, x0=self.cl_beta_max*180.0/2.0/np.pi**2)
        return rta[0]

    @property
    def cd_min(self):
        cl0 = 0.0
        f = lambda cl: self.cd_cl(cl)
        rta = fmin(f, x0=cl0, full_output=True, disp=False)
        return rta[1]

    def _b_max_func(self):
        cl0 = 1.0
        f = lambda cl: - self.b(cl)
        return fmin(f, x0=cl0, full_output=True, disp=False)

    @property
    def b_max(self):
        return -self._b_max_func()[1]

    @property
    def cl_b_max(self):
        return self._b_max_func()[0]

    @property
    def dbeta_dalpha(self):
        """"Return the slope of beta=f(alpha) at beta max"""
        f = lambda alpha: self.beta(self.cl_aoa(alpha))
        return derivative(f, x0=self.alpha_beta_max, dx=1e-4, n=1)

    @property
    def cuspide(self):
        """ Return the radius if curvature of cl=f(alpha) at cl_max.
        The biggest this value is, the smother the stall is """
        dcl_dalpha1 = derivative(self.cl_aoa, x0=self.alpha_cl_max, dx=0.001, n=1)
        dcl_dalpha2 = derivative(self.cl_aoa, x0=self.alpha_cl_max, dx=0.001, n=2)
        if dcl_dalpha2 == 0:
            return float('Inf')
        return float((1+dcl_dalpha1**2)**1.5/abs(dcl_dalpha2))

    @property
    def cm0(self):
        """ Return the moment coefficient around the leading edge of the airfoil
        for the linear part, assuming it is constant."""
        cm0_l = [self.cm_aoa(aoa, self.reynold_number) for aoa in np.linspace(-5, 10, 15)]
        return np.mean(cm0_l)

    def get_airfoil_characteristics(self):
        for re in self.reynolds:
            self.reynold_key = re
            self.reynold_number = self.reynolds[re]
            try:
                self.AIRFOIL_DATA[re]['CL_max'] = self.cl_max
                self.AIRFOIL_DATA[re]['cuspide'] = self.cuspide
            except self.DataNotAvailableError:
                logger.warning("Not cl vs alpha data available for reynolds %r", re)
                for ch in self.characteristics: self.AIRFOIL_DATA[re][ch] = float('Nan')
            else:
                try:
                    self.AIRFOIL_DATA[re]['beta_max'] = self.beta_max
                    self.AIRFOIL_DATA[re]['CL_betamax'] = self.cl_beta_max
                    self.AIRFOIL_DATA[re]['alpha_betamax'] = self.alpha_beta_max
                    self.AIRFOIL_DATA[re]['dbeta_dalpha'] = self.dbeta_dalpha
                    self.AIRFOIL_DATA[re]['b_max'] = self.b_max
                    self.AIRFOIL_DATA[re]['CD_min'] = self.cd_min
                    self.AIRFOIL_DATA[re]['CL_CD_max'] = self.cl_max / self.cd_min

                except self.DataNotAvailableError:
                    logger.warning("Not cd vs cl data available for reynolds %r", re)
                    gen = (ch for ch in self.characteristics if ch not in ['Cl_max', 'cuspide'])
                    for ch in gen: self.AIRFOIL_DATA[re][ch]=float('Nan')

            try:
                self.AIRFOIL_DATA[re]['CM0'] = self.cm0
            except self.DataNotAvailableError:
                logger.warning("Not cm vs aoa data available for reynolds %r", re)
                self.AIRFOIL_DATA[re]['CM0'] = float('Nan')

    def build_data(self, block, data_set, polar_data):
        x_data = float(polar_data.split()[0])
        y_data = float(polar_data.split()[1])
        if block == 1:
            if data_set == 1:
                self.AIRFOIL_DATA['Re3']['AoA_Cl'].append((x_data, y_data))
            elif data_set == 2:
                self.AIRFOIL_DATA['Re6']['AoA_Cl'].append((x_data, y_data))
            elif data_set == 3:
                self.AIRFOIL_DATA['Re9']['AoA_Cl'].append((x_data, y_data))
            elif data_set == 4:
                self.AIRFOIL_DATA['std']['AoA_Cl'].append((x_data, y_data))
            elif data_set == 7:
                self.AIRFOIL_DATA['Re3']['AoA_Cm'].append((x_data, y_data))
            elif data_set == 8:
                self.AIRFOIL_DATA['Re6']['AoA_Cm'].append((x_data, y_data))
            elif data_set == 9:
                self.AIRFOIL_DATA['Re9']['AoA_Cm'].append((x_data, y_data))
            elif data_set == 10:
                self.AIRFOIL_DATA['std']['AoA_Cm'].append((x_data, y_data))

        elif block == 2:
            if data_set == 1:
                self.AIRFOIL_DATA['Re3']['Cl_Cd'].append((x_data, y_data))
            elif data_set == 2:
                self.AIRFOIL_DATA['Re6']['Cl_Cd'].append((x_data, y_data))
            elif data_set == 3:
                self.AIRFOIL_DATA['Re9']['Cl_Cd'].append((x_data, y_data))
            elif data_set == 4:
                self.AIRFOIL_DATA['std']['Cl_Cd'].append((x_data, y_data))

    def get_reynold_key(self, reynold):
        try:
            Re = float(reynold)
            if Re <= 3e6:
                return 'Re3'
            elif Re <= 6e6:
                return 'Re6'
            else:
                return 'Re9'
        except ValueError:
            return 'std'
