import unittest
from perfiles.airfoil_characteristics import Airfoil
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
# logger.setLevel(logging.INFO)

class AirfoilTesting(unittest.TestCase):

    def setUp(self):
        self.airfoil_name = 'TR824-Digitized/63(420)-422.txt'
        self.clmax = 1.302
        self.cdmin = 0.00622851
        self.beta_max = 144.0134
        self.cl_beta_max = 0.97959964
        self.alpha_beta_max = 5.586111690
        self.b_max = 142.9628
        self.cl_b_max = 0.9911587
        self.dbeta_dalpha = -0.020153690
        self.airfoil = Airfoil(self.airfoil_name, reynold=3e6)


    def testAirfoilInit(self):
        self.assertEqual(self.airfoil.reynold, 'Re_3')

    def testAoACl(self):
        for alpha, cl in self.airfoil.AIRFOIL_DATA[self.airfoil.reynold]['AoA_Cl']:
            self.assertAlmostEqual(self.airfoil.cl_aoa(alpha), cl)

    def testClCd(self):
        for cl, cd in self.airfoil.AIRFOIL_DATA[self.airfoil.reynold]['Cl_Cd']:
            self.assertAlmostEqual(self.airfoil.cd_cl(cl), cd)

    def testAoACm0(self):
        for alpha, cm0 in self.airfoil.AIRFOIL_DATA[self.airfoil.reynold]['AoA_Cm']:
            self.assertAlmostEqual(self.airfoil.cm_aoa(alpha), cm0)

    def test_clmax(self):
        logger.info("Cl max = %r", self.airfoil.cl_max)
        self.assertAlmostEqual(self.airfoil.cl_max, self.clmax, places=2)

    def test_beta_max(self):
        logger.info("Beta max = %r",self.airfoil.beta_max)
        self.assertAlmostEqual(self.airfoil.beta_max, self.beta_max, places=2)

    def test_cl_beta_max(self):
        logger.info("Cl corresponding to beta max = %r", self.airfoil.cl_beta_max)
        self.assertAlmostEqual(self.airfoil.cl_beta_max, self.cl_beta_max, places=2)

    def test_alpha_beta_max(self):
        logger.info("alpha corresponding to beta max = %r", self.airfoil.alpha_beta_max)
        self.assertAlmostEqual(self.airfoil.alpha_beta_max, self.alpha_beta_max, places=2)

    def test_cd_min(self):
        logger.info("Cd min = %r", self.airfoil.cd_min)
        self.assertAlmostEqual(self.airfoil.cd_min, self.cdmin, places=2)

    def test_b_max(self):
        logger.info("b max = %r", self.airfoil.b_max)
        self.assertAlmostEqual(self.airfoil.b_max, self.b_max, places=2)


if __name__ == '__main__':
    unittest.main()