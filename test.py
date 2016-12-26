import unittest
from profile_characteristics import Airfoil

class AirfoilTesting(unittest.TestCase):

    def setUp(self):
        self.airfoil_name = 'TR824-Digitized/63(420)-422.txt'
        self.airfoil = Airfoil(self.airfoil_name, reynold=3e6)

    def testAirfoilInit(self):
        self.assertEqual(self.airfoil.reynold, 'Re_3')
        self.assertAlmostEqual(self.airfoil.AIRFOIL_DATA[self.airfoil.reynold]['CL_max'], 1.331, places=2)

    def testAoACl(self):
        for alpha, cl in self.airfoil.AIRFOIL_DATA[self.airfoil.reynold]['AoA_Cl']:
            self.assertAlmostEqual(self.airfoil.cl_aoa(alpha), cl)

    def testClCd(self):
        for cl, cd in self.airfoil.AIRFOIL_DATA[self.airfoil.reynold]['Cl_Cd']:
            self.assertAlmostEqual(self.airfoil.cd_cl(cl), cd)

    def testAoACm0(self):
        for alpha, cm0 in self.airfoil.AIRFOIL_DATA[self.airfoil.reynold]['AoA_Cm']:
            self.assertAlmostEqual(self.airfoil.cm_aoa(alpha), cm0)


if __name__ == '__main__':
    unittest.main()