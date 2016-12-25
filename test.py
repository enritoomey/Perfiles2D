import unittest
from profile_characteristics import Airfoil

class AirfoilTesting(unittest.TestCase):

    def setUp(self):
        self.airfoil_name = 'TR824-Digitized/63(420)-422.txt'
        self.airfoil = Airfoil(self.airfoil_name, reynold=3e6)

    def testAirfoilInit(self):
        self.assertEqual(self.airfoil.reynold, 'Re_3')
        self.assertAlmostEqual(self.airfoil.AIRFOIL_DATA[self.airfoil.reynold]['CL_max'], 1.331, places=2)

if __name__ == '__main__':
    unittest.main()