__author__ = 'Enriquito'
import atmosfera_estandar
import layout_atmosfera_estandar as layout
from PySide.QtCore import *
from PySide.QtGui import *
import sys

class MainDialog(QDialog, layout.Ui_Dialog):

    def __init__(self,parent=None):
        super(MainDialog,self).__init__(parent)
        self.setupUi(self)

        self.atmosfera = {'h':[],'deltaT':[],'p':[],'t':[],'rho':[],'mu':[],'Vson':[],'calculo':[]}
        self.R = 287.00
        self.gamma = 1.4

        self.atmosfera['h'] = float(self.lineEdit_h.text())
        self.atmosfera['deltaT'] = float(self.lineEdit_deltaT.text())
        self.atmosfera['p'] = float(self.lineEdit_p.text())
        self.atmosfera['t'] = float(self.lineEdit_t.text())
        self.atmosfera['rho'] = float(self.lineEdit_rho.text())
        self.atmosfera['mu'] = float(self.lineEdit_mu.text())
        self.atmosfera['Vson'] = float(self.lineEdit_Vson.text())

        self.connect(self.lineEdit_h,SIGNAL("editingFinished()"),self.actualizarH)
        self.connect(self.lineEdit_deltaT,SIGNAL("editingFinished()"),self.actualizarH)
        self.connect(self.lineEdit_t,SIGNAL("editingFinished()"),self.actualizarT)
        self.connect(self.lineEdit_p,SIGNAL("editingFinished()"),self.actualizarP)
        self.connect(self.lineEdit_rho,SIGNAL("editingFinished()"),self.actualizarRho)


    def actualizarH(self):
        calculo = 'altura'
        results = atmosfera_estandar.atmosfera_estandar(calculo,float(self.lineEdit_h.text()),float(self.lineEdit_deltaT.text()))
        self.atmosfera['h'] = results[0]
        self.lineEdit_h.setText(str(round(results[0],1)))
        self.atmosfera['deltaT'] = results[1]
        self.lineEdit_deltaT.setText(str(round(results[1],2)))
        self.atmosfera['p'] = results[2]
        self.lineEdit_p.setText(str(round(results[2],1)))
        self.atmosfera['t'] = results[3]
        self.lineEdit_t.setText(str(round(results[3],2)))
        self.atmosfera['rho'] = results[4]
        self.lineEdit_rho.setText(str(round(results[4],3)))
        self.atmosfera['mu'] = results[5]
        self.lineEdit_mu.setText(str(round(results[5],7)))
        self.atmosfera['Vson'] = results[6]
        self.lineEdit_Vson.setText(str(round(results[6],2)))

    def actualizarP(self):
        calculo = 'presion'
        results = atmosfera_estandar.atmosfera_estandar(calculo,float(self.lineEdit_p.text()),float(self.lineEdit_deltaT.text()))
        self.atmosfera['h'] = results[0]
        self.lineEdit_h.setText(str(round(results[0],1)))
        self.atmosfera['deltaT'] = results[1]
        self.lineEdit_deltaT.setText(str(round(results[1],2)))
        self.atmosfera['p'] = results[2]
        self.lineEdit_p.setText(str(round(results[2],1)))
        self.atmosfera['t'] = results[3]
        self.lineEdit_t.setText(str(round(results[3],2)))
        self.atmosfera['rho'] = results[4]
        self.lineEdit_rho.setText(str(round(results[4],3)))
        self.atmosfera['mu'] = results[5]
        self.lineEdit_mu.setText(str(round(results[5],7)))
        self.atmosfera['Vson'] = results[6]
        self.lineEdit_Vson.setText(str(round(results[6],2)))

    def actualizarT(self):
        calculo = 'temperatura'
        temp = float(self.lineEdit_t.text())
        self.atmosfera['deltaT']=self.atmosfera['deltaT']+temp-self.atmosfera['t']
        self.lineEdit_deltaT.setText(str(round(self.atmosfera['deltaT'],2)))
        self.atmosfera['t']=temp
        self.atmosfera['rho']=self.atmosfera['p']/self.atmosfera['t']/self.R
        self.lineEdit_rho.setText(str(round(self.atmosfera['rho'],3)))

    def actualizarRho(self):
        calculo = 'densidad'
        self.atmosfera['rho']=float(round(self.lineEdit_rho.text(),3))
        temp = self.atmosfera['p']/self.atmosfera['rho']/self.R
        self.atmosfera['deltaT']=self.atmosfera['deltaT']+temp-self.atmosfera['t']
        self.lineEdit_deltaT.setText(str(round(self.atmosfera['deltaT'],2)))
        self.atmosfera['t']=temp
        self.lineEdit_t.setText(str(round(temp,2)))


app = QApplication(sys.argv)
form = MainDialog()
form.show()
app.exec_()
