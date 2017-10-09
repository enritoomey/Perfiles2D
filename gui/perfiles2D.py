__author__ = 'Enriquito'
import gui.layout_lecturaPerfil as layout
import perfiles.airfoil_characteristics as profile
from PySide.QtCore import *
from PySide.QtGui import *
import sys
import matplotlib #Para los graficos
matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4']='PySide'
# Estas lineas son un poco misteriosas, pero son las que me permite
# vincular los Widget de Qt/Pyside con matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

class MainDialog(QDialog, layout.Ui_Dialog):

    def __init__(self, parent=None):
        super(MainDialog, self).__init__(parent)
        self.setupUi(self)

        # Generamos dos figuras, cada una luego asociada a un canvas, que a su vez tiene como padre una de las pestanias
        # self.tab -> contiene la pestania titulada "Diagrama P-S"
        # self.tab_2 -> contiene la pestania titulada "Diagrama T-S"
        self.fig1 = Figure(dpi=72, facecolor=(1, 1, 1), edgecolor=(0, 0, 0))#figsize=(4.8, 3.4),
        self.axes1 = self.fig1.add_subplot(111)
        self.axes1.set_ylabel('Cl')
        self.axes1.set_xlabel('AoA')
        self.fig2 = Figure(dpi=72, facecolor=(1, 1, 1), edgecolor=(0, 0, 0))# figsize=(4.8, 3.4),
        self.axes2 = self.fig2.add_subplot(111)
        self.axes2.set_ylabel('Cd')
        self.axes2.set_xlabel('Cl')
        self.fig3 = Figure(dpi=72, facecolor=(1, 1, 1), edgecolor=(0, 0, 0))
        self.axes3 = self.fig3.add_subplot(1,1,1)
        self.axes3.set_ylabel('Cm0')
        self.axes3.set_xlabel('AoA')
        # generate the canvas to display the plot
        self.canvas1 = FigureCanvas(self.fig1)
        self.canvas1.setParent(self.frame_ClAlpha)
        self.canvas1.resize(self.frame_ClAlpha.width(), self.frame_ClAlpha.height())
        self.canvas1.show()
        self.canvas2 = FigureCanvas(self.fig2)
        self.canvas2.setParent(self.frame_CdCl)
        self.canvas2.resize(self.frame_CdCl.width(), self.frame_ClAlpha.height())
        self.canvas2.show()
        self.canvas3 = FigureCanvas(self.fig3)
        self.canvas3.setParent(self.frame_CmAlpha)
        self.canvas3.resize(self.frame_CmAlpha.width(), self.frame_CmAlpha.height())
        self.canvas3.show()
        #Ponemos el primer checkbox en "checked"
        self.checkBox_Re3.setCheckState(Qt.Checked)
        # Cargamos los numeros de Reynolds
        self.Reynolds = ['Re3', 'Re6', 'Re9', 'std']
        # Cargamos los parametros caracteristicos del perfil
        self.vHeader = ['CL_max', 'beta_max',
                    'CL_betamax', 'alpha_betamax', 'dbeta_dalpha', 'b_max',
                    'CD_min', 'CL_CD_max', 'cuspide', 'CM0']
        # Cargamos el formato en que se deben mostrar cada parametro caracteristico
        self.HeadersFormat = ['{:.2f}', '{:.1f}', '{:.2f}', '{:.1f}', '{:.3f}', '{:.1f}',
                              '{:.5f}', '{:.1f}', '{:.2}', '{:.2}']
        self.tableDatos.setVerticalHeaderLabels(self.vHeader)
        self.tableDatos.setHorizontalHeaderLabels(self.Reynolds)


        self.connect(self.openButton, SIGNAL("clicked()"), self.Open)
        self.connect(self.okButton, SIGNAL("clicked()"), self, SLOT('accept()'))
        self.connect(self.cancelButton, SIGNAL("clicked()"), self, SLOT('reject()'))
        self.connect(self.checkBox_Re3, SIGNAL("toggled(bool)"), self.plotGraph)
        self.connect(self.checkBox_Re6, SIGNAL("toggled(bool)"), self.plotGraph)
        self.connect(self.checkBox_Re9, SIGNAL("toggled(bool)"), self.plotGraph)
        self.connect(self.checkBox_SR, SIGNAL("toggled(bool)"), self.plotGraph)

    def Open(self):
        dir = '.'
        fileAdress = QFileDialog.getOpenFileName(self,"Selecionar Perfil",dir,filter="Text File (*.txt)")
        self.airfoil = profile.Airfoil(fileAdress[0])
        #print self.airfoil.AIRFOIL_DATA

        for j in range(len(self.Reynolds)):
            for i in range(len(self.vHeader)):
                if self.vHeader[i] in self.airfoil.AIRFOIL_DATA[self.Reynolds[j]].keys():
                    item = QTableWidgetItem(self.HeadersFormat[i].format(self.airfoil.AIRFOIL_DATA[self.Reynolds[j]][self.vHeader[i]]))
                    self.tableDatos.setItem(i, j, item)
        nombrePerfil = fileAdress[0].split('/')[-1].split('.')[0]
        #print nombrePerfil
        self.labelPerfil.setText('Datos Perfil:  NACA '+nombrePerfil)
        self.plotGraph()

    def plotGraph(self):
        self.axes1.clear()
        self.axes2.clear()
        self.axes3.clear()

        if self.checkBox_Re3.isChecked():
            self.plotReynolds('Re3', 'b')

        if self.checkBox_Re6.isChecked():
            self.plotReynolds('Re6', 'r')

        if self.checkBox_Re9.isChecked():
            self.plotReynolds('Re9', 'g')

        if self.checkBox_SR.isChecked():
            self.plotReynolds('std', 'y')

        self.canvas1.draw()
        self.canvas2.draw()
        self.canvas3.draw()

    def plotReynolds(self, reynolds, color):
        x1 = [value[0] for value in self.airfoil.AIRFOIL_DATA[reynolds]['AoA_Cl']]
        y1 = [value[1] for value in self.airfoil.AIRFOIL_DATA[reynolds]['AoA_Cl']]
        x2 = [value[0] for value in self.airfoil.AIRFOIL_DATA[reynolds]['Cl_Cd']]
        y2 = [value[1] for value in self.airfoil.AIRFOIL_DATA[reynolds]['Cl_Cd']]
        x3 = [value[0] for value in self.airfoil.AIRFOIL_DATA[reynolds]['AoA_Cm']]
        y3 = [value[1] for value in self.airfoil.AIRFOIL_DATA[reynolds]['AoA_Cm']]
        self.axes1.plot(x1, y1, color)
        self.axes2.plot(x2, y2, color)
        self.axes3.plot(x3, y3, color)

app = QApplication(sys.argv)
form = MainDialog()
form.show()
app.exec_()