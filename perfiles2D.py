__author__ = 'Enriquito'
import layout_lecturaPerfil as layout
import profile_characteristics as profile
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

    def __init__(self,parent=None):
        super(MainDialog,self).__init__(parent)
        self.setupUi(self)
        self.AIRFOIL_DATA = None

        # Generamos dos figuras, cada una luego asociada a un canvas, que a su vez tiene como padre una de las pestanias
        # self.tab -> contiene la pestania titulada "Diagrama P-S"
        # self.tab_2 -> contiene la pestania titulada "Diagrama T-S"
        self.fig1 = Figure(figsize=(4.8,3.4),dpi=72, facecolor=(1,1,1), edgecolor=(0,0,0))
        self.axes1 = self.fig1.add_subplot(111)
        self.axes1.set_ylabel('Cl')
        self.axes1.set_xlabel('AoA')
        self.fig2 = Figure(figsize=(4.8,3.4),dpi=72, facecolor=(1,1,1), edgecolor=(0,0,0))
        self.axes2 = self.fig2.add_subplot(111)
        self.axes2.set_ylabel('Cd')
        self.axes2.set_xlabel('Cl')
        # generate the canvas to display the plot
        self.canvas1 = FigureCanvas(self.fig1)
        self.canvas1.setParent(self.frame_ClAlpha)
        self.canvas1.show()
        self.canvas2 = FigureCanvas(self.fig2)
        self.canvas2.setParent(self.frame_CdCl)
        self.canvas2.show()
        #Ponemos el primer checkbox en "checked"
        self.checkBox_Re3.setCheckState(Qt.Checked)
        # Cargamos los numeros de Reynolds
        self.Reynolds = ['Re_3','Re_6','Re_9','Re_std']
        # Cargamos los parametros caracteristicos del perfil
        self.vHeader = ['CL_max','beta_max',
                    'CL_betamax','alpha_betamax','dbeta_dalpha','b_max',
                    'CD_min','CL_CD_max','cuspide','CM0']
        # Cargamos el formato en que se deben mostrar cada parametro caracteristico
        self.HeadersFormat = ['{:.2f}','{:.1f}','{:.2f}','{:.1f}','{:.3f}','{:.1f}',
                              '{:.5f}','{:.1f}',str,'{:.2}']
        self.tableDatos.setVerticalHeaderLabels(self.vHeader)
        self.tableDatos.setHorizontalHeaderLabels(self.Reynolds)


        self.connect(self.openButton,SIGNAL("clicked()"),self.Open)
        self.connect(self.okButton,SIGNAL("clicked()"),self,SLOT('accept()'))
        self.connect(self.cancelButton,SIGNAL("clicked()"),self,SLOT('reject()'))
        self.connect(self.checkBox_Re3,SIGNAL("toggled(bool)"),self.plotGraph)
        self.connect(self.checkBox_Re6,SIGNAL("toggled(bool)"),self.plotGraph)
        self.connect(self.checkBox_Re9,SIGNAL("toggled(bool)"),self.plotGraph)
        self.connect(self.checkBox_SR,SIGNAL("toggled(bool)"),self.plotGraph)

    def Open(self):
        dir = '.'
        fileAdress = QFileDialog.getOpenFileName(self,"Selecionar Perfil",dir,filter="Text File (*.txt)")
        self.AIRFOIL_DATA = profile.lectura_perfiles(fileAdress[0])
        #print self.AIRFOIL_DATA

        for j in range(len(self.Reynolds)):
            for i in range(len(self.vHeader)):
                if self.vHeader[i] in self.AIRFOIL_DATA[self.Reynolds[j]].keys():
                    item = QTableWidgetItem(self.HeadersFormat[i].format(self.AIRFOIL_DATA[self.Reynolds[j]][self.vHeader[i]]))
                    self.tableDatos.setItem(i,j,item)
        nombrePerfil = fileAdress[0].split('/')[-1].split('.')[0]
        #print nombrePerfil
        self.labelPerfil.setText('Datos Perfil:  NACA '+nombrePerfil)
        self.plotGraph()

    def plotGraph(self):
        self.axes1.clear()
        self.axes2.clear()
        #print 'plot'

        if self.checkBox_Re3.isChecked():
            #print 'Re3'
            x1 = []
            x2 = []
            y1 = []
            y2 = []
            for value in self.AIRFOIL_DATA['Re_3']['AoA_Cl']:
                x1.append(value[0])
                y1.append(value[1])

            for value in self.AIRFOIL_DATA['Re_3']['Cl_Cd']:
                x2.append(value[0])
                y2.append(value[1])

            self.axes1.plot(x1,y1,'b')
            self.axes2.plot(x2,y2,'b')

        if self.checkBox_Re6.isChecked():
           #print 'Re6'
            x1 = []
            x2 = []
            y1 = []
            y2 = []
            for value in self.AIRFOIL_DATA['Re_6']['AoA_Cl']:
                x1.append(value[0])
                y1.append(value[1])

            for value in self.AIRFOIL_DATA['Re_6']['Cl_Cd']:
                x2.append(value[0])
                y2.append(value[1])

            self.axes1.plot(x1,y1,'r')
            self.axes2.plot(x2,y2,'r')

        if self.checkBox_Re9.isChecked():
            #print 'Re9'
            x1 = []
            x2 = []
            y1 = []
            y2 = []
            for value in self.AIRFOIL_DATA['Re_9']['AoA_Cl']:
                x1.append(value[0])
                y1.append(value[1])

            for value in self.AIRFOIL_DATA['Re_9']['Cl_Cd']:
                x2.append(value[0])
                y2.append(value[1])

            self.axes1.plot(x1,y1,'g')
            self.axes2.plot(x2,y2,'g')

        if self.checkBox_SR.isChecked():
            #print 'Re_std'
            x1 = []
            x2 = []
            y1 = []
            y2 = []
            for value in self.AIRFOIL_DATA['Re_std']['AoA_Cl']:
                x1.append(value[0])
                y1.append(value[1])

            for value in self.AIRFOIL_DATA['Re_std']['Cl_Cd']:
                x2.append(value[0])
                y2.append(value[1])

            self.axes1.plot(x1,y1,'y')
            self.axes2.plot(x2,y2,'y')

        self.canvas1.draw()
        self.canvas2.draw()


app = QApplication(sys.argv)
form = MainDialog()
form.show()
app.exec_()