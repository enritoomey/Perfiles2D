# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 00:18:27 2015

@author: Enriquito
"""

import os
File = open(os.getcwd()+'\\NACA.LST')
data = File.readlines()
File.close()

listperfiles = []
for line in data[7:]:
    listperfiles.append((line.split(" ")[1],line.split("  ")[-1][:-1]))

for perfil in listperfiles:
    aux = open(os.getcwd()+'\\'+perfil[0],'r')
    datos = aux.read()
    aux.close()
    aux2 = open(os.getcwd()+'\\'+perfil[1]+'.txt','w')
    aux2.write(datos)
    aux2.close()


