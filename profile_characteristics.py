# -*- coding: iso-8859-1 -*-

def build_data(AIRFOIL_DATA,block,data_set,polar_data):
    x_data = float(polar_data.split()[0])
    y_data = float(polar_data.split()[1])
    if block ==1:
        if data_set ==1:
            AIRFOIL_DATA['Re_3']['AoA_Cl'].append((x_data,y_data))
        elif data_set ==2:
            AIRFOIL_DATA['Re_6']['AoA_Cl'].append((x_data,y_data))
        elif data_set ==3:
            AIRFOIL_DATA['Re_9']['AoA_Cl'].append((x_data,y_data))
        elif data_set ==4:
            AIRFOIL_DATA['Re_std']['AoA_Cl'].append((x_data,y_data))
        elif data_set ==7:
            AIRFOIL_DATA['Re_3']['AoA_Cm'].append((x_data,y_data))
        elif data_set ==8:
            AIRFOIL_DATA['Re_6']['AoA_Cm'].append((x_data,y_data))
        elif data_set ==9:
            AIRFOIL_DATA['Re_9']['AoA_Cm'].append((x_data,y_data))
        elif data_set ==10:
            AIRFOIL_DATA['Re_std']['AoA_Cm'].append((x_data,y_data))

    elif block == 2:
        if data_set ==1:
            AIRFOIL_DATA['Re_3']['Cl_Cd'].append((x_data,y_data))
        elif data_set ==2:
            AIRFOIL_DATA['Re_6']['Cl_Cd'].append((x_data,y_data))
        elif data_set ==3:
            AIRFOIL_DATA['Re_9']['Cl_Cd'].append((x_data,y_data))
        elif data_set ==4:
            AIRFOIL_DATA['Re_std']['Cl_Cd'].append((x_data,y_data))

def interpol_max(p1,p2,p3):
    x1=p1[0]
    x2=p2[0]
    x3=p3[0]
    y1=p1[1]
    y2=p2[1]
    y3=p3[1]
    det = x1**2*x2+x2**2*x3+x3**2*x1-x3**2*x2-x3*x1**2-x1*x2**2
    Xmax = (y1*(x2**2-x3**2)+y2*(x3**2-x1**2)+y3*(x1**2-x2**2))/\
            (2*(y1*(x2-x3)+y2*(x3-x1)+y3*(x1-x2)))
    Ymax = (y1*(x2**2*x3-x3**2*x2)+y2*(x3**2*x1-x1**2*x3)+y3*(x1**2*x2-x2**2*x1)
            +(Xmax/2.0)*(y1*(x3**2-x2**2)+y2*(x1**2-x3**2)+y3*(x2**2-x1**2)))/det
    return (Xmax,Ymax)
  

def lectura_perfiles(airfoil):
    aux_file = open(airfoil)
    polar_data = aux_file.readlines()
    aux_file.close()
    reynolds = ['Re_3','Re_6','Re_9','Re_std']
    AIRFOIL_DATA = dict()
    for reynold in reynolds:
        AIRFOIL_DATA[reynold]={'AoA_Cl':[],'AoA_Cm':[],'Cl_Cd':[]}

    i=0
    block = 1
    while block < 3:
        b = 0
        if polar_data[i].split()[0] == 'Data':
            i=i+1
            data_set = float(polar_data[i])
            i=i+1
            while b == 0:
                if polar_data[i].split()[-1] == 'exist.':
                    i = i+1
                elif polar_data[i].split()[0] == 'Data':
                    b = 1
                elif polar_data[i].split()[0] == 'Found':
                    block = block + 1
                    #print(block)
                    i=i+1
                    b=1
                else:
                    build_data(AIRFOIL_DATA,block,data_set,polar_data[i])
                    i = i+1
        else:
            i = i +1
    
    #Calculo de los parametros
    
    for reynold in reynolds:
        AIRFOIL_DATA[reynold]['keys']=[]
        c2=0
        c3=0
        c4=0
        c5=0
        if  (AIRFOIL_DATA[reynold]['AoA_Cl'] != [] and AIRFOIL_DATA[reynold]['Cl_Cd'] != []) :
            AIRFOIL_DATA[reynold]['keys']=['CL_max','beta_max',
                    'CL_betamax','alpha_betamax','dbeta_dalpha','b_max',
                    'CD_min','CL_CD_max','cuspide']
            CL_max = AIRFOIL_DATA[reynold]['AoA_Cl'][0][1]
            beta_max = AIRFOIL_DATA[reynold]['Cl_Cd'][0][0]\
                    /AIRFOIL_DATA[reynold]['Cl_Cd'][0][1]
            CL_betamax = AIRFOIL_DATA[reynold]['Cl_Cd'][0][0]
            CD_betamax = AIRFOIL_DATA[reynold]['Cl_Cd'][0][1]
            b_max = AIRFOIL_DATA[reynold]['Cl_Cd'][0][0]**(1.5)\
                    /AIRFOIL_DATA[reynold]['Cl_Cd'][0][1]
            CD_min = AIRFOIL_DATA[reynold]['Cl_Cd'][0][1]

            for i in range(len(AIRFOIL_DATA[reynold]['Cl_Cd'])):
                CL = AIRFOIL_DATA[reynold]['Cl_Cd'][i][0]
                CD = AIRFOIL_DATA[reynold]['Cl_Cd'][i][1]
                beta = CL/CD
                b = abs(CL)**(1.5)/CD
                if beta > beta_max:
                    beta_max = beta
                    CL_betamax = CL
                    CD_betamax = CD
                    c2 = i
                if b > b_max and CL > 0:
                    b_max = b
                    c3 = i
                if CD < CD_min:
                    CD_min = CD
                    c4 = i

            for i in range(len(AIRFOIL_DATA[reynold]['AoA_Cl'])):
                AoA = AIRFOIL_DATA[reynold]['AoA_Cl'][i][0]
                CL = AIRFOIL_DATA[reynold]['AoA_Cl'][i][1]
                if CL > CL_max:
                    CL_max = CL
                    c1 = i
                if CL > CL_betamax and i > 3:
                    c5 = i

            if c2!=0 and c2!=len(AIRFOIL_DATA[reynold]['Cl_Cd']):
                p1_beta = (AIRFOIL_DATA[reynold]['Cl_Cd'][c2+1][0],
                        AIRFOIL_DATA[reynold]['Cl_Cd'][c2+1][0]\
                        /AIRFOIL_DATA[reynold]['Cl_Cd'][c2+1][1])
                p2_beta = (CL_betamax,beta_max)
                p3_beta = (AIRFOIL_DATA[reynold]['Cl_Cd'][c2-1][0],
                        AIRFOIL_DATA[reynold]['Cl_Cd'][c2-1][0]\
                        /AIRFOIL_DATA[reynold]['Cl_Cd'][c2-1][1])
                (CL_betamax,beta_max)=interpol_max(p1_beta,p2_beta,p3_beta)
                #print(beta_max)
                #print(CL_betamax)

            alpha_betamax = (AIRFOIL_DATA[reynold]['AoA_Cl'][c5+1][0]
                   + (AIRFOIL_DATA[reynold]['AoA_Cl'][c5][0]
                   - AIRFOIL_DATA[reynold]['AoA_Cl'][c5+1][0])
                   /(AIRFOIL_DATA[reynold]['AoA_Cl'][c5][1]
                   - AIRFOIL_DATA[reynold]['AoA_Cl'][c5+1][1])
                   *(CL_betamax-AIRFOIL_DATA[reynold]['AoA_Cl'][c5+1][1]))
            #print(alpha_betamax)

            if c3!=0 and c3!=len(AIRFOIL_DATA[reynold]['Cl_Cd']):
                p1_b = (AIRFOIL_DATA[reynold]['Cl_Cd'][c3+1][0],
                        abs(AIRFOIL_DATA[reynold]['Cl_Cd'][c3+1][0])**1.5\
                        /AIRFOIL_DATA[reynold]['Cl_Cd'][c3+1][1])
                p2_b = (AIRFOIL_DATA[reynold]['Cl_Cd'][c3][0],b_max)
                p3_b = (AIRFOIL_DATA[reynold]['Cl_Cd'][c3-1][0],
                        AIRFOIL_DATA[reynold]['Cl_Cd'][c3-1][1])
                b_max = interpol_max(p1_b,p2_b,p3_b)[1]
                #print(b_max)

            if c4!=0 and c4!=len(AIRFOIL_DATA[reynold]['Cl_Cd']):
                p1_CD = (AIRFOIL_DATA[reynold]['Cl_Cd'][c4+1][0],
                        AIRFOIL_DATA[reynold]['Cl_Cd'][c4+1][1])
                p2_CD = (AIRFOIL_DATA[reynold]['Cl_Cd'][c4][0],CD_min)
                p3_CD = (AIRFOIL_DATA[reynold]['Cl_Cd'][c4-1][0],
                         AIRFOIL_DATA[reynold]['Cl_Cd'][c4-1][1])
                CD_min =interpol_max(p1_CD,p2_CD,p3_CD)[1]
                #print(CD_min)

            if c5!=0 and c5!=len(AIRFOIL_DATA[reynold]['AoA_Cl']):
                p1_CL = (AIRFOIL_DATA[reynold]['AoA_Cl'][c5+1][0],
                        AIRFOIL_DATA[reynold]['AoA_Cl'][c5+1][1])
                p2_CL = (AIRFOIL_DATA[reynold]['AoA_Cl'][c5][0],CL_max)
                p3_CL = (AIRFOIL_DATA[reynold]['AoA_Cl'][c5-1][0],
                         AIRFOIL_DATA[reynold]['AoA_Cl'][c5-1][1])
                CL_max =interpol_max(p1_CL,p2_CL,p3_CL)[1]
                #print(CL_max)

            AIRFOIL_DATA[reynold]['CL_max'] = CL_max
            AIRFOIL_DATA[reynold]['beta_max'] = beta_max
            AIRFOIL_DATA[reynold]['CL_betamax'] = CL_betamax
            AIRFOIL_DATA[reynold]['alpha_betamax'] = alpha_betamax
            AIRFOIL_DATA[reynold]['b_max'] = b_max
            AIRFOIL_DATA[reynold]['CD_min'] = CD_min

            # valoracion de la cuspide (version simple)
            # La valoracion de la cuspide se calcula como la relacion delta_CL/delta_alpha
            # en la region alrededor de CLmax. A mayor delta_CL/delta_alpha, peor la valoracion
            # de la cuspide, ya que habra un cambio mas abrupto del CL a similar variación de
            # alpha. Los limites de las tres categorias "bueno", "medio" y "malo" son arbitrarios.
            cuspide = (CL_max-min(AIRFOIL_DATA[reynold]['AoA_Cl'][c1-1][1],
                  AIRFOIL_DATA[reynold]['AoA_Cl'][c1+1][1]))\
                  /(AIRFOIL_DATA[reynold]['AoA_Cl'][c1-1][0]
                  - AIRFOIL_DATA[reynold]['AoA_Cl'][c1+1][0])
            if cuspide < 0.04:
                cuspide_name = 'bueno'
            elif cuspide > 0.08:
                cuspide_name = 'malo'
            else : cuspide_name = 'medio'
            AIRFOIL_DATA[reynold]['cuspide'] = cuspide_name

            # d\beta/d\alpha (version simple)
            dCD_dCL = (AIRFOIL_DATA[reynold]['Cl_Cd'][c2-1][1]
                   - AIRFOIL_DATA[reynold]['Cl_Cd'][c2+1][1])\
                   / (AIRFOIL_DATA[reynold]['Cl_Cd'][c2-1][0]
                   - AIRFOIL_DATA[reynold]['Cl_Cd'][c2+1][0])

            dCL_dalpha = (AIRFOIL_DATA[reynold]['AoA_Cl'][c5][1]
                -AIRFOIL_DATA[reynold]['AoA_Cl'][c5+1][1])\
                /(AIRFOIL_DATA[reynold]['AoA_Cl'][c5][0]
                -AIRFOIL_DATA[reynold]['AoA_Cl'][c5+1][0])

            dbeta_dalpha = (1-beta_max*dCD_dCL)*dCL_dalpha/CD_betamax
            AIRFOIL_DATA[reynold]['dbeta_dalpha'] = dbeta_dalpha

            CL_CD_max = CL_max/CD_min
            AIRFOIL_DATA[reynold]['CL_CD_max'] = CL_CD_max

        # calculo del Cmo promedio en la zona lineal
        if AIRFOIL_DATA[reynold]['AoA_Cm'] != []:
            AIRFOIL_DATA[reynold]['keys'].append('CM0')
            sum_CM0 = 0
            c6 = 0
            for i in range(len(AIRFOIL_DATA[reynold]['AoA_Cm'])):
                alpha = AIRFOIL_DATA[reynold]['AoA_Cm'][i][0]
                CM0 = AIRFOIL_DATA[reynold]['AoA_Cm'][i][1]
                if alpha > -5 or alpha < 10:
                    sum_CM0 = sum_CM0 + CM0
                    c6 = c6+1

            CM0 = sum_CM0/c6
            AIRFOIL_DATA[reynold]['CM0'] = CM0

    return AIRFOIL_DATA











