import numpy as np

import matplotlib.pyplot as plt
from IPython.display import display, clear_output

import package_DBR

def Lead_Lag_Discreet_RT(MV,PV,Tlead,Tlag,Ts,Kp=1,method='EBD',PVInit = 0):
    K = Ts/Tlag
    if len(PV) == 0:
           PV.append(PVInit)
    else: # MV[k+1] is MV[-1] and MV[k] is MV[-2]
        if method == 'EBD':
            PV.append((1/(1+K))*PV[-1] + (K*Kp/(1+K))*((1+(Tlead/Ts))*MV[-1]-(Tlead/Ts)*MV[-2]))
        elif method == 'EFD':
            PV.append((1-K)*PV[-1] + K*Kp*MV[-2])
        elif method == 'TRAP':
            PV.append((1/(2*T+Ts))*((2*T-Ts)*PV[-1] + Kp*Ts*(MV[-1] + MV[-2])))            
        else:
            PV.append((1/(1+K))*PV[-1] + (K*Kp/(1+K))*MV[-1])
            
def PID_RT(SP,PV,Man,MVMan,MVFF,Kc,Ti,Td,Ts,MVMin,MVMax,MV,MVP,MVI,MVD,E,alpha = 0.4,ManFF=False,PVInit = 0,method = "EBD-EBD"):
    Tfd = alpha*Td
    
    if Man:
        pass
    if len(MVI) == 0 or len(MVD)==0 or len(MVP)==0:
        if len(MVI) == 0:    
            MVI.append(PVInit)
        if len(MVD) == 0:
            MVD.append(PVInit)
        if len(MVP) == 0:
            MVP.append(PVInit)
        #print("ici")
    else:
        print("la")
        e  = (SP[-1]-PV[-1])
        E.append(e)
        temp_mvp = Kc*e
        MVP.append(temp_mvp)
        temp_mvi = MVI[-1]+Kc*(Ts/Ti)*e
        temp_mvd = (Tfd/(Tfd+Ts))*MVD[-1]+(Kc*Td/(Tfd+Ts))*(e-E[-2])
        MVD.append(temp_mvd)
        if temp_mvp+temp_mvi+temp_mvd>MVMax:  
            MVI.append(MVMax-temp_mvp-temp_mvd)
        elif temp_mvp+temp_mvi+temp_mvd<MVMin:
            MVI.append(MVMin-temp_mvp-temp_mvd)
        else:
            MVI.append(temp_mvi)
        MV.append(MVP[-1]+MVI[-1]+MVD[-1])