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
            
