import numpy as np

import matplotlib.pyplot as plt
from IPython.display import display, clear_output
from package_DBR import Delay_RT,FO_RT, Bode, Process

def Lead_Lag_Discreet_RT(MV,PV,Tlead,Tlag,Ts,Kp=1,method='EBD',PVInit = 0):
    K = Ts/Tlag
    if len(PV) == 0:
           PV.append(PVInit)
    else: # MV[k+1] is MV[-1] and MV[k] is MV[-2]
        if method == 'EBD':
            PV.append((1/(1+K))*PV[-1] + (K*Kp/(1+K))*((1+(Tlead/Ts))*MV[-1]-(Tlead/Ts)*MV[-2]))
        # elif method == 'EFD':
        #     PV.append((1-K)*PV[-1] + K*Kp*MV[-2])
        # elif method == 'TRAP':
        #     PV.append((1/(2*T+Ts))*((2*T-Ts)*PV[-1] + Kp*Ts*(MV[-1] + MV[-2])))            
        # else:
        #     PV.append((1/(1+K))*PV[-1] + (K*Kp/(1+K))*MV[-1])
            
def PID_RT(SP,PV,Man,MVMan,MVFF,Kc,Ti,Td,Ts,MVMin,MVMax,MV,MVP,MVI,MVD,E,alpha = 0.4,ManFF=False,PVInit = 0,method = "EBD-EBD"):
    Tfd = alpha*Td
    
    if len(MVI) == 0 or len(MVD)==0 or len(MVP)==0:
        if len(MVI) == 0:    
            MVI.append(PVInit)  # metter a 0 seulememnt si ca marche pas
        if len(MVD) == 0:
            MVD.append(PVInit)
        if len(MVP) == 0:
            MVP.append(PVInit)
        #print("ici")
    else:
        #print("la")
        e  = (SP[-1]-PV[-1])
        E.append(e)
        temp_mvp = Kc*e
        MVP.append(temp_mvp)
        temp_mvi = MVI[-1]+Kc*(Ts/Ti)*e
        temp_mvd = (Tfd/(Tfd+Ts))*MVD[-1]+(Kc*Td/(Tfd+Ts))*(e-E[-2])
        #print(MVD[-1])
        #print(e)
        #print(E[-1])
        MVD.append(temp_mvd)
        if temp_mvp+temp_mvi+temp_mvd>MVMax:  
            MVI.append(MVMax-temp_mvp-temp_mvd)
        elif temp_mvp+temp_mvi+temp_mvd<MVMin:
            MVI.append(MVMin-temp_mvp-temp_mvd)
        else:
            MVI.append(temp_mvi)
        if Man[-1] and ManFF:
            MV.append(MVMan[-1]-MVFF[-1])
        elif Man[-1] and not(ManFF):
            MV.append(MVMan[-1])
        elif ManFF:
            MV.append(MVP[-1]+MVI[-1]+MVD[-1]-MVFF[-1])
        else:
             MV.append(MVP[-1]+MVI[-1]+MVD[-1])
        
        
        
def FF_RT(DV,DV0,Tlead1,Tlag1,Tlead2,Tlag2,Theta1,Theta2,Kp,Kd,Ts,MVFF,PV1,PV2):
    #New_DV = [x-DV0 for x in DV]
    Lead_Lag_Discreet_RT(DV,PV1,Tlead1,Tlag1,Ts,Kp=(Kd/Kp))
    Lead_Lag_Discreet_RT(PV1,PV2,Tlead2,Tlag2,Ts)
    Delay_RT(PV2,max(0,Theta2-Theta1),Ts,MVFF)

    
def sim_tclabP(MV,PV,Ts,PVtemp1,PVtemp2):
    
    FO_RT(MV,0.395,47.83,Ts,PVtemp1)
    FO_RT(PVtemp1,1,17.39,Ts,PVtemp2)
    Delay_RT(PVtemp2,9.3138,Ts,PV)
    
def sim_tclabD(MV,PV,Ts,DV0,DVtemp1,DVtemp2):

    FO_RT(MV,0.6347,245.357,Ts,DVtemp1)
    FO_RT(DVtemp1,1,3.133,Ts,DVtemp2)
    Delay_RT(DVtemp2,0.559,Ts,PV)
    
    
def IMC_Tuning_SOPDT(Kp,T1,T2,theta,gamma):
    Tclp=gamma * T1
    Kc=Kp*((T1+T2)/(Tclp+theta))
    Ti=(T1+T2)
    Td= T1*T2/T1+T2
    return (Kc,Ti,Td)


def margins(P,C):
    omega = np.logspace(-2, 2, 10000)
    PC = Bode(P,omega,Show = False)*Bode(C,omega,Show = False)
    Bode(PC,omega)


def bodePID(C,omega,Show = True):
    s = 1j*omega
    
    PGain = C.parameters['Kc']*np.ones_like(s)
    PI = 1/(C.parameters['Ti']*s)
    PP = np.ones_like(s)
    PD = (C.parameters['Td']*s)/((C.parameters['alpha']*C.parameters['Td']*s)+1)
    
    Ps = np.multiply(PGain,(PI+PP+PD))
    if Show == True:
    
        fig, (ax_gain, ax_phase) = plt.subplots(2,1)
        fig.set_figheight(12)
        fig.set_figwidth(22)

        # Gain part
        ax_gain.semilogx(omega,20*np.log10(np.abs(Ps)),label='P(s)') 
        gain_min = np.min(20*np.log10(np.abs(Ps)/5))
        gain_max = np.max(20*np.log10(np.abs(Ps)*5))
        ax_gain.set_xlim([np.min(omega), np.max(omega)])
        ax_gain.set_ylim([gain_min, gain_max])
        ax_gain.set_ylabel('Amplitude |P| [db]')
        ax_gain.set_title('Bode plot of P')
        ax_gain.legend(loc='best')
    
        # Phase part
        ax_phase.semilogx(omega, (180/np.pi)*np.unwrap(np.angle(Ps)),label='P(s)')  
        ax_phase.set_xlim([np.min(omega), np.max(omega)])
        ph_min = np.min((180/np.pi)*np.unwrap(np.angle(Ps))) - 10
        ph_max = np.max((180/np.pi)*np.unwrap(np.angle(Ps))) + 10
        ax_phase.set_ylim([np.max([ph_min, -200]), ph_max])
        ax_phase.set_ylabel(r'Phase $\angle P$ [°]')
        ax_phase.legend(loc='best')
    else:
        return Ps
    
class PID:
    
    def __init__(self, parameters):
        
        self.parameters = parameters
        self.parameters['Kc'] = parameters['Kc'] if 'Kc' in parameters else 1.0
        self.parameters['Ti'] = parameters['Ti'] if 'Ti' in parameters else 1.0
        self.parameters['Td'] = parameters['Td'] if 'Td'in parameters else 10
        self.parameters['alpha'] = parameters['alpha'] if 'alpha' in parameters else 0.4