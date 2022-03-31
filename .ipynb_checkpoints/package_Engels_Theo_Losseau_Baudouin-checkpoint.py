import numpy as np

import matplotlib.pyplot as plt
from IPython.display import display, clear_output
from package_DBR import Delay_RT,FO_RT, Bode, Process, SelectPath_RT

def Lead_Lag_Discreet_RT(MV,PV,Tlead,Tlag,Ts,Kp=1,method='EBD',PVInit = 0):
    
    """
    The function "Lead_Lag_Discreet_RT" needs to be included in a "for or while loop".
    
    :MV: input vector
    :Kp: process gain, (optional: default value is 1)
    :Tlead: lead time constant [s]
    :Tlag: lag time constant [s]
    :Ts: sampling period [s]
    :PV: output vector
    :PVInit: (optional: default value is 0)
    :method: discretisation method (optional: default value is 'EBD')
        EBD: Euler Backward difference
        EFD: Euler Forward difference
        TRAP: Trapezoïdal method
    
    The function "Lead_Lag_Discreet_RT" appends a value to the output vector "PV".
    The appended value is obtained from a recurrent equation that depends on the discretisation method.
    """
    
    
    K = Ts/Tlag
    if len(PV) == 0:
           PV.append(PVInit)
    else: # MV[k+1] is MV[-1] and MV[k] is MV[-2]
        if method == 'EBD':
            PV.append((1/(1+K))*PV[-1] + (K*Kp/(1+K))*((1+(Tlead/Ts))*MV[-1]-(Tlead/Ts)*MV[-2]))
        elif method == 'EFD':
            PV.append((1-K)*PV[-1]+Kp*K*((Tlead/Ts)*MV[-1]+(1-(Tlead/Ts))*MV[-2]))
        elif method == 'TRAP':
            a = 1+((2*Tlead)/Ts)
            b = 1-2*Tlead/Ts
            c = 2+K
            d = -2+K
            PV.append(((Kp*K)/c)*(a*MV[-1]+b*MV[-2])-(d/c)*PV[-1])            
    return
            
def Run_LL(Tlead,Tlag,TSim,method):
    
    """
    The function "Lead_Lag_Discreet_RT" is designed to be used in an interactive widget.
    
    :Tlead: lead time constant [s]
    :Tlag: lag time constant [s]
    :TSim: Total time of simulation, set small value to zoom and see the difference between methods
    :method: discretisation method (optional: default value is 'EBD')
        EBD: Euler Backward difference
        EFD: Euler Forward difference
        TRAP: Trapezoïdal method
    
    Uses Lead_Lag_Discreet_RT to display MV and PV with default values for Ts and MV.
    
    
    :MV: input vector used in Lead_Lag_Discreet_RT = {0: 0, 2: -1,12:2,22:3}
    :Ts: sampling period [s] = 0.1
    """
    
    
    Ts = 0.1
    N = int(TSim/Ts) + 1
    MVPath = {0: 0, 2: -1,12:2,22:3}
    t = []
    MV = []
    PV = []


    for i in range(0,N):
        t.append(i*Ts)
        SelectPath_RT(MVPath,t,MV)
        Lead_Lag_Discreet_RT(MV,PV,Tlead,Tlag,Ts, method = method)
            
    plt.figure(figsize = (15,9))

    plt.step(t,MV,'b-',label='MV',where='post')
    plt.step(t,PV,'skyblue',label='PV',where='post')
    plt.ylabel('Value of MV')
    plt.legend(loc='best')
    plt.title('Path response')
    plt.xlim([0, TSim])
    return 
            
def PID_RT(SP,PV,Man,MVMan,MVFF,Kc,Ti,Td,Ts,MVMin,MVMax,MV,MVP,MVI,MVD,E,alpha = 0.4,ManFF=False,PVInit = 0,method = "EBD-EBD"):
    
    """
    The function "PID_RT" needs to be included in a "for or while loop".
    
    :MV: input vector
    :MV: input vector
    :Kp: process gain
    :Tlead: lead time constant [s]
    :Tlag: lag time constant [s]
    :Ts: sampling period [s]
    :PV: output vector
    :PVInit: (optional: default value is 0)
    :method: discretisation method (optional: default value is 'EBD')
        EBD: Euler Backward difference
        EFD: Euler Forward difference
        TRAP: Trapezoïdal method
    
    The function "PID_RT" appends a value to the output vector "PV".
    The appended value is obtained from a recurrent equation that depends on the discretisation method.
    """
    
    
    Tfd = alpha*Td
    
    if len(PV) == 0:
        PV.append(PVInit)
    e  = (SP[-1]-PV[-1])
    E.append(e)
    
    if len(MVI) == 0:
        MVI.append(Kc*(Ts/Ti)*e)
    else:
        MVI.append(MVI[-1]+Kc*(Ts/Ti)*e)
    if len(MVD) == 0:
        MVD.append(0)
    else:
        MVD.append((Tfd/(Tfd+Ts))*MVD[-1]+(Kc*Td/(Tfd+Ts))*(e-E[-2]))
    MVP.append(Kc*e)

    if MVP[-1]+MVI[-1]+MVD[-1]>MVMax:  
        MVI[-1] = MVMax-MVP[-1]-MVD[-1]
    elif MVP[-1]+MVI[-1]+MVD[-1]<MVMin:
        MVI[-1] = MVMin-MVP[-1]-MVD[-1]
    if Man[-1]:
        MVI[-1] = MVMan[-1]-MVP[-1]
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

    
def sim_tclabP(MV,PV,Ts,PVtemp1,PVtemp2,Kp,T1,T2,Theta):
    
    FO_RT(MV,Kp,T1,Ts,PVtemp1)
    FO_RT(PVtemp1,1,T2,Ts,PVtemp2)
    Delay_RT(PVtemp2,Theta,Ts,PV)
    
def sim_tclabD(MV,PV,Ts,DV0,DVtemp1,DVtemp2,Kd,T1,T2,Theta):

    FO_RT(MV,Kd,T1,Ts,DVtemp1)
    FO_RT(DVtemp1,1,T2,Ts,DVtemp2)
    Delay_RT(DVtemp2,Theta,Ts,PV)
    
    
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
        ax_gain.semilogx(omega,20*np.log10(np.abs(Ps)),label='C(s)') 
        gain_min = np.min(20*np.log10(np.abs(Ps)/5))
        gain_max = np.max(20*np.log10(np.abs(Ps)*5))
        ax_gain.set_xlim([np.min(omega), np.max(omega)])
        ax_gain.set_ylim([gain_min, gain_max])
        ax_gain.set_ylabel('Amplitude |C| [db]')
        ax_gain.set_title('Bode plot of C')
        ax_gain.legend(loc='best')
    
        # Phase part
        ax_phase.semilogx(omega, (180/np.pi)*np.unwrap(np.angle(Ps)),label='C(s)')  
        ax_phase.set_xlim([np.min(omega), np.max(omega)])
        ph_min = np.min((180/np.pi)*np.unwrap(np.angle(Ps))) - 10
        ph_max = np.max((180/np.pi)*np.unwrap(np.angle(Ps))) + 10
        ax_phase.set_ylim([np.max([ph_min, -200]), ph_max])
        ax_phase.set_ylabel(r'Phase $\angle C$ [°]')
        ax_phase.legend(loc='best')
    else:
        return Ps
    
    
def bodePC(P,C,omega,Show = True):
    Ps = np.multiply(Bode(P,omega,Show = False),bodePID(C,omega,Show = False))
    if Show == True:
    
        fig, (ax_gain, ax_phase) = plt.subplots(2,1)
        fig.set_figheight(12)
        fig.set_figwidth(22)

        # Gain part
        ax_gain.semilogx(omega,20*np.log10(np.abs(Ps)),label='PC(s)') 
        gain_min = np.min(20*np.log10(np.abs(Ps)/5))
        gain_max = np.max(20*np.log10(np.abs(Ps)*5))
        ax_gain.set_xlim([np.min(omega), np.max(omega)])
        ax_gain.set_ylim([gain_min, gain_max])
        ax_gain.set_ylabel('Amplitude |PC| [db]')
        ax_gain.set_title('Bode plot of PC')
        ax_gain.legend(loc='best')
    
        # Phase part
        ax_phase.semilogx(omega, (180/np.pi)*np.unwrap(np.angle(Ps)),label='PC(s)')  
        ax_phase.set_xlim([np.min(omega), np.max(omega)])
        ph_min = np.min((180/np.pi)*np.unwrap(np.angle(Ps))) - 10
        ph_max = np.max((180/np.pi)*np.unwrap(np.angle(Ps))) + 10
        ax_phase.set_ylim([np.max([ph_min, -200]), ph_max])
        ax_phase.set_ylabel(r'Phase $\angle PC$ [°]')
        ax_phase.legend(loc='best')
    else:
        return Ps

    
    
def find_nearest_index(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx    
    
    
def Margin(P,C,omega,Show = True):
    Ps = bodePC(P,C,omega,Show=False)
    
    
    PsGain = 20*np.log10(np.abs(Ps))
    PsPhase = (180/np.pi)*np.unwrap(np.angle(Ps))
    
    F_PhaseMargin = omega[find_nearest_index(PsGain,0)]
    F_GainMargin = omega[find_nearest_index(PsPhase,-180)]
    
    
    GainMargin = np.abs(PsGain[find_nearest_index(PsPhase,-180)])
    PhaseMargin = 180+PsPhase[find_nearest_index(PsGain,0)]
    if Show == True:
    
        fig, (ax_gain, ax_phase) = plt.subplots(2,1)
        fig.set_figheight(12)
        fig.set_figwidth(22)

        # Gain part
        ax_gain.semilogx(omega,PsGain,label='PC(s)')
        ax_gain.plot(omega,np.zeros_like(omega),'black')
        ax_gain.vlines(F_GainMargin,PsGain[find_nearest_index(PsPhase,-180)],0)
        gain_min = np.min(20*np.log10(np.abs(Ps)/5))
        gain_max = np.max(20*np.log10(np.abs(Ps)*5))
        ax_gain.set_xlim([np.min(omega), np.max(omega)])
        ax_gain.set_ylim([gain_min, gain_max])
        ax_gain.set_ylabel('Amplitude |PC| [db]')
        ax_gain.set_title('Bode plot of PC')
        ax_gain.legend(loc='best')
    
        # Phase part
        ax_phase.semilogx(omega, PsPhase,label='PC(s)')  
        ax_phase.set_xlim([np.min(omega), np.max(omega)])
        ax_phase.vlines(F_PhaseMargin,-180,PsPhase[find_nearest_index(PsGain,0)])
        ax_phase.plot(omega,-180*np.ones_like(omega),'black')
        ph_min = np.min((180/np.pi)*np.unwrap(np.angle(Ps))) - 10
        ph_max = np.max((180/np.pi)*np.unwrap(np.angle(Ps))) + 10
        ax_phase.set_ylim([np.max([ph_min, -200]), ph_max])
        ax_phase.set_ylabel(r'Phase $\angle PC$ [°]')
        ax_phase.legend(loc='best')
        
        
        
        print(GainMargin)
        print(PhaseMargin)
    else:
        return GainMargin,PhaseMargin
    

    
    
class PID:
    
    def __init__(self, parameters):
        
        self.parameters = parameters
        self.parameters['Kc'] = parameters['Kc'] if 'Kc' in parameters else 1.0
        self.parameters['Ti'] = parameters['Ti'] if 'Ti' in parameters else 1.0
        self.parameters['Td'] = parameters['Td'] if 'Td'in parameters else 10
        self.parameters['alpha'] = parameters['alpha'] if 'alpha' in parameters else 0.4