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
            
def PID_RT(SP,PV,Man,MVMan,MVFF,Kc,Ti,Td,Ts,MVMin,MVMax,MV,MVP,MVI,MVD,E,alpha = 0.4,ManFF=False,FF = True,PVInit = 0):
    
    """
    The function "PID_RT" needs to be included in a "for or while loop".
    
    Inputs:
    :SP: SetPoint vector
    :PV: Precess Value vector
    :Man: Manual vector
    :MVMan: Manual value vector
    :MVFF: FeedForward vector (ouput of FF_RT)
    
    Parameters:
    :Kc: Controller gain
    :Ti: integral time constant [s]
    :Td: derivative time constant [s]
    :Ts: sampling period [s]
    :MVMin: Minimum value for MV (used for saturation and anti wind-up)
    :MVMax: Maxomum value for MV (used for saturation and anti wind-up)
    :alpha: Tfd = alpha*Td where Tfd is the derivative filter time constant [s]
    :FF: boolean value, activates or deactivates the feedforward, default value is true
    :ManFF: boolean value, activates the feedforward in manual mode, default value is false. FF needs to be True for ManFF to have an impact.
        
    Outputs:
    :MV: Manipulated value vector
    :MVP: Proportional part of manipulated value vector
    :MVI: Integral part of manipulated value vector
    :MVD: Derivative part of manipulated value vector
    :E: Control error vector
    
    
    The function "PID_RT" appends values to MV, MVP, MVD, MVI and E.
    The values are based on the PID algorithm, the controller mode and feedforward
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
    if Man[-1] and FF and ManFF:
        if MVMan[-1]-MVFF[-1]>MVMin and MVMan[-1]-MVFF[-1]<MVMax:
            MV.append(MVMan[-1]-MVFF[-1])
        elif MVMan[-1]-MVFF[-1]<MVMin:
            MV.append(MVMin)
        else:
            MV.append(MVMax)
    elif (Man[-1] and not(FF)) or (Man[-1] and not(ManFF)):
        MV.append(MVMan[-1])
    elif FF:
        if MVP[-1]+MVI[-1]+MVD[-1]-MVFF[-1]>=MVMin and MVP[-1]+MVI[-1]+MVD[-1]-MVFF[-1]<=MVMax:
            MV.append(MVP[-1]+MVI[-1]+MVD[-1]-MVFF[-1])
        elif MVP[-1]+MVI[-1]+MVD[-1]-MVFF[-1]<MVMin:
            MV.append(MVMin)
        elif MVP[-1]+MVI[-1]+MVD[-1]-MVFF[-1]>MVMax:
            MV.append(MVMax)
    else:
        MV.append(MVP[-1]+MVI[-1]+MVD[-1])        
        
        
def FF_RT(DV,Tlead1,Tlag1,Tlead2,Tlag2,Theta1,Theta2,Kp,Kd,Ts,MVFF,PV1,PV2):
    
    """
    The function "FF_RT" needs to be included in a "for or while loop".
    
    :DV: input vector
    :Tlead1: First lead time constant [s]
    :Tlag1: First lag time constant [s]
    :Tlead2: Second lead time constant [s]
    :Tlag2: Second lag time constant [s]
    :Theta1: First Process delay
    :Theta2: Seconde Process delay
    :Kp: First process gain
    :Kd: Second process gain
    
    :Ts: sampling period [s]
    :MVFF: output vector
    :PV1: Transision vector needs to be empty at start of loop
    :PV2: Transision vector needs to be empty at start of loop
    
    The function "FF_RT" appends a value to the output vector "MVFF".
    The appended value is obtained through 2 lead_lag (Lead_Lag_Discreet_RT) and a delay (Delay_RT).
    """
    Lead_Lag_Discreet_RT(DV,PV1,Tlead1,Tlag1,Ts,Kp=(Kd/Kp))
    Lead_Lag_Discreet_RT(PV1,PV2,Tlead2,Tlag2,Ts)
    Delay_RT(PV2,max(0,Theta2-Theta1),Ts,MVFF)

    
def sim_tclabP(MV,PV,Ts,PVtemp1,PVtemp2,Kp,T1,T2,Theta):
    """
    The function "sim_tclabP" needs to be included in a "for or while loop".
    
    :MV: input vector
    :Ts: sampling period [s]
    :Kp: Process gain
    :T1: First Time constant [s]
    :T2: Second Time constant [s]
    :Theta: Delay time [s]
    
    
    :PV: output vector
    :PVtemp1: Transision vector needs to be empty at start of loop
    :PVtemp2: Transision vector needs to be empty at start of loop
    
    The function "sim_tclabP" appends a value to the output vector "PV".
    The appended value is obtained through 2 first order systems (FO_RT) and a delay (Delay_RT).
    """
    
    FO_RT(MV,Kp,T1,Ts,PVtemp1)
    FO_RT(PVtemp1,1,T2,Ts,PVtemp2)
    Delay_RT(PVtemp2,Theta,Ts,PV)
    
# def sim_tclabD(MV,PV,Ts,DVtemp1,DVtemp2,Kd,T1,T2,Theta):

#     FO_RT(MV,Kd,T1,Ts,DVtemp1)
#     FO_RT(DVtemp1,1,T2,Ts,DVtemp2)
#     Delay_RT(DVtemp2,Theta,Ts,PV)
    
    
def IMC_Tuning_SOPDT(Kp,T1,T2,theta,gamma):
    """
    The function "sim_tclabP" does NOT need to be included in a "for or while loop".
    
    :Kp: Process gain
    :T1: First Time constant [s]
    :T2: Second Time constant [s]
    :Theta: Delay time [s]
    :gamma: Tclp = gamma*T1 Tclp is the desired closed-loop time constant
    
    The function "IMC_Tuning_SOPDT" return a tuple containing values for Kc, Ti and Td obtained from the formuilas given in row I of the IMC table with T3 = 0.
    
    :Kc: Controller gain
    :Ti: Integration time constant
    :Td: Derivative time constant
    """
    Tclp=gamma * T1
    Kc=Kp*((T1+T2)/(Tclp+theta))
    Ti=(T1+T2)
    Td= T1*T2/T1+T2
    return (Kc,Ti,Td)


# def margins(P,C):
#     omega = np.logspace(-2, 2, 10000)
#     PC = Bode(P,omega,Show = False)*Bode(C,omega,Show = False)
#     Bode(PC,omega)


def bodePID(C,omega,Show = True):
    """
    :C: Controller as defined by the class "PID".     
        
    :omega: frequency vector (rad/s); generated by a command of the type "omega = np.logspace(-2, 2, 10000)".
    :Show: boolean value (optional: default value = True). If Show = True, the Bode diagram is shown. Otherwise Cs (C(j omega)) (vector of complex numbers) is returned.
    
    The function "bodePID" generates the Bode diagram of the controller C.
    
    This function is based on the "Bode" function found in package_DBR.py.
    """ 
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
    """
    :P: Process as defined by the class "Process"
    :C: Controller as defined by the class "PID".     
        
    :omega: frequency vector (rad/s); generated by a command of the type "omega = np.logspace(-2, 2, 10000)".
    :Show: boolean value (optional: default value = True). If Show = True, the Bode diagram is shown. Otherwise Ls (L(j omega)) (vector of complex numbers) is returned where L(s) = P(s)*C(s).
    
    The function "bodePID" generates the Bode diagram of the controller C.
    
    This function is based on the "Bode" function found in package_DBR.py.
    """ 
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
    """
    :array: array in wich we wish to find the value
    :value: value to find in array
    
    Find the index of the value in array closest to value
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx    
    
    
def Margin(P,C,omega,Show = True):
    """
    :P: Process as defined by the class "Process"
    :C: Controller as defined by the class "PID".     
        
    :omega: frequency vector (rad/s); generated by a command of the type "omega = np.logspace(-2, 2, 10000)".
    :Show: boolean value (optional: default value = True). If Show = True, the Bode diagram is shown. Otherwise the values of the margins are returned
    
    The function "Margin" generates the Bode diagram of L(s)=P(s)*C(s) and shows the phase and gain margins.
    
    This function is based on the "Bode" function found in package_DBR.py.
    """ 

    
    Ps = bodePC(P,C,omega,Show=False)
    
    
    PsGain = 20*np.log10(np.abs(Ps))
    PsPhase = (180/np.pi)*np.unwrap(np.angle(Ps))
    
    F_PhaseMargin = omega[find_nearest_index(PsGain,0)]
    F_GainMargin = omega[find_nearest_index(PsPhase,-180)]
    
    
    GainMargin = np.abs(PsGain[find_nearest_index(PsPhase,-180)])
    GainMargin = 10**(GainMargin/20)
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
    else:
        return GainMargin,PhaseMargin
    

    
    
def Run_PID_Interactive(gamma,alpha,Kp,Kd,Tlead1,Tlag1,Tlead2,Tlag2,Theta1,Theta2,DV0,MV0,PV0, ManFF, FF):
     """
    The Run_PID_Interactive function needs to be used in a interactive function from ipywidgets.
    
    Parameters:
    :Kp: Process gain
    :Kd: Disturbance process gain
    :Tlead1: First lead time constant [s]
    :Tlag1: First lag time constant [s]
    :Tlead2: Second lead time constant [s]
    :Tlag2: Second lag time constant [s]
    :Theta1: Process delay
    :Theta2: Disturbance process delay
    :DV0: working point of DV
    :MV0: working point of MV
    :PV0: Value of PV after sabilisation at MV = MV0 and DV = DV0
     
    The function "Run_PID_Interactive" displays graphs of a PID controller calculated using the 'PID_RT' function.
    """
    
    SPPath = {0: 70,2000:80,2600:65}
    ManPath = {0:1,1000:0}
    MVManPath={0:50}
    DVPath = {0:50,500:70,1500:60,2000:50}
    
    satMin = 0
    satMax = 100
    PVInit = 50
    
    
    TSim = 3000
    Ts = 1
    N = int(TSim/Ts) + 1
    t = []
    MV = []
    MVMan = []
    Man = []
    PV = []
    SP = []
    MVP = []
    MVI = []
    MVD = []
    E = []
    Et = []
    MVFF = []
    ODV = []
    OPV = []
    DV = []
    PVtemp1 = []
    PVtemp2 = []
    DVtemp1 = []
    DVtemp2 = []

    FFtemp1 = []
    FFtemp2 = []
    Kc,Ti,Td = IMC_Tuning_SOPDT(Kp,Tlead1,Tlead2,Theta1,gamma)
    for i in range(0,N):
        t.append(i*Ts)
        SelectPath_RT(SPPath,t,SP)
        SelectPath_RT(DVPath,t,DV)
        SelectPath_RT(ManPath,t,Man)
        SelectPath_RT(MVManPath,t,MVMan)
        DV[-1]-=DV0


        FF_RT(DV,Tlead1,Tlag1,Tlead2,Tlag2,Theta1,Theta2,Kp,Kd,Ts,MVFF,FFtemp1,FFtemp2)
        PID_RT(SP,PV,Man,MVMan,MVFF,Kc,Ti,Td,Ts,satMin,satMax,MV,MVP,MVI,MVD,E,alpha=alpha,FF = FF,ManFF = ManFF)
        sim_tclabP(MV,OPV,Ts,PVtemp1,PVtemp2,Kp,Tlead1,Tlead2,Theta1)
        sim_tclabP(DV,ODV,Ts,DVtemp1,DVtemp2,Kd,Tlag1,Tlag2,Theta2)
        PV.append(OPV[-1]+ODV[-1]+PV0-Kp*MV0)
    
    fig,(ax,bx,cx,dx) = plt.subplots(4)
    fig.set_figheight(16)
    fig.set_figwidth(22)
    
    ax.plot(t,Man,'black')
    ax.set(ylabel='Value of MAN')

    bx.plot(t,MV,'blue',label = "MV")
    bx.set(ylabel='Value of MV')

    cx.plot(t,PV[:-1],'green',label="PV")
    cx.plot(t,SP,'black',label = "SP")
    cx.set(ylabel='Value of PV and SP')
    cx.legend(loc='best')


    dx.plot(t,DV,'red')
    dx.set(ylabel='Value of DV')
    return
    
    
class PID:
    
    def __init__(self, parameters):
        
        self.parameters = parameters
        self.parameters['Kc'] = parameters['Kc'] if 'Kc' in parameters else 1.0
        self.parameters['Ti'] = parameters['Ti'] if 'Ti' in parameters else 1.0
        self.parameters['Td'] = parameters['Td'] if 'Td'in parameters else 10
        self.parameters['alpha'] = parameters['alpha'] if 'alpha' in parameters else 0.4