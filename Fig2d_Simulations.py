
from brian2 import *
import numpy as np


#function use firing rate calculation
def bin_array(array, BIN, time_array):
    N0 = int(BIN/(time_array[1]-time_array[0]))
    N1 = int((time_array[-1]-time_array[0])/BIN)
    return array[:N0*N1].reshape((N1,N0)).mean(axis=1)


for Nseed in range(100):
    for NAmp in range(10):
        for NbS in range(10):
            NbSim=NbS
            Nsim=NbS



            seed(Nseed)
            start_scope()
            DT=0.1
            defaultclock.dt = DT*ms
            N1 = 2000#2000
            N2 = 8000#8000


            TotTime=4000
            duration = TotTime*ms



       # equation of the AdEx Model with "conductance-based" model of synapses 
            eqs='''
            dv/dt = (-GsynE*(v-Ee)-GsynI*(v-Ei)-gl*(v-El)+ gl*Dt*exp((v-Vt)/Dt)-w + Is)/Cm : volt (unless refractory)
            dw/dt = (a*(v-El)-w)/tau_w:ampere
            dGsynI/dt = -GsynI/Tsyn : siemens
            dGsynE/dt = -GsynE/Tsyn : siemens
            Is:ampere
            Cm:farad
            gl:siemens
            El:volt
            a:siemens
            tau_w:second
            Dt:volt
            Vt:volt
            Ee:volt
            Ei:volt
            Tsyn:second
            '''#% neuron_params

        # Populations----------------------------------------------------------------------------------

        # Population 1 - FS
            b1 = 0.0*pA
            G1 = NeuronGroup(N1, eqs, threshold='v > -47.5*mV', reset='v = -65*mV', refractory='5*ms', method='heun')
        #init:
            G1.v = -65*mV
            G1.w = 0.0*pA
            G1.GsynI=0.0*nS
            G1.GsynE=0.0*nS
        #parameters
            G1.Cm = 200.*pF
            G1.gl = 10.*nS
            G1.El = -65.*mV
            G1.Vt = -48.*mV
            G1.Dt = 0.5*mV
            G1.tau_w = 1.0*ms
            G1.a = 0.0*nS
            G1.Is = 0.0 

            G1.Ee=0.*mV
            G1.Ei=-80.*mV
            G1.Tsyn=5.*ms

        # Population 2 - RS
            b2 = 100.*pA
            G2 = NeuronGroup(N2, eqs, threshold='v > -40.*mV', reset='v = -65*mV; w += b2', refractory='5*ms',  method='heun')
            G2.v = -65.*mV
            G2.w = 0.0*pA
            G2.GsynI=0.0*nS
            G2.GsynE=0.0*nS
            G2.Cm = 200.*pF
            G2.gl = 10.*nS
            G2.El = -65.*mV
            G2.Vt = -50.*mV
            G2.Dt = 2.*mV
            G2.tau_w = 1000.*ms
            G2.a = 0.*nS
            G2.Is = 0.0*nA 
    
            G2.Ee=0.*mV
            G2.Ei=-80.*mV
            G2.Tsyn=5.*ms
    
    
            # external drive and seizure-like perturabation----------------------------------------------
            AmpStim=NAmp*5.+60  #80. #92.
            plat = 1000
            def heaviside(x):
                return 0.5 * (1 + np.sign(x))
    
    
            def input_rate(t, t1_exc, tau1_exc, tau2_exc, ampl_exc, plateau):
                    # t1_exc=10. # time of the maximum of external stimulation
                    # tau1_exc=20. # first time constant of perturbation = rising time
                    # tau2_exc=50. # decaying time
                    # ampl_exc=20. # amplitude of excitation
                inp = ampl_exc * (np.exp(-(t - t1_exc) ** 2 / (2. * tau1_exc ** 2)) * heaviside(-(t - t1_exc)) + \
                                  heaviside(-(t - (t1_exc+plateau))) * heaviside(t - (t1_exc))+ \
                                  np.exp(-(t - (t1_exc+plateau)) ** 2 / (2. * tau2_exc ** 2)) * heaviside(t - (t1_exc+plateau)))
                return inp
    
    
            t2 = np.arange(0, TotTime, DT)
            test_input = []
            TauP=20.+8*NbS
            for ji in t2:
                test_input.append(6.+input_rate(ji, 2000., TauP, TauP, AmpStim, plat))
            stimulus=TimedArray(test_input*Hz, dt=DT*ms)
            P_ed=PoissonGroup(8000, rates='stimulus(t)') #, dt=0.01*ms)
    
            # connections-----------------------------------------------------------------------------
    
            Qi=5.0*nS
            Qe=1.5*nS
    
            prbC= 0.05 #0.05
            prbC2= 0.05#0.065
            S_12 = Synapses(G1, G2, on_pre='GsynI_post+=Qi') #'v_post -= 1.*mV')
            S_12.connect('i!=j', p=prbC)
    
            S_11 = Synapses(G1, G1, on_pre='GsynI_post+=Qi')
            S_11.connect('i!=j',p=prbC)
    
            S_21 = Synapses(G2, G1, on_pre='GsynE_post+=Qe')
            S_21.connect('i!=j',p=prbC)
    
            S_22 = Synapses(G2, G2, on_pre='GsynE_post+=Qe')
            S_22.connect('i!=j', p=prbC)




            S_ed_in = Synapses(P_ed, G1, on_pre='GsynE_post+=Qe')
            S_ed_in.connect(p=prbC2)
    
            S_ed_ex = Synapses(P_ed, G2, on_pre='GsynE_post+=Qe')
            S_ed_ex.connect(p=prbC)#0.05)

            # monitor tools to record during simulation-------------------------------------------------
            #FRG1 = PopulationRateMonitor(G1)
            FRG2 = PopulationRateMonitor(G2)
            FRPed= PopulationRateMonitor(P_ed)

            # Run the simulation ----------------------------------------------------------------------
            Sim=(Nseed+1)*(NAmp+1)*(NbS+1)
            print('Starts simulation #'+str(Sim))
            run(duration)
            print('Ends simulation #'+str(Sim))  

            # Prepare and save data---------------------------------------------------------------------
            BIN=10
            time_array = np.arange(int(TotTime/DT))*DT

            #LfrG1=np.array(FRG1.rate/Hz)
            #TimBinned,popRateG1=bin_array(time_array, BIN, time_array),bin_array(LfrG1, BIN, time_array)
            LfrG2=np.array(FRG2.rate/Hz)
            TimBinned,popRateG2=bin_array(time_array, BIN, time_array),bin_array(LfrG2, BIN, time_array)
            LfrPed=np.array(FRPed.rate/Hz)
            TimBinned,popRatePed=bin_array(time_array, BIN, time_array),bin_array(LfrPed, BIN, time_array)

            np.save('Results/AD_popRateExc_Sim_'+str(TauP)+'_Amp_'+str(NAmp)+'Nseed_'+str(Nseed)+'.npy', popRateG2)
     











