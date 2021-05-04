############################################################
# Code prepared by Damien Depannemaecker and Mallory Carlu #
# associated to the figure 2(d) of the paper entitled:     #
# Seizure-like propagation in spiking network models       #
############################################################

import matplotlib.pyplot as plt
import numpy as np
from brian2 import *


NProp1=[]
a=0
for NAmp in range(10):
    NProp1.append([])
    b=0
    for NbS in range(10):
        AmpStim=NAmp*5.+60
        TauP=20.+8*NbS
        NProp1[a].append(0)   
        for Nseed in range(100):
            NbSim=NbS
            Nsim=NbS

            FRexc1 = np.load('Results/AD_popRateExc_Sim_'+str(TauP)+'_Amp_'+str(NAmp)+'Nseed_'+str(Nseed)+'.npy')

            if max(FRexc1)>(AmpStim):
                NProp1[a][b]=NProp1[a][b]+1 
        b=b+1
    a=a+1   


NProp=np.array(NProp1[::-1])

plt.imshow(NProp[2:10])
x = np.array([int(NAmp*5.+60) for NAmp in range(10)])
y = np.array([int(20.+8*NbS) for NbS in range(8)][::-1])
print(x)
plt.xticks(range(10), x)
plt.yticks(range(8), y)

plt.xlabel('Amplitude (Hz)')
plt.ylabel(r'Slope, $\tau$ (ms)')
clb=plt.colorbar(orientation='vertical', fraction=0.038)
clb.set_label('number of Non-Propagations', labelpad=-0, y=0.5, rotation=90)
plt.show()

