###---Modulación en Frecuencia con SDR para el radar meteorológico Sophy---###

import numpy as np
#import matplotlib.pyplot as plt
#from chirp_final import repet_val,valor_n,rep_seq
# DIRECTORIO UBICACION: /usr/local/lib/python3/dist-packages/gnuradio

def rep_seq(x, rep=10):
    L = len(x) * rep
    res = np.zeros(L, dtype=x.dtype)
    idx = np.arange(len(x)) * rep
    for i in np.arange(rep):
        res[idx + i] = x
    return res


##------##
def chirpStandard(sr_rx,A,fs_Hz,rep_Hz,f0_Hz,f1_Hz,phase_rad=0,mode=0):      		##rep_Hz inversa De la duracion del CHIRP
    T_rep = 1/rep_Hz                 		# Periodo CHIRP.             		//Tiempo de duración del CHIRP
    k   = (f1_Hz-f0_Hz)/T_rep        		# Chirp rate in Hz/s.        		//CHIRP rate indepenediente de lo demás
    # CAMBIO AQUI  n   = int(fs_Hz/rep_Hz)          		# Samples per repetition.    		//Número de puntitos en la duración del CHIRP
    T_rep= round(T_rep*1e6)
    print("help",T_rep)
    n = int(sr_rx*T_rep*1e-6)
    T_rep=T_rep*1e-6
    n_rep = int(fs_Hz/sr_rx)			# Número de repeticiones
    print("n_rep",n_rep)
    if mode==0:
      t_s = np.linspace(0,T_rep,n)   		# Chirp Sample times.        		//Arreglo de tiempos para la duración del CHIRP [0 ... Trep]
    else:
      t_s = np.linspace(-T_rep/2.0,T_rep/2.0,n)


    # Especial T_S
    T_S = np.linspace(0,T_rep,n*n_rep)
    # Phase,phi_Hz, is integral of frequency, f(t)=kt/2+f0
    freq_Hz = k*t_s/2.0+f0_Hz
    freq_Hz = rep_seq(freq_Hz,n_rep)
    print("freq_Hz",freq_Hz[0:32])
    phi_Hz  = freq_Hz*T_S
    #phi_Hz  = freq_Hz*t_s  			# PHI INSTANTENEA               	//f(t)*t -> F(t)     # ORIGINAL
   
    print("-------------REPITIENDO----------")
    print("phi_Hz",len(phi_Hz))
    #print('Longitud de frecuencia (puntos): ',len(phi_Hz))
    #phi_rad = 2*np.pi*phi_Hz         		# Convert to radians         		//Argumento de la señal CHIRP
    phi_rad = [i * np.pi for i in phi_Hz]
    #print(len(phi_rad))
    #phi_rad += phase_rad             		# Offser by user-specified inital phase
    return  T_S,A*np.exp([j * 1j for j in phi_rad]),T_rep,n*n_rep
    #return t_s,A*np.sin(phi_rad),T_rep,n  	#Just real or I, component.

##------##
def chirpCentral(sr_rx,A,fs_Hz,fc_Hz,rep_Hz,bw_Hz,phase_rad=0,mode=0):
     f0_Hz = fc_Hz - bw_Hz/2.0
     f1_Hz = fc_Hz + bw_Hz/2.0
     print("Band Width: ",bw_Hz, "Hz")
     print("f0: ",f0_Hz, "Hz")
     print("f1: ",f1_Hz, "Hz")
     
     return chirpStandard(sr_rx,A,fs_Hz=fs_Hz,rep_Hz=rep_Hz,f0_Hz=f0_Hz,f1_Hz=f1_Hz,mode=mode)

"""

#CHIRP
#-Parameters-#

IPP = 100/30               ##seconds
DC = 15                    ##percent %
f_central = 150            #Hz
samp_rate = 2000           #Hz
B1 = DC*samp_rate/100      #Hz

"""

def chirp_mod(sr_rx,A,ipp,dc,sr,fc,b):
  f_s = sr
  f_c = fc
  T_rep = dc*ipp/100
  rep_Hz = 1/T_rep
  B = b 
  
  return chirpCentral(sr_rx,A,f_s,f_c,rep_Hz,B)
  
"""

CHIRP USRP


"""
def chirp_usrp(sr_rx,A,ipp,dc,sr,fc,b):
  t_c,chirp_m,T_c,N = chirp_mod(sr_rx,A,ipp,dc,sr,fc,b)
  n_3 = N*100/dc
  n_2 = n_3-len(chirp_m)
  mod_chirp_ = np.hstack((chirp_m,np.zeros(int(n_2))))

  return mod_chirp_
  
a,b,c,d = chirpStandard(2500000,1,20000000,int(100/(15*0.0004)),500000,1500000,phase_rad=0,mode=0)
a = chirp_usrp(4000000,1,0.0004,15,20000000,0,4000000)
#print(a)
#print(len(a))

#(np.fromfile("/home/soporte/Descargas/chirp_amp_1_ipp_0.0004_dc_15_sr_20000000_fc_2343750_bw_4000000.bin",dtype=np.complex64)).tolist()
#(np.fromfile("/home/soporte/Descargas/chirp_amp_1_ipp_0.0004_dc_15_sr_20000000_fc_2001000_bw_4000000.bin",dtype=np.complex64)).tolist()
