import numpy as np
import math
import matplotlib.pyplot as plt
from script_3 import Time

def hildebrand_sekhon(data, navg):
    data = data.copy()
    sortdata = np.sort(data,axis=None)
    lenOfData = len(sortdata)
    nums_min = lenOfData/10
    if (lenOfData/10) > 2:
      nums_min = lenOfData/10
    else:
      nums_min = 2
    sump = 0.
    sumq = 0.
    j = 0
    cont = 1
    while((cont==1)and(j<lenOfData)):
      sump += sortdata[j]
      sumq += sortdata[j]**2
      j += 1
      if j > nums_min:
        rtest = float(j)/(j-1) + 1.0/navg
        if ((sumq*j) > (rtest*sump**2)):
          j = j - 1
          sump  = sump - sortdata[j]
          sumq =  sumq - sortdata[j]**2
          cont = 0
    lnoise = sump /j
    stdv = np.sqrt((sumq - lnoise**2)/(j - 1))
    return lnoise

def getNoisebyHildebrand(data_spc):
    noise    = np.zeros(1,dtype='f')
    data_spc = data_spc.reshape(1,data_spc.shape[0],data_spc.shape[1])
    for channel in range(1):
        daux = data_spc[channel,:,:]
        noise[channel] = hildebrand_sekhon(daux, 1)
    return noise
    
def getNoisebyHeights(data_spc):
    noise = np.zeros((1,data_spc.shape[1]),dtype='f')
    data_spc = data_spc.reshape(1,data_spc.shape[0],data_spc.shape[1])
    for channel in range(1):
    	daux = data_spc[channel,:,:]
    	for height in range(daux.shape[1]):
    		noise[channel][height] = hildebrand_sekhon(daux[:,height], 1)

    return noise
    
def freq_vel_Range(ipp,freq,nFFTPoints):
    PRF      = 1/ipp
    fmax     = PRF/2
    C        = 3.0e8
    _lambda_ = C/(freq*1e6)
    vmax     = fmax*_lambda_/2.0
    print("Freq operation:", freq)
    print("Fmax:", fmax)
    print("Lambda:",round(_lambda_,5))
    print("Vmax:", round(vmax,2))
    deltafreq = fmax/nFFTPoints
    print("DeltaFreq",deltafreq)
    freqrange = 2*deltafreq*(np.arange(nFFTPoints)-nFFTPoints/2.0)
    print("Freq_min:",freqrange[0],"Freq_max:",freqrange[-1])
    deltavel  = vmax/nFFTPoints
    print("Deltavel",deltavel)
    velrange = 2*deltavel*(np.arange(nFFTPoints)-nFFTPoints/2.0)
    print("Vel_min:",round(velrange[0],2),"Vel_max:",round(velrange[-1],2))
    return freqrange,velrange
    
def removeDC(data_spc,mode=2):
        jspectra     = data_spc.reshape(1,data_spc.shape[0],data_spc.shape[1])
        num_chan = jspectra.shape[0]
        num_hei = jspectra.shape[2]
        freq_dc = int(jspectra.shape[1] / 2)
        ind_vel = np.array([-2, -1, 1, 2]) + freq_dc
        ind_vel = ind_vel.astype(int)

        if ind_vel[0] < 0:
            ind_vel[list(range(0, 1))] = ind_vel[list(range(0, 1))] + self.num_prof

        if mode == 1:
            jspectra[:, freq_dc, :] = (jspectra[:, ind_vel[1], :] + jspectra[:, ind_vel[2], :]) / 2  # CORRECCION

        if mode == 2:
            vel = np.array([-2, -1, 1, 2])
            xx = np.zeros([4, 4])

            for fil in range(4):
                xx[fil, :] = vel[fil]**np.asarray(list(range(4)))

            xx_inv = np.linalg.inv(xx)
            xx_aux = xx_inv[0, :]

            for ich in range(num_chan):
                yy = jspectra[ich, ind_vel, :]
                jspectra[ich, freq_dc, :] = np.dot(xx_aux, yy)
                junkid = jspectra[ich, freq_dc, :] <= 0
                cjunkid = sum(junkid)
                if cjunkid.any():
                    jspectra[ich, freq_dc, junkid.nonzero()] = (jspectra[ich, ind_vel[1], junkid] + jspectra[ich, ind_vel[2], junkid]) / 2

        return jspectra[0]

def GetMoments(data_spc,absc,n):
        '''
        Function GetMoments()

        Input:
            spec: spectrogram
            absc: abscissaList
            n : noise

        Output:
            params : -Radial velocity(DOPPLER)
                     -power
                     -spectralwidth
        '''
        #data_param = numpy.zeros((data.shape[0], 4, data.shape[2]))
        #This change has been done for HF processing
        print("\nShape",data_spc.shape)
        spec    = data_spc.reshape(1,data_spc.shape[0],data_spc.shape[1])

        data_param = np.zeros((spec.shape[0], 4, spec.shape[2])) # Jm : harcoded to show the limits frequencies
        print("Error 20",data_param.shape)
        data_param[0]= CalculateMoments(data_spc, absc, n)
        return data_param

def CalculateMoments(oldspec, oldfreq, n0, nicoh = None, graph = None, smooth = None, type1 = None, fwindow = None, snrth = None, dc = None, aliasing = None, oldfd = None, wwauto = None):

    if (nicoh == None): nicoh = 1
    if (graph == None): graph = 0
    if (smooth == None): smooth = 0
    elif (smooth < 3): smooth = 0

    if (type1 == None): type1 = 0
    if (fwindow == None): fwindow = np.zeros(oldfreq.size) + 1 #abscissaList
    if (snrth == None): snrth = -3
    if (dc == None): dc = 0
    if (aliasing == None): aliasing = 0
    if (oldfd == None): oldfd = 0
    if (wwauto == None): wwauto = 0

    if (n0 < 1.e-20):   n0 = 1.e-20

    freq = oldfreq # Doppler velocity values
    vec_power = np.zeros(oldspec.shape[1])
    vec_fd = np.zeros(oldspec.shape[1])
    vec_w = np.zeros(oldspec.shape[1])
    vec_snr = np.zeros(oldspec.shape[1])
    vec_fv = np.zeros(oldspec.shape[1])#First Valid frequency
    vec_lv = np.zeros(oldspec.shape[1])#Last Valid Frequency

    for ind in range(oldspec.shape[1]):

        spec = oldspec[:,ind]
        #TODO : if snr = (spec2.mean()-n0)/n0 SNR es menor que 0.3dB no hagas el resto.
        #TODO : hacer un noise special para el slice metodo privado de ParametersProc
        aux = spec*fwindow[0:len(spec)] #Jm:hardcoded to match with lenghts
        max_spec = aux.max()
        m = list(aux).index(max_spec)
        #Smooth
        if (smooth == 0):   spec2 = spec
        else:   spec2 = scipy.ndimage.filters.uniform_filter1d(spec,size=smooth)

        #    Calculo de Momentos
        bb = spec2[range(m,spec2.size)]
        bb = (bb<n0).nonzero()
        bb = bb[0]

        ss = spec2[range(0,m + 1)]
        ss = (ss<n0).nonzero()
        ss = ss[0]

        if (bb.size == 0):
            bb0 = spec.size - 1 - m
        else:
            bb0 = bb[0] - 1
            if (bb0 < 0):
                bb0 = 0

        if (ss.size == 0):   ss1 = 1
        else: ss1 = max(ss) + 1

        if (ss1 > m):   ss1 = m

        valid = np.asarray(range(int(m + bb0 - ss1 + 1))) + ss1
        #print 'valid[0]:',freq[valid[0]]
        #print 'valid[-1]:',freq[valid[-1]]
        power = ((spec2[valid] - n0)*fwindow[valid]).sum() # m_0 = first moments
        #Todo probar la estimacion de fd con el calculo de ruido por perfil.
        fd = ((spec2[valid]- n0)*freq[valid]*fwindow[valid]).sum()/power # m_1=radial velocity = frequecy doppler?
        w = math.sqrt((  (spec2[valid] - n0)*fwindow[valid]  *(freq[valid]- fd)**2   ).sum()/power)
        snr = (spec2.mean()-n0)/n0

        if (snr < 1.e-20) :
            snr = 1.e-20

        vec_power[ind] = power
        vec_fd[ind] = fd
        vec_w[ind] = w
        vec_snr[ind] = snr
        vec_fv[ind]=freq[valid[0]]
        vec_lv[ind]=freq[valid[-1]]
        #vec_sw[ind] = sw

        #else : vec_power[ind] = un numero x, fd , w y snr igual.
    #moments = np.vstack((vec_snr, vec_power, vec_fd, vec_w,vec_fv,vec_lv))
    moments = np.vstack((vec_snr, vec_power, vec_fd, vec_w))

    return moments



####################--Parametros--#########################
print('##--Parameters--###')
ipp        = 0.0004    #IPP equivalente a 0.1 seg
freq       = 9345      #FREQ 0 radar 9345
nFFTPoints = 240.0     #600
station=['JRO_A']
option='freq'
frequency = 9345       #9345 [2.72,3.64]
DC = 1 		       #[0,1]
axis = 10 	       #list(np.arange(0.5,275.5,0.5))

###########################################################

def spectraPlot(j,station=station,frequency=frequency,option=option,DC=DC,axis=axis):
    fig, axes = plt.subplots(nrows=1,ncols=1,figsize=(20,6))
    power     = Time()
    power     = np.fft.fftshift(power,axes=0)
    noise     = getNoisebyHildebrand(power)
    nFFTPoints=power.shape[0]
    if DC==1:
        power   = removeDC(data_spc=power,mode=2)
    print("IPP",ipp)
    print("Frequency",frequency)
    print("nFFTPoints",nFFTPoints)
    freq_range,vel_range = freq_vel_Range(ipp,frequency,nFFTPoints)
    RangeMin, RangeMax = 0, 60
    if option=='freq':
      if axis>5:
        axis=5
      param=GetMoments(power,freq_range,noise)
      data_param=param[:,1:,:]
      data_snr  =param[:,0,:]
      MinFreq, MaxFreq = freq_range[0],freq_range[-1]
      plt.imshow(10*np.log10(power.T), origin='lower', cmap='jet',vmin = 45, vmax = 110, interpolation='spline36', \
            extent=[MinFreq, MaxFreq, RangeMin, RangeMax], aspect='auto')
      #y=np.arange(0,60000,7.5) #?
      #radialvelocity = data_param[0, 1, :]
      #SpectralWidth = data_param[0, 2, :]
      #plt.plot(SpectralWidth*3.0+radialvelocity, y, color="red", linestyle="solid", lw=1)
      #plt.plot(radialvelocity-SpectralWidth*3.0, y, color="red", linestyle="solid", lw=1)
      #plt.plot(radialvelocity, y, color='black', linestyle="solid", lw=1)
      plt.xlabel(r'$\mathrm{Frequency \,(Hz)}$')
      plt.ylabel(r'$\mathrm{Distance \,(Km)}$')
    else:
      if frequency>3:
        if axis>205:
          axis=205
      else:
        if axis>275:
          axis=275
      param=GetMoments(power,vel_range,noise)
      data_param=param[:,1:,:]
      data_snr  =param[:,0,:]
      MinVel, MaxVel = vel_range[0],vel_range[-1]
      plt.imshow(10*np.log10(power.T), origin='lower', cmap='jet',vmin = 40, vmax = 55, interpolation='spline36', \
            extent=[MinVel, MaxVel, RangeMin, RangeMax], aspect='auto')
      id_line_counter = 0
      #y=np.arange(0,60000,7.5)
      #radialvelocity = data_param[0, 1, :]
      #plt.plot(radialvelocity, y, color='black', linestyle="solid", lw=1)
      #SpectralWidth = data_param[0, 2, :]
      #plt.plot(SpectralWidth*3.0+radialvelocity, y, color="red", linestyle="solid", lw=1)
      #plt.plot(radialvelocity-SpectralWidth*3.0, y, color="red", linestyle="solid", lw=1)
      plt.xlabel(r'$\mathrm{Velocity \,(Hz)}$')
      plt.ylabel(r'$\mathrm{Distance \,(Km)}$')
      
    cbar=plt.colorbar()
    cbar.set_label(r'$\log_{10}\max_t P(f,r,t)$', size=10)
    axes =plt.gca()

    #fig.suptitle('Spectra ', fontsize=15)

    plt.show()

spectraPlot(0)

