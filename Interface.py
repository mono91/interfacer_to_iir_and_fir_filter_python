# -*- coding: utf-8 -*-
"""
Created on Sat Apr 18 00:15:23 2020

@author: jdmar
"""

from tkinter import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("TKAgg")
import math
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from numpy import cos, sin, pi, absolute, arange
from scipy.signal import kaiserord, lfilter, firwin, freqz
from pylab import figure, clf, plot, xlabel, ylabel, xlim, ylim, title, grid, axes, show
from scipy import signal as sp
import scipy.signal as signal
from matplotlib.pyplot import  close as close
from time import sleep



raiz=Tk()
raiz.title("Filtros")
raiz.config(bg="red")
frame=Frame(raiz)
frame.pack()
#☺frame.config(bg="green")
#frame.config(width="800",height="800")
#frame.config(bd=30)
#frame.config(relief="sunken")
frame.config(cursor="pirate")
#Funciones
texto0=Label(frame, text="FUNCIÓN = ")
texto0.grid(row=0,column=0)
ventana0=StringVar(frame)
ventana0.set('Ruido')
opciones0=['Ruido','Seno','Tren de pulsos','Triangular','Impulso','Real','Real_2','Voz']
menu0=OptionMenu(frame,ventana0,*opciones0)
menu0.config(width=11)
menu0.grid(row=0,column=2)
#Filtro
texto1=Label(frame, text="FILTRO      = ")
texto1.grid(row=1,column=0)
ventana1=StringVar(frame)
ventana1.set('FIR')
opciones1=['FIR','IIR']
menu1=OptionMenu(frame,ventana1,*opciones1)
menu1.config(width=11)
menu1.grid(row=1,column=2)
#Tipo
ventana11=StringVar(frame)
ventana11.set('Pasa-Bajas')
opciones11=['Pasa-Bajas','Pasa-Altas','Pasa-Banda']
menu11=OptionMenu(frame,ventana11,*opciones11)
menu11.config(width=11)
menu11.grid(row=1,column=3)
#Ventanas
texto2=Label(frame, text="VENTANA = ")
texto2.grid(row=2,column=0)
ventana2=StringVar(frame)
ventana2.set('kaiser')
opciones2=['kaiser','hamming','hanning','blackmanharris']
menu2=OptionMenu(frame,ventana2,*opciones2)
#menu2.config(width=15)
menu2.grid(row=2,column=2)
#Frecuencia de corte 1
texto3=Label(frame, text="Fc1 = ")
texto3.grid(row=0,column=4)
nombre0=StringVar()
cuadrotexto0=Entry(frame, textvariable=nombre0)
cuadrotexto0.grid(row=0,column=5)
#Frecuencia de corte 2
texto4=Label(frame, text="Fc2  = ")
texto4.grid(row=1,column=4)
nombre1=StringVar()
cuadrotexto1=Entry(frame, textvariable=nombre1)
cuadrotexto1.grid(row=1,column=5)
#Orden
texto5=Label(frame, text="Orden = ")
texto5.grid(row=2,column=4)
nombre2=StringVar()
cuadrotexto2=Entry(frame, textvariable=nombre2)
cuadrotexto2.grid(row=2,column=5)
#Ripple
texto6=Label(frame, text="Ripple(dB) = ")
texto6.grid(row=0,column=6)
nombre3=StringVar()
cuadrotexto3=Entry(frame, textvariable=nombre3)
cuadrotexto3.grid(row=0,column=7)
#Roll-of
texto7=Label(frame, text="Roll-off = ")
texto7.grid(row=1,column=6)
nombre4=StringVar()
cuadrotexto4=Entry(frame, textvariable=nombre4)
cuadrotexto4.grid(row=1,column=7)

#image1=PhotoImage(file="giphy.gif")
#imagen2=Label(frame, image=image1)
#imagen2.grid(row=1,column=2)
#texto=Label(frame, text="Entradas: ",padx=10)
#texto.grid(row=0,column=0, sticky="e")
#texto.config(bg="blue")
#nombre=StringVar()
#cuadrotexto=Entry(frame, textvariable=nombre)
#cuadrotexto.grid(row=0,column=1)

def codigoboton():
    
    func=ventana0.get()
    if func=='Ruido':
        sample_rate = 100.0
        nsamples = 400
        y = np.random.normal(0,1,nsamples)
        t =np.arange(len(y)) / sample_rate

        
    if func=='Seno':
        sample_rate = 100.0
        nsamples = 400
        t = arange(nsamples) / sample_rate
        y = cos(2*pi*0.5*t) + 0.2*sin(2*pi*2.5*t+0.1) + \
        0.2*sin(2*pi*15.3*t) + 0.1*sin(2*pi*16.7*t + 0.1) + \
            0.1*sin(2*pi*23.45*t+.8)
        
    if func=='Tren de pulsos':
        sample_rate = 100.0
        nsamples = 400
        t = arange(nsamples) / sample_rate
        amplitud=1
        y=((sp.square(2 * t)) * (amplitud / 2.0)) + (amplitud / 2.0)

    if func=='Triangular':
        sample_rate = 100.0
        nsamples = 400
        t = np.linspace(0, 1, 2000)
        simetria=0.5
        y=signal.sawtooth(2 * np.pi * 5 * t, simetria)
        
    if func=='Impulso':
        sample_rate = 100.0
        nsamples = 400
        y = signal.unit_impulse(100, 'mid')
        t =np.arange(len(y))
    if func=='Real':
        sample_rate = 250.0
        y= [-0.26122448, -0.26464245, -0.21141979, -0.24559948, -0.27782604, -0.04393932, 
            -0.16259167, -0.20018932, -0.22704479, -0.20214245, -0.33251354, -0.29589245, 
            -0.29296276, -0.16259167, -0.14257214, -0.08446667, -0.16210339, -0.05077526, 
            -0.03954479, -0.20849011, -0.1747987,  -0.17821667, -0.19432995, -0.18944714, 
            -0.17284557, -0.16259167, -0.11278698, -0.23632214, -0.2216737,  -0.19725964, 
            -0.15477917, -0.24804089, -0.21434948, -0.23143932, -0.24169323, -0.15722057, 
            -0.07860729, -0.17528698, -0.22850964, -0.19140026, -0.19432995, -0.14599011,
            -0.11083386, -0.19335339, -0.19481823, -0.2138612,  -0.29296276, -0.27099011,
            -0.19091198, -0.12206432, -0.14550182, -0.11718151, -0.12987682, -0.16991589,
            -0.24950573, -0.30370495, -0.37645886, -0.21727917, -0.16405651, -0.16112682,
            -0.18651745, -0.06884167, -0.22118542, -0.18700573, -0.20702526, -0.19481823,
            -0.19921276, -0.18065807, -0.0732362,  -0.13817761, -0.11962292, -0.04100964,
            -0.10155651, -0.13036511, -0.09179089, -0.02001354, -0.10643932, -0.12450573,
            -0.13768932, -0.16063854, -0.10106823, -0.04149792, -0.09960339, -0.06102917,
            -0.04638073,  0.01953724,  0.04981068, -0.03026745, -0.01513073,  0.02734974,
             0.04102161,  0.03906849,  0.06836536,  0.00195911, -0.00048229, -0.00829479,
             0.00293568, -0.02831432,  0.02783802, -0.07079479]
        t =np.arange(len(y)) / sample_rate
    if func=='Real_2':
        sample_rate = 250.0
        y=[-2.29025143e-01, -2.31466549e-01, -2.08029049e-01, -1.67501706e-01,
           -1.87521237e-01,  5.80842320e-02, -8.84001430e-02, -1.37716549e-01,
           -9.13298305e-02, -8.00993618e-02, -2.05587643e-01, -1.74825924e-01,
           -1.30392331e-01, -4.20134243e-02, -7.52165493e-02,  1.16975132e-02,
           -2.12367717e-05,  6.78498570e-02,  6.98029820e-02, -8.25407680e-02,
           -6.74040493e-02, -1.05489987e-01, -8.54704555e-02, -7.22868618e-02,
           -4.15251430e-02, -6.98454555e-02, -3.81071743e-02, -8.79118618e-02,
           -1.03536862e-01, -6.69157680e-02, -4.10368618e-02, -1.30392331e-01,
           -1.06466549e-01, -1.28439206e-01, -1.19161862e-01, -1.36931118e-02,
            8.76782573e-03, -6.39860805e-02, -1.11837643e-01, -8.98649868e-02,
           -3.81071743e-02, -2.63884243e-02,  1.51154820e-02, -7.47282680e-02,
           -9.18181118e-02, -7.37517055e-02, -1.39669674e-01, -8.69352993e-02,
            2.87873570e-02,  9.66584507e-02,  1.85334507e-02,  3.88501323e-03,
            1.41389195e-02,  9.55325728e-04, -6.15446743e-02, -1.02560299e-01,
           -1.70431393e-01, -3.90837368e-02, -2.78532680e-02, -4.49431118e-02,
           -6.88688930e-02,  6.49201695e-02, -1.33810299e-01, -6.10563930e-02,
           -1.31857174e-01, -1.15255612e-01, -9.42595180e-02, -1.08419674e-01,
           -5.12907680e-02, -6.74040493e-02, -8.88884243e-02, -4.54313930e-02,
           -7.76579555e-02, -1.05001706e-01, -4.25017055e-02,  1.26740757e-02,
           -7.22868618e-02, -7.91227993e-02, -9.76774868e-02, -9.81657680e-02,
           -6.54509243e-02,  3.51350132e-02, -1.71110805e-02,  1.21857945e-02,
            3.12287632e-02,  5.07600132e-02,  2.04865757e-02, -2.68767055e-02,
            2.34162632e-02,  6.78498570e-02,  4.73420445e-02,  6.54084507e-02,
            8.68928257e-02,  3.02522007e-02,  5.46662632e-02,  7.66389195e-02,
            6.68732945e-02,  6.10139195e-02,  1.03982669e-01,  8.76782573e-03]
        t =np.arange(len(y)) / sample_rate
    
    if func=='Voz':
        sample_rate = 44100
        nsamples = 400              
        y=P
        t =np.arange(len(y)) / sample_rate
    
    # The Nyquist rate of the signal.
    nyq_rate = sample_rate / 2.0
            
        
    filtro = ventana1.get()
    
    if filtro=='FIR':
        print(1)
        # The desired width of the transition from pass to stop,
        # relative to the Nyquist rate.  We'll design the filter
        # with a 5 Hz transition width.        
        width = int(nombre4.get())/pi        
        
        # The desired attenuation in the stop band, in dB.
        ripple_db= int(nombre3.get())        
        
        ventana=ventana2.get()
        
        filtro_tipo=ventana11.get()
        if filtro_tipo=='Pasa-Bajas':
                     
            # The cutoff frequency of the filter.
            cutoff_hz = float(nombre0.get())
            
            if ventana=='kaiser':
                # Compute the order and Kaiser parameter for the FIR filter.
                N, beta = kaiserord(ripple_db, width)
                # Use firwin with a Kaiser window to create a lowpass FIR filter.
                taps = firwin(N, cutoff_hz/pi, window= (ventana, beta))  
            else:
                N=int(nombre2.get())
                taps = firwin(N, cutoff_hz/pi, window= ventana)
                
            # Use lfilter to filter x with the FIR filter.
            filtered_x = lfilter(taps, 1.0, y)
            
            w, h = freqz(taps, worN=8000)
            
        if filtro_tipo=='Pasa-Altas':

            # The cutoff frequency of the filter.
            cutoff_hz = float(nombre0.get())
            
            if ventana=='kaiser':
                N, beta = kaiserord(ripple_db, width)
                M=N%2
                if M==0:
                    N=N+1
                # Use firwin with a Kaiser window to create a highpass FIR filter.
                taps = signal.firwin(N, cutoff_hz/pi , window =(ventana, beta), pass_zero=False)
            else:
                N=int(nombre2.get())
                taps = signal.firwin(N, cutoff_hz/pi , window =ventana, pass_zero=False)

            # Use lfilter to filter x with the FIR filter.
            filtered_x = lfilter(taps, 1.0, y)
            
            w, h = freqz(taps, worN=8000)
            
        if filtro_tipo=='Pasa-Banda':

            # The cutoff frequency of the filter.
            cutoff_hz = float(nombre0.get())
            cutoff_hz2= float(nombre1.get())
            
            if ventana=='kaiser':  
                N, beta = kaiserord(ripple_db, width)
                M=N%2
                if M==0:
                    N=N+1
                # Use firwin with a Kaiser window to create a bandpass FIR filter.
                taps = signal.firwin(N, cutoff = [cutoff_hz/pi, cutoff_hz2/pi], window = (ventana, beta), pass_zero = False)
            else:
                N=int(nombre2.get())
                taps = signal.firwin(N, cutoff = [cutoff_hz/pi, cutoff_hz2/pi], window = ventana, pass_zero = False)

            # Use lfilter to filter x with the FIR filter.
            filtered_x = lfilter(taps, 1.0, y)
            
            w, h = freqz(taps, worN=8000)
        

    if filtro=='IIR':
        print(2)
        filtro_tipo=ventana11.get()
        if filtro_tipo=='Pasa-Bajas':
            cutoff_hz = float(nombre0.get())
            b, a = signal.iirfilter(17, cutoff_hz/pi, rs=60,btype='low', analog=False)
            taps=b
            filtered_x = signal.filtfilt(b, a, y)
            w, h = signal.freqs(b, a, 1000)            
        
        if filtro_tipo=='Pasa-Altas':
            cutoff_hz = float(nombre0.get())
            b, a = signal.iirfilter(17, cutoff_hz/pi, rs=60,btype='high', analog=False)
            taps=b
            filtered_x = signal.filtfilt(b, a, y)            
            w, h = signal.freqs(b, a, 1000)
            
        if filtro_tipo=='Pasa-Banda':
            cutoff_hz = float(nombre0.get())
            cutoff_hz2= float(nombre1.get())
            b, a = signal.iirfilter(17, [cutoff_hz/pi, cutoff_hz2/pi], rs=60,btype='band', analog=False)
            taps=b
            filtered_x = signal.filtfilt(b, a, y)
            w, h = signal.freqs(b, a, 1000)



    #------------------------------------------------
    # Plot the FIR filter coefficients.
    #------------------------------------------------
    
    if filtro=='FIR':
        print(3)
        fig = Figure()
        ax3 = fig.add_subplot(221)
        ax3.set_title('Coeficientes del filtro')
        ax3.set_ylabel('Amplitud')
        ax3.grid(True)
        ax3.plot(taps, 'bo-', linewidth=2)
        #------------------------------------------------
        # Plot the magnitude response of the filter en dB.
        #------------------------------------------------

        ax2 = fig.add_subplot(222)
        #ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel('Gain (dB)')
        ax2.set_title('Respuesta en Frecuencia')
        #ax2.set_xlim((0,50))  #Para visualizar bien la respuesta en frecuencia del filtro IIR
        ax2.grid(True)
        y=absolute(h)
        db=[]
        for i in range(len(y)):
            db.append(10*math.log(y[i],10))
        ax2.plot(w, db, linewidth=2)
    
        #ax6 = fig.add_subplot(325)
        #ax6.set_xlabel('Frequencia normalizada (Hz)')
        #ax6.set_ylabel('Fase (rad)')
        #ax6.set_title('Respuesta en fase')
        #ax6.grid(True)
        #Fase_filtro=np.unwrap(np.arctan2(np.imag(h),np.real(h)))
        #ax6.plot(w/max(w),Fase_filtro)

    
    


        #------------------------------------------------
        # Plot the original and filtered signals.
        #------------------------------------------------

        # The phase delay of the filtered signal.
        delay = 0.5 * (N-1) / sample_rate

        # Plot the original signal.
        #fig.add_subplot(224).plot(t, y)
        # Plot the filtered signal, shifted to compensate for the phase delay.
        ax4 = fig.add_subplot(223)
        ax4.set_title('Señal filtrada')
        ax4.set_ylabel('Amplitud')
        ax4.set_xlabel('Muestras [n]')
        ax4.grid(True)
        ax4.plot(t-delay, filtered_x, 'r-')
        # Plot just the "good" part of the filtered signal.  The first N-1
        # samples are "corrupted" by the initial conditions.
        #fig.add_subplot(224).plot(t[N-1:]-delay, filtered_x[N-1:], 'g', linewidth=4)
    
    
        ax6 = fig.add_subplot(224)
        ax6.set_xlabel('Frequencia(Hz)')
        ax6.set_ylabel('Amplitud')
        ax6.set_title('Espectro de la señal filtrada')
        ax6.grid(True)
        w= np.linspace(-np.pi, np.pi, len(filtered_x), endpoint=False)
        Transf2=np.zeros(len(w),'complex')
        for n in np.arange(len(filtered_x)):
            Transf2 += filtered_x[n]*np.exp(-1j*w*n)        
        ax6.plot(w, abs(Transf2))  #Espectro de la señal filtrada   
        
    
    if filtro=='IIR':
        print(4)
        fig = Figure()
        ax3 = fig.add_subplot(221)
        ax3.set_title('Coeficientes del filtro')
        ax3.set_ylabel('Amplitud')
        ax3.grid(True)
        ax3.plot(taps, 'bo-', linewidth=2)
        #------------------------------------------------
        # Plot the magnitude response of the filter en dB.
        #------------------------------------------------

        ax2 = fig.add_subplot(222)
        #ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel('Gain (dB)')
        ax2.set_title('Respuesta en Frecuencia')  #Para visualizar bien la respuesta en frecuencia del filtro IIR
        #ax2.axis((10, 1000, -100, 10))
        ax2.grid(True)            
        ax2.plot(w, 20 * np.log10(abs(h)))
        #ax2.set_xlim(0,pi)
        
    
        #------------------------------------------------
        # Plot the original and filtered signals.
        #------------------------------------------------

        # Plot the filtered signal, shifted to compensate for the phase delay.
        ax4 = fig.add_subplot(223)
        ax4.set_title('Señal filtrada')
        ax4.set_ylabel('Amplitud')
        ax4.set_xlabel('Muestras [n]')
        #ax4.set_xlim(20,len(filtered_x))
        #ax4.set_ylim(-10,10)
        ax4.grid(True)
        ax4.plot(50*filtered_x, 'r-')
          
    
        ax6 = fig.add_subplot(224)
        ax6.set_xlabel('Frequencia(rad)')
        ax6.set_ylabel('Amplitud')
        ax6.set_title('Espectro de la señal filtrada')
        ax6.grid(True)
        #ax6.set_ylim(0,100)
        #ax6.set_xlim(-0.2,0.2)
        w= np.linspace(-np.pi, np.pi, len(filtered_x), endpoint=False)
        Transf2=np.zeros(len(w),'complex')
        for n in np.arange(len(filtered_x)):
            Transf2 += filtered_x[n]*np.exp(-1j*w*n)        
        ax6.plot(w,100*abs(Transf2))  #Espectrograma de la señal filtrada        
        
    global canvas 
    canvas = FigureCanvasTkAgg(fig, raiz)
    canvas.draw() 
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
    global toolbar2
    toolbar2 = NavigationToolbar2Tk(canvas, raiz)
    toolbar2.update()
        

def codigoboton2():
    
    func=ventana0.get()
    if func=='Ruido':
        sample_rate = 100.0
        nsamples = 400
        y = np.random.normal(0,1,nsamples)
        t =np.arange(len(y)) / sample_rate

        
    if func=='Seno':
        sample_rate = 100.0
        nsamples = 400
        t = arange(nsamples) / sample_rate
        y = cos(2*pi*0.5*t) + 0.2*sin(2*pi*2.5*t+0.1) + \
        0.2*sin(2*pi*15.3*t) + 0.1*sin(2*pi*16.7*t + 0.1) + \
            0.1*sin(2*pi*23.45*t+.8)
        
    if func=='Tren de pulsos':
        sample_rate = 100.0
        nsamples = 400
        t = arange(nsamples) / sample_rate
        amplitud=1
        y=4*((sp.square(2 * t)) * (amplitud / 2.0)) + (amplitud / 2.0)

    if func=='Triangular':
        sample_rate = 100.0
        nsamples = 400
        t = np.linspace(0, 1, 2000)
        simetria=0.5
        y=signal.sawtooth(2 * np.pi * 5 * t, simetria)
        
    if func=='Impulso':
        sample_rate = 100.0
        nsamples = 400
        y = signal.unit_impulse(100, 'mid')
        t =np.arange(len(y))
    if func=='Real':
        sample_rate = 250.0
        y= [-0.26122448, -0.26464245, -0.21141979, -0.24559948, -0.27782604, -0.04393932, 
            -0.16259167, -0.20018932, -0.22704479, -0.20214245, -0.33251354, -0.29589245, 
            -0.29296276, -0.16259167, -0.14257214, -0.08446667, -0.16210339, -0.05077526, 
            -0.03954479, -0.20849011, -0.1747987,  -0.17821667, -0.19432995, -0.18944714, 
            -0.17284557, -0.16259167, -0.11278698, -0.23632214, -0.2216737,  -0.19725964, 
            -0.15477917, -0.24804089, -0.21434948, -0.23143932, -0.24169323, -0.15722057, 
            -0.07860729, -0.17528698, -0.22850964, -0.19140026, -0.19432995, -0.14599011,
            -0.11083386, -0.19335339, -0.19481823, -0.2138612,  -0.29296276, -0.27099011,
            -0.19091198, -0.12206432, -0.14550182, -0.11718151, -0.12987682, -0.16991589,
            -0.24950573, -0.30370495, -0.37645886, -0.21727917, -0.16405651, -0.16112682,
            -0.18651745, -0.06884167, -0.22118542, -0.18700573, -0.20702526, -0.19481823,
            -0.19921276, -0.18065807, -0.0732362,  -0.13817761, -0.11962292, -0.04100964,
            -0.10155651, -0.13036511, -0.09179089, -0.02001354, -0.10643932, -0.12450573,
            -0.13768932, -0.16063854, -0.10106823, -0.04149792, -0.09960339, -0.06102917,
            -0.04638073,  0.01953724,  0.04981068, -0.03026745, -0.01513073,  0.02734974,
             0.04102161,  0.03906849,  0.06836536,  0.00195911, -0.00048229, -0.00829479,
             0.00293568, -0.02831432,  0.02783802, -0.07079479]
        t =np.arange(len(y)) / sample_rate
        
    if func=='Real_2':
        sample_rate = 250.0
        y=[-2.29025143e-01, -2.31466549e-01, -2.08029049e-01, -1.67501706e-01,
           -1.87521237e-01,  5.80842320e-02, -8.84001430e-02, -1.37716549e-01,
           -9.13298305e-02, -8.00993618e-02, -2.05587643e-01, -1.74825924e-01,
           -1.30392331e-01, -4.20134243e-02, -7.52165493e-02,  1.16975132e-02,
           -2.12367717e-05,  6.78498570e-02,  6.98029820e-02, -8.25407680e-02,
           -6.74040493e-02, -1.05489987e-01, -8.54704555e-02, -7.22868618e-02,
           -4.15251430e-02, -6.98454555e-02, -3.81071743e-02, -8.79118618e-02,
           -1.03536862e-01, -6.69157680e-02, -4.10368618e-02, -1.30392331e-01,
           -1.06466549e-01, -1.28439206e-01, -1.19161862e-01, -1.36931118e-02,
            8.76782573e-03, -6.39860805e-02, -1.11837643e-01, -8.98649868e-02,
           -3.81071743e-02, -2.63884243e-02,  1.51154820e-02, -7.47282680e-02,
           -9.18181118e-02, -7.37517055e-02, -1.39669674e-01, -8.69352993e-02,
            2.87873570e-02,  9.66584507e-02,  1.85334507e-02,  3.88501323e-03,
            1.41389195e-02,  9.55325728e-04, -6.15446743e-02, -1.02560299e-01,
           -1.70431393e-01, -3.90837368e-02, -2.78532680e-02, -4.49431118e-02,
           -6.88688930e-02,  6.49201695e-02, -1.33810299e-01, -6.10563930e-02,
           -1.31857174e-01, -1.15255612e-01, -9.42595180e-02, -1.08419674e-01,
           -5.12907680e-02, -6.74040493e-02, -8.88884243e-02, -4.54313930e-02,
           -7.76579555e-02, -1.05001706e-01, -4.25017055e-02,  1.26740757e-02,
           -7.22868618e-02, -7.91227993e-02, -9.76774868e-02, -9.81657680e-02,
           -6.54509243e-02,  3.51350132e-02, -1.71110805e-02,  1.21857945e-02,
            3.12287632e-02,  5.07600132e-02,  2.04865757e-02, -2.68767055e-02,
            2.34162632e-02,  6.78498570e-02,  4.73420445e-02,  6.54084507e-02,
            8.68928257e-02,  3.02522007e-02,  5.46662632e-02,  7.66389195e-02,
            6.68732945e-02,  6.10139195e-02,  1.03982669e-01,  8.76782573e-03]
        t =np.arange(len(y)) / sample_rate    
        
    if func=='Voz':
        sample_rate = 44100
        nsamples = 400  
        global P            
        P=data_int[0:nsamples]
        y=P
        t =np.arange(len(y)) / sample_rate
        
    # The Nyquist rate of the signal.
    nyq_rate = sample_rate / 2.0

    
    fig = Figure()
    ax1=fig.add_subplot(121)
    #ax1.set_xlabel('x')
    ax1.set_ylabel('Amplitd')
    ax1.grid(True)
    ax1.set_title('Señal')
    ax1.plot(t, y)
    
    ax5 = fig.add_subplot(122)
    #ax5.set_xlabel('Frequency (Hz)')
    ax5.set_ylabel('Amplitud')
    ax5.set_title('Espectro de la señal')
    ax5.grid(True)
    w= np.linspace(-np.pi, np.pi, len(y), endpoint=False)
    Transf=np.zeros(len(w),'complex')
    for n in np.arange(len(y)):
        Transf += y[n]*np.exp(-1j*w*n)        
    ax5.plot(w, abs(Transf))

        
        
    global canvas2 
    canvas2 = FigureCanvasTkAgg(fig, raiz)
    canvas2.draw() 
    canvas2.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
    global toolbar
    toolbar = NavigationToolbar2Tk(canvas2, raiz)
    toolbar.update()
    
    
    
#GRAFICAR   
boton=Button(frame, text="GRAFICAR ESPECTRO", command=codigoboton2)
boton.grid(row=0,column=8)

def quitargrafica2(canvas2):
    canvas2.get_tk_widget().destroy()
    toolbar.destroy()

#QUITAR GRÁFICA
boton=Button(frame, text="QUITAR ESPECTRO", command= lambda:  quitargrafica2(canvas2))
boton.grid(row=1,column=8 )

boton=Button(frame, text="FILTRAR SEÑAL", command=codigoboton)
boton.grid(row=0,column=9)

def quitargrafica(canvas):
    canvas.get_tk_widget().destroy()
    toolbar2.destroy()
    

#QUITAR GRÁFICA
boton=Button(frame, text="QUITAR FILTRO", command= lambda:  quitargrafica(canvas))
boton.grid(row=1,column=9 )

def GravarAudio():
    import pyaudio
    import struct
    CHUNK = 1024 * 4
    FORMAT = pyaudio.paInt16
    CHANNELS = 1
    RATE = 44100
    
    p = pyaudio.PyAudio()
    stream = p.open(
            format=FORMAT,
            channels=CHANNELS,
            rate=RATE,
            input=True,
            output=True,
            frames_per_buffer=CHUNK            
            )
    
    
    fig, ax = plt.subplots()
    x=np.arange(0,2*CHUNK,2)
    line, = ax.plot(x,np.random.rand(CHUNK))
    ax.set_ylim(-255,255)
    ax.set_xlim(0,CHUNK)
    
    while True:
        data = stream.read(CHUNK)
        global data_int
        data_int = np.array(struct.unpack(str(2 * CHUNK) + 'B', data),dtype='b')[::2]    
        line.set_ydata(data_int)
        fig.canvas.draw()
        fig.canvas.flush_events()
    

    
    
    
#ESCUCHAR SONIDO  
boton=Button(frame, text="    ESCUCHAR    ", command=GravarAudio)
boton.grid(row=0,column=10)    
      

raiz.mainloop()    