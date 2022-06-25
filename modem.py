import numpy as np
import matplotlib.pyplot as plt
import scipy.signal

class Modem:
    def __init__(self, fs, bufsz, ans=False):
        self.bit_rate = 300
        self.fi = 0
        self.bits =[]
        self.s = []
        self.fs = fs  # taxa de amostragem
        self.bufsz = bufsz  # quantidade de amostas que devem ser moduladas por vez
        # frequências de modulação (upload)
        self.tx_omega0 = 2*np.pi*(1080 + 100)
        self.tx_omega1 = 2*np.pi*(1080 - 100)
        # frequências de demodulação (download)
        self.rx_omega0 = 2*np.pi*(1750 + 100)
        self.rx_omega1 = 2*np.pi*(1750 - 100)
        # se o modem estiver atendendo uma ligação
        if ans:
            # inverte as frequências
            self.tx_omega0, self.rx_omega0 = self.rx_omega0, self.tx_omega0
            self.tx_omega1, self.rx_omega1 = self.rx_omega1, self.tx_omega1

    # Modulação

    def put_bits(self, bits):
        self.bits.extend(bits)


    def get_samples(self):
    	#MODULADOR
        y = np.zeros(self.bufsz)
        t = 0.0

        if len(self.bits) == 0:
            self.bits.append(1)

        w = (self.tx_omega1 if self.bits[0] else self.tx_omega0)
        self.fi = self.fi-w*t
        for j in range(self.bufsz):
            y[j] = np.sin(w*t+self.fi)
            t = t+1/self.fs

        self.fi = w*t + self.fi
        self.bits.pop(0)
        return y
    # Demodulação

    def put_samples(self, data):
        self.s.extend(data)

    def get_bits(self):
        plt.plot(self.s)
        plt.show()
        s = self.s
        fs=48000
        T=1/fs
        L=fs//300
        pi = np.pi
        omega0=2*pi*1180
        omega1=2*pi*980

        v0i=np.zeros(len(s))
        v0r=np.zeros(len(s))
        v1i=np.zeros(len(s))
        v1r=np.zeros(len(s))

        r=0.99

        for n in range(1,len(s)):
                v0r[n] = s[n] - r**L * np.cos(omega0*L*T)*s[n-L] + r*np.cos(omega0*T)*v0r[n-1] - r*np.sin(omega0*T)*v0i[n-1]
                v0i[n] = -r**L*np.sin(omega0*L*T)*s[n-L] + r*np.cos(omega0*T)*v0i[n-1] + r*np.sin(omega0*T)*v0r[n-1]
                v1r[n] = s[n] - r**L * np.cos(omega1*L*T)*s[n-L] + r*np.cos(omega1*T)*v1r[n-1] - r*np.sin(omega1*T)*v1i[n-1]
                v1i[n] = -r**L*np.sin(omega1*L*T)*s[n-L] + r*np.cos(omega1*T)*v1i[n-1] + r*np.sin(omega1*T)*v1r[n-1]

        rho = v1r**2+v1i**2+v0r**2+v0i**2   # carrier detection

        c = abs(v1r**2+v1i**2-v0r**2-v0i**2)
        v = np.zeros(len(c))
        y = np.zeros(len(c))
        r=0.9999
        for n in range(1,len(s)):
                v[n] = (1-r)*c[n] + 2*r*np.cos(2*pi*300/fs)*v[n-1] - r**2*v[n-2]
                y[n] = v[n] - v[n-2]

        plt.plot(v1r**2+v1i**2-v0r**2-v0i**2, 'r')
        #plt.plot(y,'g')

        filt=scipy.signal.firwin(40, 300, pass_zero='lowpass', fs=48000)
        #plt.plot(1500*((y[1:]>0)&(y[:-1]<0)), 'g')
        plt.plot(np.convolve(v1r**2+v1i**2-v0r**2-v0i**2, filt), 'r')  # 46 amostras de delay

        plt.show()
        return []
