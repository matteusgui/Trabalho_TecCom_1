import numpy as np
import scipy.signal

class Modem:
    def __init__(self, fs, bufsz, ans=False):
        self.bit_rate = 300
        self.fi = 0
        self.bits =[]
        self.s = []
        self.sbuffer= np.zeros(fs//self.bit_rate).tolist()
        self.v0i= 0
        self.v0r= 0
        self.v1i= 0
        self.v1r= 0
        self.y = None
        self.v = [0, 0, 0]
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

        self.fi = (w*t + self.fi)%(2*np.pi)
        self.bits.pop(0)
        return y
    # Demodulação

    def put_samples(self, data):
        self.s = self.s.copy() + data.tolist()

    def get_bits(self):
        

        fs= self.fs
        T = 1/fs
        L = fs//self.bit_rate #self.bit_rate = 300
        pi = np.pi
        omega0=self.rx_omega0 
        omega1=self.rx_omega1

        #print("Tam = ",tam)

        cutoff_index = len(self.s) - len(self.s)%L
        if(cutoff_index < L):
            return []
        
        s = self.sbuffer + self.s[:cutoff_index]
        leftover_s = self.s[cutoff_index:]

        v0i= self.v0i
        v0r= self.v0r
        v1i= self.v1i
        v1r= self.v1r

        V0i= np.zeros(len(s) - L)
        V0r= np.zeros(len(s) - L)
        V1i= np.zeros(len(s) - L)
        V1r= np.zeros(len(s) - L)

        r=0.99

        #print("len s=", len(s))
        #print("len v0=", len(v0i),len(v0r),len(v1i),len(v1r))

        for n in range(L,len(s)):
            V0r[n-L] = s[n] - r**L * np.cos(omega0*L*T)*s[n-L] + r*np.cos(omega0*T)*v0r - r*np.sin(omega0*T)*v0i
            V0i[n-L] = -r**L*np.sin(omega0*L*T)*s[n-L] + r*np.cos(omega0*T)*v0i + r*np.sin(omega0*T)*v0r
            V1r[n-L] = s[n] - r**L * np.cos(omega1*L*T)*s[n-L] + r*np.cos(omega1*T)*v1r - r*np.sin(omega1*T)*v1i
            V1i[n-L] = -r**L*np.sin(omega1*L*T)*s[n-L] + r*np.cos(omega1*T)*v1i + r*np.sin(omega1*T)*v1r

            v0r = V0r[n-L]
            v0i = V0i[n-L]
            v1r = V1r[n-L]
            v1i = V1i[n-L]
        
        self.v0i = v0i
        self.v0r = v0r
        self.v1i = v1i
        self.v1r = v1r

        rho = V1r**2+V1i**2+V0r**2+V0i**2   # carrier detection

        c = abs(V1r**2+V1i**2-V0r**2-V0i**2)
        v = self.v.copy()

        

        y = np.zeros(len(c))
        
        r = 0.9999
        for n in range(len(c)):
                v[2] = (1-r)*c[n] + 2*r*np.cos(2*pi*300/fs)*v[1] - r**2*v[0]
                y[n] = v[2] - v[0]

                v[0] = v[1]
                v[1] = v[2]

        self.v = v.copy()

        #print(f"help me = {self.y}")
        if self.y != None:
            y = np.concatenate((self.y, y))
        
        delta = V1r**2+V1i**2-V0r**2-V0i**2

        ponto_amostra = 1500*((y[1:]>0)&(y[:-1]<0))

        bits = []
        for i in range (len(ponto_amostra)):
            if ponto_amostra[i] != 0:
                bits.append(1 if delta[i] > 0.5 else 0)
        #for i in range(L//2 + 46, len(delta), L):
        #   bits.append(1 if delta[i] > 0.5 else 0)

        self.s = leftover_s.copy()
        self.sbuffer = s[len(s)-L:] 
        self.y = np.array([y[-1]])

        return bits
