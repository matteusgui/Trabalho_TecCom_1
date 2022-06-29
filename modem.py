import numpy as np
import matplotlib.pyplot as plt
import scipy.signal

class Modem:
    def __init__(self, fs, bufsz, ans=False):
        self.bit_rate = 300
        self.lastn = 0
        self.lastnlist = []
        self.fi = 0
        self.bits =[]
        self.s = []
        self.sbuffer=[]
        self.v0i=np.zeros(bufsz)
        self.v0r=np.zeros(bufsz)
        self.v1i=np.zeros(bufsz)
        self.v1r=np.zeros(bufsz)
        self.y = []
        self.v = []
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
        self.s = np.append(self.s,data)

    def get_bits(self):
        s = []
        s.extend(self.sbuffer)
        if len(s) == 0:
            s.extend(np.zeros(self.bufsz).tolist())
        s.extend(self.s)
        self.sbuffer= s[len(s)-self.bufsz:]

        fs= self.fs
        T = 1/fs
        L = fs//300
        pi = np.pi
        omega0=self.rx_omega0
        omega1=self.rx_omega1

        plotar = len(self.s)<20*160
        tam = len(self.v0i)
        #print("Tam = ",tam)

        self.v0i = np.append(self.v0i,np.zeros(len(self.s)))
        self.v0r = np.append(self.v0r,np.zeros(len(self.s)))
        self.v1i = np.append(self.v1i,np.zeros(len(self.s)))
        self.v1r = np.append(self.v1r,np.zeros(len(self.s)))

        v0i=self.v0i
        v0r=self.v0r
        v1i=self.v1i
        v1r=self.v1r

        r=0.99

        #print("len s=", len(s))
        #print("len v0=", len(v0i),len(v0r),len(v1i),len(v1r))

        for n in range(tam,len(s)):
                v0r[n] = s[n] - r**L * np.cos(omega0*L*T)*s[n-L] + r*np.cos(omega0*T)*v0r[n-1] - r*np.sin(omega0*T)*v0i[n-1]
                v0i[n] = -r**L*np.sin(omega0*L*T)*s[n-L] + r*np.cos(omega0*T)*v0i[n-1] + r*np.sin(omega0*T)*v0r[n-1]
                v1r[n] = s[n] - r**L * np.cos(omega1*L*T)*s[n-L] + r*np.cos(omega1*T)*v1r[n-1] - r*np.sin(omega1*T)*v1i[n-1]
                v1i[n] = -r**L*np.sin(omega1*L*T)*s[n-L] + r*np.cos(omega1*T)*v1i[n-1] + r*np.sin(omega1*T)*v1r[n-1]

        rho = v1r**2+v1i**2+v0r**2+v0i**2   # carrier detection

        c = abs(v1r**2+v1i**2-v0r**2-v0i**2)
        v = np.zeros(len(c))

        
        #self.v = np.append(self.v,np.zeros(len(self.s)))

        if len(self.v)>0:
            v[0] = self.v[0]
            v[1] = self.v[1]
        y = np.zeros(len(c))
        if len(self.y)>0:
            y[0] = self.y[0]
            y[1] = self.y[1]
        
        r = 0.9999
        for n in range(2,len(c)):
                v[n] = (1-r)*c[n] + 2*r*np.cos(2*pi*300/fs)*v[n-1] - r**2*v[n-2]
                y[n] = v[n] - v[n-2]

        self.v = v[len(v)-2:]
        #print("len de v = ", len(self.v))
        self.y = y[len(y)-2:]
        #if plotar:
            #plt.plot(y,'r')
            #plt.plot(v,'g')
        delta = v1r**2+v1i**2-v0r**2-v0i**2
        #print(delta)
        #plt.plot(delta,'r')


        filt=scipy.signal.firwin(40, 300, pass_zero='lowpass', fs=self.fs)
        ponto_amostra = 1500*((y[1:]>0)&(y[:-1]<0))
        #ponto_amostra = np.append(ponto_amostra,[0])
        bits = []
        #plt.plot(ponto_amostra)
        #delta = np.convolve(delta, filt)
        #print("len amostra=",len(ponto_amostra))
        #print ("len delta =",len(delta))
        #if(sum(ponto_amostra)!= 3000):
            #print("sum = ", sum(ponto_amostra))


        for i in range (len(ponto_amostra)):
            if ponto_amostra[i] != 0 and not self.lastnlist.count(i):# and not plotar or plotar and i == 2*self.bufsz-50:
                bits.append(1 if delta[i] > 0 else 0)
                for j in range (-self.bufsz//4,self.bufsz//4+1):
                    self.lastnlist.append(i+j)
            #if(i==2*self.bufsz-50):
                #print(i,"=>",y[i])
                #print(" ") 
            #if(i==self.bufsz-50):
                #print(i,"=>",y[i])

        for i in range (len(self.lastnlist)):
            self.lastnlist[i] = self.lastnlist[i]-self.bufsz+2

        #print("lastnlist = ", self.lastnlist)

        if(len(self.lastnlist)>0):
            while(len(self.lastnlist)>0 and self.lastnlist[0]<0):
                self.lastnlist.pop(0)


                #self.lastn = i
                #self.lastn = delta[i]
                #print("lastn", delta[i-1],delta[i],delta[i+1])


        #print("bits =",bits)
        #print("len(self.s)=",len(self.s))
        #bits = bits[:]
        self.s = []
        #print("bits =",bits)
        self.v0i = v0i[len(v0i)-self.bufsz:]
        self.v0r = v0r[len(v0r)-self.bufsz:]
        self.v1i = v1i[len(v1i)-self.bufsz:]
        self.v1r = v1r[len(v1r)-self.bufsz:]
        #self.s = s[:-self.bufsz]

        #bits = bits[self.lastn:]
        #self.lastn = 0 if len(bits) == 0 else 1

        #print(bits)
        #print("bits = ", bits)
        #if(plotar):
            #plt.plot(delta,'r')
            #plt.plot(ponto_amostra,'g')
            #plt.show()

        #plt.show()

        return bits
