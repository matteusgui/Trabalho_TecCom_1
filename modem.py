import numpy as np



class Modem:
    def __init__(self, fs, bufsz, ans=False):
        self.bit_rate = 300
        self.bits =[]
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
        if len(bits) == 0:
            bits.extend([1,1,1,1,1,1,1,1,])
        self.bits.extend(bits)


    def get_samples(self):
    	#MODULADOR
        #x = []
        y = np.zeros(self.bufsz*len(self.bits))
        t = 0.0
        fi = 0
        #print("bits=",self.bits)
        for i in range(len(self.bits)):
            w = (self.tx_omega0 if self.bits[i] == 0 else self.tx_omega1)
            fi = fi-w*t
            for j in range(self.bufsz):
                #x.append(t)
                #y.append(np.sin(w*t+fi))
                y[j+self.bufsz*i] = np.sin(w*t+fi)
                t = t+1/self.fs
            fi = w*t + fi
        self.bits=[]
        return y

    # Demodulação

    def put_samples(self, data):
        pass

    def get_bits(self):
        return []
