import numpy as np



class Modem:
    def __init__(self, fs, bufsz, ans=False):
        self.bit_rate = 300
        self.fi = 0
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
        pass

    def get_bits(self):
        return []
