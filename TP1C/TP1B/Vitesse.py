import math

class Vitesse:
    "Classe pour calculer les proprietes de la vitesse"

    def __init__(self,atm,poids,spd,typeSpd,surf_a,long_re,Nz,typeFlaps=0):
        self.poids = poids
        self.spd = spd
        self.typeSpd = typeSpd
        
        #typeSpd : 0 ----> input : Mach
        #typeSpd : 1 ----> input : Vitesse vraie (kts)
        #typeSpd : 2 ----> input : Vitesse equivalente (kts)
        #typeSpd : 3 ----> input : Vitesse calibrée (kts)
        #typeSpd : 4 ----> input : Rapport de vitesse (V/Vsr)

        self.typeFlaps = typeFlaps
        #0 ---> Flap 00 / GU
        #1 ---> Flap 00 / GD 
        #2 ---> Flap 20 / GU 
        #3 ---> Flap 20 / GD 
        #4 ---> Flap 45 / GD 

        self.surf_a = surf_a
        self.long_re = 8.286
        self.atm = atm


        # Constantes du probleme -------------------------
        self.sigma_air= 1.4
        self.R_air = 287
        self.a0 = 661.48*1.6878

        self.nz = Nz

        self.CL_maxDict = {"0" : 1.65,
                    "1" : 1.60,
                    "2" : 1.85,
                    "3" : 1.80,
                    "4" : 2.10}

    def solve(self):
        self.a = self.a0*math.sqrt(self.atm.RapportTemp())

        if self.typeSpd == 0:
            self.qc = self.atm.Pression() * ((1 + ((self.spd)**2)/5)**(1/0.2857) - 1)
            
        elif self.typeSpd == 1:
            self.qc = ( math.exp(math.log(0.2*((self.spd*1.6878)/self.a)**2 +1)/0.2857) - 1 ) * self.atm.Pression()

        elif self.typeSpd == 2:
            self.qc = (math.exp(math.log( ( (self.atm.rho0*((self.spd*1.6878)**2)) /(7*self.atm.Pression()) ) +1 ) / 0.2857) -1 ) * self.atm.Pression()

        elif self.typeSpd == 3:
            self.qc = (math.exp((math.log(0.2*((self.spd*1.6878)/self.a0)**2 + 1)) / 0.2857 ) - 1) *self.atm.p0

        elif self.typeSpd == 4:
            self.cl = self.CL_maxDict[str(self.typeFlaps)]/(self.spd)**2
            self.spd = math.sqrt((295.37*self.poids)/(self.cl*self.surf_a))
            self.qc = (math.exp(math.log( ( (self.atm.rho0*((self.spd*1.6878)**2)) /(7*self.atm.Pression()) ) +1 ) / 0.2857) -1 ) * self.atm.Pression()


        self.Vc = self.a0 * (5* (((self.qc/self.atm.p0)+1)**0.2857 - 1)) ** 0.5
        self.Ve = (7 * (self.atm.Pression()/self.atm.rho0) * (((self.qc/self.atm.Pression()) + 1 ) ** 0.2857  -1))**0.5
        self.Vv = self.a * (5 *( ((self.qc/self.atm.Pression())+1)**0.2857 -1 ) )**0.5
        #self.Ve = self.Vv * self.atm.RapportDensite()**0.5
        self.M = (5 *( ((self.qc/self.atm.Pression())+1)**0.2857 -1 ) )**0.5

        self.q = (self.atm.rho * (self.Vv)**2 ) / 2 # Pression dynamique
        if self.typeSpd in [0,1,2]:
            self.p_t = self.atm.Pression()*(1+0.2*self.M**2)**3.5 #Pression totale 
        else :
            self.p_t = self.q + self.atm.Pression() # Pression totale pour vitesse calibrée
        self.u = (0.3125e-7*(self.atm.Temperature()**1.5)) / (self.atm.Temperature()+120) # Viscosite dynamique
        self.re = self.atm.rho * (self.Vv) * self.long_re / self.u  # reynolds
        
        self.cl = self.poids*self.nz / (0.5*self.atm.rho*(self.Vv)**2*self.surf_a)
        self.Tt = self.atm.Temperature() * (1 + 0.2*self.M**2)

    def VitesseVraie(self):
        return self.Vv

    def VitesseEq(self):
        return self.Ve

    def VitesseCal(self):
        return self.Vc

    def Mach(self):
        return self.M

    def PressionDynamique(self):
        return self.q
    
    def PressionTotale(self):
        return self.p_t

    def PressionImpact(self):
        return self.qc

    def ViscositeDynamique(self):
        return self.u

    def Reynolds(self):
        return self.re

    def Portance(self):
        return self.cl

    def TemperatureTotale(self):
        return self.Tt

    def ReturnAll(self):
        return [self.a, self.Mach(), self.VitesseVraie(), self.VitesseEq(), self.VitesseCal(), \
            self.PressionTotale(),
            self.PressionDynamique(),
            self.PressionImpact(),
            self.TemperatureTotale(),
            self.ViscositeDynamique(),
            self.Reynolds(),
            self.Portance()]

