import math
import scipy.interpolate


class Force:
    def __init__(self,vitesse,regime_moteur,pourcentage_cg,pousse_moteur=0.0):
        self.pousse_moteur = pousse_moteur
        self.vitesse = vitesse
        self.regime_moteur = regime_moteur
        #regime_moteur : 0 ----> MTO AEO

        #regime_moteur : 1 ----> MTO OEI

        #regime_moteur : 2 ----> GA AEO 

        #regime_moteur : 3 ----> GA OEI

        #regime_moteur : 4 ---->MCL AEO

        #regime_moteur : 5 ---->MCR AEO

        #regime_moteur : 6 ---->MCT OEI 

        #regime_moteur : 7 ---->Idle AEO

        #regime_moteur : 8 ---->Idle OEI
        #regime_moteur : 9 ---->Poussee totale

        self.pourcentage_cg = pourcentage_cg

        self.CdpDict = {
            "0" : 0.0206,
            "1" : 0.0206,
            "2" : 0.0465,
            "3" : 0.0465,
            "4" : 0.1386
        }

        self.KDict = {
            "0" : 0.0364,
            "1" : 0.0364,
            "2" : 0.0334,
            "3" : 0.0334,
            "4" : 0.0301
        }

        self.AoAswDict = {
            "0" : 14.7,
            "1" : 14.7,
            "2" : 14.6,
            "3" : 14.6,
            "4" : 14.4

        }
        self.MachBuffetList = [0.2750,0.3000,0.3250,0.3500,0.3750,0.4000,0.4250,0.4500,
        0.4750,0.5000,0.5250,0.5500,0.5750,0.6000,0.6250,0.6500,
        0.6750,0.7000,0.7250,0.7500,0.7750,0.8000,0.8250,0.8500,
        0.8750,0.9000]
        self.clBuffetList = [1.3424,1.3199,1.2974,1.2667,1.2310,1.1930,1.1551,1.1191,
        1.0863,1.0577,1.0337,1.0142,0.9989,0.9868,0.9764,0.9659,
        0.9530,0.9349,0.9085,0.8698,0.8149,0.7391,0.6373,0.5039,
        0.3330,0.118]

    def solve(self):
        
        self.Cdp = self.CdpDict[str(self.vitesse.typeFlaps)]
        self.K = self.KDict[str(self.vitesse.typeFlaps)]
        self.Cd = self.Cdp + self.K*self.vitesse.Portance()**2
        #------------------

        self.L = self.vitesse.poids*self.vitesse.nz
        
        #---------------------------
        if (self.vitesse.typeFlaps == 0 or self.vitesse.typeFlaps == 1):
            if 0 <= self.vitesse.Mach() <= 0.60:
                self.Cdcomp = 0
            elif 0.60 < self.vitesse.Mach() <= 0.78:
                self.Cdcomp = (0.0508-0.1748*self.vitesse.Mach()+0.1504*self.vitesse.Mach()**2)*self.vitesse.Portance()**2
            elif 0.78 < self.vitesse.Mach() <= 0.85:
                self.Cdcomp = (-99.3434+380.888*self.vitesse.Mach()-486.8*self.vitesse.Mach()**2
                +207.408*self.vitesse.Mach()**3)*self.vitesse.Portance()**2
        else:
            self.Cdcomp = 0
        # ------------------------------------
        if self.regime_moteur in [1,3,6,8]:
            self.Cdwm = 0.0030
        else:
            self.Cdwm = 0

        self.Dwm = self.Cdwm*self.vitesse.surf_a*self.vitesse.PressionDynamique()
        #-------------------------------------
        if self.pousse_moteur != 0.0:
            self.T = self.pousse_moteur
        else:
            if self.regime_moteur in [0,1,2,3,6,10]:
                if self.regime_moteur in [1,3,6,10]:
                    self.T = (8775-0.1915*self.vitesse.atm.Hp)-(8505-0.195*self.vitesse.atm.Hp)*self.vitesse.Mach()
                    if (self.vitesse.atm.typeTemp == "delISA" and self.vitesse.atm.DeltaISA() > 15):
                        self.T = self.T*(1-(self.vitesse.atm.DeltaISA()-15)/100)
                    if self.regime_moteur in [6]:
                        self.T = self.T*0.9
                    elif self.regime_moteur in [10]:
                        self.T = 2*0.9*self.T
                else:
                    self.T = 2*((8775-0.1915*self.vitesse.atm.Hp)-(8505-0.195*self.vitesse.atm.Hp)*self.vitesse.Mach())
                    if (self.vitesse.atm.typeTemp == "delISA" and self.vitesse.atm.DeltaISA() > 15):
                        self.T = self.T*(1-(self.vitesse.atm.DeltaISA()-15)/100)
            if self.regime_moteur in [4,5]:
                self.T = (5690-0.0968*self.vitesse.atm.Hp) - (1813-0.0333*self.vitesse.atm.Hp)*self.vitesse.Mach()
                if (self.vitesse.atm.typeTemp == "delISA" and self.vitesse.atm.DeltaISA() > 10):
                    self.T = self.T*(1-(self.vitesse.atm.DeltaISA()-10)/100)
                if self.regime_moteur in [4]:
                    self.T = self.T*2
                elif self.regime_moteur in [5]:
                    self.T = self.T*0.98*2
            if self.regime_moteur in [7,8]:
                self.T = 600-1000*self.vitesse.Mach()
                if self.regime_moteur in [7]:
                    self.T = self.T*2
        #-------------------------------------
        if self.regime_moteur in [1,3,6,8]:
            self.kasyma = 0.10
        else:
            self.kasyma = 0

        self.Cdcntl = self.kasyma*(self.T/(self.vitesse.surf_a*self.vitesse.PressionDynamique()))**2
        self.Dcntl = self.Cdcntl*self.vitesse.PressionDynamique()*self.vitesse.surf_a
        #------------------------------------
        
        self.Cd += self.Cdcntl + self.Cdwm + self.Cdcomp
        self.cdi = self.Cd - self.Cdp - self.Cdcomp - self.Cdcntl - self.Cdwm

        if self.vitesse.typeFlaps in [1,3,4]:
            self.Cd += 0.02

        self.D = self.Cd*self.vitesse.PressionDynamique()*self.vitesse.surf_a

        

        self.rapportLD = self.L/self.D 
        #------------------------------------
       
        #-------------------------------------
        self.clAt9 = self.vitesse.Portance()*(1+(self.vitesse.long_re/40.56)*(0.09-self.pourcentage_cg/100))
        if self.vitesse.typeFlaps in [0,2]:
            if self.vitesse.typeFlaps == 0:
                self.AoAF = (self.clAt9-0.05)/0.10
            elif self.vitesse.typeFlaps == 2:
                self.AoAF = (self.clAt9-0.25)/0.10
        else:
            if self.vitesse.typeFlaps == 1:
                self.AoAF = self.clAt9/0.10
            elif self.vitesse.typeFlaps == 3:
                self.AoAF = (self.clAt9-0.20)/0.10
            elif self.vitesse.typeFlaps == 4:
                self.AoAF = (self.clAt9-0.50)/0.10

        self.AoAsw = self.AoAswDict[str(self.vitesse.typeFlaps)]

        if self.vitesse.typeFlaps in [0,2]:
            if self.vitesse.typeFlaps == 0:
                self.clAt9_pour_nzsw = 0.05+0.10*self.AoAsw
            elif self.vitesse.typeFlaps == 2:
                self.clAt9_pour_nzsw = 0.25+0.10*self.AoAsw
        else:
            if self.vitesse.typeFlaps == 1:
                self.clAt9_pour_nzsw = 0.10*self.AoAsw
            elif self.vitesse.typeFlaps == 3:
                self.clAt9_pour_nzsw = 0.20+0.10*self.AoAsw
            elif self.vitesse.typeFlaps == 4:
                self.clAt9_pour_nzsw = 0.50+0.10*self.AoAsw

        self.cl_pour_nzsw = self.clAt9_pour_nzsw*(1+(self.vitesse.long_re/40.56)*((self.pourcentage_cg/100)-0.09))
        self.nzsw = (self.cl_pour_nzsw*self.vitesse.PressionDynamique()*self.vitesse.surf_a)/self.vitesse.poids
        #self.phisw = math.degrees(math.acos(1/self.nzsw))

        if self.vitesse.typeFlaps in [2,3,4]:
            self.nzbuffet = "N/A"
        else:
            y_interp = scipy.interpolate.interp1d(self.MachBuffetList, self.clBuffetList, fill_value="extrapolate")
            self.clBuffetAt9 = y_interp(self.vitesse.Mach())
            self.clBuffet = self.clBuffetAt9*(1+(self.vitesse.long_re/40.56)*((self.pourcentage_cg/100)-0.09))
            self.nzbuffet = self.clBuffet/self.vitesse.Portance()

        self.sfc =  0.58 +  (0.035 *self.vitesse.atm.Hp/ 10000) 

        
    def ReturnAll(self):
        return [self.vitesse.Portance(),self.L, self.Cd,self.D,self.rapportLD,
        self.Cdp,self.cdi,self.Cdcomp,self.Cdwm,self.Dwm,self.Cdcntl,self.Dcntl,self.T,self.AoAF,
        self.nzsw,self.nzbuffet]


        


        
        
            