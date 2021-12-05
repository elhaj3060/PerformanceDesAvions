import math

class Atmosphere():
    """Classe pour representer les donnees de l'atmosphere """

    def __init__(self,Hp,temp,typeTemp):
        self.temp0 = 15
        self.p0 = 2116.22
        self.rho0 = 0.002377
        self.lmbda = 0.0019812

        self.Hp = Hp
        self.temp = temp
        self.typeTemp = typeTemp

        self.theta = 0
        self.delta = 0
        self.sigma = 0
        self.delISA = 0 
        self.P = 0
        self.rho = 0

        if self.Hp < -2000 or self.Hp > 50000:
            raise Exception("Altitude non acceptée") 

    def SolveTroposphere(self):

        if self.typeTemp == "delISA":         # 1 pour delISA comme entree
            self.tempISA = self.temp0 - self.lmbda * self.Hp
            self.delISA = self.temp
            self.temp = self.delISA + self.tempISA

            
        self.theta = (self.temp + 273.15) / (self.temp0 + 273.15)
        self.delta = (1-6.87535e-6*self.Hp)**5.2559
        self.P = self.p0*self.delta
        
        self.sigma = self.delta / self.theta
        self.rho = self.sigma*self.rho0
        self.tempISA =  self.temp0 - self.lmbda*self.Hp
        self.delISA = self.temp - self.tempISA

    def SolveStratosphere(self):

        if self.typeTemp == "delISA":         # 1 pour delISA comme entree
            self.tempISA = -56.5
            self.delISA = self.temp
            self.temp = self.delISA + self.tempISA

        self.theta = (self.temp +273.15) / (self.temp0 + 273.15)
        self.delta = 0.22336 * math.exp(-(self.Hp-36089)/20806)
        self.sigma = self.delta/self.theta
        self.rho = self.rho0 * self.sigma
        self.P = self.delta * self.p0
        self.tempISA = -56.5
        self.delISA = self.temp - self.tempISA

        

    def RapportTemp(self):
        return self.theta

    def RapportPression(self):
        return self.delta

    def RapportDensite(self):
        return self.sigma
    
    def Temperature(self):
        return self.temp + 273.15
    
    def DeltaISA(self):
        return self.delISA
    
    def Pression(self):
        return self.P

    def Densite(self):
        return self.rho

    def ReturnAll(self):
        return [self.RapportTemp(), self.RapportPression(), self.RapportDensite(), self.Temperature(),
                self.DeltaISA(), self.Pression(), self.Densite()]

if __name__ == "__main__":
    atmosphere = Atmosphere(36089,-60,"T")
    atmosphere.SolveTroposphere()
    allAtmosphereResults = atmosphere.ReturnAll()
    print(allAtmosphereResults)


