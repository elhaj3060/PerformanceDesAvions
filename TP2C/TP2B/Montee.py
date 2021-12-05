import math

class Montee: 
    def __init__(self,Force,Vitesse,Atm,montee_type):
        self.force = Force
        self.vitesse = Vitesse
        self.atm = Atm
        self.montee_type = montee_type
        
        #montee_type : 0 ----> M constant 
        #montee_type : 1 ----> CAS constant

    def solve(self):
        phi_m=(1/(0.7*self.vitesse.Mach()**2))*((1+0.2*self.vitesse.Mach()**2)**3.5-1)/((1+0.2*self.vitesse.Mach()**2)**2.5)
    #conditions standards 
        if self.atm.delISA == 0:
            if self.montee_type == 0:
                if self.atm.Hp<36089: 
                #equation Af
                    self.AF=-0.133184*(self.vitesse.Mach()**2)
                elif self.atm.Hp>36089: 
                    #equation Af
                    self.AF=0
                
            elif self.montee_type == 1:
                if self.atm.Hp<36089: 
                #equation Af
                    self.AF=0.7*(self.vitesse.Mach()**2)*(phi_m-0.190263)
                elif self.atm.Hp>36089: 
                    #equation Af
                    self.AF=0.7*(self.vitesse.Mach()**2) * phi_m

        elif self.atm.delISA != 0:

            if self.montee_type == 0:
                if self.atm.Hp<36089: 
                    self.AF=-0.133184*(self.vitesse.Mach()**2)*((self.atm.tempISA+273.15)/(self.atm.temp+273.15))
                elif self.atm.Hp>36089: 
                    self.AF=0

            elif self.montee_type == 1:
                if self.atm.Hp<36089: 
                    self.AF=0.7*(self.vitesse.Mach()**2)*(phi_m-0.190263*((self.atm.tempISA+273.15)/(self.atm.temp+273.15)))
                elif self.atm.Hp>36089: 
                    self.AF=0.7*(self.vitesse.Mach()**2) * phi_m


        self.ROC = ((self.vitesse.VitesseVraie()/1.6878)*101.2686*(self.force.T-self.force.D)/self.vitesse.poids)*(1/(1+self.AF))
        #---------------------------------
        self.ROCp = self.ROC*((self.atm.tempISA+273.15)/(self.atm.temp+273.15))
         #---------------------------------
        self.gradient_angle = self.ROC / ((self.vitesse.VitesseVraie()/1.6878)*101.2686)
         #---------------------------------
        self.gradient = self.gradient_angle*100
         #---------------------------------
        self.gamma_acc = math.atan(self.gradient_angle)
         #---------------------------------
        self.acc = ((self.force.T/self.vitesse.poids) - (self.force.Cd/self.vitesse.Portance())) - self.gradient_angle
        #conditions non standards 
        #elif self.atm.delISA != 0:

        
        
    def ReturnAll(self):
        return [self.AF,self.ROC, self.ROCp,self.gradient,self.acc]