import numpy as np
import math
import sys

from Atmosphere import Atmosphere
from Vitesse import Vitesse
from Force import Force
from Montee import Montee


class Montee_mission():
    def __init__(self,M_cst,Vc_post_transition1,Hpi,Hpf,Vwind,Wi,deltaISA,demo=0):
        self.demo=demo
        self.deltaISA = deltaISA
        self.M_cst = M_cst
        self.Vc_post_trans1 = Vc_post_transition1
        self.Vc_pre_trans1 = 250
        self.Hpi = Hpi
        self.Hpf = Hpf
        self.Vwind = Vwind*1.6878
        self.Wi = Wi
        self.deltaISA = deltaISA
        self.gravite = 32.174

    def solveAltitudeTransition1(self):
        self.Hp_trans1 = 0
        
        delta_ISA= 10
        M=0
        Hp_array=np.arange(start=10000, stop= 51010, step=1)
        #print(Hp)
        self.P0 = 2116.22
        T0K = 288.15
        rho0 = 0.002377
        a0 = 661.48
        i=1

        while  self.M_cst-M > 0.00000001: 
            self.Hp_trans1=Hp_array[i]
            if self.Hp_trans1 <= 36089: 
                delta=(1-6.87535e-6*self.Hp_trans1)**(5.2559)
                T_std=T0K-0.0019812*self.Hp_trans1
                TK=T_std+delta_ISA
            elif self.Hp_trans1 > 36089:
                delta=0.22336*math.exp((-self.Hp_trans1+36089)/20806)
                T_std=T0K-0.0019812*self.Hp_trans1
                TK=T_std+delta_ISA


            P=delta*self.P0
            theta=TK/T0K
            sigma=delta/theta

            #vitesse de son

            qc = self.P0*((1 + 0.2*(self.Vc_post_trans1/a0)**2)**3.5 -1)
            M = math.sqrt(5*((qc/P + 1)**0.2857-1))

            i=i+1

    def solvePasHpMontee(self):
        
        self.Pas_Hp = []
        self.Pas_Hp.append(self.Hpi)
        if self.Hpi%1000 != 0:
            self.Hpi += self.Hpi%1000
            self.Pas_Hp.append(self.Hpi)
        Hp = self.Hpi
        while Hp + 500 < self.Hpf:
            # Add transition altitude
            if Hp <= self.Hp_trans1 <= Hp+500:
                self.Pas_Hp.append(self.Hp_trans1)
            #Add 10000 altitude
            if Hp < 10000 < Hp + 500:
                if Hp != 10000 and (Hp+500) != 10000:
                    self.Pas_Hp.append(10000)
            Hp += 500
            self.Pas_Hp.append(Hp)
        self.Pas_Hp.append(self.Hpf)


    def solvePasHpDescente(self):
        self.Pas_Hp = []
        self.Pas_Hp.append(self.Hpi)
        if self.Hpi%1000 != 0:
            self.Hpi -= self.Hpi%1000
            self.Pas_Hp.append(self.Hpi)
        Hp = self.Hpi
        while Hp - 500 > self.Hpf:
            # Add transition altitude
            if Hp-500 <= self.Hp_trans1 <= Hp:
                self.Pas_Hp.append(self.Hp_trans1)
            #Add 10000 altitude
            if Hp - 500 < 10000 < Hp :
                if Hp != 10000 and (Hp-500) != 10000:
                    self.Pas_Hp.append(10000)
            Hp -= 500
            self.Pas_Hp.append(Hp)
        self.Pas_Hp.append(self.Hpf)



    def solve(self):
        #valeur initiale qui sera verifié iterativement jusqu'a obeir au critere err < 20   
        self.Wmoy = self.Wi
        self.Wfinal = self.Wi
        self.time = 0
        self.dist = 0
        self.fuel = 0
        self.ListAltitudeCroisiere = []


        for i_altitude in range(0,len(self.Pas_Hp)-1):
            if self.Pas_Hp[i_altitude+1] == 10000:
                self.Winit10000 = wTemp2
            self.Iterate(i_altitude)
            wTemp1 = self.Wfinal
            wTemp2 = wTemp1 - self.delta_fuel
            Wmoy2 = (wTemp1 + wTemp2)/2

            while(abs(Wmoy2-self.Wmoy)> 0.0001):
                self.Wmoy = Wmoy2
                self.Iterate(i_altitude)
                wTemp2 = wTemp1-self.delta_fuel
                Wmoy2 = (wTemp1 + wTemp2)/2
            if self.Hpf > self.Hpi:
                if self.roc < 100 or (self.roc < 300 and self.demo==0):
                    pass
                else:
                    self.time += (self.delta_t + self.deltaTAcc/60)
                    self.dist += self.deltaDistAcc + self.delta_dist
                    self.fuel += self.delta_fuel + self.deltaFuelAcc
                    self.Wfinal -= (self.delta_fuel + self.deltaFuelAcc)

            else:
                self.time += (self.delta_t + self.deltaTAcc/60)
                self.dist += self.deltaDistAcc + self.delta_dist
                self.fuel += self.delta_fuel + self.deltaFuelAcc
                self.Wfinal -= (self.delta_fuel + self.deltaFuelAcc)

            if self.Hpf > self.Hpi: #montee
                if self.roc< 300 and self.Pas_Hp[i_altitude+1] <= 41000:
                    self.ListAltitudeCroisiere.append(self.Pas_Hp[i_altitude])

                    
                if self.roc < 100:
                    print(f"ROC plus petit que 100ft/s")
                    break
      

        print("Valeurs à 10000ft")
        print(f"acc : {self.acc10000}")
        print(f"Vtrue_avg : {self.Vv10000_moy }")
        print(f"delta t : {self.deltaT10000}")
        print(f"delta dist : {self.deltaDist10000}")
        print(f"delta fuel : {self.deltaFuel10000}")
        print(f"W initial : {self.Winit10000}")
        print("-------------------------------")
        print(f"Temps total (s) : {self.time*60}")
        print(f"Distance Total (NM) : {self.dist}")
        print(f"Carburant total consommé (lb) : {self.fuel}")
        print(f"altitude de transition : {self.Hp_trans1}")
        if self.Hpf > self.Hpi: #montee  
            if len(self.ListAltitudeCroisiere) ==0:
                pass
            else:
                if min(self.ListAltitudeCroisiere)%1000 != 0:
                    self.altitudeCroisiere = min(self.ListAltitudeCroisiere) -  min(self.ListAltitudeCroisiere)%1000
                else:
                    self.altitudeCroisiere = min(self.ListAltitudeCroisiere)
                print(f"Altitude de croisiere : {self.altitudeCroisiere}")

        
        
           



    def Iterate(self,i_altitude):
        Hp1 = self.Pas_Hp[i_altitude]
        Hp2 = self.Pas_Hp[i_altitude+1]

        Hp_moy = (Hp2 + Hp1)/2
        delta_Hp = Hp2-Hp1

        atm = Atmosphere(Hp_moy,self.deltaISA,"delISA")

        if Hp_moy <= 36089 :
            atm.SolveTroposphere()
        else:
            atm.SolveStratosphere()
        
        if Hp2 < 10000:

            v = Vitesse(atm, self.Wmoy, self.Vc_pre_trans1, 3,520, 8.286, 1,0)
            v.solve()


            if Hp1 > Hp2:
                f = Force(v,7,25)
            else:
                f = Force(v,4,25)
            f.solve()

            
            m = Montee(f,v,atm,1)
            m.solve()

        elif 10000 <= Hp2 < self.Hp_trans1:

            if Hp2 == 10000:
                v = Vitesse(atm, self.Wmoy,(self.Vc_post_trans1+self.Vc_pre_trans1)/2, 3,520, 8.286, 1,0)
            else:
                v = Vitesse(atm, self.Wmoy,self.Vc_post_trans1, 3,520, 8.286, 1,0)
            v.solve()

            if Hp1 > Hp2:
                f = Force(v,7,25)
            else:
                f = Force(v,4,25)
            f.solve()




            m = Montee(f,v,atm,1)
            m.solve()
            
        elif Hp2 >= self.Hp_trans1:

            v = Vitesse(atm, self.Wmoy, self.M_cst, 0 ,520, 8.286, 1,0)
            v.solve()

            if Hp1 > Hp2:
                f = Force(v,7,25)
            else:
                f = Force(v,4,25)
            f.solve()

            m = Montee(f,v,atm,0)
            m.solve()

        if Hp2 == 10000:
            
            atm10000 =  Atmosphere(10000,self.deltaISA,"delISA")
            atm10000.SolveTroposphere()
            v250 = Vitesse(atm10000, self.Wmoy, 250, 3 ,520, 8.286, 1,0)
            v250.solve()
            v275 = Vitesse(atm10000, self.Wmoy, 275, 3 ,520, 8.286, 1,0)
            v275.solve()
            v262 = Vitesse(atm10000, self.Wmoy, 262.5, 3 ,520, 8.286, 1,0)
            v262.solve() 
            
            if Hp1 > Hp2:
                f10000 = Force(v262,7,25)
            else:
                f10000 = Force(v262,4,25)
            f10000.solve()
            self.acc10000 = self.gravite*(f10000.T-f10000.D)/self.Wmoy
            self.deltaT10000 = abs((v250.VitesseVraie()-v275.VitesseVraie())/self.acc10000)
            self.deltaDist10000 = abs((self.Vwind + v.VitesseVraie())*self.deltaT10000)/6076
            if f10000.T>0:
                self.deltaFuel10000 = (f10000.sfc*f10000.T/3600)*self.deltaT10000
            else:
                self.deltaFuel10000 = (f10000.sfc*1200/3600)*self.deltaT10000
            if Hp1 > Hp2:
                self.Vv10000_moy =  (v.VitesseVraie() + v250.VitesseVraie())/2
            else:
                self.Vv10000_moy = v.VitesseVraie()

            self.deltaTAcc = self.deltaT10000
            self.deltaDistAcc = self.deltaDist10000
            self.deltaFuelAcc =  self.deltaFuel10000

        else:
            self.deltaTAcc = 0
            self.deltaDistAcc = 0
            self.deltaFuelAcc = 0


        delta_Hgeo = delta_Hp*(atm.temp+273.15)/(atm.tempISA+273.15)
        self.delta_t = delta_Hgeo/m.ROC
        self.delta_dist = (self.delta_t*((self.Vwind + v.VitesseVraie()))*60)/6076
        if f.T > 0:
            self.delta_fuel = (f.sfc*f.T/60)*self.delta_t
        else:
            self.delta_fuel = (f.sfc*1200/60)*self.delta_t
        self.roc = m.ROC


if __name__ == "__main__":

    M_constant = 0.755
    VKCAS = 277
    Hpi = 1500
    Hpf = 41000
    VWIND = 0
    Wi = 48525
    T_ISA = 0

    m_mission = Montee_mission(M_constant, VKCAS, Hpi, Hpf, VWIND, Wi, T_ISA)
    m_mission.solveAltitudeTransition1()
    if Hpi > Hpf:
        m_mission.solvePasHpDescente()  
        m_mission.solve()
    else:
        m_mission.solvePasHpMontee()
        m_mission.solve()
