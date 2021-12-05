import numpy as np
import math
import sys
from Atmosphere import Atmosphere
from Vitesse import Vitesse
from Force import Force
from Montee import Montee


class Croisiere_mission():
    def __init__(self,Hp,delISA,VWind,Wf_montee, Wi_descente,choix_vitesse_croisiere,Mach=0):
        self.Mach = Mach
        self.W =( Wf_montee + Wi_descente)/2
        self.delISA = delISA
        self.Hp = Hp
        self.choix_vitesse_croisiere = choix_vitesse_croisiere
        self.VWind = VWind

    def solve(self):
        M_mo=0.85       #vitesses limites
        V_mo=330*1.688  #vitesses limites
        ROC_min=300     #ft/min 
        #initialisation 
        M_md=0.3
        D_last=100000
        V_wind = self.VWind
        choix_vitesse = self.choix_vitesse_croisiere

        Hp = self.Hp

        UsrtempType = 1
        temp = self.delISA
        poids = self.W
        if UsrtempType == 0:
            tempType = "T"
        elif UsrtempType == 1:
            tempType = "delISA"
            modele_atmosphere = Atmosphere(Hp,temp,tempType)

        if -2000 <= Hp <= 36089:
            modele_atmosphere.SolveTroposphere()
        elif 36089 < Hp <= 50000:
            modele_atmosphere.SolveStratosphere()

        results = modele_atmosphere.ReturnAll()



        choix_flaps = 0

        choix_facteur_de_charge = 0
        if choix_facteur_de_charge == 0:
            Nz = 1

            
        v = Vitesse(modele_atmosphere, poids, M_md, 0, 520, 8.286, 1, 0)
        v.solve()
        v_values = v.ReturnAll()

        #lien module vitesse et determiner les props avec M=M_md
        #données nécéssaires D, MCR thrust, Vc, ROC_p

        choix_regime_moteur = 5

        position_cg = 25

        f = Force(v,choix_regime_moteur,position_cg)

        f.solve()
        f_values = f.ReturnAll()


        choix_type_montee = 0

        m = Montee(f,v,modele_atmosphere,choix_type_montee)
        m.solve()
        Montee_values = m.ReturnAll()

        if choix_vitesse==0:
            while m.ROC> ROC_min and f.D < f.T and M_md<M_mo and v.Vc<V_mo and f.D<D_last: 
                D_last=f.D
                M_md=M_md+0.0001

            #module vitesse avec input le nouveau M_md: Vc,D, MCR, ROC_p,
                v = Vitesse(modele_atmosphere, poids, M_md, 0, 520, 8.286, 1, 0)
                v.solve()
                v_values = v.ReturnAll()

                f = Force(v,choix_regime_moteur,position_cg)
                f.solve()

                m = Montee(f,v,modele_atmosphere,choix_type_montee)
                m.solve()
                Montee_values = m.ReturnAll()

            V_sol=v.Vv/1.6878+V_wind  # en kts

            SFC=0.58+(0.035*Hp/10000)
            W_f=SFC*f.D #Drag est le thrust requis par les deux moteurs

            SAR=(v.Vv/1.6878)/W_f
            SR=SAR*(V_sol/(v.Vv/1.6878))
            print(f"V_MD calibrée : {v.Vc/1.6878:.5G} Kts\n"
            f"V_MD Sol : {V_sol:.5G} Kts\n"
            f"Mach MD : {M_md:.5G} \n"
            f"SAR MD: {SAR:.5G} nm/lb\n"
            f"SR MD: {SR:.5G} nm/lb\n"
            f"W_f MD: {W_f:.5G} lb/hr\n")

            self.V_cal = v.Vc/1.6878
            self.V_sol = V_sol
            self.Mach = M_md
            self.SAR = SAR
            self.SR = SR
            self.W_f = W_f

        elif choix_vitesse==1:
            
            M = self.Mach # ajouter les autres options de vitess/ lien main #imporatnt!!!!!!!!!!!!!!!!!???
            v = Vitesse(modele_atmosphere, poids, M, 0, 520, 8.286, 1, 0)
            v.solve()

            f = Force(v,choix_regime_moteur,position_cg)
            f.solve()

            V_sol=v.Vv/1.6878+V_wind  # en kts
            SFC=0.58+(0.035*Hp/10000)
            W_f=SFC*f.D #Drag est le thrust requis par les deux moteurs

            SAR=(v.Vv/1.6878)/W_f
            SR=SAR*(V_sol/(v.Vv/1.6878))
            print(f"V calibrée : {v.Vc/1.6878:.5G} Kts\n"
            f"V  Sol : {V_sol:.5G} Kts\n"
            f"Mach : {M:.5G} \n"
            f"SAR : {SAR:.5G} nm/lb\n"
            f"SR : {SR:.5G} nm/lb\n"
            f"W_f : {W_f:.5G} lb/hr\n")

            self.V_cal = v.Vc/1.6878
            self.V_sol = V_sol
            self.Mach = M
            self.SAR = SAR
            self.SR = SR
            self.W_f = W_f


        elif choix_vitesse==2:
            
            M=np.arange(start=0.001, stop= M_mo, step=0.0001)
            SAR_i=[0]
            #M_SARmax=[]
            SFC=0.58+(0.035*Hp/10000)
            #Boucle pour trouver SAR MAx at MCR et M correspondant à SAR Max
            for M_i in M:
                v = Vitesse(modele_atmosphere, poids, M_i, 0, 520, 8.286, 1, 0)
                v.solve()

                f = Force(v,choix_regime_moteur,position_cg)
                f.solve()

                W_f=SFC*f.D #Drag est le thrust requis par les deux moteurs

                SAR=(v.Vv/1.6878)/W_f
        
                if SAR >= SAR_i[-1]:
                    SAR_max=SAR
                    M_SARmax=M_i
                    #print(SAR_max)
                
                SAR_i.append(SAR)


            LRC=0.99*SAR_max
            #print(LRC)
            #Boucle pour M_LRC: 1ère itération avec le M SAR Max jusqu'à M_Mo pour trouver M_LRC
            list_sar=[]
            for M_i in np.arange(start=M_SARmax, stop= M_mo, step=0.00001):

                v = Vitesse(modele_atmosphere, poids, M_i, 0, 520, 8.286, 1, 0)
                v.solve()

                f = Force(v,choix_regime_moteur,position_cg)
                f.solve()

                W_f=SFC*f.D #Drag est le thrust requis par les deux moteurs

                SAR=(v.Vv/1.6878)/W_f
        
                #print(SAR)
                if LRC-0.0000001 <= SAR <= LRC+0.0000001:
                    M_LRC=M_i
                    #print(SAR)
                    #print(M_LRC)
            
            v = Vitesse(modele_atmosphere, poids, M_LRC, 0, 520, 8.286, 1, 0)
            v.solve()

            f = Force(v,choix_regime_moteur,position_cg)
            f.solve()

            V_sol=v.Vv/1.6878+V_wind  # en kts
            W_f=SFC*f.D #Drag est le thrust requis par les deux moteurs
            #SAR=LRC  ### V.v/Wf
            SR=LRC*(V_sol/(v.Vv/1.6878))
            print(f"V calibrée LRC : {v.Vc/1.6878:.5G} Kts\n"
            f"V  Sol LRC: {V_sol:.5G} Kts\n"
            f"Mach LRC: {M_LRC:.5G} \n"
            f"SAR LRC: {LRC:.5G} nm/lb\n"
            f"SR LRC: {SR:.5G} nm/lb\n"
            f"W_f LRC: {W_f:.5G} lb/hr\n")


            self.V_cal = v.Vc/1.6878
            self.V_sol = V_sol
            self.Mach = M_LRC
            self.SAR = LRC
            self.SR = SR
            self.W_f = W_f


if __name__ == "__main__":

    croisiere  = Croisiere_mission(15000,0,0,50000,30000,0)
    croisiere.solve()