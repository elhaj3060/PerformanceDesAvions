import numpy as np
import matplotlib
import math
import sys
from Atmosphere import Atmosphere
import matplotlib.pyplot as plt

##############atmosphere
print("Choix du modèle d'atmosphere : \n")
Hp = float(input("Altitude pression Hp (pi) : "))
UsrtempType = int(input("\nType de température en entrée ( 0 - Temperature / 1 - Delta ISA) : ") )
temp = float(input("\nValeur de la température (C) : "))
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
print(f"\n\nRapport de température : {results[0]:.4G}\n"
    f"Rapport de pression : {results[1]:.4G} \n"
    f"Rapport de densité : {results[2]:.4G}\n"
    f"Temperature : {results[3]:.4G} K / {results[3]-273.15:.4G} C\n"
    f"Déviation par rapport à ISA : {results[4]:.4G} \n"
    f"Pression : {results[5]:.4G} lb/ft^2\n"
    f"Densité : {results[6]:.4G} slugs/ft^3\n")

#poids = float(input("Poids (lb): "))
Piste_disponible=float(input('longueur de piste: '))
Piste_disponible-=200 #Enlever RAD
poids= np.arange(start=33500, stop= 53020, step=10)

############
Clmax=1.85 ##config flap 20 et GU
S=520
a0=661.48
V_mca=95.0  #KCAS 
V_mcl=92.0  #KCAS
V_mcg=90.0  #KCAS 
V_1mcg=95.0 #KCAS 
Kasyma=0.10
deltaCddwm=0.0030
g=32.174
pente=0
RAD=200 
kemax=16.7*1e6
tiremax=210 #mph

twdv=   [0.000, 0.010 , 0.020 , 0.100 , 0.120  ,0.200 , 0.400, 0.60]
dvrvl=  [1.80 ,  2.30 , 2.75  , 6.50  , 7.35 , 10.80 , 19.80 , 30.0]
dvlo35= [0.00 ,  0.00 ,  1.00 ,  9.00 , 11.00 , 16.70 , 31.00 ,45.0]

twdt   = [ 0.0  ,  0.02 , 0.03 , 0.05 , 0.065 , 0.075, 0.10 , 0.18,   0.2 ,  0.4 , 0.6]
dtvrvl = [ 2.35 ,  2.25 , 2.19 , 2.09 , 2.01  ,1.96  ,1.83  ,1.41 , 1.30  ,1.30 , 1.30]
dtvlo35= [17.00 , 10.00 , 6.60 , 5.50 , 5.10  , 4.90  ,4.50  ,3.80 ,3.80 , 3.80 , 3.80]

V1_min=V_1mcg/math.sqrt(modele_atmosphere.sigma)*1.6878 #TAS en ft/s
#print("V1 min est: ", V1_min)
vecDisMinW=[]
VectVlof=[]
Vectgrad=[]
print("Calcul en cours.....")
for i in range(0,len(poids)):
    cond1=1.1*V_mca/math.sqrt(modele_atmosphere.sigma)*1.6878  #v2>1.1*Vmca
    cond2=1.13*math.sqrt(295.37*poids[i]/(Clmax*S))/math.sqrt(modele_atmosphere.sigma)*1.6878 #v2>1.13Vsr
    V2=max(cond1,cond2) #TAS ft/s
    V2_KEAS=V2*math.sqrt(modele_atmosphere.sigma)/1.6878
    VR=V1_min #supposition pour lancer la boucle
    while VR<1.05*V_mca/math.sqrt(modele_atmosphere.sigma)*1.6878 or VR<V_1mcg/math.sqrt(modele_atmosphere.sigma)*1.6878:
        ############### Calcul gradient pour AEO et OEI ###############
        M=(V2/1.6878)/(a0*math.sqrt(modele_atmosphere.theta))
        #print(M)
        if (tempType == "delISA" and temp > 15):
            MTOFN=((8775-0.1915*Hp)-(8505-0.195*Hp)*M)*(1-(temp-15)/100)
        else:
            MTOFN=(8775-0.1915*Hp)-(8505-0.195*Hp)*M
        #print(MTOFN)
        #### Cd pour AEO
        Cl_V2=(295.37*poids[i])/(V2_KEAS**2*S)
        
        Cd_AEO=0.0465+0.0334*Cl_V2**2 #bas sur flap 20 et gd
        #print(Cd_AEO)
        ### Cd pour OEI
        q=V2_KEAS**2/295.37
        #print(q)
        CT=MTOFN/(q*S)
        #print(CT)
        deltaCdcntl=Kasyma*CT**2
        #print(deltaCdcntl)
        Cd_OEI=Cd_AEO+deltaCddwm+deltaCdcntl
        #print(Cd_OEI)
        ##### calcul grad ####
        grad_AEO= 2*MTOFN/poids[i] - Cd_AEO/Cl_V2
        grad_OEI=MTOFN/poids[i] - Cd_OEI/Cl_V2
        #increment 1 est V35- Vlof OEI
        #increment 2 est Vlof-Vr   OEI
        #increment 3 est Vlof-Vr   AEO
        #increment 4 est V35- Vlof AEO
        inc1=np.interp(grad_OEI,twdv,dvlo35)
        inc2=np.interp(grad_OEI,twdv,dvrvl)
        inc3=np.interp(grad_AEO,twdv,dvrvl)
        inc4=np.interp(grad_AEO,twdv,dvlo35)

        VLOF_OEI=V2-inc1
        VR=VLOF_OEI-inc2
        VLOF_AEO=VR+inc3
        V35_AEO=VLOF_AEO+inc4
        V1_max=VR

        if VR<1.05*V_mca/math.sqrt(modele_atmosphere.sigma)*1.6878 or VR<V_1mcg/math.sqrt(modele_atmosphere.sigma)*1.6878:
            print("Warning! VR ne respecte pas la condition (processus itératif en cours)")
        V2 += 0.05
        V2_KEAS=V2*math.sqrt(modele_atmosphere.sigma)/1.6878
    # print('Cl est a v2', Cl_V2) 
    #print('gradient AEO est', grad_AEO)
    #print('gradient OEI est', grad_OEI)

    # print("V1 max est: ", V1_max)
    # print("V2 est:", V2)
    # print("VLOF OEI est:",VLOF_OEI)
    # print("VR est:",VR)
    #print("VLOF AEO est:",VLOF_AEO)
    # print("V35 AEO est:",V35_AEO)	
    ## break avec cette condition et affiche rien sinon affiche les Vi et continue le calcul
    if grad_OEI<0.024:
       print('Le gradient calculé en OEI est inférieur à 2.4% requis par FAR 25.121 (b)')
       continue 
    #else:
    ### incr/ments de temps ##
    #increment t1 est  delta temps Vlof-Vr   OEI
    #increment t2 est delta temps  V35- Vlof OEI
    #increment t3 est  delta temps Vlof-Vr   AEO
    #increment t4 est  delta temps V35- Vlof AEO
    inct1=np.interp(grad_OEI,twdt,dtvrvl)
    #print("Delta t (VLOF-VR) OEI",inct1)
    inct2=np.interp(grad_OEI,twdt,dtvlo35)
    #print("Delta t (V35-VLOF) OEI",inct2)

    inct3=np.interp(grad_AEO,twdt,dtvrvl)
    #print("Delta t (VLOF-VR) AEO",inct3)
    inct4=np.interp(grad_AEO,twdt,dtvlo35)
    #print("Delta t (V35-VLOF) AEO",inct4)

    #### calcul dstances 
    dvlovr_OEI=inct1*(VLOF_OEI+VR)/2
    dv35vlo_OEI=inct2*(V2+ VLOF_OEI)/2
    dvlovr_AEO=inct3*(VLOF_AEO+VR)/2
    dv35vlo_AEO=inct4*(V35_AEO+VLOF_AEO)/2
    # print('la Distance VLOF - VR en OEI est: ',dvlovr_OEI)
    # print('la Distance V35 - VLOF en OEI est: ',dv35vlo_OEI)
    # print('la Distance VLOF - VR en AEO est: ',dvlovr_AEO)
    # print('la Distance V35 - VLOF en AEO est: ',dv35vlo_AEO)

    #rapportv = float(input("Rapport V1/VR: "))
    #V1=rapportv*VR

    # if V1<V_1mcg/math.sqrt(modele_atmosphere.sigma)*1.6878:
    #     print('Warning! le rapport spécifié résulte en une valeur de V1 qui est inférieure à V1MCG')
    #     print("V1 changé à V1MCG ")
    #     V1=V_1mcg/math.sqrt(modele_atmosphere.sigma)*1.6878
    #     #VR=V1/rapportv
    #     print("La nouvelle vitesse V1 est:",V1)
    #     ######## recalcule prop parce que VR a changé
    #     ##### Boucle while jusque à V atteint la valeur de VR voulue
    #     while VR<=V1/rapportv:
    #         ############### Calcul gradient pour AEO et OEI ###############
    #         M=(V2/1.6878)/(a0*math.sqrt(modele_atmosphere.theta))
    #         #print(M)
    #         if (tempType == "delISA" and temp > 15):
    #             MTOFN=((8775-0.1915*Hp)-(8505-0.195*Hp)*M)*(1-(temp-15)/100)
    #         else:
    #             MTOFN=(8775-0.1915*Hp)-(8505-0.195*Hp)*M
    #         Cl_V2=(295.37*poids)/(V2_KEAS**2*S)
    #         Cd_AEO=0.0465+0.0334*Cl_V2**2 #bas sur flap 20 et gu
    #         q=V2_KEAS**2/295.37
    #         CT=MTOFN/(q*S)
    #         deltaCdcntl=Kasyma*CT**2
    #         Cd_OEI=Cd_AEO+deltaCddwm+deltaCdcntl
    #         ##### calcul grad ####
    #         grad_AEO= 2*MTOFN/poids - Cd_AEO/Cl_V2
    #         grad_OEI=MTOFN/poids - Cd_OEI/Cl_V2
    #         #increment 1 est V35- Vlof OEI
    #         #increment 2 est Vlof-Vr   OEI
    #         #increment 3 est Vlof-Vr   AEO
    #         #increment 4 est V35- Vlof AEO
    #         inc1=np.interp(grad_OEI,twdv,dvlo35)
    #         inc2=np.interp(grad_OEI,twdv,dvrvl)
    #         inc3=np.interp(grad_AEO,twdv,dvrvl)
    #         inc4=np.interp(grad_AEO,twdv,dvlo35)

    #         VLOF_OEI=V2-inc1
    #         VR=VLOF_OEI-inc2
    #         VLOF_AEO=VR+inc3
    #         V35_AEO=VLOF_AEO+inc4
    #         V1_max=VR
    #         V2 += 0.2
    #         V2_KEAS=V2*math.sqrt(modele_atmosphere.sigma)/1.6878
    #     print('gradient AEO est', grad_AEO)
    #     print('gradient OEI est', grad_OEI)
    #     print("V1 max est: ", V1_max)
    #     print("V2 est:", V2)
    #     print("VLOF OEI est:",VLOF_OEI)
    #     print("VR est:",VR)
    #     print("VLOF AEO est:",VLOF_AEO)
    #     print("V35 AEO est:",V35_AEO)
    #     ## break avec cette condition et affiche rien sinon affiche les Vi et continue le calcul
    #     if grad_OEI<0.024:
    #         print('Le gradient calculé en OEI est inférieur à 2.4% requis par FAR 25.121 (b)')
    #     else:
    #         ### incr/ments de temps ##
    #         #increment t1 est  delta temps Vlof-Vr   OEI
    #         #increment t2 est delta temps  V35- Vlof OEI
    #         #increment t3 est  delta temps Vlof-Vr   AEO
    #         #increment t4 est  delta temps V35- Vlof AEO
    #         inct1=np.interp(grad_OEI,twdt,dtvrvl)
    #         print("Delta t (VLOF-VR) OEI",inct1)
    #         inct2=np.interp(grad_OEI,twdt,dtvlo35)
    #         print("Delta t (V35-VLOF) OEI",inct2)

    #         inct3=np.interp(grad_AEO,twdt,dtvrvl)
    #         print("Delta t (VLOF-VR) AEO",inct3)
    #         inct4=np.interp(grad_AEO,twdt,dtvlo35)
    #         print("Delta t (V35-VLOF) AEO",inct4)

    #         #### calcul dstances 
    #         dvlovr_OEI=inct1*(VLOF_OEI+VR)/2
    #         dv35vlo_OEI=inct2*(V2+ VLOF_OEI)/2
    #         dvlovr_AEO=inct3*(VLOF_AEO+VR)/2
    #         dv35vlo_AEO=inct4*(V35_AEO+VLOF_AEO)/2
    #         print('la Distance VLOF - VR en OEI est: ',dvlovr_OEI)
    #         print('la Distance V35 - VLOF en OEI est: ',dv35vlo_OEI)
    #         print('la Distance VLOF - VR en AEO est: ',dvlovr_AEO)
    #         print('la Distance V35 - VLOF en AEO est: ',dv35vlo_AEO)

    #     print("La nouvelle vitesse V1 est:",V1)
    #     print("La nouvelle vitesse VR est:",VR)
    #     print("La nouvelle vitesse V2 est:",V2)

    #intervalle 1: V0 à VR AEO
    #intervalle 1: V0 à V1 OEI
    #intervalle 1: V1 à VR OEI
    #intervalle 1: V1 à V0 AEO
    #intervalle 1: VR à V0 AEO
    #V1=0.9*VR
    vecDisMinV=[]
    vecE=[]
    for V1 in range(round(V1_min+0.5), int(round(VR+0.5))):
        vectV=[0, VR, 0, V1, V1, VR, V1, 0, VR, 0]
        #print(vectV)
        vectD=[]
        for K in range(5):
            #print('Calcul segment numéro:',i+1)
            Vi=vectV[2*K]
            Vf=vectV[2*K+1]
            V_RMS=(math.sqrt(2)/2)*math.sqrt(Vi**2+Vf**2)
            #print("la vitesse RMS est:", V_RMS)
            M_RMS=V_RMS/(1.6878*a0*math.sqrt(modele_atmosphere.theta))
            #print("M RMS est: ",M_RMS)
            if (tempType == "delISA" and temp > 15):
                MTOFN=((8775-0.1915*Hp)-(8505-0.195*Hp)*M_RMS)*(1-(temp-15)/100)
            else:
                MTOFN=(8775-0.1915*Hp)-(8505-0.195*Hp)*M_RMS
            #print(MTOFN)
            if K in [0,1]:
                T=MTOFN*2
            elif K==2:
                T=MTOFN
            elif K in [3,4]:
                T=(600-1000*M_RMS)*2
            #print('T est : ', T)
            if K in [0,1]:
                Cl=0.829
                Cd=0.0750
                mu=0.02
            elif K==2:
                Cl=0.829
                Cd=0.0750+0.0030+0.0020 # add control and windmilling drag
                mu=0.02
            elif K in [3,4]:
                Cl=0.209
                Cd=0.1171
                mu=0.4
            #print('Cl est : ', Cl)
            #print('Cd est : ', Cd)
            #print('mu est : ', mu)
            q_RMS=(V_RMS*math.sqrt(modele_atmosphere.sigma)/1.6878)**2/295.37
            #print('q RMS est:',q_RMS)
            a=(g/poids[i])*((T-mu*poids[i])-(Cd-mu*Cl)*q_RMS*S-poids[i]*pente)
            #print('a est:', a)
            dt=(Vf-Vi)/a
            #print('le delta t est:',dt)
            ds=dt*(Vf+Vi)/2
            #print('la distance est:',ds)
            vectD.append(ds)
            if K==3:
                ASD=2*V1
                vectD.append(ASD)
                delta_brake= mu*(poids[i] - Cl*q_RMS*S)*ds
            elif K==4:
                ASD=2*VR
                vectD.append(ASD)
        
        FTOD=1.15*(vectD[0]+ dvlovr_AEO+dv35vlo_AEO)
        #print('La distance FTOD AEO est:', FTOD)
        TOD_OEI= vectD[1]+vectD[2]+dvlovr_OEI+dv35vlo_OEI
        #print('La distance TOD OEI (V1) est:', TOD_OEI)
        #TOD_OEI2=vectD[0]+dvlovr_OEI+dv35vlo_OEI
        #print('La distance TOD OEI (VR) est:', TOD_OEI2)
        ASD_AEO=vectD[1]+vectD[3]+vectD[4] #l'ASD aeo la plus limtante
        #print('La distance ASD AEO (V1) est:', ASD_AEO)
        #ASD_AEO2=vectD[0]+vectD[5]+vectD[6]
        #print('La distance ASD AEO (VR) est:', ASD_AEO2)
        
        if delta_brake/4 > kemax or VLOF_AEO> tiremax*1.6878:
            print('l''Énergie dépasse le maximum permis sur un frein, cette valeur de poids ne sera pas considéré!' )
            print('la vitesse dépasse le maximum permis sur le pneus, cette valeur ne sera pas cosidérer!')
        else:            
            L_min=max(FTOD,TOD_OEI,ASD_AEO)
            #print('La longueur minimale de piste requise est:',L_min)
            #print("Énergie totale absorbée par les freins (millions ft-lb): ",delta_brake/1e6)
            #print("La nouvelle vitesse V1 est:",V1)
            #print("La nouvelle vitesse VR est:",VR)
            vecDisMinV.append(L_min)
            vecE.append(delta_brake/4)

    for j in range(0,len(vecDisMinV)):
        if vecDisMinV[j]<Piste_disponible+5 and vecDisMinV[j]> Piste_disponible-5:
            vecDisMinW.append(poids[i])
            energie_max=vecE[j]
            VectVlof.append(VLOF_AEO)
            Vectgrad.append(grad_OEI)


    

#print(vecDisMinW)
poids_max=max(vecDisMinW)
print('Le gradient minimal OEI est (%) : ',Vectgrad[-1]*100)
print('L''énergie maximale sur chaque frein est (ft-lb): ',energie_max)
print('La vitesse maximale des pneus AEO est (mph): ', VectVlof[-1] /1.6878)
print('Le poids maximal pour la distance donnée (lb): ',poids_max)


# plt.plot(vecDisMin,poids)
# plt.title("Longueur minimale de la piste en fonction du poids")
# plt.ylabel("poids de l'avion (lb)")
# plt.xlabel("Longueur minimale de la piste (pi)")
# plt.show()