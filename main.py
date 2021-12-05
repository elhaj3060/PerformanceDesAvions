from TP1C.TP1A.Atmosphere import Atmosphere
from TP1C.TP1B.Vitesse import Vitesse
from TP2C.TP2A.Force import Force
from TP2C.TP2B.Montee import Montee

import math
import numpy as np

import sys
sys.path.insert(0,"./TP3C/TP3B")
sys.path.insert(0,"./TP3C/TP3A")
import Croisiere_mission
import Montee_mission

# ------------------------- choix de l'atmosphere -----------------------------
UsrSoftwareType = int(input("\nChoix du programme à rouler \n \
    0 - TP1/TP2 \n \
    1 - TP3 \n \
    : ") )
if UsrSoftwareType == 0:
    print("\n\n")
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


    # ------------------------- choix de la vitesse  -----------------------------

    print("\n\n")
    print("Choix du modèle de vitesse : \n")
    poids = float(input("Poids (lb): "))

    UsrspdType = int(input("\nType de vitesse en entrée \n \
        0 - Mach \n \
        1 - Vitesse vraie \n \
        2 - Vitesse equivalente \n \
        3 - Vitesse calibree \n \
        4 - Rapport de vitesse (V/Vsr) \n \
        : ") )

    if UsrspdType==4:
        spd = float(input("\nValeur du rapport  : "))
    else:
        spd = float(input("\nValeur de la vitesse (kts) : "))

    choix_surf_a = input("\nChoix de surface alaire autre que la valeur par defaut (520pi^2) ? (oui/non) : ")
    if choix_surf_a == "oui":
        surf_a = float(input("\n Valeur de la surface alaire : "))
    else:
        surf_a = 520

    choix_long_re = input("\nChoix de longueur de reference pour le nombre de Reynolds autre que la valeur par defaut (8.2862pi) ? (oui/non):")
    if choix_long_re == "oui":
        long_re = float(input("\n Valeur de la longueur de reference : "))
    else:
        long_re = 8.286

    choix_facteur_de_charge = int(input("\nType d'entrée pour le facteur de charge:\n \
        0 - Valeur de Nz \n \
        1 - angle en degrées \n \
        :  "))
    if choix_facteur_de_charge == 0:
        Nz = float(input("Valeur du facteur de charge (Nz) :  "))
    elif choix_facteur_de_charge == 1:
        angle_facteur = float(input("Angle pour le calcul du facteur de charge (Deg) : "))
        Nz = 1/math.cos(math.radians(angle_facteur))


    choix_flaps = int(input("\nParametres du vol :  \n \
        0 - Flap 00 / GU \n \
        1 - Flap 00 / GD \n \
        2 - Flap 20 / GU \n \
        3 - Flap 20 / GD \n \
        4 - Flap 45 / GD \n \
        : ") )

    v = Vitesse(modele_atmosphere, poids, spd, UsrspdType,surf_a, long_re, Nz,choix_flaps)
    v.solve()
    v_values = v.ReturnAll()
    print(f"\n\nVitesse du son : {v_values[0]/1.6878:.4G} kts / {v_values[0]:.4G} ft/s\n"
        f"Vitesse Mach : {v_values[1]:.4G} \n"
        f"Vitesse Vrai: {v_values[2]/1.6878:.4G} kts / {v_values[2]:.4G} ft/s \n"
        f"Vitesse Equivalente: {v_values[3]/1.6878:.4G} kts / {v_values[3]:.4G} ft/s\n"
        f"Vitesse Calibrée: {v_values[4]/1.6878:.4G} kts / {v_values[4]:.4G} ft/s\n"
        f"Pression Totale : {v_values[5]:.4G} lb/ft^2\n"
        f"Pression Dynamique : {v_values[6]:.4G} lb/ft^2\n"
        f"Pression d'impact : {v_values[7]:.4G} lb/ft^2\n"
        f"Temperature Totale : {v_values[8]:.4G} k/ {v_values[8]-273.15:.4G} C\n"
        f"Viscosité dynamique : {v_values[9]:.4G} lb-sec/ft^2\n"
        f"Nombre de Reynolds (basée sur MAC = {long_re:.4G}) : {v_values[10]:.4G} \n"
        f"Coefficient de Portance (Basée sur Nz = 1) : {v_values[11]:.4G} \n")

    # --------------------------------  Parametres Force ----------------------------------
    choix_regime_moteur = int(input("\nChoix du régime du moteur \n \
        0 - MTO AEO \n \
        1 - MTO OEI \n \
        2 - GA AEO \n \
        3 - GA OEI \n \
        4 - MCL AEO \n \
        5 - MCR AEO \n \
        6 - MCT OEI \n \
        7 - Idle AEO \n \
        8 - Idle OEI \n \
        9 - Poussee totale \n \
        10 - MCT AEO \n \
        : ") )

    position_cg = float(input("\n Position du CG (%MAC) : "))

    if choix_regime_moteur == 9:
        pousse_totale = float(input("\nValeur de la poussée totale : "))
        f = Force(v,choix_regime_moteur,position_cg,pousse_totale)

    else:
        f = Force(v,choix_regime_moteur,position_cg)

    f.solve()
    f_values = f.ReturnAll()
    print(f"\n\nCoefficient de portance / Portance : {f_values[0]:.4G} / {f_values[1]:.5G} lbs\n"
    f"Coefficient de trainée / Trainée : {f_values[2]:.4G} / {f_values[3]:.4G} lbs\n"
    f"Finesse : {f_values[4]:.4G}\n"
    f"Trainée parasite (Cdp) : {f_values[5]:.4G}\n"
    f"Trainée induite (Cdi) :  {f_values[6]:.4G} \n"
    f"Trainée de compressibilité (Cdcomp) : {f_values[7]:.4G} \n"
    f"Trainée windmill (Cdwm/Dwm) : {f_values[8]:.4G} / {f_values[9]:.4G} lbs \n"
    f"Trainée de controle (Cdcntl / Dcntl): {f_values[10]:.4G} / {f_values[11]:.4G} lbs\n"
    f"Poussée totale : {f_values[12]:.5G} lbs \n"
    f"Angle d'attaque : {f_values[13]:.4G} deg\n"
    f"Facteur de charge à l'avertissement de décrochage (Nz sw) : {f_values[14]:.4G} \n"
    f"Angle de roulis à l'avertissement de décrochage (phi sw) : {f_values[15]:.4G} deg\n"
    f"Facteur de charge à l'initiation du buffet (Nz buffet) : {(f_values[16])} \n")

    #---------------------- Parametres Montee --------------------------------------


    choix_type_montee = int(input("\nChoix du régime du vol \n \
        0 - M constant \n \
        1 - CAS constant \n \
        : ") )

    m = Montee(f,v,modele_atmosphere,choix_type_montee)
    m.solve()
    Montee_values = m.ReturnAll()
    print(f"\n\nGradient(%) : {Montee_values[3]:.4G} %\n"
    f"Taux de montée (ROC) : {Montee_values[1]:.4G} ft/min\n"
    f"Taux de montée pression (ROCp) : {Montee_values[2]:.4G} ft/min\n"
    f"Facteur d'accélération (AF) : {Montee_values[0]:.4G}\n"
    f"Accélération selon l'axe de la trajectoire de vol :  {Montee_values[4]:.4G} g\n")

    #----------------------------- TP3 ---------------------------------------------------
    #--------------------------------------------------------------------------------------
if UsrSoftwareType == 1:
    UsrMissionType = int(input("\nType de mission à analyser \n \
    0 - Montée/Descente \n \
    1 - Croisière \n \
    : ") )
    if UsrMissionType == 0:
        M_constant = float(input("\n[1/7] Valeur de M constant :"))
        VKCAS = float(input("\n[2/7] Valeur de VKCAS (kts) : "))
        Hpi = float(input("\n[3/7] Altitude initiale (fts) : "))
        Hpf = float(input("\n[4/7] Altitude finale (fts) : "))
        VWIND = float(input("\n[5/7] Valeur de VWIND (kts) : "))
        Wi = float(input("\n[6/7] Poids initial (lb) : "))
        T_ISA = float(input("\n [7/7] T_ISA (C) : "))
        print("\n")

        m_mission = Montee_mission.Montee_mission(M_constant, VKCAS, Hpi, Hpf, VWIND, Wi, T_ISA,demo=1)
        m_mission.solveAltitudeTransition1()
        if Hpi > Hpf:
            m_mission.solvePasHpDescente()  
            m_mission.solve()
        else:
            m_mission.solvePasHpMontee()
            m_mission.solve()

        
    elif UsrMissionType == 1:
        Alt_cruise =  float(input("\n[1/5] Altitude Croisière (ft) :"))
        T_ISA = float(input("\n [2/5] T_ISA (C) : "))
        VWIND = float(input("\n[3/5] Valeur de VWIND (kts) : "))
        Wi = float(input("\n[4/5] Poids moyen (lb) : "))
        choix_vitesse_croisiere = int(input("\n[5/5] Choix du vitesse : \
            \n0 - V_MD \
            \n1 - Valeur (Mach par exemple) \
            \n2 - LRC \
            \n: "))
        if choix_vitesse_croisiere == 1:
            valeur = float(input("\nValeur de la vitesse : "))
        print("\n")
        Wf_montee = Wi+1
        wf_descente = Wi-1
        if choix_vitesse_croisiere == 1:
            croisiere  = Croisiere_mission.Croisiere_mission(Alt_cruise,T_ISA,VWIND,Wf_montee,wf_descente,choix_vitesse_croisiere,valeur)

        else:
            croisiere  = Croisiere_mission.Croisiere_mission(Alt_cruise,T_ISA,VWIND,Wf_montee,wf_descente,choix_vitesse_croisiere)
        croisiere.solve()



