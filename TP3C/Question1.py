	import sys
sys.path.insert(0,"./TP3B")
sys.path.insert(0,"./TP3A")


from Montee_mission import Montee_mission
from Croisiere_mission import Croisiere_mission


import matplotlib.pyplot as plt
import numpy as np


# Distance descente environ 90/100 NM
# //// Carburant consommé : 90 / 100 lb


# ------------------- Données ---------------------------------------
W_BAL = 36175
W_ETO = 48525
M_constant = 0.755
VKCAS = 277
Hpi_montee = 1500
Hpf_montee = 41000
VWIND = 0
T_ISA = 0
Wf_TAXI = 200
Wf_TAKEOFF = 250
Wf_Approche = 200
Wf_TAXI2 = 100


# ----------------  Calcul Montee ---------------------------------
Segment_montee = Montee_mission(M_constant, VKCAS, Hpi_montee, Hpf_montee, VWIND, W_ETO, T_ISA)

Segment_montee.solveAltitudeTransition1()

if Hpi_montee > Hpf_montee:
    Segment_montee.solvePasHpDescente()  
    Segment_montee.solve()
else:
    Segment_montee.solvePasHpMontee()
    Segment_montee.solve()
W_TOC = Segment_montee.Wfinal

Dist_montee = Segment_montee.dist
fuel_montee = Segment_montee.fuel
alt_cruise  = Segment_montee.altitudeCroisiere

# ----------------  Calcul Descente  ---------------------------------

Wi_descente = W_BAL + 400  # Poids initial de descente mis égal à W_BAL + 400

Wf_descente = 100000 # Valeur initiale arbitraire (>W_BAL)

# On cherche poids final apres la descente == W_BAL donnée
while (Wf_descente -W_BAL) > 10:
    
    segment_descente = Montee_mission(M_constant, VKCAS, alt_cruise, Hpi_montee, VWIND, Wi_descente, T_ISA)
    segment_descente.solveAltitudeTransition1()
    segment_descente.solvePasHpDescente()  
    segment_descente.solve()
    Wf_descente = segment_descente.Wfinal
    
    # On diminiue notre poids initial de descente de 5 a chaque iteration
    Wi_descente -= 5
    print(Wi_descente)

#--------------- Calcul Croisiere ----------------------------------

segment_cruise = Croisiere_mission(alt_cruise,T_ISA,VWIND,W_TOC,Wi_descente,1,M_constant)
segment_cruise.solve()

print(f"\n\n\n\n\nValeurs du TP3C :\n \
        Distance de montée (nm) : {Segment_montee.dist:.5G}\n \
        Altitude de croisière (pi) : {alt_cruise:.5G}\n \
        Distance de croisière (nm) : {segment_cruise.SAR*(W_TOC-Wi_descente)}\n \
        Distance descente (nm) : {segment_descente.dist:.5G}\n \
        Distance total (nm) : {Segment_montee.dist + segment_cruise.SAR*(W_TOC-Wi_descente) + segment_descente.dist:.5G}\n \
        Carburant total consommé (lb) : {Segment_montee.fuel + segment_descente.fuel + (W_TOC-Wi_descente) + Wf_TAXI + Wf_TAKEOFF +Wf_Approche + Wf_TAXI2:.5G} \n \
        ")