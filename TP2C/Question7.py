
import sys
sys.path.insert(0, '../../')

from AER8375.TP1C.TP1A.Atmosphere import Atmosphere
from AER8375.TP1C.TP1B.Vitesse import Vitesse
from AER8375.TP2C.TP2A.Force import Force
from AER8375.TP2C.TP2B.Montee import Montee
import numpy as np
import matplotlib.pyplot as plt

# ----------------  Question5  --------------------------
#------------------- parametres -------------------------
altitude = 30000
type_temp = "delISA"
valeur_temp = 15+8+6
#valeur_spd = 1.19
type_spd = 0
surface_alaire = 520
long_re = 8.286
Nz = 1
choix_flaps = 0
regime_moteur = 5
pos_cg = 25
type_montee = 1
W = 40000

modele_atmosphere = Atmosphere(altitude,valeur_temp,type_temp)
modele_atmosphere.SolveTroposphere()

# Liste a remplir pour afficher le graphique
liste_Mach = []
liste_gradient = []
for M in range(40,85):
    M /= 100
    print(M)
    liste_Mach.append(M)
    # ------------ Configuration Vitesse ---------------------
    v = Vitesse(modele_atmosphere, W, M, type_spd,surface_alaire,long_re, Nz,choix_flaps)
    v.solve()

    #------------- Configuation Forces ---------------------
    f = Force(v,regime_moteur,pos_cg)
    f.solve()

    # ------------- Configuration Montée 
    m = Montee(f,v,modele_atmosphere,type_montee)
    m.solve()
    liste_gradient.append(m.ROC)

droite = np.polyfit(liste_Mach[35:40],liste_gradient[35:40],1)
grad_cherchee = 0
Mach_cherche = (grad_cherchee - droite[1])/droite[0]
print(f"Nombre de Mach pour un vol en palier : {Mach_cherche}")
plt.plot(Mach_cherche,grad_cherchee,"o")
plt.plot(liste_Mach,liste_gradient)
plt.legend(["Taux de montée nul"])
plt.xlabel("Mach")
plt.ylabel("Taux de montée (ft/min)")
plt.title("Taux de montée en fonction du nombre de Mach de l'avion ")
plt.show()
