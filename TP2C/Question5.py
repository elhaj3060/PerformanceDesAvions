
import sys
sys.path.insert(0, '../../')

from AER8375.TP1C.TP1A.Atmosphere import Atmosphere
from AER8375.TP1C.TP1B.Vitesse import Vitesse
from AER8375.TP2C.TP2A.Force import Force
from AER8375.TP2C.TP2B.Montee import Montee
import matplotlib.pyplot as plt
import numpy as np

# ----------------  Question5  --------------------------
#------------------- parametres -------------------------
altitude = 5000
type_temp = "delISA"
valeur_temp = 32.91
valeur_spd = 1.19
type_spd = 4
surface_alaire = 520
long_re = 8.2862
Nz = 1
choix_flaps = 2
regime_moteur = 1
pos_cg = 25
type_montee = 1

# ------- Configuration atmosphere -----------------------=
modele_atmosphere = Atmosphere(altitude,valeur_temp,type_temp)
modele_atmosphere.SolveTroposphere()

# Liste a remplir pour afficher le graphique
liste_poids = []
liste_gradients = []
for W in range(30000,50000,1000):
    liste_poids.append(W)
    print(W)
    # ------------ Configuration Vitesse ---------------------
    v = Vitesse(modele_atmosphere, W, valeur_spd, type_spd,surface_alaire,long_re, Nz,choix_flaps)
    v.solve()

    #------------- Configuation Forces ---------------------
    f = Force(v,regime_moteur,pos_cg)
    f.solve()

    # ------------- Configuration Montée 
    m = Montee(f,v,modele_atmosphere,type_montee)
    m.solve()
    liste_gradients.append(m.gradient)


droite = np.polyfit(liste_poids[10:15],liste_gradients[10:15],1)
grad_cherchee = 3.1
poids_cherche = (grad_cherchee - droite[1])/droite[0]
print(f"Poids maximal (lbs) : {poids_cherche}")
plt.plot(poids_cherche,grad_cherchee,"o")
plt.plot(liste_poids,liste_gradients)

plt.legend(["gradient cible (3.1%)"])
plt.xlabel("Poids (lbs)")
plt.ylabel("Gradient (%)")

plt.title("Gradient en fonction du poids de l'avion ")


plt.show()
