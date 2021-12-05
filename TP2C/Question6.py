
import sys
sys.path.insert(0, '../../')

from AER8375.TP1C.TP1A.Atmosphere import Atmosphere
from AER8375.TP1C.TP1B.Vitesse import Vitesse
from AER8375.TP2C.TP2A.Force import Force
from AER8375.TP2C.TP2B.Montee import Montee

import matplotlib.pyplot as plt


# ----------------  Question5  --------------------------
#------------------- parametres -------------------------
altitude = 8000
type_temp = "delISA"
valeur_temp = 10
#valeur_spd = 1.19
type_spd = 0
surface_alaire = 520
long_re = 8.286
Nz = 1
choix_flaps = 0
regime_moteur = 6
pos_cg = 25
type_montee = 1
W = 36000

# ------- Configuration atmosphere -----------------------=
modele_atmosphere = Atmosphere(altitude,valeur_temp,type_temp)
modele_atmosphere.SolveTroposphere()

# Liste a remplir pour afficher le graphique
liste_Mach = []
liste_gradient = []
for M in range(21,50):
    M /= 100
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
    liste_gradient.append(m.gamma_acc)


Mach_max = liste_Mach[liste_gradient.index(max(liste_gradient))]
grad_max = max(liste_gradient)


print(f"Gradient de montée maximal : {grad_max}")
print(f"Nombre de Mach correspondant : {Mach_max}")
plt.plot(Mach_max, grad_max, "o")
plt.plot(liste_Mach,liste_gradient)
plt.legend(["gradient de montée maximal"])
plt.xlabel("Mach")
plt.ylabel("Gradient de montée")
plt.title("Gradient de montée en fonction du nombre de Mach de l'avion ")

plt.show()
