

from TP1A.Atmosphere import Atmosphere
from TP1B.Vitesse import Vitesse
import matplotlib.pyplot as plt

Hp = 3100 # Altitude pression

tempType = 0 #temperature ambiante

temp = -273.15 # Valeur de la temperature en celsius
V_c = 278
V_v=0


while V_v < V_c-0.5 or V_v > V_c + 0.5:
    temp+=(0.5) 
    modele_atmosphere = Atmosphere(Hp,temp,tempType)

    modele_atmosphere.SolveTroposphere()
    v = Vitesse(modele_atmosphere, 40000, V_c, 3,520,  8.286,1)
    v.solve()
    V_v = v.VitesseVraie()
    print(f"Temperature : {temp:.4G}C")
    print(f"Vitesse vraie : {V_v:.4G}kts")


temp_array = [x for x in range(-273,300)]
V_v_array = []

for temp in temp_array:
    modele_atmosphere = Atmosphere(Hp,temp,tempType)

    modele_atmosphere.SolveTroposphere()
    v = Vitesse(modele_atmosphere, 40000, V_c, 3,520,  8.286,1)
    v.solve()
    V_v = v.VitesseVraie()
    V_v_array.append(V_v)

plt.plot(temp_array,V_v_array)
plt.title("Vitesse vraie en fonction de la température pour une vitesse calibrée (278kts) et altitude pression (3100pi) constants")
plt.xlabel("Temperature (C)")
plt.ylabel("Vitesse vraie (kts)")
plt.show()
