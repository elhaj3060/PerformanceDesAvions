from TP1A.Atmosphere import Atmosphere
from TP1B.Vitesse import Vitesse
import matplotlib.pyplot as plt

Hp = 0

tempType = 0 #temperature ambiante

temp = 20  # arbitraire
V_c = 273
M=0

while 0.75 < M-0.001 or 0.75 > M + 0.001:
    Hp += 10
    modele_atmosphere = Atmosphere(Hp,temp,tempType)

    if -2000 <= Hp <= 36089:
        modele_atmosphere.SolveTroposphere()
    elif 36089 < Hp <= 50000:
        modele_atmosphere.SolveStratosphere()
    v = Vitesse(modele_atmosphere, 40000, V_c, 3,520,  8.286)
    v.solve()
    M = v.Mach()
    print(f"Altitude pression : {Hp:.4G}pi")
    print(f"Nombre de mach : {M:.4G}")


Hp_array = [x for x in range(0,50000,100)]
M_array = []
for Hp in Hp_array:
    modele_atmosphere = Atmosphere(Hp,temp,tempType)
    if -2000 <= Hp <= 36089:
        modele_atmosphere.SolveTroposphere()
    elif 36089 < Hp <= 50000:
        modele_atmosphere.SolveStratosphere()

    v = Vitesse(modele_atmosphere, 40000, V_c, 3,520,  8.286)
    v.solve()
    M = v.Mach()

    M_array.append(M)

plt.plot(Hp_array,M_array)
plt.title("Nombre de Mach en fonction de l'altitude pression pour une vitesse calibrée constante (273kts)")
plt.xlabel("Altitude pression (pi)")
plt.ylabel("Nombre de Mach")
plt.show()