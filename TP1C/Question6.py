from TP1A.Atmosphere import Atmosphere
from TP1B.Vitesse import Vitesse
import matplotlib.pyplot as plt

Hp = 8000 # Altitude pression

tempType = 0 #temperature ambiante

temp = 273.15 # Valeur de la temperature



modele_atmosphere = Atmosphere(Hp,temp,tempType)


if -2000 <= Hp <= 36089:
    modele_atmosphere.SolveTroposphere()
elif 36089 < Hp <= 50000:
    modele_atmosphere.SolveStratosphere()

V_c_array = [x for x in range(1,400)]
Mach_array = []
for i in V_c_array:
    v = Vitesse(modele_atmosphere, 40000, i, 3,520,  8.286)
    v.solve()
    Mach_array.append(v.Mach())


V_c = 0
T_t =0
while T_t< 273.15+9.9 or T_t > 273.15+10.1:
    V_c += 1
    print(f"Vitesse calibrée : {V_c}kts")
    v = Vitesse(modele_atmosphere, 40000, V_c, 3,520,  8.286)
    v.solve()
    T_t = temp*(1+0.2*v.Mach()**2)
    print(f"Temperature totale : {T_t:.4G}K")

plt.plot(V_c_array,[temp*(1+0.2*mach**2)-273.15 for mach in Mach_array])
plt.title("Température totale en fonction de la vitesse calibrée")
plt.ylabel("Temperature totale (C)")
plt.xlabel("Vitesse calibrée (kts)")
plt.show()
