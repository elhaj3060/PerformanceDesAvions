from Atmosphere import Atmosphere

print("\n\n")
print("Choix du modèle d'atmosphere : \n")
Hp = float(input("Altitude pression Hp : "))
UsrtempType = int(input("\nType de température en entrée ( 0 - Temperature / 1 - Delta ISA) : ") )
temp = float(input("\nValeur de la température : "))

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
print(f"\n\n\nRapport de température : {results[0]:.4G}\n"
      f"Rapport de pression : {results[1]:.4G} \n"
      f"Rapport de densité : {results[2]:.4G}\n"
      f"Temperature : {results[3]:.4G} \n"
      f"Déviation par rapport à ISA : {results[4]:.4G} \n"
      f"Pression : {results[5]:.4G} \n"
      f"Densité : {results[6]:.4G} \n")