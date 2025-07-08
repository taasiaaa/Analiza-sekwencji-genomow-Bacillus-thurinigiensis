import pandas as pd

# Wczytanie pliku wejściowego (plik wynikowy BtToxin_digger w formacie xlsx, po filtrowaniu z kolumną zawierającą nowe nazwy toksyn (w tym wypadku 'New_Hit_id')
int = pd.read_excel("int.xlsx")
data = pd.read_excel("Data.xlsx")

# Genomy z Data.xlsx (Assembly Accession)
genomy_data = set(data['Assembly Accession'].dropna().unique())
print(f"Całkowita liczba genomów w Data.xlsx: {len(genomy_data)}")

# Szczepy z int.xlsx (Organism Qualifier)
szczepy_int = set(int['Organism Qualifier'].dropna().unique())
print(f"Całkowita liczba szczepów w int.xlsx: {len(szczepy_int)}")

# Genomy z int.xlsx (Accession) — to są te z toksynami
genomy_int = set(int['Accession'].dropna().unique())

# Genomy bez toksyn — występują w Data.xlsx, ale nie w int.xlsx
genomy_bez_toksyn = genomy_data - genomy_int
print(f"Liczba genomów bez toksyn: {len(genomy_bez_toksyn)}")
print("Genomy bez toksyn:")
for g in sorted(genomy_bez_toksyn):
    print(g)

# Szczepy z toksynami: w int mamy też kolumnę New_Hit_id, szukamy szczepów bez toksyn (New_Hit_id pusty lub NaN)
szczepy_z_toksyn = set(int.loc[int['New_Hit_id'].notna(), 'Organism Qualifier'].unique())
szczepy_bez_toksyn = szczepy_int - szczepy_z_toksyn

print(f"\nLiczba szczepów bez toksyn: {len(szczepy_bez_toksyn)}")
print("Szczepy bez toksyn:")
for s in sorted(szczepy_bez_toksyn):
    print(s)


#Wyniki dla pracy magisterskiej:
#Całkowita liczba genomów w Data.xlsx: 1197
#Całkowita liczba szczepów w int.xlsx: 1160
#Liczba genomów bez toksyn: 1
#Genomy bez toksyn:
#GCF_028401745.1
#Liczba szczepów bez toksyn: 0
#Szczepy bez toksyn:
