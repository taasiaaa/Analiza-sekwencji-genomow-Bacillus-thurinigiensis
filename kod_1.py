import pandas as pd

# Wczytanie pliku xlsx (wynik programu BtToxin_digger - Alltoxins.txt otworzony w programie excel, w którym w kolumnie Hit_id pola 'NA' uzupełniono wartościami HMM, oraz który połączono z metadanymi z NCBI)
plik = 'int.xlsx'
df = pd.read_excel(plik)

df['Hit_id'] = df['Hit_id'].astype(str).str.strip()
df['Organism Qualifier'] = df['Organism Qualifier'].astype(str).str.strip()
df['Accession'] = df['Accession'].astype(str).str.strip()

filtered_rows = []
for accession, acc_group in df.groupby('Accession'):
    for organism, org_group in acc_group.groupby('Organism Qualifier'):
        for new_hit_id, sub_group in org_group.groupby('Hit_id'):
            if sub_group['Identity'].notna().any():
                # Wybieramy najlepszy rekord po Identity, jeśli istnieje
                najlepszy_rekord = sub_group.loc[sub_group['Identity'].idxmax()]
            elif sub_group['Protein_len/Hit_len(blast)'].notna().any():
                # Jeśli Identity to NaN, sprawdzamy Protein_len/Hit_len(blast)
                najlepszy_rekord = sub_group.loc[sub_group['Protein_len/Hit_len(blast)'].idxmax()]
            else:
                # Jeśli wszystko jest NaN, bierzemy losowy
                najlepszy_rekord = sub_group.iloc[0]
            filtered_rows.append(najlepszy_rekord)

# Tworzymy nowy DataFrame z wynikami
filtered_df = pd.DataFrame(filtered_rows)
# Zapisujemy plik do Excela
filtered_df.to_excel('plik_bez_duplikatow.xlsx', index=False)
