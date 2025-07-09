import pandas as pd

# 1. Wczytanie pliku wejściowego (plik wynikowy BtToxin_digger w formacie xlsx, po filtrowaniu z kolumną zawierającą nowe nazwy toksyn (w tym wypadku 'New_Hit_id', z dodanymi ręcznie informacjami o zanieczyszczeniu genomów)
path = "int.xlsx"
df = pd.read_excel(path)

# 2. Usuń ewentualne duplikaty
df_unique = df.drop_duplicates(subset=['Accession', 'New_Hit_id', 'Protein_sequence'])

# 3. Zdefiniuj kolumny flagowe i odpowiadające im etykiety
flag_columns = {
    'contaminated': 'contaminated',
    'unverified source organism': 'unverified source organism',
    'genome length too large': 'genome length too large',
    'low quality sequence': 'low quality sequence'
}

# 4. Dla każdego accession utwórz listę flag (po przecinku)
def extract_flags_for_accession(df_acc):
    flags_set = set()
    for col, label in flag_columns.items():
        if 'yes' in df_acc[col].astype(str).str.lower().values:
            flags_set.add(label)
    return ', '.join(sorted(flags_set)) if flags_set else ''

# 5. Utwórz DataFrame z informacjami per accession
acc_info = (
    df_unique
    .groupby('Accession')
    .apply(extract_flags_for_accession)
    .reset_index(name='Contamination_flags')
)

# Dołącz Submission Date i level do acc_info
meta_cols = df_unique[['Accession', 'Submission Date', 'level']].drop_duplicates(subset='Accession')
acc_info = acc_info.merge(meta_cols, on='Accession', how='left').set_index('Accession')

# 6. Funkcja pomocnicza do agregacji po accession
def join_by_accession(acc_list, col):
    accs = sorted(set(acc_list))
    vals = []
    for a in accs:
        val = str(acc_info.loc[a, col]) if a in acc_info.index and pd.notna(acc_info.loc[a, col]) else ''
        if val:
            vals.append(val)
    return '; '.join(sorted(set(vals)))

# 7. Grupowanie po szczepie
grouped = (
    df_unique
    .groupby('Organism Qualifier')
    .agg(
        toxin_count    = ('New_Hit_id', 'nunique'),
        toxins         = ('New_Hit_id', lambda x: '; '.join(sorted(set(x)))),
        Accessions      = ('Accession',   lambda x: '; '.join(sorted(set(x)))),
        contaminated   = ('Accession',   lambda x: join_by_accession(x, 'Contamination_flags')),
        Submission_Date= ('Accession',   lambda x: join_by_accession(x, 'Submission Date')),
        level          = ('Accession',   lambda x: join_by_accession(x, 'level')),
    )
    .reset_index()
)

# 7.1 Statystyki ogólne dla toksyn
summary_stats = pd.DataFrame({
    'metric': ['mean', 'std', 'min', 'max', 'median'],
    'value': [
        grouped['toxin_count'].mean(),
        grouped['toxin_count'].std(),
        grouped['toxin_count'].min(),
        grouped['toxin_count'].max(),
        grouped['toxin_count'].median()
    ]
})

# 8. Zapis do Excela
output = "analiza_składu_toksyn_pestycydowych_szczepów.xlsx"
with pd.ExcelWriter(output) as writer:
    grouped.to_excel(writer, sheet_name="By Strain", index=False)
    summary_stats.to_excel(writer, sheet_name="Stats", index=False)

print(f"Zapisano wynik do: {output}")
