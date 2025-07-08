import pandas as pd
import re

# Wczytaj plik wynikowy z kodu 7
df_rank = pd.read_excel("plik_wynikowy_z_kodu_7.xlsx", sheet_name=0)

# Dodaj brakującą kolumnę 'mnemonic' – wyodrębnioną z 'toxin'
df_rank['mnemonic'] = df_rank['toxin'].astype(str).str.replace('like', '', regex=False)\
                                .str.extract(r'^([A-Za-z]+)', expand=False)

# Klasyfikacja rank_label
def classify_rank(label):
    if 'Known rank 3' in label:
        return 'known'
    elif 'Novel rank 1' in label:
        return 'novel rank1'
    elif 'Novel rank 2' in label:
        return 'novel rank2'
    elif 'Novel rank 3' in label:
        return 'novel rank3'
    else:
        return 'other'

df_rank['rank_group'] = df_rank['rank_label'].apply(classify_rank)
df_filtered = df_rank[df_rank['rank_group'] != 'other']

# Wczytaj oryginalne dane (te przefiltrowane z BtToxin_digger)
df_orig = pd.read_excel("input.xlsx")
df_orig['Organism Qualifier'] = df_orig['Organism Qualifier'].fillna('').astype(str)
df_orig['Accession'] = df_orig['Accession'].fillna('').astype(str)
df_orig['New_Hit_id'] = df_orig['New_Hit_id'].fillna('').astype(str)
df_orig['mnemonic'] = df_orig['New_Hit_id'].str.replace('like','', regex=False)\
                                         .str.extract(r'^([A-Za-z]+)', expand=False)

# Ustalenie zbiorów unikalnych szczepów i genomów

def get_union_sets(toxin_label, mnemonic_val):
    if toxin_label == mnemonic_val:
        subset = df_orig[df_orig['mnemonic'] == toxin_label]
    else:
        subset = df_orig[df_orig['New_Hit_id'].str.fullmatch(re.escape(toxin_label), na=False)]
        if subset.empty:
            subset = df_orig[df_orig['New_Hit_id'].str.startswith(toxin_label)]
    strain_set = set(subset['Organism Qualifier'].unique())
    genome_set = set(subset['Accession'].unique())
    strain_set.discard('')
    genome_set.discard('')
    return strain_set, genome_set

def get_unique_frequency(toxin_label, mnemonic_val):
    if toxin_label == mnemonic_val:
        subset = df_orig[df_orig['mnemonic'] == toxin_label]
    else:
        subset = df_orig[df_orig['New_Hit_id'].str.fullmatch(re.escape(toxin_label), na=False)]
        if subset.empty:
            subset = df_orig[df_orig['New_Hit_id'].str.startswith(toxin_label)]
    return subset['New_Hit_id'].nunique()

df_filtered['strain_set'] = df_filtered.apply(lambda row: get_union_sets(row['toxin'], row['mnemonic'])[0], axis=1)
df_filtered['genome_set'] = df_filtered.apply(lambda row: get_union_sets(row['toxin'], row['mnemonic'])[1], axis=1)
df_filtered['frequency_recalc'] = df_filtered.apply(lambda row: get_unique_frequency(row['toxin'], row['mnemonic']), axis=1)

# Grupowanie

grouped = df_filtered.groupby(['mnemonic', 'rank_group']).agg({
    'frequency_recalc': 'sum',
    'strain_set': lambda sets: set().union(*sets),
    'genome_set': lambda sets: set().union(*sets)
}).reset_index()

grouped['strain_count_union'] = grouped['strain_set'].apply(len)
grouped['genome_count_union'] = grouped['genome_set'].apply(len)

# Tworzenie pivotów i zapis

summary_frequency = grouped.pivot(index='mnemonic', columns='rank_group', values='frequency_recalc').fillna(0)
summary_strain = grouped.pivot(index='mnemonic', columns='rank_group', values='strain_count_union').fillna(0)
summary_genome = grouped.pivot(index='mnemonic', columns='rank_group', values='genome_count_union').fillna(0)

for summary in [summary_frequency, summary_strain, summary_genome]:
    for col in ['known', 'novel rank1', 'novel rank2', 'novel rank3']:
        if col not in summary.columns:
            summary[col] = 0
    summary = summary[['known', 'novel rank1', 'novel rank2', 'novel rank3']]

# Zapisz do excela
with pd.ExcelWriter("podsumowanie_toksyn_szczepów_genomów.xlsx") as writer:
    summary_frequency.to_excel(writer, sheet_name='toxins')
    summary_strain.to_excel(writer, sheet_name='strain_count')
    summary_genome.to_excel(writer, sheet_name='genome_count')


