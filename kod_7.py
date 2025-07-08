import pandas as pd
import re

# Wczytanie danych wejsciowych (Ścieżka do pliku BtToxin_digger po filtracji)
df = pd.read_excel("int.xlsx")

df['Protein_sequence'] = (
    df['Protein_sequence']
    .fillna('')
    .astype(str)
    .str.strip()
    .str.replace(r'\s+', '', regex=True)
    .str.upper()
)
df['Identity'] = pd.to_numeric(df['Identity'], errors='coerce')

def assign_ranks(row):
    name = row['New_Hit_id']
    identity = row.get('Identity', 0)
    if re.match(r'^[A-Za-z]+\d+[A-Z][a-z]$', name):
        return 4 if identity >= 95 else 3
    if re.match(r'^[A-Za-z]+\d+[A-Z]like\d+$', name):
        return 3
    if re.match(r'^[A-Za-z]+\d+like\d+$', name):
        return 2
    if re.match(r'^[A-Za-z]+like\d+$', name):
        return 1
    if re.match(r'^[A-Za-z]+\d+[A-Z]$', name):
        return 2
    if re.match(r'^[A-Za-z]+\d+$', name):
        return 1
    if re.match(r'^[A-Za-z]+$', name):
        return 'mnemonic'
    return 1

df['rank'] = df.apply(assign_ranks, axis=1)
df['toxin'] = df['New_Hit_id']
df['mnemonic'] = df['New_Hit_id'].str.replace('like', '', regex=False)\
                                 .str.extract(r'^([A-Za-z]+)', expand=False)

total_toxins = df['New_Hit_id'].nunique()
total_strains = df['Organism Qualifier'].nunique()
total_genomes = df['Accession'].nunique()

summary = df.groupby(['mnemonic', 'rank', 'toxin']).agg(
    frequency=('New_Hit_id', 'nunique'),
    strain_count=('Organism Qualifier', lambda x: x.dropna().nunique()),
    genome_count=('Accession', lambda x: x.dropna().nunique()),
    Protein_sequence=('Protein_sequence', lambda x: '; '.join(sorted(set(x))))
).reset_index()

summary['frequency_pct'] = (summary['frequency'] / total_toxins * 100).round(2).astype(str) + '%'
summary['strain_pct'] = (summary['strain_count'] / total_strains * 100).round(2).astype(str) + '%'
summary['genome_pct'] = (summary['genome_count'] / total_genomes * 100).round(2).astype(str) + '%'

mnemonic_totals = df.groupby('mnemonic').agg(
    frequency=('New_Hit_id', 'nunique'),
    strain_count=('Organism Qualifier', lambda x: x.dropna().nunique()),
    genome_count=('Accession', lambda x: x.dropna().nunique())
).reset_index()

mnemonic_totals['frequency_pct'] = (mnemonic_totals['frequency'] / total_toxins * 100).round(2).astype(str) + '%'
mnemonic_totals['strain_pct'] = (mnemonic_totals['strain_count'] / total_strains * 100).round(2).astype(str) + '%'
mnemonic_totals['genome_pct'] = (mnemonic_totals['genome_count'] / total_genomes * 100).round(2).astype(str) + '%'
mnemonic_totals['rank'] = 'mnemonic'
mnemonic_totals['toxin'] = mnemonic_totals['mnemonic']
mnemonic_totals['Protein_sequence'] = '-'
mnemonic_totals['rank_label'] = 'Protein family'
mnemonic_totals = mnemonic_totals[summary.columns.tolist()]

def rank_label(r):
    if r == 'mnemonic':
        return 'Protein family'
    elif r == 1:
        return 'Novel rank 1'
    elif r == 2:
        return 'Novel rank 2'
    elif r == 3:
        return 'Novel rank 3'
    elif r == 4:
        return 'Known rank 4'
    else:
        return 'Other'

summary['rank_label'] = summary['rank'].apply(rank_label)

combined = pd.concat([summary, mnemonic_totals], ignore_index=True)
combined['sort_order'] = combined['rank'].apply(lambda x: -1 if x == 'mnemonic' else int(x))
combined = combined.sort_values(by=['mnemonic', 'sort_order', 'toxin'])\
                   .drop(columns='sort_order').reset_index(drop=True)

# Dodanie aliasów toksyn (np. Cry19A → Cry19, Cry19A)

def get_aliases(toxin):
    toxin = str(toxin)
    aliases = []
    match = re.match(r"^([A-Za-z]+)(\d+)([A-Za-z]*)$", toxin)
    if match:
        base, number, letters = match.groups()
        if number:
            aliases.append(f"{base}{number}")
        if letters:
            for i in range(1, len(letters) + 1):
                aliases.append(f"{base}{number}{letters[:i]}")
    return aliases

def get_like_prefix(toxin):
    toxin = str(toxin)
    match = re.match(r"^([A-Za-z]+\d*[A-Za-z]*?(like|Alike))\d*$", toxin)
    if match:
        return match.group(1)
    return None

def get_sort_key(toxin):
    match_main = re.fullmatch(r"[A-Za-z]+", toxin)
    if match_main:
        return (0, toxin, 0, '', 0, 0)
    match = re.match(r"^([A-Za-z]+)(\d+)([A-Za-z]*)(like|Alike)?(\d*)$", toxin)
    if match:
        prefix, number, letters, like_type, like_num = match.groups()
        number = int(number)
        like_num = int(like_num) if like_num else 0
        like_order = 1 if like_type else 0
        return (1, prefix, number, letters, like_order, like_num)
    return (2, toxin, 0, '', 0, 0)

new_rows = []
for mnemonic, group in combined.groupby('mnemonic'):
    seen_aliases = set()
    group = group.sort_values(by='toxin', key=lambda col: col.map(get_sort_key))

    if not ((group['toxin'] == mnemonic).any()):
        new_rows.append({
            'mnemonic': mnemonic,
            'rank_label': '',
            'toxin': mnemonic,
            **{col: '' for col in combined.columns if col not in ['mnemonic', 'rank_label', 'toxin']}
        })
        seen_aliases.add(mnemonic)

    for _, row in group.iterrows():
        toxin = row['toxin']

        like_prefix = get_like_prefix(toxin)
        if like_prefix and like_prefix not in seen_aliases:
            new_rows.append({
                'mnemonic': mnemonic,
                'rank_label': '',
                'toxin': like_prefix,
                **{col: '' for col in combined.columns if col not in ['mnemonic', 'rank_label', 'toxin']}
            })
            seen_aliases.add(like_prefix)

        for alias in get_aliases(toxin)[:-1]:
            if alias not in seen_aliases:
                new_rows.append({
                    'mnemonic': mnemonic,
                    'rank_label': '',
                    'toxin': alias,
                    **{col: '' for col in combined.columns if col not in ['mnemonic', 'rank_label', 'toxin']}
                })
                seen_aliases.add(alias)

        new_rows.append(row.to_dict())

# Finalna tabela z aliasami i toksynami
final_df = pd.DataFrame(new_rows)
final_df = final_df.sort_values(by=['mnemonic', 'toxin'], key=lambda col: col.map(get_sort_key)).reset_index(drop=True)

# Tabela objaśnień
description_data = {
    'Kolumna': [
        'frequency_pct', 'strain_pct', 'genome_pct',
        'frequency', 'strain_count', 'genome_count', 'Protein_sequence'
    ],
    'Opis': [
        '(unikalne toksyny / wszystkie toksyny) * 100',
        '(szczepy z toksyną / wszystkie szczepy) * 100',
        '(genomy z toksyną / wszystkie genomy) * 100',
        'Liczba unikalnych nazw toksyn',
        'Liczba unikalnych szczepów z daną toksyną',
        'Liczba unikalnych genomów z daną toksyną',
        'Lista unikalnych sekwencji białkowych (oddzielona średnikiem)'
    ]
}
description_df = pd.DataFrame(description_data)

# Eksport do Excela
output_path = "analiza_unikatowych_toksyn_i_czestosci_ich_wystepowania.xlsx"
with pd.ExcelWriter(output_path, engine='xlsxwriter') as writer:
    final_df.to_excel(writer, index=False, sheet_name='Hierarchia toksyn')
    description_df.to_excel(writer, index=False, sheet_name='Objaśnienia')


