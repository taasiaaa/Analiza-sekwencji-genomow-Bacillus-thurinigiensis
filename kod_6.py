import pandas as pd
import re
import scipy.stats as stats
import itertools

# Wczytanie danych wejściowych
df = pd.read_excel("plik_z_kodu_5.xlsx")

# Uproszczenie nazw typu 'like'
def simplify_toxin_name(name):
    return re.sub(r'like\d+', 'like', name)

df['simplified_toxins'] = df['toxins'].fillna('').apply(
    lambda x: [simplify_toxin_name(tok.strip()) for tok in x.split(';') if tok.strip()]
)

# Lista unikalnych toksyn
all_toxins = sorted(set(tok for toks in df['simplified_toxins'] for tok in toks))

# Tworzenie binarnej macierzy obecności toksyn
rows = []
for idx, row in df.iterrows():
    strain = row['Organism Qualifier']
    toxin_presence = {toxin: 1 if toxin in row['simplified_toxins'] else 0 for toxin in all_toxins}
    toxin_presence['Organism Qualifier'] = strain
    rows.append(toxin_presence)

binary_matrix = pd.DataFrame(rows)

# Grupowanie duplikatów i agregacja maksymalna
binary_matrix = binary_matrix.groupby('Organism Qualifier', as_index=False).max()

# Zapisz macierz obecności toksyn
binary_matrix.to_excel("matrix.xlsx", index=False)

# Analiza współwystępowania toksyn

# Organism Qualifier jako indeks
df_matrix = binary_matrix.set_index('Organism Qualifier')

results = []
for t1, t2 in itertools.combinations(df_matrix.columns, 2):
    n = len(df_matrix)
    x1 = df_matrix[t1].sum()
    x2 = df_matrix[t2].sum()
    both = ((df_matrix[t1] == 1) & (df_matrix[t2] == 1)).sum()
    p = (x1 / n) * (x2 / n)

    pval = stats.binomtest(both, n=n, p=p, alternative='greater').pvalue

    results.append({
        'Toxin 1': t1,
        'Toxin 2': t2,
        'Observed': both,
        'Expected_prob': p,
        'empirical_prob': both / n,
        'p-value': pval
    })

res_df = pd.DataFrame(results)
res_df = res_df.sort_values(by='p-value')

# Zapis wyników analizy korelacji
res_df.to_excel("analiza_korelacji_wspolwystepowania_toksyn.xlsx", index=False)


