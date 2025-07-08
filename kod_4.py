import pandas as pd
import matplotlib.pyplot as plt

# Wczytaj dane (wybrany plik.xlsx zawierający kolumnę Submission Date)
df = pd.read_excel('Dane.xlsx')

# Oczyszczenie daty (w wypadku złego formatu)
df['Submission_Date'] = pd.to_datetime(
    df['Submission_Date'].astype(str).str.extract(r'(\d{4}-\d{2}-\d{2})')[0],
    errors='coerce'
)
df['Year'] = df['Submission_Date'].dt.year

# Usuń rok 2025, ponieważ się jeszcze nie skończył i wyniki nie będą obiektywne
df = df[df['Year'].notna() & (df['Year'] <= 2024)]

# Zachowaj tylko unikalne genomy (Accession)
df_unique = df.drop_duplicates(subset='Accession')

# Poziomy assembly (contigs, scaffolds...)
df_unique['level'] = df_unique['level'].astype(str).str.lower().str.strip()

# Grupowanie po roku i poziomie
grouped = df_unique.groupby(['Year', 'level']).size().unstack(fill_value=0)

# Skumulowanie danych
cumulative = grouped.cumsum()

# Wykres skumulowany
cumulative.plot(kind='bar', stacked=True, colormap='Set3', figsize=(10, 6))
plt.title('Skumulowana liczba genomów Bacillus thuringiensis według roku i poziomu assembly')
plt.xlabel('Rok')
plt.ylabel('Skumulowana liczba genomów')
plt.legend(title='Poziom (level)', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()
