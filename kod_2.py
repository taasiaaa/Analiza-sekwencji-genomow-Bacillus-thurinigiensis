import pandas as pd
import re
from collections import defaultdict

# Ścieżki do plików (input to wynik poprzedniego kodu (kod_1.py)
input_path = "plik_z_kodu_1.xlsx"
output_path = "nowe_nazwy_toksyn.xlsx"

# Wczytaj dane
df = pd.read_excel(input_path)

# Inicjalizuj słowniki
name_counters = defaultdict(lambda: {"counter": 0})
seq_to_name = {}

# Funkcja do parsowania bazowych składników nazw toksyn
def parse_base_parts(name):
    name = name.replace("-other", "")
    patterns = [
        r"^([A-Za-z]+?)(\d{1,2})([A-Z])([a-z])(\d+)?$",  # np. Cry1Aa2
        r"^([A-Za-z]+?)(\d{1,2})([A-Z])$",               # np. Zwa5A
        r"^([A-Za-z]+?)(\d{1,2})$",                      # np. Bmp1
        r"^([A-Za-z]+)$"                                 # np. Enhancin
    ]
    for pattern in patterns:
        match = re.match(pattern, name)
        if match:
            return match.groups() + (None,) * (5 - len(match.groups()))
    return None, None, None, None, None

# Funkcja do generowania nowych nazw toksyn
def generate_new_name(row):
    Hit_id = str(row["Hit_id"])
    nomen = str(row["Nomenclature"]).lower()
    seq = row["Protein_sequence"]

    if nomen == "nd":
        nomen = "rank1"

    if pd.isna(seq) or seq.strip() == "":
        return "UNKNOWN"

    if seq in seq_to_name:
        return seq_to_name[seq]

    if Hit_id.endswith(".hmm"):
        base = Hit_id.replace(".hmm", "")
        base_name = f"{base}like"
    else:
        mnemonic, digit, up, low, _ = parse_base_parts(Hit_id)
        if mnemonic is None:
            return "UNKNOWN"

        if nomen == "rank1":
            base_name = f"{mnemonic}like"
        elif nomen == "rank2":
            if digit:
                base_name = f"{mnemonic}{digit}like"
            else:
                base_name = f"{mnemonic}like"
        elif nomen == "rank3":
            if digit and up:
                base_name = f"{mnemonic}{digit}{up}like"
            elif digit:
                base_name = f"{mnemonic}{digit}like"
            else:
                base_name = f"{mnemonic}like"
        elif nomen == "rank4":
            if digit and up and low:
                return f"{mnemonic}{digit}{up}{low}"
            elif digit and up:
                return f"{mnemonic}{digit}{up}"
            elif digit:
                return f"{mnemonic}{digit}"
            else:
                return f"{mnemonic}"
        else:
            return "UNKNOWN"

    name_counters[base_name]["counter"] += 1
    new_id = f"{base_name}{name_counters[base_name]['counter']}"
    seq_to_name[seq] = new_id
    return new_id

# Zastosuj funkcję do danych
df["New_Hit_id"] = df.apply(generate_new_name, axis=1)

# Zapisz do pliku
df.to_excel(output_path, index=False)

