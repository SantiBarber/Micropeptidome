import pandas as pd

df = pd.read_csv("data/microproteins.csv", sep="\t", encoding="utf-8-sig")

# Diagnóstico: ver los nombres reales de las columnas
print("Columnas encontradas:")
print(df.columns.tolist())
print("\nPrimeras filas:")
print(df.head())