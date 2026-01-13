#!/usr/bin/env python3

import pandas as pd
import sys

print(r"""
          _____           _   _     _                      
 _   _   |  __ \         | | (_)   | |                     
| | | |  | |__) |__ _ __ | |_ _  __| | ___  _ __ ___   ___ 
| | | |  |  ___/ _ \ '_ \| __| |/ _` |/ _ \| '_ ` _ \ / _ \
| |_| |  | |  |  __/ |_) | |_| | (_| | (_) | | | | | |  __/
| ___/   |_|   \___| .__/ \__|_|\__,_|\___/|_| |_| |_|\___|
| |                | |                                       
|_|                |_|                  almoraco 03/11/2025
""")

print(
    "FASTA creator for Micropeptidome pipeline\n"
    "Using Riboseq data from Barnes 2023.xlsx extract as FASTA the microproteins\n"
)

# Verificar que se pasaron los argumentos correctos
if len(sys.argv) != 3:
    print("Use: python3 fastear.py <archivo.xlsx> <salida.fasta>")
    print("For example: python3 fastear.py GSE197909_mouse.xlsx Barnes2023_smORFs.fasta")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

# Leer el archivo Excel
df = pd.read_excel(input_file)

with open(output_file, 'w') as fastafile:
    for index, row in df.iterrows():
        # Obtener las coordenadas directamente de la primera columna
        coordenadas = row['hg19 coordinates and strand']
        
        # Obtener la secuencia de la columna "Microprotein Sequence"
        secuencia = row['Microprotein Sequence']
        
        # Escribir en formato FASTA
        fastafile.write(f">Barnes 2023 {coordenadas}\n")
        fastafile.write(f"{secuencia}\n")

print(f"\n✅ FASTA file created: {output_file}")