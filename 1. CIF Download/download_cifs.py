import requests
import threading
import pandas as pd

def download_file(url, path):
    try:
        r = requests.get(url)
        r.raise_for_status() # Check of de download succesvol is
        with open(path, 'wb') as f:
            f.write(r.content)
        return True
    except Exception as e:
        print(f"Fout bij downloaden van {url}: {e}")
        return False

# 1. Lees het CSV bestand in
# Pas het pad aan als het bestand ergens anders staat
csv_pad = "/perovskite_data.csv"
try:
    df = pd.read_csv(csv_pad, sep=';', header=None)
except FileNotFoundError:
    print(f"Kan het bestand {csv_pad} niet vinden. Zorg dat het in dezelfde map staat of geef het volledige pad op.")
    exit()

# 2. Definieer de elementen om te controleren
elementen_vereist = ['Pb', 'I']
elementen_uitsluiten = ['Br', 'Cl', 'Sn', 'Ge', 'Bi']

# We gaan er vanuit dat de index (CIF nummer) in kolom 0 staat, 
# en de chemische formule in kolom 1.
for i in range(1, 850):
    # Zoek de rij op voor dit CIF nummer
    # Zorg dat het datatype overeenkomt (meestal int)
    rij = df[df[0] == i]
    
    if rij.empty:
        print(f"CIF {i} niet gevonden in de data. Overslaan.")
        continue
        
    # Haal de formule op. Neem de eerste waarde als er meerdere zijn.
    formule = str(rij.iloc[0][1])
    
    # Check of de vereiste elementen erin zitten
    heeft_pb_en_i = all(elem in formule for elem in elementen_vereist)
    
    # Check of de uitgesloten elementen NIET in de formule zitten
    heeft_geen_andere = not any(elem in formule for elem in elementen_uitsluiten)
    
    if heeft_pb_en_i and heeft_geen_andere:
        url = f"http://pdb.nmse-lab.ru/collection/cif/{i}.cif"
        path = f"/cifs/{i}.cif"
        
        gelukt = download_file(url, path)
        if gelukt:
            print(f"{i} is gedownload (Formule: {formule})")
            threading.Event().wait(0.5)
    else:
        print(f"{i} is overgeslagen (Formule: {formule})")