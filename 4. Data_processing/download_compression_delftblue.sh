#!/bin/bash

# Zorg dat de hoofdmap 'data' lokaal bestaat
mkdir -p ./data

# Het lokale error-bestand
error_log="./cifs_with_lammps_error.txt"
# Zorg dat het bestand bestaat, zodat we er later veilig in kunnen zoeken
touch "$error_log" 

# Definieer de basis-map op DelftBlue 
# (Let op: ik gebruik hier ~/scratch/ op basis van je eerdere foutmelding. 
# Verander dit naar ~/home/mmpool/ als dat toch de juiste locatie is).
user="nmelikian"
remote_base="~/MIDAS_BEP/final_nicoleblue"
#remote_base="~/final_thirdbatch"

# Loop over de cif_names van 1 tot en met 100
for cif_name in {518..527}; do
    
    # Definieer de paden voor deze specifieke cif_name
    remote_dir="${remote_base}/${cif_name}"
    remote_file="${remote_dir}/${cif_name}_box_compression_continuous.txt"
    remote_input="${remote_dir}/input.out"
    local_dir="./data/${cif_name}"
    
    # Stap 1: Controleer of het txt-bestand (succesvolle run) op DelftBlue bestaat
    if ssh ${user}@login.delftblue.tudelft.nl "test -f ${remote_file}" > /dev/null 2>&1; then
        
        echo "[INFO] Txt-bestand gevonden voor CIF ${cif_name}. Downloaden..."
        mkdir -p "${local_dir}"
        
        # Download het resultaat
        scp -p "${user}@login.delftblue.tudelft.nl:${remote_file}" "${local_dir}/"
        
        if [ $? -eq 0 ]; then
            echo "  -> [SUCCES] CIF ${cif_name} is gedownload."
            
            # --- OPRUIMEN BIJ SUCCES ---
            
            # 1. Verwijder lokaal input.out als het (nog) in het mapje staat
            if [ -f "${local_dir}/input.out" ]; then
                rm "${local_dir}/input.out"
                echo "  -> Oude input.out lokaal verwijderd."
            fi
            
            # 2. Haal de cif_name uit de error_log (als deze erin staat)
            # sed zoekt naar regels die beginnen met "cif_name:" en verwijdert deze ('d' voor delete)
            sed -i "/^${cif_name}:/d" "$error_log"
            
        else
            echo "  -> [FOUT] Downloaden mislukt voor CIF ${cif_name}."
        fi
        
    else
        # --- BESTAND BESTAAT NIET: CHECK INPUT.OUT ---
        echo "[INFO] Txt-bestand NIET gevonden voor CIF ${cif_name}. Checken op input.out..."
        
        if ssh mmpool@login.delftblue.tudelft.nl "test -f ${remote_input}" > /dev/null 2>&1; then
            echo "  -> input.out gevonden op DelftBlue. Downloaden..."
            mkdir -p "${local_dir}"
            
            # Download de input.out
            scp -p "mmpool@login.delftblue.tudelft.nl:${remote_input}" "${local_dir}/"
            
            if [ $? -eq 0 ]; then
                # Lees de eennalaatste regel (tail pakt de laatste 2, head daarvan de bovenste)
                error_msg=$(tail -n 2 "${local_dir}/input.out" | head -n 1)
                
                # Check of deze cif_name al in het logbestand staat
                # De -q zorgt dat grep stil is, hij geeft alleen een 'true' of 'false' terug
                if ! grep -q "^${cif_name}:" "$error_log"; then
                    # Schrijf de naam en de error weg
                    echo "${cif_name}: ${error_msg}" >> "$error_log"
                    echo "  -> Error toegevoegd aan cifs_with_lammps_error.txt"
                else
                    echo "  -> CIF ${cif_name} stond al in de error log, is overgeslagen."
                fi
            else
                echo "  -> [FOUT] Kon input.out niet downloaden voor CIF ${cif_name}."
            fi
        else
            echo "[OVERGESLAGEN] Geen txt én geen input.out gevonden voor CIF ${cif_name}."
        fi
    fi
    
done

echo "Klaar met het verwerken van de dataset!"
done

echo "Klaar met het verwerken van alle CIFs!"