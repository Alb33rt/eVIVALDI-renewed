def LoadSpeciesStats(setup_file):
    bacteria=[];phage=[];matrix=[];superinfection=[];host_interactions={};superinfection_probabilities={}
    
    with open(setup_file, "r") as species_file:
        currently_analysing=None
        
        for line in species_file:
            if line=="#\n":continue
            if line.split("\t")[0]=="#BacteriaName":currently_analysing="Bacteria"; continue
            if line.split("\t")[0]=="#PhageName":currently_analysing="Phage"; continue
            if line.split("\t")[0]=="#HostMatrix\n":currently_analysing="BacPhageInteractions"; continue
            if line.split("\t")[0]=="#SuperinfectionProbability\n":currently_analysing="PhageSuperinfection"; continue
            
            if currently_analysing=="Bacteria":
                name=line.split("\t")[0]
                color=line.split("\t")[1]
                dna_id=line.split("\t")[2]
                inital_freq=float(line.split("\t")[3])
                growth_rate=float(line.split("\t")[4])
                ant_res=line.split("\t")[5].split(";")
                
                #Check if any mistake was made
                for a in ant_res:
                    if a not in ["Rif", "Str", "Quin", "NA"]: 
                        print("Unrecognized antibiotic resistance" + ant_res)
                        exit()
                    
                phage_res=line.split("\t")[6].split(";")                
                prophage=line.split("\t")[7].split(";")
                if (prophage[0]!="NA\n") and (prophage[0]!="NA"):
                    prophage=[p_tuple.split("\n")[0].replace("(","").replace(")","") for p_tuple in prophage]
                    prophage=[(p_tuple.split(",")[0],int(p_tuple.split(",")[1])) for p_tuple in prophage]
                else: prophage=None         
                bacteria.append({"Name":name, "Color":color, "DNAIdentifier":dna_id, "Freq":inital_freq, "Growth": growth_rate, "AntRes":ant_res, "PhageRes":phage_res, "Prophage":prophage})
            
            if currently_analysing=="Phage":
                name=line.split("\t")[0]
                color_g=line.split("\t")[1]
                phage_type=line.split("\t")[2]
                receptor=int(line.split("\t")[3])
                burst_size=int(line.split("\t")[4])
                dna_id=line.split("\t")[5]
                inital_numbers_env=int(line.split("\t")[6])
                
                lysogeny_alpha=float(line.split("\t")[7])#New in V9
                lysogeny_kappa=float(line.split("\t")[8])#New in V9
                induction_alpha=float(line.split("\t")[9])#New in V9
                induction_kappa=float(line.split("\t")[10])#New in V9
                gentransd_prob=float(line.split("\t")[11])#New in V9
                spectransd_prob=float(line.split("\t")[12])#New in V9
                
                max_cargo_size=int(line.split("\t")[13]) #New in V8.2
                phage.append({"Name":name, 
                              "color_genome":color_g, 
                              "Type":phage_type, 
                              "BurstSize": burst_size, 
                              "DNAIdentifier":dna_id, 
                              "FreeNumbers":inital_numbers_env, 
                              "Receptor":receptor,
                              "Lysogeny Alpha":lysogeny_alpha,
                              "Lysogeny Kappa":lysogeny_kappa,
                              "Induction Alpha":induction_alpha,
                              "Induction Kappa":induction_kappa,
                              "Generalized Transduction Probability": gentransd_prob,
                              "Specialized Transduction Probability": spectransd_prob, 
                              "MaxCargoSize":max_cargo_size})
                
            if currently_analysing=="BacPhageInteractions":matrix.append(line.split("\t"))
            if currently_analysing=="PhageSuperinfection":superinfection.append(line.split("\t"))
    
    
    print(host_interactions)
    if len(phage)>0:
        m_bac_species=[sp.replace("\n","") for sp in matrix[0] if sp!=""]
        m_phg_species=[]
        
        for line in matrix[1:]:
            m_phg_species.append(line[0])
            host_interactions[line[0]]={}
            for pos, entry in enumerate(line[1:]):
                host_interactions[line[0]][m_bac_species[pos]]=float(entry.split("\n")[0])
             
        phg_variants_incoming=[sp.replace("\n","") for sp in superinfection[0] if sp!=""]
        phg_variants_resident=[]
        
        for line in superinfection[1:]:
            phg_variants_resident.append(line[0])
            superinfection_probabilities[line[0]]={}
            for pos, entry in enumerate(line[1:]):
                superinfection_probabilities[line[0]][phg_variants_incoming[pos]]=float(entry.split("\n")[0])
    
    
    #Check if any mistake was made in phage names
    phage_names=[p["Name"] for p in phage]
    for b in bacteria:
        if b["Prophage"]!=None:
            for prop in b["Prophage"]:
                if prop[0] not in phage_names: 
                    print("Unrecognized prophage" + prop + phage_names)
                    exit()

    return bacteria, phage, host_interactions, superinfection_probabilities
