import os
from Parameters import model_parameters as param

dictio_entries=["Bacteria Population Size",
           "Bacteria Subpopulation Size",
           "Bacteria ARG Population Size",
           "Bacteria ARG Subpopulation Size",           
           "Bacteria OnLysis Population Size",
           "Bacteria OnLysis Subpopulation Size",                      
           "Bacteria Temperate Population Size",
           "Bacteria Temperate Subpopulation Size",           
           "Bacteria Defective Population Size",
           "Bacteria Defective Subpopulation Size",
           "Bacteria Phage Resistant Population Size",
           "Bacteria Phage Resistant Subpopulation Size",
           "Bacteria With Phage With ARG Population Size",
           "Bacteria With Phage With ARG Subpopulation Size",           
           "Phage Population Size",                    
           "Phage Subpopulation Size",
           "Phage With ARG Population Size",           
           "Phage With ARG Subpopulation Size",
           "Temperate Phage Population Size",
           "Virulent Phage Population Size",
           "Defective Phage Population Size",
           ]

extra_stats_dictio_entries=["Bacteria With Receptor","Phage With Receptor"]

log_stats={
           "Bacteria Population Size":0,
           "Bacteria Subpopulation Size":{},
           "Bacteria ARG Population Size":0,
           "Bacteria ARG Subpopulation Size":{},           
           "Bacteria OnLysis Population Size":0,
           "Bacteria OnLysis Subpopulation Size":{},                      
           "Bacteria Temperate Population Size":0,
           "Bacteria Temperate Subpopulation Size":{},           
           "Bacteria Defective Population Size":0,
           "Bacteria Defective Subpopulation Size":{},
           "Bacteria Phage Resistant Population Size":0,
           "Bacteria Phage Resistant Subpopulation Size":{},
           "Bacteria With Phage With ARG Population Size":0,
           "Bacteria With Phage With ARG Subpopulation Size":{},
           
           "Phage Population Size":0,                      
           "Phage Subpopulation Size":{},
           "Phage With ARG Population Size":0,           
           "Phage With ARG Subpopulation Size":{},
           "Temperate Phage Population Size":0,
           "Virulent Phage Population Size":0,
           "Defective Phage Population Size":0,
           }

extra_stats_log_stats={"All Bacteria With Receptor":{}, "Bacteria With Receptor":{},"Phage With Receptor":{} }

def PrintStats(iteration, log_stats, output_file):
    #If file does not exist yet    
    if (not(os.path.isfile(output_file))) or (iteration==0):
        with open(output_file, "w") as out_file:
            #Write header
            for entry in dictio_entries:                
                if isinstance(log_stats[entry], dict):
                    if entry.split(" ")[0]=="Bacteria":                   
                        for species in param["IndividualBacteriaParameters"]["Basal Reproduction Probability"]: out_file.write("\t"+species)
                    if entry.split(" ")[0]=="Phage":  
                        for species in param["IndividualPhageParameters"]["Burst Size"]: out_file.write("\t"+species+"\t"+species+"_Defective")
                else:
                    out_file.write("\t"+entry)
                    
            out_file.write("\n"+str(iteration))
            
            for entry in dictio_entries:                
                if isinstance(log_stats[entry], dict):
                    if entry.split(" ")[0]=="Bacteria":                   
                        for species in param["IndividualBacteriaParameters"]["Basal Reproduction Probability"]: out_file.write("\t"+str(log_stats[entry][species]))
                    if entry.split(" ")[0]=="Phage":  
                        for species in param["IndividualPhageParameters"]["Burst Size"]: out_file.write("\t"+str(log_stats[entry][species])+ "\t"+str(log_stats[entry][species+"_Defective"]))
                else:
                    out_file.write("\t"+str(log_stats[entry]))
    else:
        with open(output_file, "a") as out_file:
            out_file.write("\n"+str(iteration))
            for entry in dictio_entries:                
                if isinstance(log_stats[entry], dict):
                    if entry.split(" ")[0]=="Bacteria":                   
                        for species in param["IndividualBacteriaParameters"]["Basal Reproduction Probability"]: out_file.write("\t"+str(log_stats[entry][species]))
                    if entry.split(" ")[0]=="Phage":  
                        for species in param["IndividualPhageParameters"]["Burst Size"]: out_file.write("\t"+str(log_stats[entry][species])+ "\t"+str(log_stats[entry][species+"_Defective"]))
                else:
                    out_file.write("\t"+str(log_stats[entry]))
                    
                    
                    
def AppendStats_OLD(simulation, iteration, log_stats, extra_stats, output_file):
    #If file does not exist yet    
    if not(os.path.isfile(output_file)):
        with open(output_file, "w") as out_file:
            #Write header
            out_file.write("Simulation\tRandom Seed\tIteration")
            
            for entry in dictio_entries:                
                if isinstance(log_stats[entry], dict):
                    if entry.split(" ")[0]=="Bacteria":                   
                        for species in param["IndividualBacteriaParameters"]["Basal Reproduction Probability"]: out_file.write("\t"+species)
                    if entry.split(" ")[0]=="Phage":  
                        for species in param["IndividualPhageParameters"]["Burst Size"]: out_file.write("\t"+species+"\t"+species+"_Defective")
                else:
                    out_file.write("\t"+entry)
            
            for rec in range(0, param["Possible Receptors"]):
                for species in param["IndividualBacteriaParameters"]["Basal Reproduction Probability"]:out_file.write("\t"+species+"_Receptor"+str(rec))                                   
                out_file.write("\tPhage With Receptor "+str(rec))
                                   
            out_file.write("\n"+str(simulation)+"\t"+str(param["Random Seed"])+"\t"+str(iteration))
            
            for entry in dictio_entries:                
                if isinstance(log_stats[entry], dict):
                    if entry.split(" ")[0]=="Bacteria":                   
                        for species in param["IndividualBacteriaParameters"]["Basal Reproduction Probability"]: out_file.write("\t"+str(log_stats[entry][species]))
                    if entry.split(" ")[0]=="Phage":  
                        for species in param["IndividualPhageParameters"]["Burst Size"]: out_file.write("\t"+str(log_stats[entry][species])+ "\t"+str(log_stats[entry][species+"_Defective"]))
                else:
                    out_file.write("\t"+str(log_stats[entry]))
                                
            
            for rec in range(0, param["Possible Receptors"]):
                out_file.write("\t"+str(extra_stats["Bacteria With Receptor"][rec])+"\t"+str(extra_stats["Phage With Receptor"][rec]))

    else:
        with open(output_file, "a") as out_file:
            out_file.write("\n"+str(simulation)+"\t"+str(param["Random Seed"])+"\t"+str(iteration))
            for entry in dictio_entries:                
                if isinstance(log_stats[entry], dict):
                    if entry.split(" ")[0]=="Bacteria":                   
                        for species in param["IndividualBacteriaParameters"]["Basal Reproduction Probability"]: out_file.write("\t"+str(log_stats[entry][species]))
                    if entry.split(" ")[0]=="Phage":  
                        for species in param["IndividualPhageParameters"]["Burst Size"]: out_file.write("\t"+str(log_stats[entry][species])+ "\t"+str(log_stats[entry][species+"_Defective"]))
                else:
                    out_file.write("\t"+str(log_stats[entry]))
            
            for rec in range(0, param["Possible Receptors"]):
                for species in param["IndividualBacteriaParameters"]["Basal Reproduction Probability"]:                
                    out_file.write("\t"+str(extra_stats["Bacteria With Receptor"][rec][species]))
                                   
                out_file.write("\t"+str(extra_stats["Phage With Receptor"][rec]))

                
def AppendStats(simulation, iteration, log_stats, extra_stats, output_file):
    #If file does not exist yet, write headers    
    if not(os.path.isfile(output_file)) or ((iteration==-param["Iterations Outgrowth"]) and (simulation==0)):
        with open(output_file, "w") as out_file:            
            out_file.write("Simulation\tRandom Seed\tIteration")
            
            for entry in dictio_entries:                
                if isinstance(log_stats[entry], dict):
                    if entry.split(" ")[0]=="Bacteria":                   
                        for species in param["IndividualBacteriaParameters"]["Basal Reproduction Probability"]: out_file.write("\t"+species)
                    if entry.split(" ")[0]=="Phage":  
                        for species in param["IndividualPhageParameters"]["Burst Size"]: out_file.write("\t"+species+"\t"+species+"_Defective")
                else:
                    out_file.write("\t"+entry)            
                        
            for rec in range(0, param["Possible Receptors"]):
                for species in param["IndividualBacteriaParameters"]["Basal Reproduction Probability"]:out_file.write("\t"+species+"_Receptor"+str(rec))                                   
                out_file.write("\tPhage With Receptor "+str(rec))
               
    #Always write entry
    with open(output_file, "a") as out_file:
        out_file.write("\n"+str(simulation)+"\t"+str(param["Random Seed"])+"\t"+str(iteration))
        for entry in dictio_entries:                
            if isinstance(log_stats[entry], dict):
                if entry.split(" ")[0]=="Bacteria":                   
                    for species in param["IndividualBacteriaParameters"]["Basal Reproduction Probability"]: out_file.write("\t"+str(log_stats[entry][species]))
                if entry.split(" ")[0]=="Phage":  
                    for species in param["IndividualPhageParameters"]["Burst Size"]: out_file.write("\t"+str(log_stats[entry][species])+ "\t"+str(log_stats[entry][species+"_Defective"]))
            else:
                out_file.write("\t"+str(log_stats[entry]))            
        
        for rec in range(0, param["Possible Receptors"]):            
            for species in param["IndividualBacteriaParameters"]["Basal Reproduction Probability"]:                 
                out_file.write("\t"+str(extra_stats["Bacteria With Receptor"][rec][species]))
                                               
            out_file.write("\t"+str(extra_stats["Phage With Receptor"][rec]))

        

def GetStats(bacteria_list, world):
    
    #Bacteria
    for bac_species in param["IndividualBacteriaParameters"]["Basal Reproduction Probability"]:
        log_stats["Bacteria Subpopulation Size"][bac_species]=0
        log_stats["Bacteria ARG Subpopulation Size"][bac_species]=0
        log_stats["Bacteria OnLysis Subpopulation Size"][bac_species]=0
        log_stats["Bacteria Temperate Subpopulation Size"][bac_species]=0
        log_stats["Bacteria Defective Subpopulation Size"][bac_species]=0
        log_stats["Bacteria Phage Resistant Subpopulation Size"][bac_species]=0
        log_stats["Bacteria With Phage With ARG Subpopulation Size"][bac_species]=0
        
    for b in bacteria_list:
        log_stats["Bacteria Subpopulation Size"][b.name]+=1
        if b.resistances["AR_Rif"]: log_stats["Bacteria ARG Subpopulation Size"][b.name]+=1
        if b.phage_resistance==True: log_stats["Bacteria Phage Resistant Subpopulation Size"][b.name]+=1
        if len(b.phages)>0:            
            phg_types=[phg["Type"] for phg in b.phages]
            if ("Virulent" in phg_types) or ("Induced_Temperate" in phg_types): log_stats["Bacteria OnLysis Subpopulation Size"][b.name]+=1
            if "Temperate" in phg_types: log_stats["Bacteria Temperate Subpopulation Size"][b.name]+=1
            if "Defective" in phg_types: log_stats["Bacteria Defective Subpopulation Size"][b.name]+=1
            for phg in b.phages:
                if (not(phg["Type"]=="Virulent")) and (not(phg["Type"]=="Induced_Temperate")):
                    #if (b.color=="green") and (b.resistances["AR_Rif"]): print phg["Cargo"]                              
                    idxs = [phg["Cargo"][i][3] for i, x in enumerate(zip(*phg["Cargo"])[0]) if "AR_Rif" in x]                
                    if "Resistant" in idxs: log_stats["Bacteria With Phage With ARG Subpopulation Size"][b.name]+=1
                
        if (b.resistances["AR_Rif"]):
            for gene in b.genome:                
                if not(isinstance(gene, dict)) and ("AR" in gene[0]) and (gene[3]=="Resistant") and (len(gene)==5) and ("_rec" in gene[4]):
                    log_stats["Bacteria With Phage With ARG Subpopulation Size"][b.name]+=1
                    break #Only counts once!!

                                
    log_stats["Bacteria Population Size"]=sum(log_stats["Bacteria Subpopulation Size"].values())
    log_stats["Bacteria ARG Population Size"]=sum(log_stats["Bacteria ARG Subpopulation Size"].values())
    log_stats["Bacteria OnLysis Population Size"]=sum(log_stats["Bacteria OnLysis Subpopulation Size"].values())
    log_stats["Bacteria Temperate Population Size"]=sum(log_stats["Bacteria Temperate Subpopulation Size"].values())
    log_stats["Bacteria Defective Population Size"]=sum(log_stats["Bacteria Defective Subpopulation Size"].values())
    log_stats["Bacteria Phage Resistant Population Size"]=sum(log_stats["Bacteria Phage Resistant Subpopulation Size"].values())
    log_stats["Bacteria With Phage With ARG Population Size"]=sum(log_stats["Bacteria With Phage With ARG Subpopulation Size"].values())
    
    #Phage
    for phg_species in param["IndividualPhageParameters"]["Burst Size"]:
        log_stats["Phage Subpopulation Size"][phg_species]=0
        log_stats["Phage Subpopulation Size"][phg_species+"_Defective"]=0
        log_stats["Phage With ARG Subpopulation Size"][phg_species]=0
        log_stats["Phage With ARG Subpopulation Size"][phg_species+"_Defective"]=0
    
    log_stats["Temperate Phage Population Size"]=0
    log_stats["Virulent Phage Population Size"]=0
    log_stats["Defective Phage Population Size"]=0
    
    for loc in world.world:
        for phg in world.world[loc].free_phages:
            phage_fam=(phg[1]["Family"]+"_Defective" if phg[1]["Type"]=="Defective" else phg[1]["Family"])
            log_stats["Phage Subpopulation Size"][phage_fam]+=1            
            idxs = [phg[1]["Cargo"][i][3] for i, x in enumerate(zip(*phg[1]["Cargo"])[0]) if "AR_Rif" in x]
            if "Resistant" in idxs:log_stats["Phage With ARG Subpopulation Size"][phage_fam]+=1
            
            if phg[1]["Type"]=="Temperate": log_stats["Temperate Phage Population Size"]+=1
            if phg[1]["Type"]=="Virulent": log_stats["Virulent Phage Population Size"]+=1
            if phg[1]["Type"]=="Defective": log_stats["Defective Phage Population Size"]+=1
            
    log_stats["Phage Population Size"]=sum(log_stats["Phage Subpopulation Size"].values())
    log_stats["Phage With ARG Population Size"]=sum(log_stats["Phage With ARG Subpopulation Size"].values())
        
    return log_stats


def GetExtraStats(bacteria_list, world):
    for rec in range(0, param["Possible Receptors"]):
        extra_stats_log_stats["All Bacteria With Receptor"][rec]=0
        extra_stats_log_stats["Bacteria With Receptor"][rec]={}
        for species in param["IndividualBacteriaParameters"]["Basal Reproduction Probability"]:extra_stats_log_stats["Bacteria With Receptor"][rec][species]=0
        extra_stats_log_stats["Phage With Receptor"][rec]=0
    
    for b in bacteria_list:
        for rec in range(0, param["Possible Receptors"]):                      
            if rec in b.phage_receptor_resistances: 
                extra_stats_log_stats["Bacteria With Receptor"][rec][b.name]+=1
                extra_stats_log_stats["All Bacteria With Receptor"][rec]+=1

    for loc in world.world:
        for phg in world.world[loc].free_phages: extra_stats_log_stats["Phage With Receptor"][phg[1]["Receptor"]]+=1

    return extra_stats_log_stats

# this file has been completely translated