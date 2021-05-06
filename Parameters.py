
model_parameters={
                  "Main Directory": None,
                  
                  #SIMULATIONS#
                  "Random Seed": None,
                  "Number Simulations": None,
                  "Iterations": None,
                  "Iterations Outgrowth": None,
                  "Log Interval": None,
                  "Log Sample Genomes": [],
                  
                  #WORLD#
                  "World Size": None,
                  "World Type":None,
                  "Environmental Diffusion": None,
                  "Antibiotic Degradation Lambda RIF": None,
                  "Antibiotic Degradation Lambda STR": None,
                  "Antibiotic Degradation Lambda QUIN": None,
                  "Free Phage Degradation Lambda": None,
                                    
                  #ANTIBIOTICS#
                  "Max Antibiotic Concentration RIF": None,
                  "Max Antibiotic Concentration STR": None,
                  "Max Antibiotic Concentration QUIN": None,                  
                  "Antibiotic Times RIF": [],
                  "Antibiotic Times STR": [],
                  "Antibiotic Times QUIN": [],
                  
                  "Antibiotic Exposure Structured Environments": None,
                  "Death Curve A": None,
                  "Death Curve B": None,
                  "Death Curve M": None,
                  
                  #BACTERIA#
                  "Number Genes": None,
                  "Minimum Gene Size": None,
                  "Maximum Gene Size": None,
                  "Possible Receptors": None,
                  "Resistance Phage Receptor Mutation Probability": None,
                  "Resistance Phage CRISPR Mutation Probability": None,
                  "Resistance Phage Cost": None, 
                  "Gene Loss Probability": None,

                  #ARGS#
                  "ARG Mutation Probability": None,
                  "ARG Rif Cost": None,
                  "ARG Str Cost": None,
                  "ARG Quin Cost": None,
                  
                  #PHAGES#
                  "Phage Carrying Capacity": None,
                  "Abortive Infection Probability": None,                  
                  "Phage Receptor Mutation Probability": None,
                  "Phage CRISPR Mutation Probability": None,
                  "Phage Host Range Mutation Probability": None,
                  "Phage Recombination Probability": None,
                  "Burst Probability": None,                          
                  "Curing Probability": None,
                  "Specialized Transduction Genomic Distance": None,
                  "Infection Distance": None,
                  
                  "Failed Phage Infection DNA Prophage Integration":None,
                  "Failed Phage Infection DNA Chromosome Integration":None,
                  }

##Load parameters into main dictionary
def LoadParameters(pfile="Parameters.txt"):
    with open(pfile, 'r') as parameters_file:
        for line in parameters_file:
            if line=='\n': continue
            if line[0]=='#': continue
            
            param, value=line.replace("\n", "").split("=")
            
            if param=="Main Directory": model_parameters[param]=value
            
            #SIMULATIONS#
            if param=="Random Seed": 
                if value!=None: model_parameters[param]=int(value)
            if param=="Number Simulations": model_parameters[param]=int(value)        
            if param=="Iterations": model_parameters[param]=int(value)
            if param=="Iterations Outgrowth": model_parameters[param]=int(value)
            if param=="Log Interval": model_parameters[param]=int(value)
            if param=="Log Sample Genomes": model_parameters[param]=[int(t) for t in value.split(",")]
            
            #WORLD#
            if param=="World Size": model_parameters[param]=int(value)
            if param=="World Type": model_parameters[param]=value
            if param=="Environmental Diffusion": model_parameters[param]=value
            if param=="Antibiotic Degradation Lambda RIF": model_parameters[param]=float(value)            
            if param=="Antibiotic Degradation Lambda STR": model_parameters[param]=float(value)
            if param=="Antibiotic Degradation Lambda QUIN": model_parameters[param]=float(value)
            if param=="Free Phage Degradation Lambda": model_parameters[param]=float(value)
            
            #ANTIBIOTICS#    
            if param=="Max Antibiotic Concentration RIF": model_parameters[param]=int(value)
            if param=="Max Antibiotic Concentration STR": model_parameters[param]=int(value)
            if param=="Max Antibiotic Concentration QUIN": model_parameters[param]=int(value)
            
            if param=="Antibiotic Times RIF": model_parameters[param]=[int(t) for t in value.split(",")]
            if param=="Antibiotic Times STR": model_parameters[param]=[int(t) for t in value.split(",")]
            if param=="Antibiotic Times QUIN": model_parameters[param]=[int(t) for t in value.split(",")]
            
            if param=="Antibiotic Exposure Structured Environments": model_parameters[param]=value
            
            if param=="Death Curve (A, B, M)": 
                model_parameters["Death Curve A"]=float(value.split(",")[0].strip())
                model_parameters["Death Curve B"]=float(value.split(",")[1].strip())
                model_parameters["Death Curve M"]=float(value.split(",")[2].strip())
            
            #BACTERIA#
            if param=="Number Genes": model_parameters[param]=int(value)
            if param=="Minimum Gene Size": model_parameters[param]=int(value)
            if param=="Maximum Gene Size": model_parameters[param]=int(value)
            if param=="Possible Receptors": model_parameters[param]=int(value)
            if param=="Resistance Phage Receptor Mutation Probability": model_parameters[param]=float(value)
            if param=="Resistance Phage CRISPR Mutation Probability": model_parameters[param]=float(value)
            if param=="Resistance Phage Cost": model_parameters[param]=float(value)
            if param=="Gene Loss Probability": model_parameters[param]=float(value)
            
            #ARGS#
            if param=="ARG Mutation Probability": model_parameters[param]=float(value)
            if param=="ARG Rif Cost": model_parameters[param]=float(value)
            if param=="ARG Str Cost": model_parameters[param]=float(value)
            if param=="ARG Quin Cost": model_parameters[param]=float(value)
            
            #PHAGES#
            if param=="Phage Carrying Capacity": model_parameters[param]=int(value)
            if param=="Abortive Infection Probability": model_parameters[param]=float(value)
            if param=="Phage Receptor Mutation Probability": model_parameters[param]=float(value)
            if param=="Phage CRISPR Mutation Probability": model_parameters[param]=float(value)
            if param=="Phage Host Range Mutation Probability": model_parameters[param]=float(value)
            if param=="Phage Recombination Probability": model_parameters[param]=float(value)
            if param=="Burst Probability": model_parameters[param]=float(value)
                        
            if param=="Curing Probability": model_parameters[param]=float(value)
            if param=="Specialized Transduction Genomic Distance": model_parameters[param]=int(value)
            if param=="Infection Distance": model_parameters[param]=int(value)
            
            if param=="Failed Phage Infection DNA Prophage Integration": model_parameters[param]=float(value)
            if param=="Failed Phage Infection DNA Chromosome Integration": model_parameters[param]=float(value)
            
            
    for p in model_parameters:
        if (model_parameters[p]==None) or (model_parameters[p]==[]):
            print(p + "Missing!")
            exit()
    

##Print conditions used for this simulation
#TODO