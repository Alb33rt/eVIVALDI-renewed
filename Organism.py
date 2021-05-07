import random, math, itertools, copy, numpy, string
from Parameters import model_parameters as param
import string

class Bacteria():
    def __init__(self, name=None, genome=None, color=None, ar_locus=[0,0,0], movement=1):
        
        self.name=name
        self.genome=copy.copy(genome)
        
        self.color=color
        
        self.ar_locus=ar_locus #Bacteria start as sensitive to all antibiotics
        self.location=None
        
        self.just_reproduced=False
        
        self.movement=movement
        
        self.resistances={"AR_Rif":False, "AR_Str":False, "AR_Quin":False} #Bacteria start as sensitive to all antibiotics
        
        self.phage_receptor_resistances=[]        
        self.phage_resistance=False
        self.phage_resistance_cassettes=[]
        
        self.transformable=True        
        self.is_ARtransformed=False
        
        self.phages=[] 
        
        
    def Set_Location(self, location, set_un=True):
        self.location=location
        self.location.occupants.append(self)
        if set_un: self.location.Add_Unavailable()
        
    def Leave_Current_Location(self, set_a=True):
        self.location.occupants.remove(self)
        if set_a: 
            self.location.Add_Available()
        self.location=None
    
    def Get_Death_Probability(self):
        
        death_rate_modifiers=0 #will have the sum of the factors that might increase the probability of death (antibiotics in the environment, phages...)
        
        a=param["Death Curve A"]
        b=param["Death Curve B"]
        m=param["Death Curve M"]
        
        ant_concentration=0
        ant_concentration+=(0 if self.resistances["AR_Rif"] else self.location.rif_conc)
        ant_concentration+=(0 if self.resistances["AR_Str"] else self.location.str_conc)
        ant_concentration+=(0 if self.resistances["AR_Quin"] else self.location.quin_conc)
        
        death_rate_modifiers+=a+((1-a)/(1+math.exp(-b*(ant_concentration-m))))
        
        #Upon phage excision, death is certain        
        if self.InnerPhageLifeCycle(): 
            death_rate_modifiers+=1
        
        return death_rate_modifiers
    
    def Death(self):
        
        coord=copy.deepcopy(self.location.coordinates)
        #self.GenerateTransformableDNA() #Future function
        self.Leave_Current_Location(set_a=False)
        
        ##Do other stuff: DNA release, phage, etc...
        return coord
    
    def AssertCandidacy(self):        
        base_prob_reproduction=param["IndividualBacteriaParameters"]["Basal Reproduction Probability"][self.name]
    
        #Add a reproductive cost for each ARG
        cumulative_cost=base_prob_reproduction
        for gene in self.genome:
            if isinstance(gene, tuple): #If it is a gene originally in the bacteria
                if (("AR_Rif" in gene[0]) and (gene[3]=="Resistant")):
                    cumulative_cost*=(1-param["ARG Rif Cost"])
                if (("AR_Str" in gene[0]) and (gene[3]=="Resistant")):
                    cumulative_cost*=(1-param["ARG Str Cost"])
                if (("AR_Quin" in gene[0]) and (gene[3]=="Resistant")):
                    cumulative_cost*=(1-param["ARG Quin Cost"])
            if isinstance(gene, dict): 
                pass 
                #It is a mobile element, should it also count as a cost
        
        #Add reproductive cost for changes in the receptor conferring phage resistance
        # decide if it should depend on the number of receptor resistances
        if len(self.phage_receptor_resistances)>0:
            cumulative_cost*=(1-param["Resistance Phage Cost"])        
        
        #So that it does not go below zero (if everyone is zero, nobody is chosen)
        prob_reproduction=max(0.000001, cumulative_cost)
        
        if self.just_reproduced==True: 
            prob_reproduction=0 #Added v4: Only case when reproduction is zero is when the organism has already reproduced
        for p in self.phages: 
            if (p["Type"]=="Virulent") or (p["Type"]=="Induced_Temperate"): 
                prob_reproduction=0 
                #Added v5: Another case when reproduction is zero: there is a lytic phage within the cell
        
        if random.random()<prob_reproduction: 
            return prob_reproduction
        else: 
            return None        
        
    def Mutation(self):

        #SNP Mutations in loci
        for idx, gene in enumerate(self.genome):
            if isinstance(gene, dict):
                for idx_cargo, cargo_gene in enumerate(gene["Cargo"]):
                    #ARG mutations
                    if (("AR_Rif" in cargo_gene[0]) or ("AR_Str" in cargo_gene[0]) or ("AR_Quin" in cargo_gene[0])): 
                        if random.random()<param["ARG Mutation Probability"]:
                            if cargo_gene[3]=="Sensitive": 
                                self.genome[idx]["Cargo"][idx_cargo]=(cargo_gene[0], cargo_gene[1], cargo_gene[2], "Resistant")
                                self.resistances[cargo_gene[0].split(":")[1]]=True #Acquisition of resistance
                            else: 
                                self.genome[idx]["Cargo"][idx_cargo]=(cargo_gene[0], cargo_gene[1], cargo_gene[2], "Sensitive")
                                self.resistances[cargo_gene[0].split(":")[1]]=False #Loss of resistance
                                            
                    #Everything else
                    else: 
                        pass 
                        # Define mutations in non-ARG genes                        
            else:
                #ARG mutations            
                if (("AR_Rif" in gene[0]) or ("AR_Str" in gene[0]) or ("AR_Quin" in gene[0])) and (random.random()<param["ARG Mutation Probability"]):
                    if gene[3]=="Sensitive": 
                        self.genome[idx]=(gene[0], gene[1], gene[2], "Resistant")
                        self.resistances[gene[0].split(":")[1]]=True #Acquisition of resistance
                    else: 
                        self.genome[idx]=(gene[0], gene[1], gene[2], "Sensitive")
                        self.resistances[gene[0].split(":")[1]]=False #Loss of resistance
                                    
                #Everything else
                else: 
                    pass 
                    # Define mutations in non-ARG genes
                
        if ((random.random()<param["Resistance Phage Receptor Mutation Probability"]) and (len(self.phage_receptor_resistances)<param["Possible Receptors"])):

            #Which ones are missing?
            existing_receptors=set(self.phage_receptor_resistances)
            all_possible_receptors = set(range(0,param["Possible Receptors"]))
            
            #Choose a random one from those        
            if len(all_possible_receptors - existing_receptors)!=0:            
                self.phage_receptor_resistances.append(random.choice(list(all_possible_receptors - existing_receptors)))
                
            # Add in genome? Also subject to gene loss                            
        
        #Gene loss
        def GeneLoss():            
            print("GeneLoss FUNCTION is deactivated in (v9.2). IT NEEDS TO BE REFACTORED TO USE MULTIPLE ANTIBIOTIC RESISTANCE GENES")
            exit()
                
            if isinstance(gene, dict): #mobile genes, good to delete

                if len(gene["Cargo"])==1: 
                    return False #Always keeps at least one cargo gene!!!
                else: 
                    save_gene=gene["Cargo"][0]
                                
                gene["Cargo"][:]=[cargo_gene for cargo_gene in gene["Cargo"] if random.random()>param["Gene Loss Probability"]]
                            
                if len(gene["Cargo"])==0: 
                    gene["Cargo"][:]=[save_gene] #Always keeps at least one cargo gene!!! (mroe than one can be deleted in the line above if gene loss probability is high enough...)
                
                return False #Always returns false because the mobile element is still there (at least the "structure")
#             else:
#                 if (gene[0]=="N") and (random.random()<param["Gene Loss Probability"]):return True #Non-coding genes, good to delete            
#                 else: return False            
#             print "Should not be here!"; exit()   
        
        
        #TEMPORARY: REMOVE WHEN GENE LOSS IS FUNCTIONAL AGAIN!    
        return 
    
        length_genome_before=len(self.genome)
        lengths_mobile_before=[len(phage["Cargo"]) for phage in self.phages]
        self.genome[:] = [gene for gene in self.genome if not GeneLoss()]
                
        if len(self.genome)!=length_genome_before: #Check if there were deletions
            #Update mobile elements position accordingly
            for phage in self.phages:phage["Position"]=self.genome.index(phage)
        
        if lengths_mobile_before!=[len(phage["Cargo"]) for phage in self.phages]:
            #Check if ARG was deleted
            arg_res_found=False
            for gene in self.genome:                
                if isinstance(gene, tuple): #If it is a gene originally in the bacteria
                    if (("AR_Rif" in gene[0]) and (gene[3]=="Resistant")): 
                        arg_res_found=True
                        break
                if isinstance(gene, dict):
                    for sub_gene in gene["Cargo"]: #If it is a gene within a mobile element
                        if (("AR_Rif" in sub_gene[0]) and (sub_gene[3]=="Resistant")): 
                            arg_res_found=True
                            break
                        
            if not(arg_res_found):self.resistances["AR_Rif"]=False
                
                            

    def Reproduce(self):
        #Creates new offspring. Let's leave the genome empty for now
        offspring=Bacteria(name=self.name, genome=None, color=self.color)
        
        #Let's add the "accessory" mobile elements
        offspring.phages=[copy.deepcopy(ph) for ph in self.phages] 
        #Try to make this faster...
        
        #Let's create the genome with the correct references for mobile elements
        new_genome=[]
        for gene in self.genome:
            if isinstance(gene, dict):
                pass 
                #Add only the other genes for now                
            else: 
                new_genome.append(gene) #"normal" gene
        
        #Finally insert the mobile elements in the correct positions
        for phg in offspring.phages:
            new_genome.insert(phg["Position"], phg)
        
        #And set the new genome to the offspring
        offspring.genome=new_genome
        
        #Resistance profiles should be the same
        offspring.resistances=copy.deepcopy(self.resistances)
        offspring.phage_resistance=self.phage_resistance #Used for "Simple" phage resistance
        offspring.phage_resistance_cassettes=copy.deepcopy(self.phage_resistance_cassettes) #Used for CRISPR based phage resistance
        offspring.phage_receptor_resistances=copy.deepcopy(self.phage_receptor_resistances) #Used for receptor based phage resistance
                
        #Mutate the new offspring (maybe)
        offspring.Mutation()                   
           
        self.just_reproduced=True #Set this so that it doesn't reproduce multiple times per iteration
        return offspring
    
    def GenerateTransformableDNA(self):
        print("Alpha function, not usable")
        #Add exogenous DNA resulting from cell death
        eDNA_distance=3
        eDNA_pieces=numpy.random.randint(0,3)
        eDNA_AR_prob=0.1
        
        deposit_locations=self.location.Get_Reachable(eDNA_distance, desired_type="LocationObj")
        deposit_locations.append(self.location)
        
        for _ in range(eDNA_pieces):
            if random.random()<eDNA_AR_prob: 
                random.choice(deposit_locations).eDNA.append(("AR", self.genome))
            else: 
                random.choice(deposit_locations).eDNA.append(("Trash", self.genome))
            
    def Transformation(self):
        print("Alpha function, not usable")
        #Need to decide how to to this. Get one at random or have the same probability for each available DNA piece
        transf_prob=0.01
        eDNA_trans_distance=2
        
        if random.random()>transf_prob: 
            return 
        
        import_locations=self.location.Get_Reachable(eDNA_trans_distance, desired_type="LocationObj")
        import_locations.append(self.location)
        
        for loc in import_locations:
            if len(loc.eDNA)>0: 
                if random.random()<transf_prob: 
                    imported_dna=random.choice(loc.eDNA)
                    if imported_dna[0]=="AR": 
                        self.ar_locus[0]=1
                        self.is_ARtransformed=True
                    else: 
                        pass 
                        #randomly insert in other genes
                    loc.eDNA.remove(imported_dna) 
                    #remove imported DNA from environment
                    
    def GetPhageARG(self):
        ARG_inPhage=[]
        for phg in self.phages:
            for gene in phg["Cargo"]: 
                if "AR" in gene[0]: 
                    ARG_inPhage.append(gene)
        return ARG_inPhage
    
    
    def EnvironmentalInfection(self):
        #Get surrounding areas from which infection can arise
        phage_origins=self.location.Get_Reachable(param["Infection Distance"], desired_type="LocationObj")        
        phage_origins.append(self.location)#Add self location as possible source of infection
        
        for location in phage_origins:
            if len(location.free_phages)>0:

                to_remove=-1
                for idx, phage in enumerate(location.free_phages):
                    
                    #Successful infection removes phage from the environment       
                    if self.PhageInfection(phage_type=phage[1]["Type"],
                                           family=phage[1]["Family"], 
                                           donor_genome=phage[0], 
                                           cargo=phage[1]["Cargo"], 
                                           crispr_sequence=phage[1]["crispr_seq"], 
                                           phage_receptor=phage[1]["Receptor"],
                                           maxcargo_size=phage[1]["MaxCargoSize"]):                        
                        to_remove=idx
                        break
                
                if to_remove!=-1: 
                    del location.free_phages[to_remove] 
                    #Should the phage be removed if infection is unsuccessful?
                    break 
                    #After infection, stop looking 
                    #Does this make sense? MOI?
                      
    
    def IntegratePhageDNA(self, infecting_phage_cargo):
                
        if len(infecting_phage_cargo)==0: 
            print("Phage has less than 1 gene. This shouldn't happen")
        
        #Choose random gene from phage cargo
        incoming_gene=random.choice(infecting_phage_cargo)
        
        #If there is a resident prophage (active or defective) calculate probability to integrate within that prophage (per phage)
        if len(self.phages)>0:
            for resident_phage in self.phages:
                if ((resident_phage["Type"]=="Temperate") or (resident_phage["Type"]=="Defective")) and random.random()<param["Failed Phage Infection DNA Prophage Integration"]:
                     
                    #Insert DNA at random position and leave function, indicating the phage has been chopped
                    insertion_pos=random.choice(range(0, len(resident_phage["Cargo"])))
                    resident_phage["Cargo"].insert(insertion_pos, incoming_gene)                    
                    resident_phage["Cargo"][insertion_pos]=resident_phage["Cargo"][insertion_pos][0:-1]+(resident_phage["Cargo"][insertion_pos][-1]+"_rec",)                    
                    
                    #Becomes defective if too much DNA within phage
                    if resident_phage["Cargo"]>resident_phage["MaxCargoSize"]: 
                        resident_phage["Type"]="Defective"
                                                        
                    #Check if it is an antibiotic resistance gene
                    if ("AR" in incoming_gene[0]) and (incoming_gene[3]=="Resistant"):
                        self.resistances[incoming_gene[0].split(":")[1]]=True                          
                    return True       
        
        #Probability to integrate within the bacterial chromosome
        if random.random()<param["Failed Phage Infection DNA Chromosome Integration"]:
            insertion_pos=random.choice(range(0, len(self.genome)-1))
            self.genome.insert(insertion_pos, incoming_gene);     
            
            self.genome[insertion_pos]=self.genome[insertion_pos][0:-1]+(self.genome[insertion_pos][-1]+"_rec",)
        
            #Check if it is an antibiotic resistance gene 
            if ("AR" in incoming_gene[0]) and (incoming_gene[3]=="Resistant"):
                self.resistances[incoming_gene[0].split(":")[1]]=True
            return True
        
        return False
                
    def PhageInfection(self, phage_type, family, donor_genome, cargo, crispr_sequence, phage_receptor, maxcargo_size):
              
        #Infection occurs with a certain probability (before defences and superinfection protection), accounting for host range        
        infection_probability=param["InfectionMap"][family][self.name]
                
        if random.random()<infection_probability:     
                        
            #Check superinfection protection
            if len(self.phages)>0:
                incoming_phage=(family+"_Defective" if phage_type=="Defective" else family)
                all_superinf_probs=[]
                for phg in self.phages:
                    resident_phage=(phg["Family"]+"_Defective" if phg["Type"]=="Defective" else phg["Family"])                     
                    all_superinf_probs.append(param["Superinfection"][resident_phage][incoming_phage]) #Collect the probabilities from all integrated phage, then choose the highest
                                
                if random.random()<max(all_superinf_probs):
                    return self.IntegratePhageDNA(cargo) #; return False 

                
            #Check resistance by receptor
            if phage_receptor in self.phage_receptor_resistances:
                return self.IntegratePhageDNA(cargo)#; return False
            
            #Check resistance by CRISPR (existing cassettes)
            if crispr_sequence in self.phage_resistance_cassettes: 
                return self.IntegratePhageDNA(cargo)
        
            #Check resistance by CRISPR: probability of acquisition of new crispr cassette (automatically fails infection)
            if random.random()<param["Resistance Phage CRISPR Mutation Probability"]: 
                self.phage_resistance_cassettes.append(crispr_sequence) 
                self.IntegratePhageDNA(cargo)
                return True #CRISPR acquisition always titers phage from the environment
                                                                    
            ###If this point is reached, phage will definitely infect###
            
            #Check if there was cargo
            if cargo==None: 
                cargo=[]
                                
            #Choose location of insertion in genome
            position=random.choice([pos for pos in range(0, len(self.genome))])
            self.phages.append({"Family":family, "Position":position, "Cargo": cargo, "Time":0, 
                                "Type":phage_type,
                                "Receptor":phage_receptor,
                                "crispr_seq": crispr_sequence,
                                "MaxCargoSize": maxcargo_size})
            
            #Only integrates if it is not Virulent...
            if (phage_type!="Virulent"):
                #... or a temperate that does not go into lysogeny and instead it enters the lytic pathway
                
                #Calculate probability of lysogeny based on MOI/surrounding phage
                surrounding_phage=self.location.GetNumberOfNearbyPhage()
                
                lys_alpha=param["IndividualPhageParameters"]["Lysogeny Alpha"][family]
                lys_kappa=param["IndividualPhageParameters"]["Lysogeny Kappa"][family]
                
                lysogeny_probability=(float(1)/(1+lys_alpha*math.exp(-lys_kappa*surrounding_phage)))                                
                
                    
                if (phage_type=="Temperate") and (random.random()>lysogeny_probability):                    
                    #print "Lyticized with", surrounding_phage, lysogeny_probability                                      
                    self.phages[-1]["Type"]="Induced_Temperate"
                
                else:
                    #It is either defective or temperate going into lysogeny           
                                        
                    #Change genome so that it includes the phage and the genes it carries
                    self.genome.insert(self.phages[-1]["Position"], self.phages[-1])
                                        
                    for gene in cargo:
                        #Check if ARG comes along, add to resistance profile...                
                        if ("AR" in gene[0]) and (gene[3]=="Resistant"):
                            self.resistances[gene[0].split(":")[1]]=True

                        #Check if phage resistance comes along, add to resistance profile...
                        if ("Phage_Res" in gene[0]) and (gene[3]=="Resistant"): 
                            self.phage_resistance=True
                               
                    #Effects of insertion (KO of gene, machinery, plasmids...)  
            return True #Will take phage out of the environment
        else: 
            return False   
        
    def InnerPhageLifeCycle(self, force_it=False):#Returns true if host is going to die, false otherwise
                
        if len(self.phages)==0: 
            return False #No phage, no cry (no death)
        
        #Excision occurs with certain probability.
        for phg in self.phages:
            if phg["Type"]=="Defective":phg["Time"]+=1
            continue 
            #Increase time that phage lived in this genome
            
            if phg["Type"]=="Temperate":
                if random.random()<param["Curing Probability"]:
                    phg["Type"]="Defective"
                    continue #Calculate probability of being "cured" (i.e, deactivated)
                            
                #Calculate probability of entering the lytic cycle
                alpha=param["IndividualPhageParameters"]["Induction Alpha"][phg["Family"]];
                kappa=param["IndividualPhageParameters"]["Induction Kappa"][phg["Family"]] #These parameters control the shape of the induction curve
                                                
                #Antibiotic stress does not apply if bacteria is resistant                
                if not(self.resistances["AR_Rif"]):
                    lysogen_to_lysis_prob_rif=float(1)/(1+alpha*math.exp((-kappa*self.location.rif_conc)))
                else: 
                    lysogen_to_lysis_prob_rif=float(1)/(1+alpha*math.exp((-kappa*0)))
                                                
                if not(self.resistances["AR_Str"]):
                    lysogen_to_lysis_prob_str=float(1)/(1+alpha*math.exp((-kappa*self.location.str_conc)))
                else: 
                    lysogen_to_lysis_prob_str=float(1)/(1+alpha*math.exp((-kappa*0)))
                
                if not(self.resistances["AR_Quin"]):
                    lysogen_to_lysis_prob_quin=float(1)/(1+alpha*math.exp((-kappa*self.location.quin_conc)))
                else: 
                    lysogen_to_lysis_prob_quin=float(1)/(1+alpha*math.exp((-kappa*0)))
                
                #The maximum of all 3 antibiotics is taken (NOT CUMULATIVE!)
                lysogen_to_lysis_prob=max([lysogen_to_lysis_prob_rif, lysogen_to_lysis_prob_str, lysogen_to_lysis_prob_quin])
                
                #Draw the dice...                                
                if random.random()>lysogen_to_lysis_prob:
                    phg["Time"]+=1
                    continue 
                    #Exits anyway without entering the lytic cycle
                else:                                                                            
                    phg["Type"]="Induced_Temperate"
                    #phg["Time"]=0 #So it does not burst in this same generation 
                    # See how to deal with this. If I activate this, when antibiotics are used there is almost no time for bacteria to survive long enough to spread phage 
                    
                    #Check probability of packaging some nearby host genes with the phage DNA
                    
                    spec_transd_prob=param["IndividualPhageParameters"]["Specialized Transduction Probability"][phg["Family"]]
                    if (random.random()<spec_transd_prob): 
                        self.SpecializedTransduction(phg)
                    
                    #It will now enter the lytic cycle (will not enter here in the next iteration)
            
            #If it gets to here is either lytic or induced temperate, it will behave according to the probability of lysis, which is assumed to be high
            if (random.random()<param["Abortive Infection Probability"]) and not(force_it): 
                return True #abortive infection
                         
            if ((random.random()<param["Burst Probability"]) and (phg["Time"]>0)) or force_it:              
                self.PhageBurst(phg)
                return True #Will activate death and generate eDNA
            else: 
                phg["Time"]+=1#Increase time that phage lived in this genome
                   
        return False #If none of these bad things happens, host lives to fight another day 
    
    def SpecializedTransduction(self, lytic_phage=None):
        genomic_distance=param["Specialized Transduction Genomic Distance"]
          
        neighbouring_genes=[pos for pos in range(lytic_phage["Position"]-genomic_distance, lytic_phage["Position"]+genomic_distance+1) if ((pos>=0) and (pos<len(self.genome)) and (pos!=lytic_phage["Position"]) )] #Check boundary 0 of genome
        
        gene=self.genome[random.choice(neighbouring_genes)]#Randomly select one nearby gene
        
        if isinstance(gene, dict):#It means that it is a mobile element            
            gene=random.choice(gene["Cargo"])#Randomly selects one of the genes from the cargo. Here position is not taken into account
                
        if len(gene)>3: #Means it is a resistance gene (ARG or phage resistance)
            phage_new_cargo=(gene[0], gene[1], "Mobile:Ph", gene[3])
        else: 
            phage_new_cargo=(gene[0], gene[1], "Mobile:Ph")
        
        #print "Transduced, Specialized", self.color, gene, neighbouring_genes
        lytic_phage["Cargo"].append(phage_new_cargo) #Append it to phage cargo
        
        
    def GeneralizedTransduction(self, lytic_phage=None):
        genes=random.sample(self.genome, lytic_phage["MaxCargoSize"])
        phage_new_cargo=[]
        for gene in genes:
            if isinstance(gene, dict): #Means it is a mobile element
                gene=random.choice(gene["Cargo"]) #Take a gene at random from the MGE
            if len(gene)>3: #Means it is a resistance gene (ARG or phage resistance)
                phage_new_cargo.append((gene[0], gene[1], "Mobile:Ph", gene[3]))
            else:
                phage_new_cargo.append((gene[0], gene[1], "Mobile:Ph"))
        return phage_new_cargo
    
    
    def PhageBurst(self, lysed_phage=None):

        def CRISPR_Mutation():
            if random.random()<param["Phage CRISPR Mutation Probability"]: 
                return ''.join(random.choice(string.ascii_lowercase) for _ in range(len(lysed_phage["crispr_seq"])))
            else: 
                return lysed_phage["crispr_seq"]
        
        #New in v7.1
        def HostRangeMutation(): #ONLY USE IF THERE IS ONE UNAFFECTED SPECIES
            
            possible_new_hosts=[bac for bac in param["InfectionMap"][lysed_phage["Family"]] if param["InfectionMap"][lysed_phage["Family"]][bac]==0]

            if ((len(possible_new_hosts)>0) and (random.random()<param["Phage Host Range Mutation Probability"])):
                
                print("WARNING: Host range mutation is in beta version, might not work as intended.")
                
                new_host=random.choice(possible_new_hosts)

                new_family=lysed_phage["Family"]+"_Variant"+str(len(param["InfectionMap"]))                
                param["IndividualPhageParameters"]["Burst Size"][new_family]=param["IndividualPhageParameters"]["Burst Size"][lysed_phage["Family"]]
                                
                param["InfectionMap"][new_family]={}
                
                for bac in param["InfectionMap"][lysed_phage["Family"]]:
                    if bac==new_host: param["InfectionMap"][new_family][bac]=0.5
                    else:
                        param["InfectionMap"][new_family][bac]=param["InfectionMap"][lysed_phage["Family"]][bac]
                
                #Update superinfection table                        
                for other_phage in param["Superinfection"]:
                    param["Superinfection"][other_phage][new_family]=param["Superinfection"][other_phage][lysed_phage["Family"]]
                    param["Superinfection"][other_phage][new_family+"_Defective"]=param["Superinfection"][other_phage][lysed_phage["Family"]+"_Defective"]
                
                param["Superinfection"][new_family]={}
                param["Superinfection"][new_family+"_Defective"]={}
                
                for other_phage in param["Superinfection"]:
                    if other_phage==new_family: 
                        param["Superinfection"][new_family][other_phage]=1 #Protection against self
                        param["Superinfection"][new_family][other_phage+"_Defective"]=0            
                        param["Superinfection"][new_family+"_Defective"][other_phage]=0
                        param["Superinfection"][new_family+"_Defective"][other_phage+"_Defective"]=0 #No Protection against defective self
                    elif other_phage!=new_family+"_Defective": 
    
                        param["Superinfection"][new_family][other_phage]=param["Superinfection"][lysed_phage["Family"]][other_phage]
                        param["Superinfection"][new_family+"_Defective"][other_phage]=param["Superinfection"][lysed_phage["Family"]+"_Defective"][other_phage]
                return new_family
            
            else: 
                return lysed_phage["Family"] #Still the same phage as its parent
        
        #New in v7.1
        def ReceptorMutation():          
            if random.random()<param["Phage Receptor Mutation Probability"]:
                all_possible_receptors = set(range(0,param["Possible Receptors"]))   
                     
                #Choose a random one from those      
                new_receptor=random.choice(list(all_possible_receptors - set([lysed_phage["Receptor"]])))
                 
                return new_receptor
            else: 
                return lysed_phage["Receptor"]
        
        #Recombination with other phages
        if len(self.phages)>1: #If there is more than one phage
            for phg1, phg2 in itertools.combinations(self.phages, 2): #Create all pairwise combinations
                 
                if random.random()<param["Phage Recombination Probability"]: # This probability should depend on the phage homology
                    
                    #Select gene positions to be swapped
                    cargo1_gene_pos=random.choice(range(0, len(phg1["Cargo"]))) 
                    cargo2_gene_pos=random.choice(range(0, len(phg2["Cargo"])))

                    new_crg1=[]
                    new_crg2=[]
                    
                    #Recreate cargo of phage1 with the new gene from phage2 (and without the swapped one)
                    for pos in range(0, len(phg1["Cargo"])):
                        new_crg1.append((phg1["Cargo"][pos] if pos!=cargo1_gene_pos else phg2["Cargo"][cargo2_gene_pos]))
                    
                    #Recreate cargo of phage2 with the new gene from phage1 (and without the swapped one)                        
                    for pos in range(0, len(phg2["Cargo"])):
                        new_crg2.append((phg2["Cargo"][pos] if pos!=cargo2_gene_pos else phg1["Cargo"][cargo1_gene_pos]))

                    #Reassign new cargoes
                    phg1["Cargo"]=new_crg1
                    phg2["Cargo"]=new_crg2
        
        
        #If it was induced, reset status to Temperate
        if lysed_phage["Type"]=="Induced_Temperate":
            lysed_phage["Type"]="Temperate"
        
        
        #Check if cargo has reached maximum capacity, due to specialized transduction. Inactivate phage if this is the case
        if len(lysed_phage["Cargo"])>=lysed_phage["MaxCargoSize"]:
            lysed_phage["Type"]="Defective"
        
        
        def MutatePhageGenes(non_mutated_cargo):
            
            def mutate_string(gene_string):
                mutated_pos=random.choice(range(len(gene_string)))
                return gene_string[:mutated_pos] + random.choice(string.letters).upper() + gene_string[mutated_pos+1:]                
            if random.random()<0.5:                
                mutated_pos=random.choice(range(len(non_mutated_cargo)))                
                non_mutated_cargo[mutated_pos]=(mutate_string(non_mutated_cargo[mutated_pos][0]),)+non_mutated_cargo[mutated_pos][1:]                       
            return non_mutated_cargo
            
        
        #Will deposit a number of phages in its current location
        for _ in range(0, param["IndividualPhageParameters"]["Burst Size"][lysed_phage["Family"]]):
            
            #Each virion has a probability to perform generalized transduction             
            gen_transd_prob=param["IndividualPhageParameters"]["Generalized Transduction Probability"][lysed_phage["Family"]]
            
            if (random.random()<gen_transd_prob):   
                                      
                self.location.AddFreePhage(family=lysed_phage["Family"], donor_genome=self, cargo=copy.copy(self.GeneralizedTransduction(lysed_phage)), 
                                   phage_type="Defective",
                                   receptor=lysed_phage["Receptor"],
                                   crispr_sequence=lysed_phage["crispr_seq"], 
                                   maxcargo=lysed_phage["MaxCargoSize"])
                                   # Check if deepcopy is not needed...
                
                
            else:
                #If the phage is intact, there is a probability to mutate its receptor and its CRISPR sequence            
                self.location.AddFreePhage(family=HostRangeMutation(), donor_genome=self, cargo=copy.copy(lysed_phage["Cargo"]),                 
                                   phage_type=lysed_phage["Type"],
                                   receptor=ReceptorMutation(),
                                   crispr_sequence=CRISPR_Mutation(),
                                   maxcargo=lysed_phage["MaxCargoSize"])
                                   # Check if deepcopy is not needed...
                

            
            # If there are memory issues, remove the donor_genome=self, because it will keep around dead bacteria         


# This file is fully translated.   