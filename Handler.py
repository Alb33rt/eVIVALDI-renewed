
import Organism
import ParseSpecies
import Parameters
import Stats
import os
import string
import random
import pprint
import numpy as np

from Environment import World
from Parameters import model_parameters as param
from Visual import WorldPainter, PrintGenomes
from Environment import World
from copy import copy, deepcopy

##Check input parameters
import argparse

parser = argparse.ArgumentParser(description='Phageology')
parser.add_argument(
    '-id', dest='id', help='The identification for this simulation', required=True)
parser.add_argument('-paramfile', dest='paramfile',
                    help='The file with the parameters for this simulation', required=True)
parser.add_argument('-setupfile', dest='setupfile',
                    help='The file with simulation setup', required=True)
parser.add_argument('--visualize', action='store_true',
                    help='Include this flag if you want plots to appear in realtime')

args = parser.parse_args()

print("Parameter file:" + args.paramfile)
print("Setup file:" + args.setupfile)

Parameters.LoadParameters(args.paramfile)


#Initialize random machine with seed
if param["Random Seed"] != None:
    random.seed(param["Random Seed"])
    np.random.seed(param["Random Seed"])


#Create subdirectories that hold the printed structures
if not os.path.exists(param["Main Directory"]+"/Genomes/"):
    print("Directory missing (" +
          param["Main Directory"]+"/Genomes/"+"), creating...")
    os.makedirs(param["Main Directory"]+"/Genomes/")
if not os.path.exists(param["Main Directory"]+"/Bacteria/"):
    print("Directory missing (" +
          param["Main Directory"]+"/Bacteria/"+"), creating...")
    os.makedirs(param["Main Directory"]+"/Bacteria/")
if not os.path.exists(param["Main Directory"]+"/World/"):
    print("Directory missing (" +
          param["Main Directory"]+"/World/"+"), creating...")
    os.makedirs(param["Main Directory"]+"/World/")
if not os.path.exists(param["Main Directory"]+"/Joined/"):
    print("Directory missing (" +
          param["Main Directory"]+"/Joined/"+"), creating...")
    os.makedirs(param["Main Directory"]+"/Joined/")
if not os.path.exists(param["Main Directory"]+"/PlotDynamics/"):
    print("Directory missing (" +
          param["Main Directory"]+"/PlotDynamics/"+"), creating...")
    os.makedirs(param["Main Directory"]+"/PlotDynamics/")

#Accessory function - Roullette wheel selection


def WeightedReproductionChoice(bacteria_list):

    def weighted_choice(choices):
        total = sum(w for c, w in choices)
        if total == 0:
            return None
        r = random.uniform(0, total)
        upto = 0
        for c, w in choices:
            if upto + w >= r:
                return (c, w)
            upto += w

    candidates = [(b, b.AssertCandidacy()) for b in bacteria_list]
    candidates = filter(lambda x: x[1] != None, candidates)

    result = weighted_choice(candidates)

    if (result != None) and (len(result[0].phages) > 0):
        for p in result[0].phages:
            if p["Type"] == "Virulent":
                print("You should not be here")
                exit()

    if result == None:
        return None
    else:
        return result[0]

#Accessory function - Logic check for bacteria survival


def IsDEAD(bacteria):
    if bacteria.location == None:
        return False
    else:
        return True


def AddPhages(future_phage, world):
    for pos in future_phage:
        phg = future_phage[pos]
        print("Adding phage", phg["Type"], phg["Family"])
        world.world[pos].AddFreePhage(
            phage_type=copy(phg["Type"]),
            family=copy(phg["Family"]),
            donor_genome=None,
            cargo=copy(phg["Cargo"]),
            crispr_sequence=copy(phg["crispr_seq"]),
            receptor=copy(phg["Receptor"]),
            maxcargo=copy(phg["MaxCargoSize"])
        )


sims = param["Number Simulations"]

for sim in range(sims):
    #First set up the world
    size = param["World Size"]
    w = World(lines=size, rows=size)

    #Second, let's create a visualizer
    vanGogh = WorldPainter(
        main_location=param["Main Directory"], world_size=size)

    #Third, Define bacterial genomes (along with which genes are potential ARG)
    def DefineGenome(block):

        #Define main DNA
        genome_type = [(random.choice([block, block, block, block, "N"])
        (random.uniform(param["Minimum Gene Size"],
        param["Maximum Gene Size"])), "Core") for _ in range(param["Number Genes"])]

        #Define antibiotic resistance genes (x3)
        for resistance_allele in ["Rif", "Str", "Quin"]:
            ar_locus = random.choice([pos for pos in range(0, len(genome_type)) if (
                (genome_type[pos][0] != "N") and (not("AR" in genome_type[pos][0])))])
            genome_type[ar_locus] = (genome_type[ar_locus][0]+":AR_" +
                                     resistance_allele, genome_type[ar_locus][1], "Core", "Sensitive")

        print(genome_type)
        return genome_type

    ##Get population composition from input file
    setup_bacs, setup_phgs, setup_inter, setup_superinfection = ParseSpecies.LoadSpeciesStats(
        args.setupfile)

    ##Setup phage
    initial_phages = []
    future_phage = {}
    param["IndividualPhageParameters"] = {}
    param["IndividualPhageParameters"]["Burst Size"] = {}
    param["IndividualPhageParameters"]["Lysogeny Alpha"] = {}
    param["IndividualPhageParameters"]["Lysogeny Kappa"] = {}
    param["IndividualPhageParameters"]["Induction Alpha"] = {}
    param["IndividualPhageParameters"]["Induction Kappa"] = {}
    param["IndividualPhageParameters"]["Generalized Transduction Probability"] = {}
    param["IndividualPhageParameters"]["Specialized Transduction Probability"] = {}

    genomes_colors_map = {}

    for p in setup_phgs:
        genomes_colors_map[p["DNAIdentifier"]] = (p["color_genome"], "Phage")

        if "ARG" in p["Name"]:
            ar_type = p["Name"].split("ARG_")[-1]
            initial_phages.append({
                "Type": p["Type"],
                "Position": None,
                "Family": p["Name"],
                "crispr_seq": ''.join(random.choice(string.ascii_lowercase) for _ in range(10)),
                "Cargo": [(p["DNAIdentifier"], 10, 'Mobile:Ph'),
                        ('A:AR_'+ar_type, 9, 'Core', 'Resistant', 'Mobile:Ph')],
                "Receptor": p["Receptor"],
                "MaxCargoSize": p["MaxCargoSize"],
                "Time": 10
            })

        else:
            initial_phages.append({
                "Type": p["Type"],
                "Position": None,
                "Family": p["Name"], "crispr_seq": ''.join(random.choice(string.ascii_lowercase) for _ in range(10)),
                "Cargo": [(p["DNAIdentifier"], 10, 'Mobile:Ph'),
                        (p["DNAIdentifier"], 20, 'Mobile:Ph')],
                "Receptor": p["Receptor"],
                "MaxCargoSize": p["MaxCargoSize"],
                "Time": 0,
            })

        param["IndividualPhageParameters"]["Burst Size"][p["Name"]] = p["BurstSize"]
        param["IndividualPhageParameters"]["Lysogeny Alpha"][p["Name"]
            ] = p["Lysogeny Alpha"]
        param["IndividualPhageParameters"]["Lysogeny Kappa"][p["Name"]
            ] = p["Lysogeny Kappa"]
        param["IndividualPhageParameters"]["Induction Alpha"][p["Name"]
            ] = p["Induction Alpha"]
        param["IndividualPhageParameters"]["Induction Kappa"][p["Name"]
            ] = p["Induction Kappa"]
        param["IndividualPhageParameters"]["Generalized Transduction Probability"][p["Name"]
            ] = p["Generalized Transduction Probability"]
        param["IndividualPhageParameters"]["Specialized Transduction Probability"][p["Name"]
            ] = p["Specialized Transduction Probability"]

        for _ in range(0, p["FreeNumbers"]):

            #Clustered in the middle of the environment
            #pos=(int(random.uniform(int(param["World Size"]/2)-5, int(param["World Size"]/2)+5)),int(random.uniform(int(param["World Size"]/2)-5, int(param["World Size"]/2)+5)))

            #OR Randomly sparsed throughout the environment
            pos = (int(random.uniform(0, param["World Size"])), int(
                random.uniform(0, param["World Size"])))

            #Save to add at a later point
            future_phage[pos] = initial_phages[-1]

    ##Setup the infection map (host range)
    param["InfectionMap"] = setup_inter

    ##Setup the superinfection protection amongst phage
    param["Superinfection"] = setup_superinfection

    ##Setup bacteria
    genome_types = []
    colors_used = []
    bacs = []

    param["Bacteria Species"] = {}
    param["IndividualBacteriaParameters"] = {}
    param["IndividualBacteriaParameters"]["Basal Reproduction Probability"] = {}

    for b in setup_bacs:
        if not(b["Color"] in colors_used):
            colors_used.append(b["Color"])
            genome_types.append((DefineGenome(b["DNAIdentifier"]), b["Color"]))
            param["Bacteria Species"][b["Name"]] = b["Color"]
            genomes_colors_map[b["DNAIdentifier"]] = (b["Color"], "Bacteria")

        genome_to_use = None
        for gen in genome_types:
            if gen[1] == b["Color"]:
                genome_to_use = gen[0]
                param["IndividualBacteriaParameters"]["Basal Reproduction Probability"][b["Name"]] = b["Growth"]

        for b_individual in range(0, int(param["World Size"]*param["World Size"]*b["Freq"])):
            if b["Color"] == None:
                print(b_individual, b)
                new_bacteria = Organism.Bacteria(b["Name"], copy(
                    genome_to_use), b["Color"], ar_locus=[0, 0, 0], movement=1)

            if b["AntRes"] != "NA":
                for res in b["AntRes"]:
                    for idx, gene in enumerate(new_bacteria.genome):
                        if "AR_"+res in gene[0]:
                            new_bacteria.genome[idx] = (
                                gene[0], gene[1], gene[2], "Resistant")
                        new_bacteria.resistances[gene[0].split(":")[1]] = True

            if b["PhageRes"][0] != "NA":
                new_bacteria.phage_receptor_resistances = [
                    int(recept) for recept in b["PhageRes"]]

            if b["Prophage"] != None:
                for phg in b["Prophage"]:
                    for phg_template in initial_phages:
                        if phg_template["Family"] == phg[0]:
                            #print arg_pos
                            new_bacteria.phages.append({
                                "Type": copy(phg_template["Type"]),
                            "Position": copy(phg[1]),
                            "Family": copy(phg_template["Family"]),
                            "crispr_seq": copy(phg_template["crispr_seq"]),
                            "Cargo": copy(phg_template["Cargo"]),
                            "Receptor": copy(phg_template["Receptor"]),
                            "MaxCargoSize": copy(phg_template["MaxCargoSize"]),
                            "Time": copy(phg_template["Time"])})
                            new_bacteria.genome.insert(phg[1], new_bacteria.phages[-1])
                            for gene in phg_template["Cargo"]:
                                #Check if ARG comes along, add to resistance profile...
                                if ("AR" in gene[0]) and (gene[3] == "Resistant"):
                                    new_bacteria.resistances[gene[0].split(":")[
                                                                           1]] = True
            bacs.append(new_bacteria)
    
    
    print("Setup done. Here are the parameters to be simulated")
    pprint.pprint(param)
    
    print("Positioning bacteria...")
    print(len(bacs))
    all_vacancies=deepcopy(w.Get_All_Free_Spaces()) #Gather empty spots (deepcopy, because we are changing it)
    random.shuffle(all_vacancies)
    random.shuffle(bacs) 
    #Shuffle everything (places, bacteria)

    for idx in range(len(bacs)):bacs[idx].Set_Location(w.world[all_vacancies[idx]], set_un=False) 
    #Assign each bacteria to a place (sequential, now that everything is shuffled) 

    w.Set_Many_Unavailable(all_vacancies[0:len(bacs)]) 
    #Mark assigned places as unavailable (the first part of the location list, since we assigned them sequentially)

    #w.Set_All_Unavailable() #Comment previous line and uncomment this if all lattice is filled (faster like this)
    
    #Fifth, and finally the life cycle
    iterations=param["Iterations"]
    iterations_outgrowth=param["Iterations Outgrowth"]
    
    #Real time plots    
    if args.visualize: 
        import VisualStats
        VisualStats.LaunchPlot()
        VisualStats.ResetPlots()
    
    phage_diversity={}    
    bac_diversity={}
    
    ####PART 1: Getting the transducing particles        
    for iteration in range(iterations_outgrowth,iterations):

        #Add phages if outgrowth is finished
        if iteration==0:
            print("Adding phages to the environment...")
            AddPhages(future_phage, w)
                
        #Log genomes (if required)
        if iteration in param["Log Sample Genomes"]:            
            
            #Select 30 random genomes to print
            to_print=random.sample(bacs, min(30, len(bacs)))
            PrintGenomes([b.genome for b in to_print], genomes_colors_map, param["Main Directory"]+"Genomes/"+args.id+"_"+str(param["Random Seed"])+"_"+str(sim)+"_Iter"+str(iteration))
        
        #Log population stats(if required)        
        if (iteration==0) or (iteration%param["Log Interval"]==0):            
            sts=Stats.GetStats(bacs, w)
            e_sts=Stats.GetExtraStats(bacs, w)     
            Stats.AppendStats(sim, iteration, sts, e_sts, param["Main Directory"]+"OutStats_"+args.id+".txt")
                                   
        #Apply Antibiotics - Do this before visual plots to see the dynamics (a lot will happen after antibiotics are applied)
        ant_to_apply=[]
        if iteration in param["Antibiotic Times RIF"]: 
            ant_to_apply.append("RIF")
        if iteration in param["Antibiotic Times STR"]: 
            ant_to_apply.append("STR")
        if iteration in param["Antibiotic Times QUIN"]: 
            ant_to_apply.append("QUIN")
                
        #if iteration in param["Antibiotic Times"]:
        for ant in ant_to_apply:
            
            #Skip if there are no antibiotics
            if param["Max Antibiotic Concentration "+ant]==0: 
                continue
            
            #Uniformly, if environment is liquid
            if (param["Environmental Diffusion"]=="Liquid") or (param["Antibiotic Exposure Structured Environments"]=="Homogeneous"):                        
                for loc_x in range(0, param["World Size"]):
                    #for loc_y in range(0, param["World Size"]):w.world[loc_x, loc_y].AddAntibiotic("RIF",param["Max Antibiotic Concentration"])
                    for loc_y in range(0, param["World Size"]):
                        w.world[loc_x, loc_y].AddAntibiotic(ant,param["Max Antibiotic Concentration "+ant])  
            else: 
                #Hetereogeneously in the environment, otherwise
                                           
                #Choose random spots to apply drug
                spread=2
                random_x=random.sample(range(0, param["World Size"]), int(param["World Size"]/spread))
                random_y=random.sample(range(0, param["World Size"]), int(param["World Size"]/spread))
                   
                for pos_index in range(0, int(param["World Size"]/spread)):
                    antibiotic_pos=(random_x[pos_index],random_y[pos_index])                                      
                    w.world[antibiotic_pos[0], antibiotic_pos[1]].AddAntibiotic(ant,random.choice(range(param["Max Antibiotic Concentration "+ant]/2,param["Max Antibiotic Concentration "+ant])))
                    w.world[antibiotic_pos[0], antibiotic_pos[1]].Diffuse(ant)
                    
                    
        print("(SIM "+str(sim)+"_"+args.id+") Iteration "+ str(iteration)+", with "+str(len(bacs))+" bacteria and "+str(w.GetAllFreePhages())+" free phage...")
        
        if args.visualize: 
            vanGogh.Paint(w.world, iteration)
        
        #Real time plots
        if args.visualize: 
            sts=Stats.GetStats(bacs, w)
        extra_sts=Stats.GetExtraStats(bacs, w)
        VisualStats.UpdateStats(iteration, sim, sts, extra_sts, print_plots=True)
            
        print("\tShuffling bacteria access...")
        random.shuffle(bacs) #Randomize order of assessing each cell
        
        #Death probability
        print("\tDeath...")
        empty_spots=[] #Will hold the locations of dead bacteria
        for b in bacs:
            if random.random()<b.Get_Death_Probability(): 
                empty_spots.append(b.Death())
        
        print("\tRemoving "+str(len(empty_spots))+" dead cells...")
        bacs=filter(IsDEAD, bacs)
        w.Set_Many_Available(empty_spots)

        print("\tEnvironmental phage infection...")
        for b in bacs: 
            b.EnvironmentalInfection()
        
        #Reproduction/movement into newly open spaces]
        for b in bacs:
            b.just_reproduced = False
        
        print("\tReproducing and repopulating...")
        new_bacteria_list=[]
        used_locations_list=[]
        empty_locations=deepcopy(w.Get_All_Free_Spaces()) #Need to do deepcopy because the world list will be changed below (messes up the loop)
        random.shuffle(empty_locations)
        
        if len(empty_locations)>0:        
            for idx, location in enumerate(empty_locations):
                if (param["Environmental Diffusion"]=="Liquid") or (param["Environmental Diffusion"]=="Solid"):
                    candidates=w.world[location].Get_Reachable(1, desired_type="Bacteria")
                else: #Semi-solid case
                    candidates=w.world[location].Get_Reachable(3, desired_type="Bacteria") 
                    
                if len(candidates)==0: continue #No cell might be reachable
                selected_candidate=WeightedReproductionChoice(candidates)
                if selected_candidate!=None:                    
                    new_bacteria_list.append((selected_candidate.Reproduce(), location))
                    used_locations_list.append(location)

        
        for idx, new_b in enumerate(new_bacteria_list):
            new_b[0].Set_Location(w.world[new_b[1]], set_un=False)
            bacs.append(new_b[0])
        w.Set_Many_Unavailable(used_locations_list)
        
        #Disperse bacteria if environment is well-mixed
        if param["Environmental Diffusion"]=="Liquid":
            #Shuffle
            w.Set_All_Available()
            if len(bacs)>0:
                all_vacancies=deepcopy(w.Get_All_Free_Spaces())
                random.shuffle(all_vacancies)
                random.shuffle(bacs)
                for idx, bac in enumerate(bacs): 
                    bac.Leave_Current_Location(set_a=False)
                    bac.Set_Location(w.world[all_vacancies[idx]], set_un=False)        
                if idx<len(all_vacancies)-1:
                    w.Set_Many_Unavailable(all_vacancies[0:idx+1])
                else: 
                    w.Set_All_Unavailable()
        
        #Degrade antibiotics, phages...
        print("\tEnvironmental changes...")
        w.UpdateLocations()

        # This file has been completely translated.
