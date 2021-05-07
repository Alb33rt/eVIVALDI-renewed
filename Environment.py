import random, math, copy, itertools
from Parameters import model_parameters as param

class Place():
    def __init__(self, position, world):
        self.my_world = world

        self.coordinates = position
        self.occupants = []
        
        self.neighbours = {}
        
        #(block, original_genome)
        #block is ARG or "trash"
        self.eDNA = [] #Holds DNA of dead cells (more than one block of DNA). Is degraded at a rate. 
        self.free_phages = []
        
        self.rif_conc=0
        self.rif_time=0 
        self.rif_inoculum=0

        self.str_conc=0
        self.str_time=0
        self.str_inoculum=0
        self.quin_conc=0
        self.quin_time=0
        self.quin_inoculum=0
        
    def Add_Available(self):
        self.my_world.Set_Available(self.coordinates)
    
    def Add_Unavailable(self):
        self.my_world.Set_Unavailable(self.coordinates)
    
    def Get_Distance_To(self, other_location):
        return max(abs(self.coordinates[0]-other_location[0]), abs(self.coordinates[1]- other_location[1]))
  
        
    def Get_Reachable(self, distance, desired_type="Location"):
        new_ds = []

        for dimension in self.coordinates:
            new_d = [dimension + d for d in range(-distance, distance + 1)]
            new_new_d = []
            for d in new_d: 
                if param["World Type"]=="Toroidal":
                    if d<0: 
                        new_new_d.append(self.my_world.max_height-int(math.fabs(d)))
                        continue
                    if d>=self.my_world.max_height: 
                        new_new_d.append(d-self.my_world.max_height)
                        continue
                new_new_d.append(d)
            new_ds.append(new_new_d)

        full_locations=list(itertools.product(*new_ds))
        
        possible_locations=[]
        
        for pos in full_locations:
            if (pos[0]==self.coordinates[0]) and (pos[1]==self.coordinates[1]): 
                continue 
            #skip if it is my own coordinates
            if (pos[0]<0) or (pos[0]>=self.my_world.max_height) or (pos[1]<0) or (pos[1]>=self.my_world.max_width): 
                continue 
            #skip if it is out of boundaries
            possible_locations.append(pos)
    
        if desired_type=="Location": 
            return possible_locations
        if desired_type=="LocationObj":
            return [self.my_world.world[loc] for loc in possible_locations]
        if desired_type=="Bacteria":
            reachable_bacteria=[self.my_world.world[loc].occupants[0] for loc in possible_locations if len(self.my_world.world[loc].occupants)>0]                        
            return reachable_bacteria
    
                
    def AddFreePhage(self, phage_type, family, donor_genome, cargo, crispr_sequence, receptor, maxcargo):
        if len(self.free_phages)>=param["Phage Carrying Capacity"]: 
            print("NO SPACE FOR MORE PHAGE")
            exit()
            return 
        self.free_phages.append((donor_genome, {
            "Family":family, 
            "Position":None, 
            "Cargo": cargo, 
            "Time":0, 
            "Type":phage_type,
            "Receptor":receptor,
            "crispr_seq": crispr_sequence,
            "MaxCargoSize": maxcargo,
            }))
        
    def GetNumberOfNearbyPhage(self):
        n_phage=0
        for l in self.Get_Reachable(1, desired_type="Location"): 
            n_phage+=len(self.my_world.world[l].free_phages)
        n_phage+=len(self.free_phages)
        #Add also phage in the current location
        return n_phage
    
    def AddAntibiotic(self, ant_type, quantity):
        if ant_type=="RIF":
            self.rif_time=0
            self.rif_conc+=quantity
            self.rif_inoculum=self.rif_conc
            
        if ant_type=="STR":
            self.str_time=0
            self.str_conc+=quantity
            self.str_inoculum=self.str_conc
            
        if ant_type=="QUIN":
            self.quin_time=0
            self.quin_conc+=quantity
            self.quin_inoculum=self.str_conc                    
        
    def Diffuse(self, ant_type):

        lambda_diffusion=0.17
        
        if ant_type=="RIF": 
            current_conc=self.rif_conc
        if ant_type=="STR": 
            current_conc=self.str_conc
        if ant_type=="QUIN": 
            current_conc=self.quin_conc
        
        distance=1
                
        deployed_in=[]
        while current_conc>1:
        # CURRENT PROGRESS FOR TRANSLATION
            if ant_type=="RIF": current_conc=int(self.rif_conc*math.exp(-lambda_diffusion*distance))
            if ant_type=="STR": current_conc=int(self.str_conc*math.exp(-lambda_diffusion*distance))
            if ant_type=="QUIN": current_conc=int(self.quin_conc*math.exp(-lambda_diffusion*distance))
            
            locs=self.Get_Reachable(distance, desired_type="Location")
            for l in locs:
                if not(l in deployed_in):                
                    self.my_world.world[l].AddAntibiotic(ant_type, int(random.choice(range(current_conc/2, current_conc))))                
                    deployed_in.append(l)

            distance+=1
                

        
    def Erode_Elements(self, current_iteration):
        
        #Degrade antibiotics                
        if self.rif_conc>0:
            antibiotic_degradation_rate=param["Antibiotic Degradation Lambda RIF"]            
            self.rif_conc=self.rif_inoculum*math.exp(-antibiotic_degradation_rate*self.rif_time)
            self.rif_time+=1
            if self.rif_conc<pow(10, -5): self.rif_conc=0 #To save some time
            
        if self.str_conc>0:
            antibiotic_degradation_rate=param["Antibiotic Degradation Lambda STR"]            
            self.str_conc=self.str_inoculum*math.exp(-antibiotic_degradation_rate*self.str_time)
            self.str_time+=1
            if self.str_conc<pow(10, -5): self.str_conc=0 #To save some time
                        
        if self.quin_conc>0:
            antibiotic_degradation_rate=param["Antibiotic Degradation Lambda QUIN"]            
            self.quin_conc=self.quin_inoculum*math.exp(-antibiotic_degradation_rate*self.quin_time)
            self.quin_time+=1
            if self.quin_conc<pow(10, -5): self.quin_conc=0 #To save some time
                                        
        
        #Degrade phage
        def CheckPhageSurvival(phage):
            prob_phage_survival=1*math.exp(-param["Free Phage Degradation Lambda"]*(-phage["Time"]))                    
            if random.random()>prob_phage_survival: return False #Doesn't make it
            else: return True #Still survives
        
        #Check which phage survive outside
        self.free_phages[:]=[p for p in self.free_phages if CheckPhageSurvival(p[1])]
        for phage in self.free_phages: phage[1]["Time"]-=1#Then increase counter of surviving phage lifetime outside a host   
                    
        #Same for exogenous DNA...        
        if len(self.eDNA)>0:
            if random.random()>0.5: #TODO: This more efficiently (as a function of "age" of eDNA?)
                self.eDNA.pop(random.randrange(len(self.eDNA)))
                
        

class World():
    def __init__(self, lines, rows):
        self.max_width=lines
        self.max_height=rows        
        
        self.world={}
        self.available_locations=[]

        print("Setting world...")
        for i in range(0, lines):
            for j in range(0, rows): 
                self.world[(i,j)]=Place((i,j), self)
                self.available_locations.append((i,j))
      
     
    def GetMooreNeighborhood(self, origin, distance):
        adjs=[]
         
        for i in range(-distance,distance+1,1):
            for j in range(-distance,distance+1,1):
                if i != 0 or j != 0: #Do not add self as neighbor
                    if self.world.has_key((origin[0]+i, origin[1]+j)):
                        adjs.append((origin[0]+i, origin[1]+j))
        return adjs
    
         
    
    def Set_Available(self, coordinates):
        if coordinates in self.available_locations: 
            print("Bronca SET")
            exit()
        self.available_locations.append(coordinates)
        
    def Set_Unavailable(self, coordinates):
        if not(coordinates in self.available_locations): 
            print("Bronca UNSET")
            exit()
        self.available_locations.remove(coordinates)

        
    def Get_Free_Space(self):return self.world[random.choice(self.available_locations)]
    def Get_All_Free_Spaces(self):return self.available_locations

    def Set_All_Unavailable(self):
        self.available_locations=[]
    
    def Set_Many_Unavailable(self, coordinate_list):         
        self.available_locations=list(set(self.available_locations)-set(coordinate_list)) #Careful with this. It changes order of list and removes duplicates. Might work for how we manage this list
        
    def Set_Many_Available(self, coordinate_list):self.available_locations.extend(coordinate_list)
    def Set_All_Available(self): self.available_locations=self.world.keys()
    
    def GetAllFreePhages(self):
        count_all_phages=0;
        for loc in self.world:count_all_phages+=len(self.world[loc].free_phages)               
        return count_all_phages
    
    def UpdateLocations(self, current_iteration=None):
        #Degrade elements (eDNA, antibiotics...)
        for _, loc in enumerate(self.world): self.world[loc].Erode_Elements(current_iteration)
        
        #Shuffle according to environmental Structure
        if param["Environmental Diffusion"]=="Liquid":   
                       
            #Redistribute antibiotic randomly across locations
            all_rifs=[self.world[loc].rif_conc for loc in self.world]#Gather all individual concentrations in array
            random_positions=random.sample(self.world.keys(), len(self.world.keys()))#Randomize locations            
            for loc in random_positions:self.world[loc].rif_conc=all_rifs.pop()#Assign new antibiotic values
            
            
            #PHAGE SHUFFLING - randomize each viral particle
            all_phages=[] #Remove all phages from positions and collect all in a single list
            for loc in self.world:all_phages.extend(copy.copy(self.world[loc].free_phages));self.world[loc].free_phages=[]
                        
            factor_chunk=int(round(float(len(all_phages))/(self.max_height*self.max_width)))#This is how many phages each location is going to "receive"
    
    
            if factor_chunk>=1:#When there are too many phages, divide evenly (faster)                                    
                random.shuffle(all_phages)
                chunks=[all_phages[x:x+factor_chunk] for x in range(0, len(all_phages), factor_chunk)]
                
                for loc in random_positions:
                    if len(chunks)<=0:break
                    else:self.world[loc].free_phages=chunks.pop()
                    
                while len(chunks)>0: #If there are remainders, distribute randomly                                                        
                    random_pos=random.choice(random_positions)                    
                    self.world[random_pos].free_phages.extend(chunks.pop())

    
            else:#If not enough phages, distribute at random by the locations                
                for loc in random_positions:
                    if len(all_phages)<=0: break
                    self.world[loc].free_phages.append(all_phages.pop())                                 
            
            
        elif param["Environmental Diffusion"]=="Semi-Solid": 
            distance_to_use=3
            
            #Create dictionary that holds the future (i.e., after diffusion) values (so that it does not chain into more changes than it should per generation)
            new_values={}
            for loc in self.world:new_values[loc]={"rif":-1,"phage":[]}
            
            #Exchange value of antibiotic with a random neighbour
            for loc in self.world:
                
                partner_loc=random.choice(self.world[loc].Get_Reachable(distance=distance_to_use, desired_type="LocationObj"))                
                new_values[loc]["rif"]=partner_loc.rif_conc;
                new_values[partner_loc.coordinates]["rif"]=self.world[loc].rif_conc


                if len(self.world[loc].free_phages)>0:
                    partner_loc=self.world[loc].Get_Reachable(distance=distance_to_use, desired_type="LocationObj")
                    partner_loc.append(self.world[loc])#Add self as possible recipient
                    
                    while len(self.world[loc].free_phages)>0:
                        new_loc=random.choice(partner_loc)
                        new_values[new_loc.coordinates]["phage"].append(self.world[loc].free_phages.pop())                        
                        
            #After everything is calculated, update values
            for loc in self.world:
                self.world[loc].rif_conc=new_values[loc]["rif"]             
                self.world[loc].free_phages=copy.copy(new_values[loc]["phage"])#TODO: Check if deepcopy is needed
