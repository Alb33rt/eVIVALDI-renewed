import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib import gridspec
import numpy as np
import os

from Parameters import model_parameters as param

sns.set(font_scale=1.5)
sns.set_style("white")

#fig = plt.figure(figsize=(23,12))
fig = plt.figure(figsize=(18,9))

#gs = gridspec.GridSpec(4, 4)
gs = gridspec.GridSpec(4, 6)

#Bacteria stats
bacteria_zone = fig.add_subplot(gs[0,0:2])
total_bac_pop, = bacteria_zone.plot([], '-', color='black', label="Bacteria Pop. Size")
AR_bac_pop, = bacteria_zone.plot([], '-', color='yellow', label="Antibiotic Resistant")


#Phage stats
phage_zone = fig.add_subplot(gs[0,2:4])
total_phage_pop, = phage_zone.plot([], '-', color='black', label="Phage Pop. Size")
AR_phage_pop, = phage_zone.plot([], '-', color='red', label="Phage Carrying ARG")
defective_phage_pop, = phage_zone.plot([], '--', color='green', label="Defective Phage")

#Spatio-visual representation  
spatial_zone = fig.add_subplot(gs[1:3,0:4])

#Bacterial subspecies
sub_species_bacteria_zone = fig.add_subplot(gs[3,0:2])
all_bac_species=[]
for index, bac_species in enumerate(param["Bacteria Species"]):
    temp_line, = sub_species_bacteria_zone.plot([], '-', color=param["Bacteria Species"][bac_species], label=bac_species)
    all_bac_species.append(temp_line)

#Phage receptors and resistance
phage_receptors_zone = fig.add_subplot(gs[3,2:4])
all_receptors_bac=[]
all_receptors_phage=[]
colors_maps={0: "grey", 1:"green", 2:"yellow", 3:"blue", 4:"purple", 5:"orange", 6:"red", 7:"brown", 8:"black", 9:"pink"}

for index in range(0,param["Possible Receptors"]):
    temp_line_bac, = phage_receptors_zone.plot([], '-', color=colors_maps[index], label="Rec "+str(index))
    temp_line_phage, = phage_receptors_zone.plot([], '--', color=colors_maps[index], label="Rec "+str(index))
    all_receptors_bac.append(temp_line_bac)
    all_receptors_phage.append(temp_line_phage)



#Genomes
genomes_zone = fig.add_subplot(gs[0:5,4:6])

gs.update(wspace=0.35, hspace=0.35)


def LaunchPlot():
    
    max_bacteria=(param["World Size"]*param["World Size"])*1.5
    max_phage=(param["World Size"]*param["World Size"]*param["Phage Carrying Capacity"])*1.5
    
    bacteria_zone.set_yscale('log');bacteria_zone.set_ylim(pow(10, 0), 1.1*max_bacteria); bacteria_zone.set_xlim(-param["Iterations Outgrowth"], param["Iterations"])    
    phage_zone.set_yscale('log');phage_zone.set_ylim(pow(10, 0), 1.1*max_phage); phage_zone.set_xlim(-param["Iterations Outgrowth"], param["Iterations"])
    sub_species_bacteria_zone.set_yscale('log');sub_species_bacteria_zone.set_ylim(pow(10, 0), 1.1*max_bacteria); sub_species_bacteria_zone.set_xlim(-param["Iterations Outgrowth"], param["Iterations"])
    phage_receptors_zone.set_yscale('log');phage_receptors_zone.set_ylim(pow(10, 0), 1.1*max_phage); phage_receptors_zone.set_xlim(-param["Iterations Outgrowth"], param["Iterations"])
    
    img=mpimg.imread("Empties/EmptyWorld200.bmp")
    img_genomes=mpimg.imread("Empties/EmptyWorld200.bmp")
    
    global img_artist, img_artist2 

    img_artist = spatial_zone.imshow(img, interpolation="nearest", aspect='auto')
    spatial_zone.axis('off')
    
    img_artist2 = genomes_zone.imshow(img_genomes, aspect='auto')
    genomes_zone.axis('off')

    plt.ion()
    plt.show()
    plt.pause(0.01)


def ResetPlots():
    
    #Bacteria panel
    total_bac_pop.set_xdata([])
    AR_bac_pop.set_xdata([])
    
    for index, bac_species in enumerate(param["IndividualBacteriaParameters"]["Basal Reproduction Probability"]):        
        all_bac_species[index].set_xdata([])
    
    total_bac_pop.set_ydata([])
    AR_bac_pop.set_ydata([])
    
    #Bacteria-species panel
    for index, bac_species in enumerate(param["IndividualBacteriaParameters"]["Basal Reproduction Probability"]):
        all_bac_species[index].set_ydata([])
    
    #Phage panel
    total_phage_pop.set_xdata([])
    AR_phage_pop.set_xdata([])
    defective_phage_pop.set_xdata([])
    
    total_phage_pop.set_ydata([])                              
    AR_phage_pop.set_ydata([]) 
    defective_phage_pop.set_ydata([])
        

    for index in range(0,param["Possible Receptors"]):
        all_receptors_bac[index].set_xdata([])
        all_receptors_phage[index].set_xdata([])
        all_receptors_bac[index].set_ydata([])
        all_receptors_phage[index].set_ydata([])
    

    bacteria_zone.figure.canvas.draw();phage_zone.figure.canvas.draw(); spatial_zone.figure.canvas.draw(); sub_species_bacteria_zone.figure.canvas.draw(); phage_receptors_zone.figure.canvas.draw()
        
    plt.pause(0.001)

    

def UpdateStats(iteration, simulation, all_stats, extra_stats, print_plots=False):    
    
    #Bacteria panel
    total_bac_pop.set_xdata(np.append(total_bac_pop.get_xdata(), iteration))
    AR_bac_pop.set_xdata(np.append(AR_bac_pop.get_xdata(), iteration))
    
    for index, bac_species in enumerate(param["IndividualBacteriaParameters"]["Basal Reproduction Probability"]):        
        all_bac_species[index].set_xdata(np.append(all_bac_species[index].get_xdata(), iteration))
    
    total_bac_pop.set_ydata(np.append(total_bac_pop.get_ydata(), all_stats["Bacteria Population Size"]))
    AR_bac_pop.set_ydata(np.append(AR_bac_pop.get_ydata(), all_stats["Bacteria ARG Population Size"]))
    
    #Bacteria-species panel
    for index, bac_species in enumerate(param["IndividualBacteriaParameters"]["Basal Reproduction Probability"]):
        all_bac_species[index].set_ydata(np.append(all_bac_species[index].get_ydata(), all_stats["Bacteria Subpopulation Size"][bac_species]))
    
    #Phage panel
    total_phage_pop.set_xdata(np.append(total_phage_pop.get_xdata(), iteration))
    AR_phage_pop.set_xdata(np.append(AR_phage_pop.get_xdata(), iteration))
    defective_phage_pop.set_xdata(np.append(defective_phage_pop.get_xdata(), iteration))
    
    total_phage_pop.set_ydata(np.append(total_phage_pop.get_ydata(), all_stats["Phage Population Size"]))
    AR_phage_pop.set_ydata(np.append(AR_phage_pop.get_ydata(), all_stats["Phage With ARG Population Size"]))  
    defective_phage_pop.set_ydata(np.append(defective_phage_pop.get_ydata(), sum([all_stats["Phage Subpopulation Size"][p] for p in all_stats["Phage Subpopulation Size"] if "_Defective" in p])))
        

    for index in range(0,param["Possible Receptors"]):
        all_receptors_bac[index].set_xdata(np.append(all_receptors_bac[index].get_xdata(),iteration))
        all_receptors_phage[index].set_xdata(np.append(all_receptors_phage[index].get_xdata(),iteration))
        all_receptors_bac[index].set_ydata(np.append(all_receptors_bac[index].get_ydata(), extra_stats["All Bacteria With Receptor"][index]))
        all_receptors_phage[index].set_ydata(np.append(all_receptors_phage[index].get_ydata(), extra_stats["Phage With Receptor"][index]))
    
    #Spatial zone update    
    if iteration>-param["Iterations Outgrowth"]:        
        img=mpimg.imread(param["Main Directory"]+"Joined/Generation"+str(iteration).zfill(3)+".bmp")
        img_artist.set_data(img)
    
    #Genomes zone update    
    if os.path.isfile(param["Main Directory"]+"Genomes/"+str(param["Random Seed"])+"_"+str(simulation)+"_Iter"+str(iteration)+"_Loci.png"):
        img_artist2.set_data(mpimg.imread(param["Main Directory"]+"Genomes/"+str(param["Random Seed"])+"_"+str(simulation)+"_Iter"+str(iteration)+"_Loci.png"))

    
    #Define legends
    h0, l0 = bacteria_zone.get_legend_handles_labels()
    legend=bacteria_zone.legend(h0, l0, bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, borderaxespad=0.)
    for label in legend.get_texts():label.set_fontsize(12)
    for label in legend.get_lines():label.set_linewidth(0.75)
    
    
    h0, l0 = phage_zone.get_legend_handles_labels()
    legend=phage_zone.legend(h0, l0, bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, borderaxespad=0.)
    for label in legend.get_texts():label.set_fontsize(12)
    for label in legend.get_lines():label.set_linewidth(0.75)
    
    
    h0, l0 = sub_species_bacteria_zone.get_legend_handles_labels()
    legend=sub_species_bacteria_zone.legend(h0, l0, bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, borderaxespad=0.)
    for label in legend.get_texts():label.set_fontsize(12)
    for label in legend.get_lines():label.set_linewidth(0.75)
    
        
    #h0, l0 = phage_receptors_zone.get_legend_handles_labels()
    #legend=phage_receptors_zone.legend(h0, l0, bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=10, borderaxespad=0.)
    #for label in legend.get_texts():label.set_fontsize(6)
    #for label in legend.get_lines():label.set_linewidth(0.3)
    
    
    bacteria_zone.figure.canvas.draw();phage_zone.figure.canvas.draw(); spatial_zone.figure.canvas.draw(); sub_species_bacteria_zone.figure.canvas.draw(); phage_receptors_zone.figure.canvas.draw()
        
    plt.pause(0.001)
    
    if print_plots: plt.savefig(param["Main Directory"]+"PlotDynamics/Dynamics"+str(iteration).zfill(3)+'.png')
        

#plt.show()

