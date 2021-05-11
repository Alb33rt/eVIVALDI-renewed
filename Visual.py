import PIL
from PIL import Image, ImageDraw
import matplotlib
import numpy
import random
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from Parameters import model_parameters as param

color_map = {"black": (0, 0, 0),
           "white": (255, 255, 255),
           "red": (125, 0, 0),
           'purple': (153, 50, 204),
           "green": (0, 125, 0),
           "lightgreen": (162, 205, 162),
           "yellow": (255, 255, 0),
           "blue": (0, 0, 255),
           "light blue": (0, 0, 255),
           "grey": (128, 128, 128),
           "lightgrey": (215, 215, 215),
           "aqua": matplotlib.colors.cnames["aqua"],
           "rosybrown": matplotlib.colors.cnames["rosybrown"]}

max_ant_conc = 15

crispr_colors = {}


class WorldPainter():
    def __init__(self, main_location, world_size):
        self.world_size = world_size
        self.main_location = main_location
        self.empty_world = "Empties/EmptyWorld"+str(self.world_size)+".bmp"
        self.colors_generated = list(plt.prism(numpy.linspace(0, 1, 1000)))  # @UndefinedVariable
        random.shuffle(self.colors_generated)

    def Paint(self, world, generation):

        world_snapshot = Image.open(self.empty_world)
        visual_environment = ImageDraw.Draw(world_snapshot)

        for coordinate, space in world.iteritems():
            if len(space.occupants) > 0:
                visual_environment.point(
                    coordinate, fill=color_map[space.occupants[0].color])
        del visual_environment
        file_name = param["Main Directory"] + \
            "/Bacteria/Generation"+str(generation).zfill(3)+".bmp"
        world_snapshot.save(file_name, "BMP")

        world_snapshot = Image.open(self.empty_world)
        visual_environment = ImageDraw.Draw(world_snapshot)

        for coordinate, space in world.iteritems():

            total_ant_conc = space.rif_conc+space.str_conc+space.quin_conc

            if total_ant_conc > 0:
                if abs(total_ant_conc) > 0:
                    new_color = (240, 240, 240)
                if abs(total_ant_conc) > (max_ant_conc *
                       0.1): new_color = (225, 225, 225)
                if abs(total_ant_conc) > (max_ant_conc*0.2):
                    new_color = (200, 200, 200)
                if abs(total_ant_conc) > (max_ant_conc *
                       0.3): new_color = (170, 170, 170)
                if abs(total_ant_conc) > (max_ant_conc*0.4):
                    new_color = (145, 145, 145)
                if abs(total_ant_conc) > (max_ant_conc *
                       0.5): new_color = (125, 125, 125)
                if abs(total_ant_conc) > (max_ant_conc*0.6):
                    new_color = (100, 100, 100)
                if abs(total_ant_conc) > (
                    max_ant_conc*0.7): new_color = (80, 80, 80)
                if abs(total_ant_conc) > (max_ant_conc*0.8):
                    new_color = (60, 60, 60)
                if abs(total_ant_conc) > (
                    max_ant_conc*0.9): new_color = (40, 40, 40)
                if abs(total_ant_conc) >= max_ant_conc:
                    new_color = (20, 20, 20)

                visual_environment.point(coordinate, fill=new_color)
        del visual_environment
        file_name = param["Main Directory"] + \
            "/World/WorldGeneration"+str(generation).zfill(3)+".bmp"
        world_snapshot.save(file_name, "BMP")

        """
        ##AR RESISTANCE PANEL##
        """

        world_snapshot = Image.open(self.empty_world)
        visual_environment = ImageDraw.Draw(world_snapshot)

        for coordinate, space in world.iteritems():
            if len(space.occupants) > 0:
                if space.occupants[0].resistances["AR_Rif"] == True:
                    visual_environment.point(
                        coordinate, fill=color_map["yellow"])
                else: visual_environment.point(coordinate, fill=color_map["lightgrey"])
        del visual_environment

        file_name = param["Main Directory"] + \
            "/World/Resistance"+str(generation).zfill(3)+".bmp"
        world_snapshot.save(file_name, "BMP")

        """
        ##FREE PHAGE PANEL##
        """

        world_snapshot = Image.open(self.empty_world)
        visual_environment = ImageDraw.Draw(world_snapshot)

        for coordinate, space in world.iteritems():
            if len(space.free_phages) > 0:
                new_color = (240, 240, 240)

                if len(space.free_phages) > (100*0.01):
                    new_color = (170, 170, 170)
                if len(space.free_phages) > (
                    100*0.03): new_color = (125, 125, 125)
                if len(space.free_phages) > (100*0.06):
                    new_color = (100, 100, 100)
                if len(space.free_phages) > (
                    100*0.09): new_color = (60, 60, 60)
                if len(space.free_phages) >= 9:
                    new_color = (20, 20, 20)

                visual_environment.point(coordinate, fill=new_color)

        del visual_environment
        file_name = param["Main Directory"] + \
            "/World/FreePhage"+str(generation).zfill(3)+".bmp"
        world_snapshot.save(file_name, "BMP")

        """
        ##PHAGE INFECTION (AND AR TRANSDUCTION) PANEL##
        """

        world_snapshot = Image.open(self.empty_world)
        visual_environment = ImageDraw.Draw(world_snapshot)

        for coordinate, space in world.iteritems():
            if len(space.occupants) > 0:
                if len(space.occupants[0].phages) > 0:

                    found_mobile_resistant = False

                    for phg in space.occupants[0].phages:

                        if (phg["Type"] == "Virulent") or (phg["Type"] == "Induced_Temperate"):
                            visual_environment.point(
                                coordinate, fill=color_map["blue"])

                        if (phg["Type"] == "Temperate"):
                            visual_environment.point(
                                coordinate, fill=color_map["rosybrown"])
                            for gene in phg["Cargo"]:
                                if ("AR_Rif" in gene[0]) and (gene[3] == "Resistant"):
                                    found_mobile_resistant = "red"

                        if (phg["Type"] == "Defective"):
                            visual_environment.point(
                                coordinate, fill=color_map["lightgreen"])
                            for gene in phg["Cargo"]:
                                if ("AR_Rif" in gene[0]) and (gene[3] == "Resistant"):
                                    found_mobile_resistant = "green"
                    if found_mobile_resistant: visual_environment.point(
                        coordinate, fill=color_map["red"])

                    if found_mobile_resistant:
                        visual_environment.point(
                            coordinate, fill=color_map[found_mobile_resistant])

        del visual_environment
        file_name = param["Main Directory"] + \
            "/World/PhageInfection"+str(generation).zfill(3)+".bmp"
        world_snapshot.save(file_name, "BMP")

        """
        ##PHAGE CRISPR PANEL##
        """

        world_snapshot = Image.open(self.empty_world)
        visual_environment = ImageDraw.Draw(world_snapshot)

        for coordinate, space in world.iteritems():
            if len(space.free_phages) > 0:

                last_phage = space.free_phages[-1][1]
                visual_environment.point(
                    coordinate, fill=color_map["lightgrey"])
                if not(last_phage["crispr_seq"] in crispr_colors.keys()):
                    #new_color=self.colors_generated.next()
                    new_color = self.colors_generated.pop()
                    crispr_colors[last_phage["crispr_seq"]] = (
                        int(new_color[0]*255), int(new_color[1]*255), int(new_color[2]*255))
                visual_environment.point(
                    coordinate, fill=crispr_colors[last_phage["crispr_seq"]])

        del visual_environment
        file_name = param["Main Directory"] + \
            "/World/PhageCrispr"+str(generation).zfill(3)+".bmp"
        world_snapshot.save(file_name, "BMP")

        """
        ##PHAGE CRISPR RESISTANCE PANEL##
        """

        world_snapshot = Image.open(self.empty_world)
        visual_environment = ImageDraw.Draw(world_snapshot)

        for coordinate, space in world.iteritems():
            if len(space.occupants) > 0:

                if len(space.occupants[0].phage_resistance_cassettes) > 0:
                       if not(space.occupants[0].phage_resistance_cassettes[-1] in crispr_colors.keys()):
                            new_color = self.colors_generated.pop()
                            crispr_colors[space.occupants[0].phage_resistance_cassettes[-1]] = (
                                int(new_color[0]*255), int(new_color[1]*255), int(new_color[2]*255))
                            visual_environment.point(coordinate, fill=crispr_colors[space.occupants[0].phage_resistance_cassettes[-1]])
                else:
                    visual_environment.point(coordinate, fill=color_map["lightgrey"])                            
                     
        del visual_environment
        file_name = param["Main Directory"]+"/World/PhageResistanceCrispr"+str(generation).zfill(3)+".bmp"
        world_snapshot.save(file_name, "BMP")


        """
        ##PHAGE RECEPTOR RESISTANCE PANEL##
        """

        world_snapshot = Image.open(self.empty_world)
        visual_environment = ImageDraw.Draw(world_snapshot)

        for coordinate, space in world.iteritems():
            if len(space.occupants) >0:

                if len(space.occupants[0].phage_receptor_resistances)>0:
                    visual_environment.point(coordinate, fill=color_map["red"])  
                else:visual_environment.point(coordinate, fill=color_map["lightgrey"])                            

        del visual_environment
        file_name = param["Main Directory"]+"/World/PhageResistanceReceptor"+str(generation).zfill(3)+".bmp"
        world_snapshot.save(file_name, "BMP")

        self.ConcatenateImages8(generation)

    def ConcatenateImages8(self, generation):
        imageBacteria = Image.open(
            param["Main Directory"]+"Bacteria/Generation"+str(generation).zfill(3)+".bmp")
        imageWorld = Image.open(
            param["Main Directory"]+"World/WorldGeneration"+str(generation).zfill(3)+".bmp")
        imageResistance = Image.open(
            param["Main Directory"]+"World/Resistance"+str(generation).zfill(3)+".bmp")
        imageFreePhage = Image.open(
            param["Main Directory"]+"World/FreePhage"+str(generation).zfill(3)+".bmp")
        imagePhageCrispr = Image.open(
            param["Main Directory"]+"World/PhageCrispr"+str(generation).zfill(3)+".bmp")
        imagePhageResistanceCrispr = Image.open(
            param["Main Directory"]+"World/PhageResistanceCrispr"+str(generation).zfill(3)+".bmp")

        imagePhageInfection = Image.open(
            param["Main Directory"]+"World/PhageInfection"+str(generation).zfill(3)+".bmp")
        imagePhageResistanceReceptor = Image.open(
            param["Main Directory"]+"World/PhageResistanceReceptor"+str(generation).zfill(3)+".bmp")

        blank_image = Image.new(
            "RGB", ((imageBacteria.size[1]*4)+4, (imageBacteria.size[1]*2)+1), "white")

        #Print borders
        plate_draw = ImageDraw.Draw(blank_image)

        #Vertical lines
        for y in range(0, blank_image.size[0]):
            plate_draw.point((imageBacteria.size[1], y), "black")
        for y in range(0, blank_image.size[0]):
            plate_draw.point((imageBacteria.size[1]*2+1, y), "black")
        for y in range(0, blank_image.size[0]):
            plate_draw.point((imageBacteria.size[1]*3+2, y), "black")

            #Middle horizontal black line
        for x in range(0, blank_image.size[0]):
            plate_draw.point((x, imageBacteria.size[1]), "black")
        for x in range(0, blank_image.size[0]):
            plate_draw.point((x, imageBacteria.size[1]+1), "black")

        del plate_draw

            #Top
        location = (0, 0,imageBacteria.size[0],imageBacteria.size[1])
        blank_image.paste(imageBacteria, location)
        location = (imageBacteria.size[1]+1, 0)
        blank_image.paste(imageWorld, location)
        location = (2*imageBacteria.size[1]+2, 0)
        blank_image.paste(imageResistance, location)
        location = (3*imageBacteria.size[1]+3, 0)
        blank_image.paste(imagePhageInfection, location)

            #Bottom
            #location=(0, imageBacteria.size[1]+2);blank_image.paste(imageeDna, location)
        location = (0, imageBacteria.size[1]+2)
        blank_image.paste(imageFreePhage, location)
        location = (imageBacteria.size[1]+1, imageBacteria.size[1]+2)
        blank_image.paste(imagePhageCrispr, location)
        location = (2*imageBacteria.size[1]+2, imageBacteria.size[1]+2)
        blank_image.paste(imagePhageResistanceCrispr, location)
        location = (3*imageBacteria.size[1]+3, imageBacteria.size[1]+2)
        blank_image.paste(imagePhageResistanceReceptor, location)

        blank_image.save(
            param["Main Directory"]+"Joined/Generation"+str(generation).zfill(3)+".bmp")


#########
#########GENOME PRINTING CODE
#########

mobile_color_genomes_map ={"red":"lightsalmon",
                          "green": "palegreen",
                          "yellow": "lightgoldenrodyellow",
                          "blue": "lightblue",
                          "purple": "violet",
                          "grey": "lightgrey"
                          }


def PrintGenome_Draw(genome, colors_map, vertical_location, dr, width_block, height_block, space_between):
    print("Drawing..." + vertical_location)
    current_pos = 0

    for _, loci in enumerate(genome):
        if isinstance(loci, dict):  # Mobile element

            #Phage
            if (loci["Type"] =="Virulent") or (loci["Type"]=="Temperate") or (loci["Type"]=="Induced_Temperate") or (loci["Type"]=="Defective"):
                for gene in loci["Cargo"]:
                    if "AR_Rif" in gene[0]:
                        dr.polygon(((current_pos, (vertical_location+height_block/2)),
                                    (current_pos+(width_block/2), vertical_location),
                                    (current_pos+width_block, (vertical_location+height_block/2)),
                                    (current_pos+(width_block/2),vertical_location+height_block)), fill=("burlywood" if gene[3] =="Sensitive" else "darkorange"))
                        current_pos += width_block+space_between

                    elif "AR_Str" in gene[0]:
                        dr.polygon(((current_pos, (vertical_location+height_block/2)),
                                    (current_pos+(width_block/2), vertical_location),
                                    (current_pos+width_block, (vertical_location+height_block/2)),
                                    (current_pos+(width_block/2),vertical_location+height_block)), fill=("pink" if gene[3] =="Sensitive" else "hotpink"))
                        current_pos += width_block+space_between

                    elif "AR_Quin" in gene[0]:
                        dr.polygon(((current_pos, (vertical_location+height_block/2)),
                                    (current_pos+(width_block/2), vertical_location),
                                    (current_pos+width_block, (vertical_location+height_block/2)),
                                    (current_pos+(width_block/2),vertical_location+height_block)), fill=("grey" if gene[3] =="Sensitive" else "black"))
                        current_pos += width_block+space_between

                    elif gene[0] =="N":
                        dr.rectangle(((current_pos, vertical_location+(height_block/2)),(current_pos+width_block,vertical_location+(height_block/2)+1)), fill="black" )
                        current_pos += width_block+space_between
                    else:
                        color_to_use = None
                        for entry in colors_map:
                            if ((entry==gene[0]) and (colors_map[entry][1]=="Phage")):
                                color_to_use=colors_map[gene[0]][0]                            
                        if color_to_use==None: color_to_use=mobile_color_genomes_map[colors_map[gene[0]][0]]

                        dr.ellipse(((current_pos, vertical_location),(current_pos+width_block,vertical_location+height_block)), fill=color_to_use, outline = "black")
                        current_pos += width_block+space_between

        elif loci[0] =="N":
            dr.rectangle(((current_pos, vertical_location+(height_block/2)),(current_pos+width_block,vertical_location+(height_block/2)+1)), fill="black" )
            current_pos += width_block+space_between
        elif "AR_Rif" in loci[0]:
            dr.polygon(((current_pos, (vertical_location+height_block/2)),
                        (current_pos+(width_block/2), vertical_location),
                        (current_pos+width_block, (vertical_location+height_block/2)),
                        (current_pos+(width_block/2),vertical_location+height_block)), fill=("burlywood" if loci[3] =="Sensitive" else "darkorange"))
            current_pos += width_block+space_between
        elif "AR_Str" in loci[0]:
            dr.polygon(((current_pos, (vertical_location+height_block/2)),
                        (current_pos+(width_block/2), vertical_location),
                        (current_pos+width_block, (vertical_location+height_block/2)),
                        (current_pos+(width_block/2),vertical_location+height_block)), fill=("pink" if loci[3] =="Sensitive" else "hotpink"))
            current_pos += width_block+space_between
        elif "AR_Quin" in loci[0]:
            dr.polygon(((current_pos, (vertical_location+height_block/2)),
                        (current_pos+(width_block/2), vertical_location),
                        (current_pos+width_block, (vertical_location+height_block/2)),
                        (current_pos+(width_block/2),vertical_location+height_block)), fill=("grey" if loci[3] =="Sensitive" else "black"))
            current_pos += width_block+space_between
        else:
            dr.rectangle(((current_pos, vertical_location),(current_pos+width_block,vertical_location+height_block)), fill=colors_map[loci[0]][0], outline = "black")
            current_pos += width_block+space_between


def PrintGenomeExtended_Draw(genome, vertical_location, dr, width_block, height_block, space_between):
    current_pos = 0
    for _, loci in enumerate(genome):
        if isinstance(loci, dict):  # Mobile element

            #Phage
            if (loci["Type"] =="Virulent") or (loci["Type"]=="Temperate") or (loci["Type"]=="Induced Temperate") or (loci["Type"]=="Defective"):
                for gene in loci["Cargo"]:
                    if "AR_Rif" in gene[0]:
                        dr.polygon(((current_pos, (vertical_location+height_block/2)),
                                    (current_pos+(gene[1]/2), vertical_location),
                                    (current_pos+gene[1], (vertical_location+height_block/2)),
                                    (current_pos+(gene[1]/2),vertical_location+height_block)), fill=("burlywood" if gene[3] =="Sensitive" else "darkorange"))


                    if "Phage_Res" in gene[0]:
                        dr.polygon((    (current_pos,(vertical_location+height_block/2)),
                                                (current_pos+(gene[1]/2),vertical_location),
                                                (current_pos+gene[1], (vertical_location+height_block/2)),
                                                (current_pos+(gene[1]/2),vertical_location+height_block)), fill=("paleturquoise" if gene[3] =="Sensitive" else "darkcyan"))

                    else:
                        if gene[0]=="A":
                            dr.ellipse(((current_pos,vertical_location),(current_pos+gene[1],vertical_location+height_block)), fill="lightgreen", outline = "black")
                        if gene[0]=="T":dr.ellipse(((current_pos,vertical_location),(current_pos+gene[1],vertical_location+height_block)), fill="papayawhip", outline = "black")
                        if gene[0]=="C":
                            dr.ellipse(((current_pos,vertical_location),(current_pos+gene[1],vertical_location+height_block)), fill="salmon", outline = "black")
                        if gene[0]=="G":dr.ellipse(((current_pos,vertical_location),(current_pos+gene[1],vertical_location+height_block)), fill="lightcyan", outline = "black")
                        if gene[0]=="N":
                            dr.ellipse(((current_pos,vertical_location),(current_pos+gene[1],vertical_location+height_block)), fill="grey", outline = "black")
                        
                    current_pos += gene[1]+space_between
        else:
            if loci[0]=="A":
                dr.rectangle(((current_pos,vertical_location),(current_pos+loci[1],vertical_location+height_block)), fill="green", outline = "black")
            if loci[0]=="T":dr.rectangle(((current_pos,vertical_location),(current_pos+loci[1],vertical_location+height_block)), fill="yellow", outline = "black")
            if loci[0]=="C":
                dr.rectangle(((current_pos,vertical_location),(current_pos+loci[1],vertical_location+height_block)), fill="red" , outline = "black")
            if loci[0]=="G":dr.rectangle(((current_pos,vertical_location),(current_pos+loci[1],vertical_location+height_block)), fill="lightblue", outline = "black")
            if loci[0]=="N":
                dr.rectangle(((current_pos,vertical_location+(height_block/2)),(current_pos+loci[1],vertical_location+(height_block/2)+1)), fill="black" )
            
            if "AR_" in loci[0]:
                dr.polygon((    (current_pos,(vertical_location+height_block/2)),
                                                (current_pos+(loci[1]/2),vertical_location),
                                                (current_pos+loci[1], (vertical_location+height_block/2)),
                                                (current_pos+(loci[1]/2),vertical_location+height_block)), fill=("burlywood" if loci[3] =="Sensitive" else "darkorange"))

            if "Phage_Res" in loci[0]:
                dr.polygon((    (current_pos,(vertical_location+height_block/2)),
                                                (current_pos+(loci[1]/2),vertical_location),
                                                (current_pos+loci[1], (vertical_location+height_block/2)),
                                                (current_pos+(loci[1]/2),vertical_location+height_block)), fill=("paleturquoise" if loci[3] =="Sensitive" else "darkcyan"))

            current_pos += loci[1]+space_between


def PrintGenomes(genome_list, colors_map, filename):
    if len(genome_list) ==0: 
        print("No genomes")
        return

    width_block = 30
    height_block = 20
    space_between_genes = 2
    space_between_genomes = height_block

    max_loci = 0
    max_genome = 0

    for genome in genome_list:

        size_loci = 0
        size_genome = 0

        for gene in genome:
            if isinstance(gene, dict):
                for cargo_gene in gene["Cargo"]:
                    size_loci += 1; size_genome += cargo_gene[1]
            else: size_loci += 1; size_genome += gene[1]

        if size_loci>max_loci:
            max_loci=size_loci
        if size_genome>max_genome: max_genome=size_genome


    im = Image.new('RGB', ((max_loci*(width_block+space_between_genes)), height_block*(len(genome_list)*2)), (255, 255,255))
    dr = ImageDraw.Draw(im)

    order = 0
    for genome in genome_list:
        PrintGenome_Draw(genome, colors_map, order, dr,
                         width_block, height_block, space_between_genes)
        order += (space_between_genomes*2)
    im.save(filename+"_Loci.png")


# this file is compatible with python 3, so assume translated
