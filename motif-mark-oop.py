#!/usr/bin/env python
# coding: utf-8

# In[168]:


from itertools import product
import re
import cairo
import argparse 


# In[169]:

def get_args():
    parser = argparse.ArgumentParser(description="visualize motifs on sequences")
    parser.add_argument("-f", "--fastafile", help="fasta file", type= str)
    parser.add_argument("-m", "--motiffile", help="motif file", type= str)


    return parser.parse_args()

args = get_args()




sequence = 'aaaabbbaaacccdddeeefff'
query_set = {'aaa','bbb'}

def motif_position_calculator(query_set,sequence):
    motif_position_set= set()
    for query in query_set:  
        result = re.finditer(f'(?=({query}))', sequence)
    
        for i in result:
            motif_position_set.add((i.start(1), i.end(1)))
    return motif_position_set

    
motif_position_calculator(query_set,sequence)


# In[170]:


motif="YgcY"

def filler(word, from_char, to_char):
    options = [(c,) if c != from_char else (from_char, to_char) for c in word]
    return (''.join(o) for o in product(*options))

def motif_set_getter(motif):

    motif=motif.upper()
    
    motif_set=set(filler(motif, "Y", "C"))
    
    new_motif_set=set()
    for sequence in motif_set:
        new_motif_set.add(sequence.replace('Y','T'))

    return new_motif_set

motif_set_getter(motif)


# In[171]:


instring = "ABcdEFgh"

def intron_exon_split(instring):
    uplow = instring[0].isupper() # boolean to check whether the first letter is uppercase or not
    
    splitted = []
    new_word = ""
    position_counter = 0
    
    for c in instring:
        if c.isupper() == uplow:
            new_word += c
            position_counter+=1
        else:
            if c.isupper()!=True:
                splitted.append(("E",new_word,len(new_word),position_counter-len(new_word)))
            else:
                splitted.append(("I",new_word,len(new_word),position_counter-len(new_word)))
                
            new_word = c
            position_counter+=1
            uplow = not uplow
    
    if c.isupper()==True:
                splitted.append(("E",new_word,len(new_word),position_counter-len(new_word)))
    else:
                splitted.append(("I",new_word,len(new_word),position_counter-len(new_word)))
    
      # bc no change at end of string
    
    return splitted
intron_exon_split(instring)


# In[172]:


def oneline_fasta(file_input,first_output):
    first_line=True
    with open(first_output,"w") as op: ##to get the files on a single line
        with open(file_input,"r") as file:
    
            for line in file:
                line=line.strip('\n')
                if line.startswith(">"):
                    if first_line==True:
                        op.writelines([line,'\n'])
                        first_line=False
                    else:
                        op.writelines(['\n',line,'\n'])
                else:
                    op.writelines(line)


def oneline_fasta1(file_input):
    first_line=True
    list_records = []
    with open(file_input,"r") as file:
        line_seq=''
        for line in file:
            line=line.strip('\n')
            if line.startswith(">"):
                if first_line==True:
                    list_records.append(line)
                    first_line=False
                else:
                    list_records.append(line_seq)
                    list_records.append(line)
                    
                    line_seq=''
            else:
                line_seq +=line
    list_records.append(line_seq)
    return list_records

# In[173]:

inputf = args.fastafile   #"Figure_1.fasta"
motiff = args.motiffile



#print( oneline_fasta1(inputf))


#oneline_fasta(inputf,f"{inputf.split('.f')[0]}_oneline.fasta")    



#input_file = open(f"{inputf.split('.f')[0]}_oneline.fasta","r")
input_file = oneline_fasta1(inputf)
final_figure_name = f"{inputf.split('.f')[0]}.png"
motif_file = open(motiff,'r')#'Fig_1_motifs.txt'


# In[174]:


dict_DNA_or_RNA={"A":"T","G":"C","C":"G","T":"A","N":"N","A":"U","U":"A",
          "a":"t","g":"c","c":"g","t":"a","n":"n","a":"u","u":"a"}

def rev_complement(sequence: str) -> str:
    '''Takes a DNA sequence and for each base in the sequence, we complement it with the associated base i.e, A>T, G>C, C>G, T>A. We then reverse the
sequence ('hello world'[::-1]) to output finally the reverse complement'''
    revcomp=""
    for letter in sequence:
        revcomp=dict_DNA_or_RNA[letter]+revcomp
    return revcomp


# In[175]:


class Gene:
    def __init__(self,header,the_sequence):
        self.name = header.strip('\n')
        self.sequence = the_sequence.strip('\n')
        self.sequence_length = len(self.sequence)
        self.sequence_singlecase= self.sequence.upper()
        self.intron_exon_list=intron_exon_split(self.sequence)
        self.rev_comp= False
        
        if "reverse complement" in header:
            self.rev_comp= True
                            


# In[176]:


#color_list = ['#006400', '#00008b', '#b03060', '#ff0000', '#ffff00', '#deb887', '#00ff00', '#00ffff','#ff00ff','#6495ed']
color_list = [(1,0,1),(0,1,0),(0,0,1),(0,1,1),(1,1,0)] #add more

class Motif:
    motif_times_called = 0
    def __init__(self,motif):
        self.ogname=motif.strip('\n')
        if "u" in self.ogname or "U" in self.ogname:
            self.name = self.ogname.replace('u','T')
            self.name = self.ogname.replace('U','T')
        else:
            self.name = self.ogname
        self.length=len(self.name)
        self.sequences= motif_set_getter(self.name)
        self.color = color_list[Motif.motif_times_called]
        Motif.motif_times_called+=1
        
        
# In[177]:


class MotifFinderPerGene:
    motif_gene_name_set = set()
    def __init__(self,gene_object,motif):
        self.motif_name = motif.ogname
        self.gene_name = gene.name
        self.motif_positions = motif_position_calculator(motif.sequences,gene_object.sequence_singlecase)
        self.color = motif.color
        MotifFinderPerGene.motif_gene_name_set.add(self.gene_name)
        
        # if gene.rev_comp == False:
        #     self.motif_positions = motif_position_calculator(motif.sequences,gene_object.sequence_singlecase)
        # else:
        #     motif_set = set()
        #     for motif in motif.sequences:
        #         motif_set.add(rev_complement(motif))
        #     self.motif_positions = motif_position_calculator(motif_set,gene_object.sequence_singlecase)
        #     print(motif_set)


# In[178]:


motif_object_list = []
for line in motif_file:
    motif_object_list.append(Motif(line))
#print(motif_object_list)
# for motif in motif_object_list:
#     print(motif.sequences)
#     print(motif.ogname)
#     print(motif.color)


# In[179]:


counter=0
gene_list=[]
for line in input_file:
    if counter%2==0:
        gene_name= line
    else:
        gene_sequence = line
        gene_list.append(Gene(gene_name,gene_sequence))
    counter+=1

# for gene in gene_list:
#     print(gene.name)
#     print(gene.sequence_singlecase)
#     print(gene.rev_comp)


# In[180]:


motif_finder_list = []
for gene in gene_list:
    for motif in motif_object_list:
        motif_finder_list.append(MotifFinderPerGene(gene,motif))

# for i in motif_finder_list:
#     print(i.gene_name)
#     print(i.motif_name)
#     print(i.motif_positions)
#print(MotifFinderPerGene.motif_gene_name_set)


# In[192]:

fig_width = 60+150*len(MotifFinderPerGene.motif_gene_name_set) + 20*len(MotifFinderPerGene.motif_gene_name_set)+25
# code 
# importing pycairo 

surface = cairo.ImageSurface (cairo.FORMAT_ARGB32, 1100, fig_width)

context = cairo.Context(surface) 

context.set_source_rgb(1, 1, 1)
context.rectangle(0, 0, 1100, 2000)
context.fill()

context.set_source_rgb(0, 0, 0) 
context.set_font_size(20) 
  
context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL) 

counter = 50
for gene in gene_list:
    # position for the text 
    context.move_to(50, (counter)) 
    # displays the text 
    context.show_text(gene.name) 
    counter +=150
# stroke out the color and width property 
    context.stroke() 

counter = 125
gene_positions = {}
for gene in gene_list:
    for feature in gene.intron_exon_list:
        #print((50+feature[3],counter))
        #print((50+feature[3]+feature[2],counter))
        
        if feature[0] == 'I':
            context.set_line_width(4)
        else:
            context.set_line_width(12)
        context.move_to(50+feature[3],counter)        #(x,y)
        context.line_to(50+feature[3]+feature[2],counter)
        context.stroke()
    gene_positions[gene.name]=counter
    counter +=150

for motif in motif_finder_list:
    for gene in motif.motif_gene_name_set:
        if motif.gene_name == gene:
            #print(gene,gene_positions[gene])
            
            for position in motif.motif_positions :
                #print(position[0],gene_positions[gene]+5,position[1]-position[0],10)
                #print(motif.color[0], motif.color[1], motif.color[2])
                context.set_source_rgba(motif.color[0], motif.color[1], motif.color[2],0.4) 
                context.rectangle(50+position[0],gene_positions[gene]-25,position[1]-position[0],50)
                context.fill()
#print(counter)


context.set_font_size(15)
counter = counter -90
for motif in motif_object_list:
     
    context.set_line_width(15)
    context.set_source_rgba(motif.color[0], motif.color[1], motif.color[2],0.5)
    context.move_to(50,counter)        #(x,y)
    context.line_to(100,counter)
    context.stroke()
    
    context.set_source_rgb(0, 0, 0) 
    context.move_to(120, counter+5)
    context.show_text(motif.ogname)
    context.stroke()
    counter+=20
    #print(counter)
    

context.set_source_rgb(0, 0, 0) 
context.move_to(120, counter+5)
context.show_text('intron')
context.stroke()

context.set_line_width(4)
context.move_to(50,counter)       
context.line_to(100,counter)
context.stroke()

counter+=20

context.move_to(120, counter+5)
context.show_text('exon')
context.stroke()

context.set_line_width(12)
context.move_to(50,counter)       
context.line_to(100,counter)
context.stroke()

# printing message when file is saved 
print(f"File Saved to {final_figure_name}") 
surface.write_to_png(final_figure_name)







