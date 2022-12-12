#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Import relevent libraries
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
import seaborn as sns
import pandas as pd
from colorama import Back, Style, Fore
import matplotlib.pyplot as plt
import matplotlib


# # Load Covid, MERS, SARS and  Ebola virus 
# 
# - Ebola virus link: https://www.ncbi.nlm.nih.gov/nuccore/MZ854251.1
# - SARS virus link : https://www.ncbi.nlm.nih.gov/nuccore/NC_004718.3
# - MERS virus link: https://www.ncbi.nlm.nih.gov/nuccore/NC_019843
# - Covid-19 Virus link: https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2

# In[2]:


# load the fasta file of the viruses
# using SeqIO.read we can load the entire fasta file and the output will be  the complete records of the virus genome

# list of  fasta file genomes that will be importetd
# Filename format: virusName_sequence.fasta
filenames = ['corona_sequence.fasta', 'MERS_sequence.fasta','SARS_sequence.fasta','EBOLA_sequence.fasta']

#empty dict to store all the information of the fasta files 
genome_dict = {}

#for through list of genomes to add them to dictionary
for file in filenames:
    key=file.partition("_")[0] # extract the first part of the string "virusName_sequence.fasta" until the "_" to make the dict key
    record=SeqIO.read(file, "fasta") ##function that reads in the DNA/RNA fasta file
    genome_dict[key] = record # add the record to genome_ dict


# In[3]:


#example of genome dict entry
print(genome_dict["corona"])


# In[4]:


class Genome():
    """
    Genome class that contains all the attibutes and methods that pertain to the genome

    Attributes
    -------------
    sequnece: str
        Nucleotide sequence of genome 
        example: ATGCCCCTTT.....
        
    name: str
        name of virus sequence 
        example: "EBOLA" or  MERS
    """
    
    # constructor to initialize the object
    def __init__(self,sequence,name):
        self.sequence= str(sequence) # One strand sequence of nucleotides 
        self.name=str(name) # custom short name of virus
        
    # instance method   
    def base_composition(self,plot=False):
        """
        Calculates the nucleotide base compostion of the genome as a proportion of the whole genome
        
        input
        -------
        
        plot: bool
             The function will output a barplot of the total base compostion of the genome
            
        Returns:
            dict: nucleotide base compositon
            
        if plot== True:
            Barplot of dict representing base composition
        """
        bases = {'A': 0, 'C': 0, 'G': 0, 'T': 0} # empty  dict that will store the number bases
        
        #for loop to that calculates the number of base in the genome
        for nucleotide in bases:
            bases[nucleotide] = self.sequence.count(nucleotide)/len(self.sequence)*100 # store the base 
            
        # conditional that outputs the bar graph
        if plot == True:
            
            ndf = pd.DataFrame.from_dict(bases, orient ='index') # convert base dict to dataframe 
            ndf = ndf.reset_index() # reset the index of the dataframe
            ndf = ndf.rename(columns={"index": "Nucleotide", 0: "Composition"}) # rename columnms of the dataframa
            
            #create the plot opf the base composition
            ax = sns.barplot(x="Nucleotide", y="Composition",
                             data=ndf, palette="Set2")#.set(title=f"Base compostion of {self.name}")
            plt.title(f"Base compostion of {self.name}")
            
            plt.grid()# gridlines on the plot
            
            # set font details on figure
            font = {'family' : 'sans-serif',
            'weight' : 'bold',
            'size'   : 15} 
            matplotlib.rc('font', **font)
            
            
            return plt.plot()
        
        else:
            return print(f"{self.name} has a base composition of A:{round(bases['A'],2) }%, C:{round(bases['C'],2) }%, G:{round(bases['G'],2) }%,  T:{round(bases['T'],2) }%")
    

    def gc_content(self):
    
        """
        Calculate the GC content of DNA/RNA sequence
        
        --------
        Returns:
            float: gc content percentage of base sequence

        """
        num_G_and_C=float(self.sequence.count('C') + self.sequence.count('G')) # extract the number of G and C of in sequence
        res=num_G_and_C/len(self.sequence)*100 # calculate as a proportion of nucleotide length
        
        return round(res,2)
    
    # for transcription
    def transcribe(self):
        
        """
        Transcribe a DNA sequence to an RNA sequence
        ---------
        Returns:
            string: RNA base string sequence
        """

        # convert string into list
        RNA_seq_list = list(self.sequence)  

        for i in range(len(self.sequence)):

          if(RNA_seq_list[i]=='G'):
              RNA_seq_list[i]='C'

          elif(RNA_seq_list[i]=='C'):
              RNA_seq_list[i]='G'

          elif (RNA_seq_list[i] == 'T'):
              RNA_seq_list[i] = 'A'

          elif (RNA_seq_list[i] == 'A'):
              RNA_seq_list[i] = 'U'

          else:
              print('Invalid Input')  

        return ''.join(RNA_seq_list)# converts the list of single characters into one long neucletotide sequqence
    
    # for transcription
    def transcript(self):
        """
        Create mRNA transcript that will be used during the translation process to creat proteins        
        ---------
        Returns:
            string: mRNA base sequence that will be used to translation process
        """
        # convert string into list
        RNA_seq_list = list(self.sequence)  

        for i in range(len(self.sequence)):

          if (RNA_seq_list[i] == 'T'):
              RNA_seq_list[i] = 'U'

        return ''.join(RNA_seq_list)# converts the list of single characters into one long neucletotide sequqence

    def translate(self):
        """
        Translate the genome to a string sequence of amino acids
        "*" represent the stop condons
        "-" represent nulls in the  string
        ---------
        Returns:
            string: list of amino acids
        """
        rna=Genome(self.sequence,self.name).transcript()
        
        # RNA codon table
        rna_codon = {"UUU" : "F", "CUU" : "L", "AUU" : "I", "GUU" : "V",
                   "UUC" : "F", "CUC" : "L", "AUC" : "I", "GUC" : "V",
                   "UUA" : "L", "CUA" : "L", "AUA" : "I", "GUA" : "V",
                   "UUG" : "L", "CUG" : "L", "AUG" : "M", "GUG" : "V",
                   "UCU" : "S", "CCU" : "P", "ACU" : "T", "GCU" : "A",
                   "UCC" : "S", "CCC" : "P", "ACC" : "T", "GCC" : "A",
                   "UCA" : "S", "CCA" : "P", "ACA" : "T", "GCA" : "A",
                   "UCG" : "S", "CCG" : "P", "ACG" : "T", "GCG" : "A",
                   "UAU" : "Y", "CAU" : "H", "AAU" : "N", "GAU" : "D",
                   "UAC" : "Y", "CAC" : "H", "AAC" : "N", "GAC" : "D",
                   "UAA" : "*", "CAA" : "Q", "AAA" : "K", "GAA" : "E",
                   "UAG" : "*", "CAG" : "Q", "AAG" : "K", "GAG" : "E",
                   "UGU" : "C", "CGU" : "R", "AGU" : "S", "GGU" : "G",
                   "UGC" : "C", "CGC" : "R", "AGC" : "S", "GGC" : "G",
                   "UGA" : "*", "CGA" : "R", "AGA" : "R", "GGA" : "G",
                   "UGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G" ,
                   "NNN" : "_"
                   }
       

        protein_string = ""

        # Generate protein string
        for i in range(0, len(rna)-(3+len(rna)%3), 3): 
        #goes through sequence of nucleotides and groups them into sets of 3
             protein_string += rna_codon[rna[i:i+3]]

        return protein_string
    
    def amino_acid_composition(self, plot=False):
        """
        Calculate the amino the 
        ---------
        Returns:
            string: list of amino acids
        """
        protein_seq=Genome(self.sequence,self.name).translate()# translaate the genome to a protein sequence
        
        # dictionary to store the amino acids numbers that are genrated during translation
        amino_acids_dict = {'Y': 0, 'A': 0, 'Q': 0, 'D': 0, 'C': 0, 'G': 0, 'V': 0, 'T': 0, 'E': 0, 'N': 0, 
                       'K': 0, 'R': 0, 'S': 0, 'I': 0, 'H': 0, 'M': 0, 'F': 0, 'L': 0, 'W': 0, 'P': 0}
        
        #count the amino acids
        for amino_acid in amino_acids_dict:
            amino_acids_dict[amino_acid] = protein_seq.count(amino_acid)/len(protein_seq)*100
        
        else:
            return amino_acids_dict
    
    def three_nucleotide_seq(self):
    # function to split the nuclotitie seq into a 3 nucleotide bundle 
        nucloetides = []
        for pos, char in enumerate(self.sequence):
            if pos != 0 and pos%3 == 0:
                nucloetides.append(' ')
            nucloetides.append(char)
        return ''.join(nucloetides)
    

    
    def seq_repr(self, strand="dna",end_of_str=-1,remove_gaps=False):
        #color codes that elements o the stings sequence
        
        if strand == 'dna': #if sequence given is DNA exceute code below
            
            # color code to represent genome sequences
            nu_clr_switcher = {
                # standard color-codes
                'A': Back.GREEN,
                'C': Back.YELLOW,
                'G': Back.RED,
                'T': Back.BLUE,
                ' ': Style.RESET_ALL
            }
            
            genome_str=Genome(self.sequence[0:end_of_str],self.name).three_nucleotide_seq()
            
            #genome_str = three_nucleotide_seq(seq=genome_str[0:end_of_str]) # use the three_nucleotide_seq to break up nuclotides
            line_break_cntr = 0

            for i in range(len(genome_str)):# for loop that truns through the entrire strand


                if genome_str[i] == ' ': # Count the number of spaces added by the three_nucleotide_seq fucntion
                    line_break_cntr += 1

                    if line_break_cntr>0 and line_break_cntr%20==0:# Use a new line if there are more than 20 tripltes in one line
                        text = "\n"
                    else:
                        text = nu_clr_switcher[genome_str[i]] + genome_str[i]

                else:
                    text = nu_clr_switcher[genome_str[i]] + genome_str[i]
                print(text, end="")

            Style.RESET_ALL

        if strand == 'protein': #if sequence given is protien exceute code below
            
            protein_clr_switcher = {
            # color-code by proteinfamily's polarity
            'A': Back.BLUE,
            'V': Back.BLUE,
            'I': Back.BLUE,
            'L': Back.BLUE,
            'M': Back.BLUE,
            'F': Back.BLUE,
            'Y': Back.CYAN,
            'W': Back.BLUE,
            'H': Back.CYAN,
            'R': Back.RED,
            'K': Back.RED,
            'N': Back.GREEN,
            'Q': Back.GREEN,
            'E': Back.MAGENTA,
            'D': Back.MAGENTA,
            'S': Back.GREEN,
            'T': Back.GREEN,
            'G': Back.YELLOW,
            'P': Back.YELLOW,
            'C': Back.BLUE,
            ' ': Style.RESET_ALL,
            '*': Style.RESET_ALL,
            "_" :Style.RESET_ALL
            
        }
            
            genome_str=Genome(self.sequence,self.name).translate() # translate nucleotide sequnce to amino acid sequence
            
            # removes the "*" when printing out the nucletide sequence
            if remove_gaps==True:
                genome_str=''.join(list(filter(None,genome_str.split("*"))))
                
            for i in range(len(genome_str)):# for loop that truns through the entrire strand

                if genome_str[i] in protein_clr_switcher: # check if the amino acid is in the the  protien color swticher dict
                        text = protein_clr_switcher[genome_str[i]] + genome_str[i]

                else:
                    Style.RESET_ALL
                    text = genome_str[i]
                print(text, end="")
                
    def __repr__(self):
        return f"The instance contains the genome of the {self.name} virus which sequence has {len(self.sequence)} pairs"
     


# # Use Genome class to extract and comapare different virus information

# In[5]:


# Create instanances of each virus in Genome class

Covid=Genome(genome_dict["corona"].seq,"Covid-19")
SARS=Genome(genome_dict["SARS"].seq,"SARS")
MERS=Genome(genome_dict["MERS"].seq,"MERS")
Ebola=Genome(genome_dict["EBOLA"].seq,"Ebola")


class_list = [Covid,SARS,MERS,Ebola] # enter variable of instance manually into a list


# # Comparing base compostion between genomes
# 
# The information in DNA is stored as a code made up of four chemical bases: adenine (A), guanine (G), cytosine (C).The order, or sequence, of these bases determines the information available for building and maintaining an organism, similar to the way in which letters of the alphabet appear in a certain order to form words and sentences.DNA bases pair up with each other, A with T and C with G, to form units called base pairs. In this section, we explore the proportion of bases that each genome has.

# In[6]:


# Visualize the base composition virus nucleotides

for count,ele in enumerate(class_list):
    plt.subplot(int(len(class_list)**0.5), #nrows is the square root of the total number of instances in class
                int(len(class_list)**0.5), #ncols is the square root of the total number of instances in class
                count+1) #index of plot
    plot=ele.base_composition(plot=True)
    plt.plot(plot)
    ele.base_composition()
    
    # set the spacing between subplots
    plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=1.9,
                    top=1.9,
                    wspace=0.5,
                    hspace=0.5)
    


# # Comparing GC-content
# 
# In molecular biology and genetics, GC-content (or guanine-cytosine content) is the percentage of nitrogenous bases in a DNA or RNA molecule that are either guanine (G) or cytosine (C). This measure indicates the proportion of G and C bases out of an implied four total bases, also including adenine and thymine in DNA and adenine and uracil in RNA.
# 
# GC-content may be given for a certain fragment of DNA or RNA or for an entire genome. When it refers to a fragment, it may denote the GC-content of an individual gene or section of a gene (domain), a group of genes or gene clusters, a non-coding region, or a synthetic oligonucleotide such as a primer.
# 
# G-C base-pairs have stronger interactions (than A-T BPs) arising from their ability to form three hydrogen bonds in water. A-T base-pairing yeilds only two hydrogen bonds. That is why the melting point of double stranded DNA is higher for high G-C content DNA as well as for longer pieces of DNA. High GC content DNA can make it difficult to perform PCR amplification, because it is difficult to design a primer that is long enough to provide great specificity, while maintaining a melting point below the optimal temperature of DNA polymerase. DMSO, Glycerol, Salt or MgCl2 can all be used to perturb the effective melting point.

# In[7]:


# print the GC content of the genomes
for ele in class_list:
    print(f"The GC content of {ele.name} is {ele.gc_content()}%")

    
#visualize the GC content on one plot
virus_name = [str(ele.name) for ele in class_list] # for loop to create names of the virus 
virus_values = [ele.gc_content() for ele in class_list] # for loop to apply the gc_content function to the all the viruss
  
fig = plt.figure(figsize = (10, 5))#

 
# creating the bar plot
plt.bar(virus_name, virus_values, color ='blue', width = 0.8)

plt.grid()
plt.xlabel("Virus")
plt.ylabel("GC content percentage")
plt.title("GC content of different viruses")
plt.show()


# # Show 1st 120 base pairs of genomes
# 
# Here, we aim to show the color coded nucleotide sequence of the genome using standard color-codes.

# In[8]:


end_of_str=120 # define the number of base pairs to visualize


# ### The first 120 base pairs of Covid-19

# In[9]:


Covid.seq_repr(strand="dna",end_of_str=end_of_str)


# ### The first 120 base pairs of MERS

# In[10]:


MERS.seq_repr(strand="dna",end_of_str=end_of_str)


# ### The first 120 base pairs of SARS

# In[11]:


SARS.seq_repr(strand="dna",end_of_str=end_of_str)


# ### The first 120 base pairs of Ebola

# In[12]:


Ebola.seq_repr(strand="dna",end_of_str=end_of_str)


# # Transcribing and translating translating the genomes to find amino acid compostion 

# In[13]:


amino_names=list(Ebola.amino_acid_composition().keys())# get the names of the colunms names of the amino acids

#creating a dataframe to make a style chart to show the different 
data = {'amino': amino_names}

#for loop to add genome values to dict
for ele in class_list:
    data[ele.name]=list(ele.amino_acid_composition().values())

#create dataframe and add it too plot a stlye bar that represents the 
amino = pd.DataFrame.from_dict(data)
r1 = amino.sort_values(by='Covid-19', ascending=False).style.bar(subset=["Covid-19"],color='#')
r1.background_gradient(cmap='Reds')


# # Visualizing the amino acid sequences  for each genome
# 
# The amino acid base sequence is the order in which amino acids are arranged in a protein after translation.The sequence of amino acids determines the protein's primiary, secondary, tertiray as well its quaternary structure whoch dictates its three-dimensional structure pertaining to its specific biological function.During the process of protein synthesis, the nucleotide sequence is first transcribed into messenger RNA (mRNA), and then the sequence of nucleotides in the mRNA is translated into the sequence of amino acids in the protein.

# # Covid-19

# In[14]:


Covid.seq_repr(strand="protein")


# # MERS

# In[15]:


MERS.seq_repr(strand="protein")


# # SARS

# In[16]:


SARS.seq_repr(strand="protein")


# # Ebola

# In[17]:


Ebola.seq_repr(strand="protein")


# In[ ]:




