import os
import sys

# Add the `genedesign` directory to sys.path
#os.chdir('/content/drive/MyDrive/VSCode/Genetic_Frame/UCB_BioE134_GeneDesign/genedesign/')
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__))))

from rbs_chooser import RBSChooser
from models.transcript import Transcript
from checkers.codon_checker import CodonChecker
from checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from checkers.hairpin_checker import hairpin_checker
from checkers.internal_promoter_checker import PromoterChecker
from checkers.GC_checker import GC_content
from checkers.RNasE_Checker import RNaseE_Checker
from checkers.Interference_Checker import Interference_Checker
from seq_utils.hairpin_counter import hairpin_counter
import numpy as np
import pandas as pd
import time
import re

class TranscriptDesigner:
    """
    Reverse translates a protein sequence into a DNA sequence and chooses an RBS using the highest CAI codon for each amino acid.
    """

    def __init__(self):
        self.aminoAcidToCodon = {}
        self.rbsChooser = None

    def initiate(self) -> None:
        """
        Initializes the codon table and the RBS chooser.
        """
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()

        self.cod_check = CodonChecker() #Requires List of strings
        self.forbid = ForbiddenSequenceChecker()
        self.promote = PromoterChecker()
        self.rnase = RNaseE_Checker() #Requires RNA Sequence
        self.interfere = Interference_Checker() #Requries RNA Sequence
        #hairpin has no initiate

        self.forbid.initiate()
        self.cod_check.initiate()
        self.promote.initiate()
        self.rnase.initiate() 
        self.interfere.initiate()

        # Codon table with highest CAI codon for each amino acid (for E. coli)
        self.aminoAcidToCodon = {
            'A': "GCG", 'C': "TGC", 'D': "GAT", 'E': "GAA", 'F': "TTC",
            'G': "GGT", 'H': "CAC", 'I': "ATC", 'K': "AAA", 'L': "CTG",
            'M': "ATG", 'N': "AAC", 'P': "CCG", 'Q': "CAG", 'R': "CGT",
            'S': "TCT", 'T': "ACC", 'V': "GTT", 'W': "TGG", 'Y': "TAC"
        }
        np.random.seed(101)

        #Optimizing Monte Carlo by just having it run all codons for the AA's here
        df = pd.read_csv('./data/codon_usage.txt', sep="\s+", names = ['Codon', 'AA', 'Freq', '1', '2', '3'])
        df = df.drop(['1','2', '3'], axis = 1)
        self.aas = {}
        all_aas = df['AA'].unique()
        np.random.seed(101)
        for aa in all_aas:
          ls = df[df['AA'] == aa][['Codon', 'Freq']].set_index('Codon')['Freq'].to_dict()
          ls = {codon: (freq / 5 if freq < 0.1 else freq) for codon, freq in ls.items()}
          norms = list(ls.values())/np.sum(list(ls.values()))
          #smoothed_probabilities = (1 - 0.4) * norms + 0.4 * norms.mean()
          #smoothed_probabilities /= smoothed_probabilities.sum()
          self.aas[aa] = np.random.choice(list(ls.keys()), p = norms, size = 3000)


    def run(self, peptide: str, ignores: set) -> Transcript:
        """
        Translates the peptide sequence to DNA and selects an RBS.
        
        Parameters:
            peptide (str): The protein sequence to translate.
            ignores (set): RBS options to ignore.
        
        Returns:
            Transcript: The transcript object with the selected RBS and translated codons.
        """
        '''
        # Translate peptide to codons
        codons = [self.aminoAcidToCodon[aa] for aa in peptide]

        # Append the stop codon (TAA in this case)
        codons.append("TAA")  

        # Build the CDS from the codons
        cds = ''.join(codons)
        '''
        def seq_to_list(seq):
          return [seq[i:i+3] for i in range(0, len(seq), 3)]
        
        peptide = peptide + '*'
        #Sliding Window
        cds = "ATG"
        for ind in range(1,len(peptide),3):
          if ind - 12 < 0:
            preamble = cds
          else:
            preamble = cds[len(cds) - 36:] #36
          oligo = peptide[ind:ind+3]
          down = peptide[ind+3:ind+15] #15
          
          '''
          sliced = tuple(self.aas[key] for key in list(oligo+down))
          potentials = tuple(preamble + "".join(group) for group in zip(*sliced))
          '''
          
          sliced = [self.aas[key] for key in list(oligo + down)]
          #potentials = tuple(preamble + "".join(np.random.choice(group) for group in sliced) for _ in range(len(sliced[0])))
          
          ideal_cds = None
          somewhat_ideal = None
          even_less_ideal = None
          #This is where I would do tests to try and figure out the ideal oligo

          for _ in range(len(sliced[0])): #len(sliced[0])
            potent = preamble + "".join(np.random.choice(group) for group in sliced)
            potent_RNA = potent.replace('T', 'U')

            CAI = self.cod_check.run(seq_to_list(potent))[0] if len(seq_to_list(potent)) > 31 else 1
            forbid = self.forbid.run(potent)[0]
            hairpin = hairpin_checker(potent)[0]
            #hairpin_2 = hairpin_counter(potent)[0] < 3
            int_promote = self.promote.run(potent)[0]
            GC = GC_content(potent)[0]
            RNaseE = self.rnase.run(potent_RNA)
            interfer = self.interfere.run(potent_RNA)

            if (CAI and forbid and hairpin and int_promote and GC and RNaseE and interfer):
              ideal_cds = potent
              break
            elif (forbid and hairpin and int_promote):
            #elif(CAI and forbid and int_promote and GC and RNaseE and interfer):
              somewhat_ideal = potent
            elif(forbid and int_promote):
              even_less_ideal = potent

          if ideal_cds == None:
            #print('Got here 1')
            if somewhat_ideal != None:
              ideal_cds = somewhat_ideal
            elif even_less_ideal != None:
              ideal_cds = even_less_ideal
            else:
              #print('Nope')
              ideal_cds = preamble + "".join(np.random.choice(group) for group in sliced)
          #print('Got here')
          cds += ideal_cds[len(preamble):len(preamble) + 9]
        codons = seq_to_list(cds)
        # Choose an RBS
        selectedRBS = self.rbsChooser.run(cds, ignores) 
        # Return the Transcript object
        #print(hairpin_checker(cds))
        #print(self.cod_check.run(codons))
        #print(self.forbid.run(cds)[0])
        #print(self.promote.run(cds)[0])
        #print(cds)
        #print(codons)
        return Transcript(selectedRBS, peptide[:-1], codons)

if __name__ == "__main__":
    # Example usage of TranscriptDesigner
    #peptide = "MYPFIRTARMTV"
    #peptide = 'MKYIVVTGGVMSGLGKGITAASIGRLFVDMGYRVIPIKIDPYINIDAGTMNPFQHGEVYVLKDGTEVDLDLGHYERFIGEEVTGDHNITTGKIYKRVIEKERKGDYLGQTVQIIPHVTDEIKSWIRRVAKESNAEICLVEIGGTVGDIEGMPFLEAIRQMHNEEKEEDFALVHVTLVPLDAGGEQKTKPTQHSVKELRELGLHPDVIVGRCSERLKPATKKKIALFCDVPEEAVISNEDAEDIYEVPLIFKREKLDEYLMRKLNLRAKESRKEWEEMVKRMKTLYEEASIAIVGKYVDVRDAYLSIKEALKHGGIEAGCKVNIVWVDSEDLENVDDFTLDIDGILVPGGFGARGAEGKIRAIEYARENGVPFLGICFGFQLAVIEFARNVVGFSEANSTELDENTPHPVIDLLPEQKGIDEMGGTMRLGDIEVTIKPGTIAHKLYGSEKVVERHRHRYEVNPEYIEKIESKGLVFSAYSDGGRRMEIAELPDHPFFFATQFHPEFKSRPYRPSPPFVGFVRAALKYRREEEI'
    peptide = 'MNPSDVFQIIEGHTKLMRDSIPLIASENLTSLSVRRCYVSDLGHRYAEGRVGERFYEGCKYVDQIESMAIELTRKIFEAEHANVQPISGVVANLAAFFALTNVGDTIMSISVPCGGHISHDRVSAAGLRGLRVIHYPFNSEEMSVDVDETRKVAERERPKLFILGSTLILFRQPVKEIREIADEIGAYVMYDASHVLGLIAGKAFQNPLKEGADVMTGSTHKTFFGPQRAIIASRKELAEKVDRAVFPGVVSNHHLNTLAGYVVAAMEMLEFGEDYAKQVVRNAKALAEELYSLGYKVLGEKRGFTETHQVAVDVREFGGGERVAKVLENAGIILNKNLLPWDSLEKTANPSGIRIGVQEVTRIGMKEEEMRAIAEIMDAAIKEKKSVDELRNEVKELKERFNVIKYSFDESEAYHFPDLR'
    #peptide = 'MRRGLVIVGHGSQLNHYREVMELHRKRIEESGAFDEVKIAFAARKRRPMPDEAIREMNCDIIYVVPLFISYGLHVTEDLPDLLGFPRGRGIKEGEFEGKKVVICEPIGEDYFVTYAILNSVFRIGRDGKGEE'
    #peptide = 'MNFEKVVERICDFIRGVVSSSGSTGVVLGLSGGVDSATVAYLCVRALGSERVFALIMPETGVTPEQDVEDAINVAESLGMEYKLIEINDIVRVFKEKAGEGSKIAEANLKPRIRMVLNYYHANSMNRLVAGTGNKSELMVGYFTKYGDGGVDFLPIGDLYKTEVFQLAAYLGVPRRIIEKKPSARLWPGQTDEEEMGISYAELDEILKLIEKGERRDDEKFRRVVQMVERSRHKREMPPVARVRDLL'
    #peptide = 'MHLERIRGRRSRKLERKKIVLGVTGSIAAVETVKLARELVRRGADVTAVMSRAARKIIHPYALEFATGKRVVTEITGSIEHVNLLGEYGDADLFLIAPCTANTISKIAQGIDDTPVTTFATTALGSGKPIIIVPAMHEAMMRNKAVLENIQRLIDMGIEFVQPRIEEGKAKFPSTETICLHVERELYPKEMKGKRVVVTSGPTYEQIDPIRFISNKSSGRMGLEIALEFWRRGADVVHVTSKPSGMSLPNYKEIRVWSVEDMMKAVLYEIGKGCDLFVSSAAAADFIVDAEAKKIKTAPELVLKLKESPKIIKEVRKIYSGHIIGFKAETGMSDDELLKVASEKMADDNLNMVVANDVLERGMGTEDTRVLILTPKRQEWVEGLKQHVAERIVEAYLEDCL'
    designer = TranscriptDesigner()
    designer.initiate()

    ignores = set()
    transcript = designer.run(peptide, ignores)
    
    # Print out the transcript information
    print(transcript)
