from models.rbs_option import RBSOption
from seq_utils.Translate import Translate
from seq_utils.hairpin_counter import hairpin_counter
from seq_utils.calc_edit_distance import calculate_edit_distance
from seq_utils.merger import merger
from typing import Set


class RBSChooser:
    """
    A class to choose the best RBS for a given CDS sequence.
    """
    def initiate(self) -> None:
        """
        Initialization method for RBSChooser.
        """
        # TODO: Implement initialization logic here
        self.rbs_options = []
        self.translator = Translate()
        self.translator.initiate()
        merged_dict = merger()
        for gene_name, sequences in merged_dict.items():
            utr = sequences['UTR']
            cds = sequences['CDS']

            first_codon_seq = cds[:18]
            first_six_aas = self.translator.run(first_codon_seq)
            self.rbs_options.append(RBSOption(utr, cds, gene_name, first_six_aas))

    def run(self, cds: str, ignores: Set[RBSOption]) -> RBSOption:
        """
        Executes the RBS selection process for the given CDS.

        Parameters:
        - cds (str): The coding sequence to pair with an RBS.
        - ignores (Set[RBSOption]): A set of RBSOption instances to ignore during selection.

        Returns:
        - RBSOption: The selected RBSOption that best pairs with the given CDS.
        """
        rbs = {}
        min_hairpins = float('inf')
        min_edit_distance = float('inf')
        best_option = None
        for rbs_option in self.rbs_options:
            if (rbs_option in ignores): #RBSOption 3 criteria: Not in ignores list, 5'UTR doesn't form secondary structure, Translation of parameter CDS resembles source gene CDS (RBSOption CDS)
                continue
            hairpin = hairpin_counter(rbs_option.utr + cds)[0]
            distance = calculate_edit_distance(rbs_option.first_six_aas, self.translator.run(cds[:18]))

            if hairpin < min_hairpins or (hairpin == min_hairpins and distance < min_edit_distance):
              best_option = rbs_option
              min_hairpins = hairpin
              min_edit_distance = distance
        return best_option  # Placeholder return value for testing purposes


class RBSChooser_simple:
    """
    A simple RBS selection algorithm that chooses an RBS from a list of options, excluding any RBS in the ignore set.
    """

    def __init__(self):
        self.rbsOptions = []

    def initiate(self) -> None:
        """
        Populates the RBS options list with predefined options.
        """
        # Add RBS options based on their sequences and properties
        opt1 = RBSOption(
            utr="aaagaggagaaatactag",
            cds="atggcttcctccgaagacgttatcaaagagttcatgcgtttcaaagttcgtatggaaggttccgttaacggtcacgagttcgaaatcgaaggtgaaggtgaaggtcgtccgtacgaaggtacccagaccgctaaactgaaagttaccaaaggtggtccgctgccgttcgcttgggacatcctgtccccgcagttccagtacggttccaaagcttacgttaaacacccggctgacatcccggactacctgaaactgtccttcccggaaggtttcaaatgggaacgtgttatgaacttcgaagacggtggtgttgttaccgttacccaggactcctccctgcaagacggtgagttcatctacaaagttaaactgcgtggtaccaacttcccgtccgacggtccggttatgcagaaaaaaaccatgggttgggaagcttccaccgaacgtatgtacccggaagacggtgctctgaaaggtgaaatcaaaatgcgtctgaaactgaaagacggtggtcactacgacgctgaagttaaaaccacctacatggctaaaaaaccggttcagctgccgggtgcttacaaaaccgacatcaaactggacatcacctcccacaacgaagactacaccatcgttgaacagtacgaacgtgctgaaggtcgtcactccaccggtgcttaa",
            gene_name="BBa_b0034",
            first_six_aas="MASSED"
        )
        opt2 = RBSOption(
            utr="tcacacaggaaagtactag",
            cds="atgactcaacgtatcgcatatgtaactggtggtatgggtggtatcggtactgcaatttgccagcgtctggcgaaagacggtttccgtgttgttgcgggctgcggtccgaactccccgcgtcgtgaaaagtggctggaacaacagaaagccctgggcttcgacttcattgcctccgagggtaatgtagctgactgggattccaccaagactgccttcgataaagttaaatctgaagtgggcgaagtagatgtactgatcaacaacgccggtattactcgtgatgtcgtattccgcaaaatgacccgtgcagactgggatgcagttatcgacaccaacctgacgtctctgttcaacgttaccaaacaggttattgatggtatggctgaccgtggctggggccgcatcgtgaacatctctagcgttaacggccaaaaaggccaatttggtcagacgaattacagcacggctaaagcaggcctgcacggtttcaccatggcactggcgcaggaagtggcgaccaaaggtgttaccgttaataccgtttctccaggttacatcgccaccgatatggttaaggctatccgccaagatgttctggacaagatcgtggctaccattccggttaaacgcctgggcctgccggaagaaattgcgtccatctgtgcgtggctgagctccgaagagtctggtttttccaccggtgcggatttctctctgaacggtggtctgcacatgggttga",
            gene_name="BBa_b0032",
            first_six_aas="MTQRIA"
        )
        opt3 = RBSOption(
            utr="CCATACCCGTTTTTTTGGGCTAACAGGAGGAATTAAcc",
            cds="atgGacacAattaacatcgctaagaacgacttctctgacatcgaactggctgctatcccgttcaacactctggctgaccattacggtgagcgtttagctcgcgaacagttggcccttgagcatgagtcttacgagatgggtgaagcacgcttccgcaagatgtttgagcgtcaacttaaagctggtgaggttgcggataacgctgccgccaagcctctcatcactaccctactccctaagatgattgcacgcatcaacgactggtttgaggaagtgaaagctaagcgcggcaagcgcccgacagccttccagttcctgcaagaaatcaagccggaagccgtagcgtacatcaccattaagaccactctggcttgcctaaccagtgctgacaatacaaccgttcaggctgtagcaagcgcaatcggtcgggccattgaggacgaggctcgcttcggtcgtatccgtgaccttgaagctaagcacttcaagaaaaacgttgaggaacaactcaacaagcgcgtagggcacgtctacaagaaagcatttatgcaagttgtcgaggctgacatgctctctaagggtctactcggtggcgaggcgtggtcttcgtggcataaggaagactctattcatgtaggagtacgctgcatcgagatgctcattgagtcaaccggaatggttagcttacaccgccaaaatgctggcgtagtaggtcaagactctgagactatcgaactcgcacctgaatacgctgaggctatcgcaacccgtgcaggtgcgctggctggcatctctccgatgttccaaccttgcgtagttcctcctaagccgtggactggcattactggtggtggctattgggctaacggtcgtcgtcctctggcgctggtgcgtactcacagtaagaaagcactgatgcgctacgaagacgtttacatgcctgaggtgtacaaagcgattaacattgcgcaaaacaccgcatggaaaatcaacaagaaagtcctagcggtcgccaacgtaatcaccaagtggaagcattgtccggtcgaggacatccctgcgattgagcgtgaagaactcccgatgaaaccggaagacatcgacatgaatcctgaggctctcaccgcgtggaaacgtgctgccgctgctgtgtaccgcaaggacaaggctcgcaagtctcgccgtatcagccttgagttcatgcttgagcaagccaataagtttgctaaccataaggccatctggttcccttacaacatggactggcgcggtcgtgtttacgctgtgtcaatgttcaacccgcaaggtaacgatatgaccaaaggactgcttacgctggcgaaaggtaaaccaatcggtaaggaaggttactactggctgaaaatccacggtgcaaactgtgcgggtgtcgataaggttccgttccctgagcgcatcaagttcattgaggaaaaccacgagaacatcatggcttgcgctaagtctccactggagaacacttggtgggctgagcaagattctccgttctgcttccttgcgttctgctttgagtacgctggggtacagcaccacggcctgagctataactgctcccttccgctggcgtttgacgggtcttgctctggcatccagcacttctccgcgatgctccgagatgaggtaggtggtcgcgcggttaacttgcttcctagtgaaaccgttcaggacatctacgggattgttgctaagaaagtcaacgagattctacaagcagacgcaatcaatgggaccgataacgaagtagttaccgtgaccgatgagaacactggtgaaatctctgagaaagtcaagctgggcactaaggcactggctggtcaatggctggcttacggtgttactcgcagtgtgactaagcgttcagtcatgacgctggcttacgggtccaaagagttcggcttccgtcaacaagtgctggaagataccattcagccagctattgattccggcaagggtctgatgttcactcagccgaatcaggctgctggatacatggctaagctgatttgggaatctgtgagcgtgacggtggtagctgcggttgaagcaatgaactggcttaagtctgctgctaagctgctggctgctgaggtcaaagataagaagactggagagattcttcgcaagcgttgcgctgtgcattgggtaactcctgatggtttccctgtgtggcaggaatacaagaagcctattcagacgcgcttgaacctgatgttcctcggtcagttccgcttacagcctaccattaacaccaacaaagatagcgagattgatgcacacaaacaggagtctggtatcgctcctaactttgtacacagccaagacggtagccaccttcgtaagactgtagtgtgggcacacgagaagtacggaatcgaatcttttgcactgattcacgactccttcggtaccattccggctgacgctgcgaacctgttcaaagcagtgcgcgaaactatggttgacacatatgagtcttgtgatgtactggctgatttctacgaccagttcgctgaccagttgcacgagtctcaattggacaaaatgccagcacttccggctaaaggtaacttgaacctccgtgacatcttagagtcggacttcgcgttcgcAtaa",
            gene_name="Pbad_rbs",
            first_six_aas="MDTINI"
        )
        self.rbsOptions.extend([opt1, opt2, opt3])

    def run(self, cds: str, ignores: set) -> RBSOption:
        """
        Selects an RBS that is not in the ignore set.
        
        Parameters:
            cds (str): The coding sequence.
            ignores (set): A set of RBS options to ignore.
        
        Returns:
            RBSOption: The selected RBS option.
        """
        for rbsopt in self.rbsOptions:
            if rbsopt not in ignores:
                return rbsopt
        raise Exception("No valid RBS option available.")

if __name__ == "__main__":
    # Example usage of RBSChooser
    cds = "ATGGTAAGAAAACAGTTGCAGAGAGTTGAATTTAA" #Originally not TAA at end but "..."", this however doesn't work since code doesn't know what "..."" is

    # Initialize the chooser
    chooser = RBSChooser()
    chooser.initiate()

    # Choose RBS with no ignores
    ignores = set()
    selected1 = chooser.run(cds, ignores)
    
    # Add the first selection to the ignore list
    ignores.add(selected1)
    
    # Choose another RBS option after ignoring the first
    selected2 = chooser.run(cds, ignores)

    # Print the selected RBS options
    print("Selected1:", selected1)
    print("Selected2:", selected2)
