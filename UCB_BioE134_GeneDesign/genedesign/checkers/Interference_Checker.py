import sys
import os

# Add the `genedesign` directory to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from seq_utils.merger import merger

class Interference_Checker:
  def initiate(self):
    self.native_rna_dict = merger()

  def run(self, designed_rna, min_duplex_length=7):
      """
      Checks for potential duplex formation between the designed RNA and native RNAs.
      
      Parameters:
          designed_rna (str): The RNA sequence that has been designed.
          native_rna_dict (dict): A dictionary where keys are gene names and values are native RNA sequences.
          min_duplex_length (int): Minimum length of complementary sequence to consider as a potential duplex.
          
      Returns:
          bool or tuple: True if no duplex is predicted; otherwise, returns a tuple (problematic sequence, gene name).
      """
      # Convert the designed RNA to its reverse complement
      complement = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
      reverse_complement = ''.join(complement[base] for base in reversed(designed_rna))
      
      # Check each native RNA for potential duplex formation
      for gene, native_rna in self.native_rna_dict.items():
          for i in range(len(native_rna) - min_duplex_length + 1):
              # Extract a segment of the native RNA to match against
              segment = native_rna[i:i + min_duplex_length]
              
              # Check if this segment is in the reverse complement of the designed RNA
              if segment in reverse_complement:
                  return False, segment  # Return problematic sequence and gene name
      
      return True, None  # No duplex formation detected

# Example usage
def main():
  designed_rna = "AUGCCGAUUCGGAUC"
  checker = Interference_Checker()
  checker.initiate()
  result = checker.run(designed_rna)
  print(result)  # Output should be (problematic sequence, gene name) if duplex is found, otherwise True.

if __name__ == "__main__":
    main()