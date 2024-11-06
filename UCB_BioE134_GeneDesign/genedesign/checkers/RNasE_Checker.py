import re
import sys
import os

# Add the `genedesign` directory to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

class RNaseE_Checker:
  def initiate(self):
    self.patterns = [
        r'AUU',       # AUU site
        r'AU..',      # AUNN pattern where N can be any nucleotide
        r'A[AU]{2,}', # AU-rich region, e.g., AAU, AUU, AAAUU, etc.
    ]

  def run(self, rna_sequence):
      """
      Check for RNase E cleavage sites in the RNA sequence.
      
      Parameters:
          rna_sequence (str): RNA sequence to check for RNase E cleavage sites.
      
      Returns:
          bool or str: True if no cleavage sites are found; otherwise, returns the problematic sequence.
      """
      # Define RNase E consensus cleavage patterns
      
      # Search for each pattern in the DNA sequence
      for pattern in self.patterns:
          match = re.search(pattern, rna_sequence)
          if match:
              # If a cleavage site is found, return the problematic sequence
              return False, match.group(0)
      
      # Return True if no cleavage sites are found
      return True, None

# Example usage
def main():
  checker = RNaseE_Checker()
  checker.initiate()
  rna_sequence = "GAUGCGAUUUACGGAUGC"  # This sequence has "AUU"
  result = checker.run(rna_sequence)
  print(result)  # Output should be the first problematic sequence or True if no cleavage site is found.

if __name__ == "__main__":
    main()
