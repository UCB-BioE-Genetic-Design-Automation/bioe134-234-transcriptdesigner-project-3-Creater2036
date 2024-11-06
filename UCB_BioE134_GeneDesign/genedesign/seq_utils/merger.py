from Bio import SeqIO
from collections import defaultdict
import numpy as np
import pandas as pd
import subprocess

def merger():
  file_urls = [
    'https://drive.google.com/uc?id=1U7AKYm2n0O1KDOdcYcJCYvrZQw5Ul_TZ',
    'https://drive.google.com/uc?id=1gALv5ZIoWCXWGk4U93pAiJUwbwTCCti6'
  ]

  # Download the files
  '''
  for url in file_urls:
      subprocess.run(['gdown', url, '-O', './seq_utils/'])
  '''
  genbank_file = './seq_utils/sequence.gb'
  gene_dict = defaultdict(dict)  # Dictionary to store gene info
  for record in SeqIO.parse(genbank_file, "genbank"):
      for feature in record.features:
          if feature.type == "gene":
              locus_tag = feature.qualifiers.get("locus_tag", [None])[0]
              gene_name = feature.qualifiers.get("gene", [None])[0]

              # CDS information
              cds_feature = None
              for cds in record.features:
                  if cds.type == "CDS" and cds.qualifiers.get("locus_tag") == [locus_tag]:
                      cds_feature = cds
                      break

              if cds_feature:
                  start, end = cds_feature.location.start, cds_feature.location.end
                  strand = cds_feature.location.strand
                  if strand == 1:  # Forward strand
                      utr_start = max(0, start - 50)
                      utr_seq = record.seq[utr_start:start]
                  else:  # Reverse strand, we need to reverse complement
                      utr_start = end
                      utr_seq = record.seq[utr_start:utr_start + 50].reverse_complement()

                  cds_seq = cds_feature.extract(record.seq)
                  # Save the gene information in the dictionary
                  gene_dict[locus_tag] = {
                      "gene": gene_name,
                      "UTR": utr_seq,
                      "CDS": cds_seq
                  }
  genes_info = gene_dict
  df = pd.read_csv('./seq_utils/511145-WHOLE_ORGANISM-integrated.txt', sep='\t', comment= '#', skiprows = 8, header = None, names = ['Locus_Tag', 'Abundance_Pairs'])
  df['Locus_Tag'] = df['Locus_Tag'].str.replace('511145.', '', regex=False)
  df = df.nlargest(int(len(df)*0.05),'Abundance_Pairs')

  prot = df.set_index('Locus_Tag')['Abundance_Pairs'].to_dict()
  merged_dict = {}

  # Merge the dictionaries
  for locus_tag, abundance_pair in prot.items():
      if locus_tag in genes_info:
          # Extract UTR and CDS from default_dict
          utr = genes_info[locus_tag]['UTR']
          cds = genes_info[locus_tag]['CDS']

          # Store UTR and CDS in the merged dictionary with locus_tag as key
          merged_dict[locus_tag] = {
              'UTR': utr,
              'CDS': cds
          }
  return merged_dict

def main():
  print(merger())
if __name__ == "__main__":
    main()