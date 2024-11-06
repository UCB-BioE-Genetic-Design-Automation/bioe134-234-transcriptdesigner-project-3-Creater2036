def GC_content(seq):
  content = (seq.count('C') + seq.count('G'))/len(seq)*100
  if content < 40 or content > 60:
    return False, content
  return True, None

if __name__ == "__main__":
    result, content = GC_content("GCATCGATCGATCGGGATGCCCATGCTAGGGCAAAAATCCTATCTAGCTAC")
    print(result, content)