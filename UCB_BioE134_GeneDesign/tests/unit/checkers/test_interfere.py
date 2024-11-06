import pytest
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))
from genedesign.checkers.Interference_Checker import Interference_Checker

@pytest.fixture
def checker():
    inter = Interference_Checker()
    inter.initiate()
    return inter

def test_fails(checker):
  fail_seqs = ["GAUCCGAUUCGGAUUAA", "UUCGGCAUUUUUUUU", "UCGGAUCUUUGA"]
  for fail in fail_seqs:
    result = checker.run(fail)
    assert result == False

def test_success(checker):
  true_seqs = ["UUUAAAACCCGGG", "CGAUGCUACGAU", "AUGCGGCCAAUUG"]
  for success in true_seqs:
    result = checker.run(success)
    assert result == True
