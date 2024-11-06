import pytest
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))
from genedesign.checkers.RNasE_Checker import RNaseE_Checker

@pytest.fixture
def checker():
    r = RNaseE_Checker()
    r.initiate()
    return r

def test_fails(checker):
  fail_seqs = ["UUGACAAUUAAUCAUCGAACUAGTAUAAU", "AUGACAAUUAAUCAUCGAACUAGTAUAAU", "AAUCGAGUGGGCA"]
  for fail in fail_seqs:
    result = checker.run(fail)
    assert result == False

def test_success(checker):
  true_seqs = ["AACAATGTATTCACCG", "AATCACAGGAGTCGCCC", "TTGAAGAGAGC"]
  for success in true_seqs:
    result = checker.run(success)
    assert result == True
