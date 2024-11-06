import pytest
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))

from genedesign.checkers.GC_checker import GC_content

@pytest.fixture
def test_working_GC():
  seq = 'ACGGACCCCCCCTCTTTCAAGGTCTGGTAACACGCTTCGGCCGAAAAGGT'
  ans, content = GC_content(seq)
  assert ans == True
  assert content == None

def test_failed_GC():
  seq = 'GAAAAATCTATACTAATTTCTGATATAATAATATAATTTAGATTAGATCA'
  ans, content = GC_content(seq)
  assert ans == False
  assert content == 16.0