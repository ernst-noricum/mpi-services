from typing import List

from pydantic import BaseModel

#
# DATA MODEL
#


class Sequence(BaseModel):
    name: str
    type: str
    seq: str


class CodonOptimizeRequest(BaseModel):
    sequences: List[Sequence]


class CodonOptimizedSequence(BaseModel):
    original: Sequence
    codon_optimized: Sequence


class CodonOptimizeResponse(BaseModel):
    sequences: List[CodonOptimizedSequence]
