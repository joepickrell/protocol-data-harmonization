"""
FHE Encryption Module for SecureGenomics Protocol.

This module handles the encryption of encoded genomic data using
the FHE public context for secure computation.
"""

import tenseal as ts
from typing import List

def encrypt_data(encoded_data: List[int], public_crypto_context: bytes) -> bytes:
    """
    Encrypt encoded genomic data using FHE public context.
    """
    
    # Add a 1 to the end of the encoded data to count the number of alleles (see its twin in decrypt.py)
    encoded_data = encoded_data + [1]
    
    return ts.bfv_vector(ts.context_from(public_crypto_context), list(encoded_data)).serialize()