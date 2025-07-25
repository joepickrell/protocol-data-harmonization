"""FHE Decryption and Result Interpretation for SecureGenomics Protocol."""

import tenseal as ts
from typing import Dict, Any

def decrypt_result(encrypted_result: bytes, private_crypto_context: bytes) -> Dict[str, Any]:
    private_crypto_context = ts.context_from(private_crypto_context)
    encrypted_result = ts.bfv_vector_from(private_crypto_context, encrypted_result)
    result = encrypted_result.decrypt()
    return result

def interpret_result(result):
    # remove the extra one at the end to count the number of alleles
    result, num_genomes = result[:-1], result[-1]
    
    print('# Genomes:', num_genomes)
    
    num_alleles = num_genomes * 2
    allele_freqs = [r / num_alleles for r in result]
    
    print('Allele frequencies:', allele_freqs)
    return {
        'allele_frequencies': allele_freqs,
        'num_genomes': num_genomes
    }