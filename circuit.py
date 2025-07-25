"""FHE Computation Circuit for SecureGenomics Protocol."""

from typing import List
import tenseal as ts

def compute(encrypted_datasets: List[bytes], public_crypto_context: bytes) -> bytes:
    # Deserialize the public context
    context = ts.context_from(public_crypto_context)
    
    # Deserialize all encrypted vectors using the context
    vectors = [ts.bfv_vector_from(context=context, data=data) for data in encrypted_datasets]
    
    # multiplication depth: 1
    
    # Homomorphic sum
    encrypted_result = vectors[0]
    for vec in vectors[1:]:
        encrypted_result += vec

    # Serialize the result for return
    return encrypted_result.serialize()
