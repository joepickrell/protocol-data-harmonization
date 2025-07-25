import tenseal as ts

def generate_keys():
    # see https://github.com/OpenMined/TenSEAL/blob/main/tutorials/Tutorial%200%20-%20Getting%20Started.ipynb
    private_context = ts.context(
        ts.SCHEME_TYPE.BFV, 
        poly_modulus_degree=4096, 
        plain_modulus=1032193
    )
    
    # Enable encryption keys
    private_context.generate_galois_keys()
    private_context.generate_relin_keys()
    
    private_context_bytes = private_context.serialize(save_public_key=True, save_secret_key=True)
    public_context_bytes = private_context.serialize(save_public_key=True, save_secret_key=False)
    
    return public_context_bytes, private_context_bytes