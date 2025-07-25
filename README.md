# Data harmonization protocol -- input any genetic data and get a VCF in a standard format


Why? Any locally run protocol needs to make some assumptions about the file it's working on. E.g.
- what genome build are coordinates on? B37 or B38? (or something else?)
- are chromosomes named like 'chr1' or '1'?
- should missing variants be treated as reference/reference calls (e.g. in WGS data)? Or as missing/untyped variants (as in genotyping array data)?

## Usage
- **Aggregated**: Secure computation across encrypted datasets
- **Local**: Direct analysis on local VCF files

## Files
- `protocol.yaml`: Configuration and parameters
- `crypto_context.py`: FHE key generation
- `encode.py`: VCF data encoding
- `encrypt.py`: Data encryption
- `circuit.py`: FHE computation circuit
- `decrypt.py`: Result decryption
- `local_analysis.py`: Local-only analysis 
