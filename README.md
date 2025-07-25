# Alzheimer's Disease Sensitive Allele Frequency Protocol

Template protocol for secure multi-party computation of Alzheimer's disease susceptibility variant frequencies with privacy-preserving genomic analysis.

## Variants Analyzed
### Primary APOE Variants (Critical for AD Risk)
- **rs429358**: APOE ε4 defining SNP (19:44908684 T>C) - 3-12x increased AD risk
- **rs7412**: APOE ε2 defining SNP (19:44908822 C>T) - protective variant

### Secondary High-Impact Variants
- **rs2075650**: TOMM40 gene (19:44892362 A>G) - 2-4x higher AD risk, linked to APOE
- **rs199768005**: APOE rare variant (19:44909057 T>A) - pathogenic significance
- **rs6857**: NECTIN2 gene (19:44888997 G>A) - age-at-onset modifier

*Note: These variants show dramatic population frequency differences (5-37%) across ethnic groups, making them sensitive for genomic research requiring privacy protection.*

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