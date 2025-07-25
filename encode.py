"""VCF Encoding Module for SecureGenomics Protocol - Alzheimer's Disease Allele Frequency Analysis."""

from typing import List, Dict, Any, TextIO
import pysam

# Target variants for Alzheimer's disease sensitive allele frequency analysis
# These are the most clinically significant and population-variable SNPs for AD research
TARGET_VARIANTS = [
    # Primary APOE variants (most critical for AD risk)
    "rs429358",    # APOE ε4 defining SNP - 3-12x increased AD risk
    "rs7412",      # APOE ε2 defining SNP - protective variant
    
    # Secondary high-impact variants
    "rs2075650",   # TOMM40 gene - 2-4x higher AD risk, linked to APOE
    "rs199768005", # APOE rare variant - pathogenic significance
    "rs6857",      # NECTIN2 gene - age-at-onset modifier
    
    # Alternative specification by genomic coordinates (GRCh38)
    # Useful if rsIDs are not available in VCF
    ("19", 44908684, "C"),   # rs429358 - APOE ε4 (C=risk allele)
    ("19", 44908822, "T"),   # rs7412 - APOE ε2 (T=protective allele)  
    ("19", 44892362, "G"),   # rs2075650 - TOMM40 (G=risk allele)
    ("19", 44909057, "A"),   # rs199768005 - APOE rare variant
    ("19", 44888997, "A"),   # rs6857 - NECTIN2 variant
]

def make_record_map(vcf_path):
    """Create mapping from variant identifiers to genotype values."""
    record_map = {}
    vcf_reader = pysam.VariantFile(vcf_path)
    for record in vcf_reader.fetch():
        # Get genotype value (0=ref/ref, 1=ref/alt, 2=alt/alt)
        # Handle missing genotypes gracefully
        try:
            gt = record.samples[0]['GT']
            if None in gt:
                alt_count = 0  # Treat missing as reference
            else:
                alt_count = sum(gt)
        except (KeyError, IndexError):
            alt_count = 0
            
        # Map by rsID if available
        if record.id and record.id != '.':
            record_map[record.id] = alt_count
            
        # Map by chromosomal coordinates
        if record.alts:
            key = (str(record.chrom), record.pos, record.alts[0])
            record_map[key] = alt_count
            
    return record_map

def encode_on_variant_list(record_map: dict, filter_list: list):
    """Encode specific variants from the target list."""
    for pos_or_id in filter_list:
        if isinstance(pos_or_id, tuple):
            # Genomic coordinate specification
            key = (str(pos_or_id[0]), pos_or_id[1], pos_or_id[2])
        elif isinstance(pos_or_id, str):
            # rsID specification
            key = pos_or_id
        else:
            raise ValueError(f"Invalid variant specification: {pos_or_id}")
            
        # Return genotype count (0, 1, or 2) or 0 if variant not found
        yield record_map.get(key, 0)

def encode_vcf(vcf_path: str) -> List[int]:
    """
    Encode VCF genomic data to integer vectors suitable for FHE.
    
    Returns list of integers representing genotype counts for each target variant:
    - 0: homozygous reference (no risk alleles)
    - 1: heterozygous (one risk allele) 
    - 2: homozygous alternate (two risk alleles)
    
    For privacy-sensitive allele frequency analysis of Alzheimer's disease variants.
    """
    record_map = make_record_map(vcf_path)
    encoded_data = list(encode_on_variant_list(record_map, TARGET_VARIANTS))
    return encoded_data