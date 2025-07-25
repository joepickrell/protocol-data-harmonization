"""Local-only Analysis for SecureGenomics Protocol - Alzheimer's Disease Allele Frequency Analysis."""

from typing import Dict, Any, List, Union, Tuple
from encode import encode_vcf, TARGET_VARIANTS

# Mapping of variants to their clinical information
VARIANT_INFO = {
    "rs429358": {
        "gene": "APOE",
        "variant_id": "rs429358",
        "position": "19:44908684",
        "risk_allele": "C",
        "clinical_significance": "APOE ε4 defining SNP - 3-12x increased AD risk"
    },
    "rs7412": {
        "gene": "APOE", 
        "variant_id": "rs7412",
        "position": "19:44908822",
        "risk_allele": "T",
        "clinical_significance": "APOE ε2 defining SNP - protective variant"
    },
    "rs2075650": {
        "gene": "TOMM40",
        "variant_id": "rs2075650", 
        "position": "19:44892362",
        "risk_allele": "G",
        "clinical_significance": "2-4x higher AD risk, linked to APOE"
    },
    "rs199768005": {
        "gene": "APOE",
        "variant_id": "rs199768005",
        "position": "19:44909057", 
        "risk_allele": "A",
        "clinical_significance": "APOE rare variant - pathogenic significance"
    },
    "rs6857": {
        "gene": "NECTIN2",
        "variant_id": "rs6857",
        "position": "19:44888997",
        "risk_allele": "A", 
        "clinical_significance": "Age-at-onset modifier"
    },
    # Coordinate-based mappings (for backup specifications)
    ("19", 44908684, "C"): {
        "gene": "APOE",
        "variant_id": "rs429358",
        "position": "19:44908684",
        "risk_allele": "C",
        "clinical_significance": "APOE ε4 defining SNP - 3-12x increased AD risk"
    },
    ("19", 44908822, "T"): {
        "gene": "APOE",
        "variant_id": "rs7412", 
        "position": "19:44908822",
        "risk_allele": "T",
        "clinical_significance": "APOE ε2 defining SNP - protective variant"
    },
    ("19", 44892362, "G"): {
        "gene": "TOMM40",
        "variant_id": "rs2075650",
        "position": "19:44892362", 
        "risk_allele": "G",
        "clinical_significance": "2-4x higher AD risk, linked to APOE"
    },
    ("19", 44909057, "A"): {
        "gene": "APOE",
        "variant_id": "rs199768005",
        "position": "19:44909057",
        "risk_allele": "A",
        "clinical_significance": "APOE rare variant - pathogenic significance"
    },
    ("19", 44888997, "A"): {
        "gene": "NECTIN2", 
        "variant_id": "rs6857",
        "position": "19:44888997",
        "risk_allele": "A",
        "clinical_significance": "Age-at-onset modifier"
    }
}

def get_variant_info(variant: Union[str, Tuple]) -> Dict[str, str]:
    """Get clinical information for a variant."""
    if variant in VARIANT_INFO:
        return VARIANT_INFO[variant]
    else:
        # Default info for unknown variants
        return {
            "gene": "Unknown",
            "variant_id": str(variant),
            "position": str(variant),
            "risk_allele": "Unknown",
            "clinical_significance": "Unknown variant"
        }

def calculate_allele_frequency(genotype_count: int, total_samples: int = 1) -> float:
    """
    Calculate allele frequency from genotype count.
    
    Args:
        genotype_count: 0 (ref/ref), 1 (ref/alt), or 2 (alt/alt)
        total_samples: Number of samples (for single sample analysis = 1)
    
    Returns:
        Allele frequency of the alternate allele
    """
    if total_samples == 0:
        return 0.0
    
    # For a single sample, convert genotype count to allele frequency
    # 0 genotypes = 0/2 alleles = 0.0 frequency
    # 1 genotype = 1/2 alleles = 0.5 frequency  
    # 2 genotypes = 2/2 alleles = 1.0 frequency
    return genotype_count / 2.0

def analyze_local(vcf_file_path: str, protocol_config: Dict[str, Any] = None) -> Dict[str, Any]:
    """
    Perform local-only analysis without encryption.
    
    This simulates the end-to-end encode/encrypt/circuit flow but without actual encryption,
    providing a clear view of what the secure protocol would compute.
    """
    # Encode VCF data (this replaces the encode step in the secure protocol)
    encoded_data = encode_vcf(vcf_file_path)
    
    # Simulate the circuit computation locally (without encryption)
    variant_count = len(TARGET_VARIANTS)
    
    if variant_count == 0:
        return {
            'analysis_type': 'Local Alzheimer Disease Allele Frequency Analysis',
            'error': 'No target variants defined',
            'variants': []
        }
    
    if len(encoded_data) != variant_count:
        return {
            'analysis_type': 'Local Alzheimer Disease Allele Frequency Analysis', 
            'error': f'Data length mismatch: expected {variant_count}, got {len(encoded_data)}',
            'variants': []
        }
    
    # Format results (this simulates the decrypt step)
    results = {
        'analysis_type': 'Local Alzheimer Disease Allele Frequency Analysis',
        'sample_count': 1,  # Single sample analysis
        'protocol_simulation': 'encode -> compute -> decrypt (without encryption)',
        'variants': []
    }
    
    # Process each variant
    for i, variant in enumerate(TARGET_VARIANTS):
        if i < len(encoded_data):
            genotype_count = encoded_data[i]
            allele_frequency = calculate_allele_frequency(genotype_count)
            
            # Get variant information
            variant_info = get_variant_info(variant)
            
            # Calculate risk assessment
            risk_level = "Unknown"
            if genotype_count == 0:
                risk_level = "No risk alleles"
            elif genotype_count == 1:
                risk_level = "One risk allele (heterozygous)"
            elif genotype_count == 2:
                risk_level = "Two risk alleles (homozygous)"
            
            results['variants'].append({
                'gene': variant_info['gene'],
                'variant_id': variant_info['variant_id'],
                'position': variant_info['position'],
                'risk_allele': variant_info['risk_allele'],
                'clinical_significance': variant_info['clinical_significance'],
                'genotype_count': genotype_count,
                'allele_frequency': allele_frequency,
                'risk_level': risk_level,
                'raw_encoding': genotype_count  # Shows what the FHE circuit would process
            })
    
    return results

def print_analysis_results(results: Dict[str, Any]) -> None:
    """Pretty print the analysis results."""
    print(f"\n=== {results['analysis_type']} ===")
    
    if 'error' in results:
        print(f"Error: {results['error']}")
        return
        
    print(f"Sample Count: {results['sample_count']}")
    print(f"Protocol Simulation: {results['protocol_simulation']}")
    print(f"\nVariant Analysis:")
    print("-" * 100)
    
    for variant in results['variants']:
        print(f"Gene: {variant['gene']}")
        print(f"Variant: {variant['variant_id']} ({variant['position']})")
        print(f"Risk Allele: {variant['risk_allele']}")
        print(f"Clinical Significance: {variant['clinical_significance']}")
        print(f"Genotype Count: {variant['genotype_count']}")
        print(f"Allele Frequency: {variant['allele_frequency']:.3f}")
        print(f"Risk Level: {variant['risk_level']}")
        print(f"Raw Encoding: {variant['raw_encoding']} (what FHE circuit processes)")
        print("-" * 100)

if __name__ == "__main__":
    # Example usage
    import sys
    if len(sys.argv) > 1:
        vcf_path = sys.argv[1]
        results = analyze_local(vcf_path)
        print_analysis_results(results)
    else:
        print("Usage: python local_analysis.py <vcf_file_path>") 