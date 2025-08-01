"""Local-only Analysis for SecureGenomics Protocol - data haromonization"""

from typing import Dict, Any
import os, mimetypes, gzip, subprocess, uuid


def is_gzipped(filepath: str) -> bool:
    try:
        with open(filepath, "rb") as f:
            magic = f.read(2)
            return magic == b'\x1f\x8b'
    except Exception:
        return False

def open_file(filepath):
    if is_gzipped(filepath):
        return gzip.open(filepath, "rt", encoding="utf-8", errors="ignore")
    return open(filepath, "r", encoding="utf-8", errors="ignore")

def is_vcf(filepath):
    """
    Check if a file is a valid VCF using bcftools view.
    """
    try:
        subprocess.run(
            ["bcftools", "view", "-Ov", filepath],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True
        )
        return True
    except subprocess.CalledProcessError:
        return False

def convert_gvcf_to_vcf(input_file, output_file):
    cmd = (
        f"bcftools norm -m - {input_file} | "
        f"bcftools view -i 'ALT!=\"<NON_REF>\" && GT!=\"0/0\"' -Oz -o {output_file}"
    )
    subprocess.run(cmd, shell=True, check=True)
    subprocess.run(f"bcftools index {output_file}", shell=True, check=True)
    return output_file

def is_gvcf(vcf_file, sample_lines=100):
    with open_file(vcf_file) as f:
        count = 0
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 5:
                continue
            alt = fields[4]
            if "<NON_REF>" in alt:
                return True
            count += 1
            if count >= sample_lines:
                break
    return False
        
def detect_23andme_v_ancestry(input_path):
    """
    Detects if a tsv is in 23andme or ancestry style
    23andme-style means 4 columns and chr names like X Y MT
    ancestry-style means 5 columns and chr names like 23, 24, 25
    Output: '23andme', 'ancestry', 'undetermined'
    """
    sharedchrs = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"}
    ttamchrs = {"X", "Y", "MT"}
    ancchrs = {"23", "24", "25"}
    sharedcount= 0
    ttcount = 0
    anccount = 0
    fourfieldcount = 0
    fivefieldcount = 0
    with open_file(input_path) as fin:
        for line in fin:
            if line.startswith("#") or line.lower().startswith("rsid"):
                continue
            fields = line.strip().split("\t")
            if len(fields) == 4:
                fourfieldcount = fourfieldcount+1
            elif len(fields) == 5:
                fivefieldcount = fivefieldcount+1
            else:
                return 'undetermined'
            chr = fields[1]
            if chr in sharedchrs:
                sharedcount = sharedcount+1
            elif chr in ttamchrs:
                ttcount = ttcount+1
            elif chr in ancchrs:
                anccount = anccount+1
    if ttcount > 0 and anccount ==0 and fivefieldcount ==0:
        print("Detected 23andMe format")
        return '23andme'
    elif anccount >0 and ttcount ==0 and fourfieldcount ==0:
        print("Detected ancestry format")
        return 'ancestry'
    return 'undetermined'

def detect_genome_build_tsv(input_path, b37_positions_set, b38_positions_set, threshold_ratio=5):

    tsv_positions = set()
    with open_file(input_path) as fin:
        for line in fin:
            if line.startswith("#") or line.lower().startswith("rsid"):
                continue
            chrom, pos = line.strip().split("\t")[1:3]
            key = f"chr{chrom}\t{pos}"
            tsv_positions.add(key)

    matches_b37 = len(tsv_positions & b37_positions_set)
    matches_b38 = len(tsv_positions & b38_positions_set)

    print(f"Matched {matches_b37} positions to b37, {matches_b38} to b38.")

    if matches_b37 > threshold_ratio * matches_b38:
        return "b37"
    elif matches_b38 > threshold_ratio * matches_b37:
        return "b38"
    else:
        return "ambiguous"

def load_positions(file_path):
    with gzip.open(file_path, 'rt') as f:
        return set(line.strip() for line in f if line.strip())

def convert_ancestry_to_23andme_format(input_path, uuid):
    """
    Converts AncestryDNA format file to 23andMe-style format.
    Output columns: rsid, chromosome, position, genotype
    """
    chrom_map = {"23": "X", "24": "Y", "25": "X"} # 25 is PAR
    output_path = f"{uuid}.converted.txt"
    with open_file(input_path) as fin, open(output_path, "w") as fout:
        fout.write("# This file was converted from AncestryDNA format to 23andMe-style format\n")
        fout.write("# rsid\tchromosome\tposition\tgenotype\n")

        for line in fin:
            if line.startswith("#") or line.lower().startswith("rsid"):
                continue
            fields = line.strip().split()
            if len(fields) != 5:
                continue
            rsid, chrom, pos, allele1, allele2 = fields
            chrom = chrom_map.get(chrom, chrom)
            if allele1 =="0" or allele2 =="0":
                continue
            genotype = allele1 + allele2
            fout.write(f"{rsid}\t{chrom}\t{pos}\t{genotype}\n")
    return output_path

def is_tsv(filepath):
    try:
        with open_file(filepath) as f:
            for line in f:
                if line.startswith("# rsid") or line.startswith("rs"):
                    parts = line.strip().split("\t")
                    if len(parts) >= 4 and parts[0].startswith("rs"):
                        return True
        return False
    except Exception:
        return False

def is_supported_genetic_file(filepath, uuid, b37_position_file, b38_position_file):
    if not os.path.exists(filepath):
        return False, "File not found.", None

    mimetype, _ = mimetypes.guess_type(filepath)
    if mimetype and not (mimetype.startswith("text") or mimetype == "application/octet-stream"):
        return False, f"Rejected MIME type: {mimetype}", None

    b37_positions = load_positions(b37_position_file)
    b38_positions = load_positions(b38_position_file)

    if is_vcf(filepath):
        if is_gvcf(filepath):
            print("Detected gVCF, converting to VCF...")
            filepath = convert_gvcf_to_vcf(filepath, uuid)
        newpath = fix_vcf_main_chromosomes_to_chr(filepath, uuid)
        build = detect_genome_build(newpath, b37_positions, b38_positions)

        if build == "b37":
            #lifted_path = liftover_b37_to_b38_crossmap(newpath, uuid)
            return True, "VCF file detected as b37.", newpath
        elif build == "b38":
            return True, "VCF file detected as b38.", newpath
        else:
            return False, "Genome build ambiguous, unable to process.", None
    if is_tsv(filepath):
        type = detect_23andme_v_ancestry(filepath)
        if type == "undetermined":
            return False, "Unable to determine file format", None
        build = detect_genome_build_tsv(filepath, b37_positions, b38_positions)
        if build == "ambiguous":
            return False, "Unable to determine genome build", None
        if type == "23andme" and build =="b37":
            #vcfpath = convert_23_to_vcf(filepath, uuid, b37_genome_file)
            #newpath = fix_vcf_main_chromosomes_to_chr(vcfpath, uuid)
            #lifted_path = liftover_b37_to_b38_crossmap(newpath, uuid)
            return True, "23andMe file detected as b37", filepath
        elif type == "23andme" and build =="b38":
            #vcfpath = convert_23_to_vcf(filepath, uuid, b38_genome_file)
            #newpath = fix_vcf_main_chromosomes_to_chr(vcfpath, uuid)
            return True, "23andMe file detected as b38.", filepath
        elif type == "ancestry" and build =="b37":
            ttstylepath = convert_ancestry_to_23andme_format(filepath, uuid)
            #vcfpath = convert_23_to_vcf(ttstylepath, uuid, b37_genome_file)
            #newpath = fix_vcf_main_chromosomes_to_chr(vcfpath, uuid)
            #lifted_path = liftover_b37_to_b38_crossmap(newpath, uuid)
            return True, "Ancestry file detected as b37", ttstylepath
        elif type == "ancestry" and build =="b38":
            ttstylepath = convert_ancestry_to_23andme_format(filepath, uuid)
            #vcfpath = convert_23_to_vcf(ttstylepath, uuid, b38_genome_file)
            #newpath = fix_vcf_main_chromosomes_to_chr(vcfpath, uuid)
            return True, "Ancestry file detected as b38.", ttstylepath
        return False, "Unsupported file format. Must be a VCF or TSV file.", None
    # if is_23andme(filepath):
    #     convert_23andme_to_vcf(filepath, uuid)

    return False, "Unsupported file format. Must be a VCF file.", None

def detect_genome_build(vcf_file, b37_positions_set, b38_positions_set, threshold_ratio=5):
    result = subprocess.run(
        ["bcftools", "query", "-f", "%CHROM\t%POS\n", vcf_file],
        stdout=subprocess.PIPE,
        text=True,
        check=True
    )

    vcf_positions = set()
    for line in result.stdout.strip().splitlines():
        chrom, pos = line.split("\t")
        key = f"{chrom}\t{pos}"
        vcf_positions.add(key)

    matches_b37 = len(vcf_positions & b37_positions_set)
    matches_b38 = len(vcf_positions & b38_positions_set)

    print(f"Matched {matches_b37} positions to b37, {matches_b38} to b38.")

    if matches_b37 > threshold_ratio * matches_b38:
        return "b37"
    elif matches_b38 > threshold_ratio * matches_b37:
        return "b38"
    else:
        return "ambiguous"

def fix_vcf_main_chromosomes_to_chr(filepath, uuid):
    """
    Converts main chromosomes (1-22, X, Y, MT) to 'chr*' style if needed.
    Raises error if mixed styles are found.
    Always outputs to {uuid}.vcf.gz when valid.
    """
    main_chr = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}
    main_chr_chr = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}
    output_filepath = f"{uuid}.vcf.gz"

    # Extract chromosome names
    result = subprocess.run(
        ["bcftools", "query", "-f", "%CHROM\n", filepath],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
        text=True
    )
    chroms = set(result.stdout.strip().splitlines())
    main_present = {c for c in chroms if c in main_chr or c in main_chr_chr}

    # Detect style
    in_chr_style = {c for c in main_present if c in main_chr_chr}
    in_nochr_style = {c for c in main_present if c in main_chr}

    if in_chr_style and in_nochr_style:
        raise ValueError(f"❌ Mixed main chromosome styles detected: chr* = {in_chr_style}, no-chr = {in_nochr_style}")

    if in_chr_style:
        print("Main chromosomes already in chr* style. Copying to output.")
        subprocess.run([
            "bcftools", "view", "-Oz",
            "-o", output_filepath,
            filepath
        ], check=True)
        subprocess.run(["bcftools", "index", "-f", output_filepath], check=True)
        return output_filepath

    if in_nochr_style:
        print(f"Converting main chromosomes to chr* style → {output_filepath}")
        mapping_file = f"{uuid}_chrom_map.txt"
        with open(mapping_file, "w") as f:
            for c in chroms:
                if c in main_chr:
                    new_c = "chrM" if c == "MT" else f"chr{c}"
                    f.write(f"{c}\t{new_c}\n")
        subprocess.run([
            "bcftools", "annotate",
            "--rename-chrs", mapping_file,
            "-Oz", "-o", output_filepath, filepath
        ], check=True)
        subprocess.run(["bcftools", "index", "-f", output_filepath], check=True)
        os.remove(mapping_file)
        return output_filepath

    print("No main chromosomes found")
    raise ValueError(f"❌ No main chromosome styles detected: chr* = {in_chr_style}, no-chr = {in_nochr_style}")


def analyze_local(input_vcf_file_path: str, output_vcf_file_path: str, protocol_config: Dict[str, Any] = None) -> Dict[str, Any]:
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
    if len(sys.argv) > 2:
        input_vcf_path = sys.argv[1]
        output_vcf_path = sys.argv[2]
        uuid = uuid.uuid4()
        ok, msg, newpath = is_supported_genetic_file(input_vcf_path, uuid, b37_position_file="data/GSA-b37.pos.gz", b38_position_file="data/GSA-b38.pos.gz")
        if ok:
            try:
                os.rename(newpath, output_vcf_path)
                print(f"File successfully renamed to: {output_vcf_path}")
            except Exception as e:
                print(f"Error renaming file: {e}")   
        print(msg)
    else:
        print("Usage: python local_analysis.py <input_vcf_file_path> <output_vcf_file_path>") 