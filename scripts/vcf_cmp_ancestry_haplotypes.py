#!/usr/bin/env python3

#

import argparse
import contextlib
import gzip
import os
import numpy as np

def parse_vcf_filepath(vcf_file_path):
    if not os.path.exists(vcf_file_path):
        raise FileNotFoundError(f"File not found: {vcf_file_path}")
    zipped = False
    path, filename = os.path.split(vcf_file_path)
    if filename.endswith('.vcf.gz'):
        zipped = True
        prefix = os.path.splitext(os.path.splitext(filename)[0])[0]
    elif filename.endswith('.vcf'):
        prefix = os.path.splitext(filename)[0]
    else:
        raise ValueError(f"Unexpected file extension: {filename}. VCF file must end in .vcf or .vcf.gz")
    return {'path': path, 'prefix': prefix, 'zipped': zipped}

def extract_tracts(vcf=str, msp=str, num_ancs=int, output_dir=None, output_vcf=None, compress_output=None):
    if not os.path.exists(vcf):
        raise ValueError(f"The path '{vcf}' does not exist.")
    if not os.path.exists(msp):
        raise ValueError(f"The path '{msp}' does not exist.")
    if not isinstance(num_ancs, int):
        raise ValueError("The 'num_ancs' must be an integer.")
    if output_dir is not None and not os.path.exists(output_dir):
        raise ValueError(f"The output directory '{output_dir}' does not exist.")
    if not isinstance(output_vcf, bool):
        raise ValueError("The 'output_vcf' must be a boolean value.")
    if not isinstance(compress_output, bool):
        raise ValueError("The 'compress_output' must be a boolean value.")

    vcf_info = parse_vcf_filepath(vcf)
    vcf_path, vcf_prefix, zipped = vcf_info['path'], vcf_info['prefix'], vcf_info['zipped']

    output_files = {}
    output_path = f"{os.path.join(output_dir, vcf_prefix) if output_dir else os.path.join(vcf_path, vcf_prefix)}"
    file_extension = f"{'.gz' if compress_output else ''}"

    for i in range(num_ancs):
        output_files[f"dos{i}"] = f"{output_path}.anc{i}.dosage.txt{file_extension}"
        output_files[f"ancdos{i}"] = f"{output_path}.anc{i}.hapcount.txt{file_extension}"
        if output_vcf:
            output_files[f"vcf{i}"] = f"{output_path}.anc{i}.vcf{file_extension}"

    with open(msp) as mspfile, gzip.open(vcf, "rt") if zipped else open(vcf) as vcf, contextlib.ExitStack() as stack:
        files = {
            fname: stack.enter_context(gzip.open(output_file, "wt") if compress_output else open(output_file, "w"))
            for fname, output_file in output_files.items()
        }
        vcf_header = ""
        window = ("", 0, 0)

        for line in vcf:
            skip_line = False
            if line.startswith("##"):
                if output_vcf:
                    vcf_header += line
                continue
            elif line.startswith("#"):
                if output_vcf:
                    vcf_header += line
                anc_header = "\t".join([line.strip("# ").split("\t", 9)[item] for item in [0, 1, 2, 3, 4, 9]])

                for i in range(num_ancs):
                    files[f"dos{i}"].write(anc_header)
                    files[f"ancdos{i}"].write(anc_header)
                    if output_vcf:
                        files[f"vcf{i}"].write(vcf_header)

            if not line.startswith("#"):
                row = line.strip().split("\t", 9)
                row[8] = "GT"
                vcf_out = "\t".join(row[:9])
                dos_anc_out = "\t".join(row[:5])

                genos = row[9].split("\t")
                pos = int(row[1])
                output_lines = {}
                pop_genos = {}

                for i in range(num_ancs):
                    output_lines[f"dos{i}"] = dos_anc_out
                    output_lines[f"ancdos{i}"] = dos_anc_out
                    if output_vcf:
                        output_lines[f"vcf{i}"] = vcf_out

                while not (row[0] == window[0] and (window[1] <= pos < window[2])):
                    if row[0] == window[0] and window[1] > pos:
                        skip_line = True
                        break
                    ancs = mspfile.readline()
                    if ancs.startswith("#"):
                        continue
                    if not ancs:
                        break
                    ancs_entry = ancs.strip().split("\t", 6)
                    calls = np.array([int(x) for x in ancs_entry[6].split("\t")])
                    window = (ancs_entry[0], int(ancs_entry[1]), int(ancs_entry[2]))
                    if row[0] == window[0] and window[1] > pos:
                        skip_line = True
                        break
                if skip_line:
                    continue

                # Arreglo para almacenar los resultados por individuo
                counts = np.zeros((len(genos), num_ancs), dtype=int)
                anc_counts = np.zeros((len(genos), num_ancs), dtype=int)

                for i, geno in enumerate(genos):
                    geno_parts = list(map(int, geno.split(":")[0].split("|")))
                    geno_a, geno_b = geno_parts[:2]
                    call_a, call_b = calls[2 * i], calls[2 * i + 1]

                    counts[i, call_a] += 1 if geno_a == 1 else 0
                    counts[i, call_b] += 1 if geno_b == 1 else 0
                    anc_counts[i, call_a] += 1
                    anc_counts[i, call_b] += 1

                for j in range(num_ancs):
                    output_lines[f"dos{j}"] += "\t" + "\t".join(map(str, counts[:, j]))
                    output_lines[f"ancdos{j}"] += "\t" + "\t".join(map(str, anc_counts[:, j]))
                    if output_vcf:
                        output_lines[f"vcf{j}"] += f"\t{pop_genos.get(j, '.|.')}"

                for j in range(num_ancs):
                    files[f"dos{j}"].write(output_lines[f"dos{j}"] + "\n")
                    files[f"ancdos{j}"].write(output_lines[f"ancdos{j}"] + "\n")
                    if output_vcf:
                        files[f"vcf{j}"].write(output_lines[f"vcf{j}"] + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True, help="Path to phased genotypes, VCF file (*.vcf or *.vcf.gz)")
    parser.add_argument("--msp", required=True, help="Path to ancestry calls, MSP file (*.msp or *.msp.tsv)")
    parser.add_argument("--num-ancs", type=int, required=True, help="Number of ancestral populations within the VCF file.")
    parser.add_argument("--output-dir", help="Directory for output files. Directory must already exist.")
    parser.add_argument("--output-vcf", help="Boolean. If ancestry-specific VCF files need to be generated", action="store_true")
    parser.add_argument("--compress-output", help="gzip (not bgzip) all output files", action="store_true")

    args = parser.parse_args()
    extract_tracts(**vars(args))

