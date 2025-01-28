import numpy as np
from pathlib import Path
import multiprocessing as mp
from tqdm import tqdm
import biotite.structure as struc
import biotite.structure.io as strucio
from biotite.structure.io.pdbx import CIFFile
import subprocess
import tempfile
import os

def get_mean_plddt(struct):
    """Get mean pLDDT from B-factor column for CA atoms"""
    try:
        print(f"Available annotation categories: {struct.get_annotation_categories()}")
        
        # Try to get B-factor annotation
        if "b_factor" in struct.get_annotation_categories():
            b_factors = struct.get_annotation('b_factor')
        elif "B_factor" in struct.get_annotation_categories():
            b_factors = struct.get_annotation('B_factor')
        else:
            print("No B-factor annotation found in structure")
            return 0.0
            
        ca_mask = struct.atom_name == "CA"
        print(f"Found {np.sum(ca_mask)} CA atoms")
        plddt_scores = b_factors[ca_mask]
        return np.mean(plddt_scores) / 100.0
    except Exception as e:
        print(f"Error calculating mean pLDDT: {str(e)}")
        return 0.0

def calculate_sasa_ss8(filepath, struct):
    """Calculate SASA using biotite and SS8 using mkdssp"""
    try:
        # Get CA atoms for SASA calculation
        ca_mask = struct.atom_name == "CA"
        ca_struct = struct[ca_mask]
        
        # Calculate SASA
        sasa = struc.sasa(ca_struct, vdw_radii="Single", probe_radius=1.4)
        
        # Calculate SS8 using mkdssp
        with tempfile.NamedTemporaryFile(suffix='.dssp', delete=False) as temp_dssp:
            try:
                result = subprocess.run(
                    ['/usr/bin/mkdssp', str(filepath), temp_dssp.name],
                    check=False,
                    capture_output=True,
                    text=True
                )
                if result.returncode != 0:
                    print(f"mkdssp error output: {result.stderr}")
                    raise Exception(f"mkdssp failed with error: {result.stderr}")
                
                # Debug: Print first 20 lines of DSSP output
                print("\nFirst 20 lines of DSSP output:")
                with open(temp_dssp.name, 'r') as f:
                    for i, line in enumerate(f):
                        if i < 20:
                            print(f"Line {i+1}: {line.rstrip()}")
                
                # Parse DSSP output
                ss8 = []
                with open(temp_dssp.name, 'r') as f:
                    for line in f:
                        # Look for the start of the secondary structure section
                        # This might need to be adjusted based on the actual output format
                        if line.startswith('  #  RESIDUE'):
                            print("Found header line")
                            break
                    
                    for line in f:
                        if len(line) > 16:
                            ss_type = line[16]
                            if ss_type == ' ' or ss_type == 'S':  # Map space and bend to coil
                                ss_type = 'C'
                            ss8.append(ss_type)
                    print(f"Found {len(ss8)} secondary structure assignments")
                
                return sasa, ss8
            finally:
                if os.path.exists(temp_dssp.name):
                    os.unlink(temp_dssp.name)
        
    except Exception as e:
        print(f"Error calculating SASA/SS8 for {filepath}: {str(e)}")
        return None, None

def check_globularity(filepath):
    """Process a single structure file"""
    try:
        # Load structure with B-factors
        cif_file = CIFFile.read(str(filepath))
        
        # Create structure from CIF data, explicitly requesting b_factor field
        struct = strucio.pdbx.get_structure(
            cif_file,
            model=1,
            extra_fields=["b_factor"]  # Explicitly request B-factors
        )
        
        # Get mean pLDDT and apply first filter
        mean_plddt = get_mean_plddt(struct)
        if mean_plddt < 0.7:
            return filepath.name, False, mean_plddt, None, None
            
        # Get CA coordinates
        ca_mask = struct.atom_name == "CA"
        ca_coords = struct.coord[ca_mask]
        L = len(ca_coords)
        
        # Skip small proteins
        if L < 30:
            return filepath.name, False, mean_plddt, None, None
            
        # Calculate long-range contacts
        dists = np.linalg.norm(ca_coords[:, None] - ca_coords, axis=2)
        seq_dists = np.abs(np.arange(L)[:, None] - np.arange(L))
        long_range_contacts = np.sum((seq_dists >= 12) & (dists <= 8.0))
        is_globular = long_range_contacts >= (L/2)
        
        # Calculate SASA and SS8 only if structure passes filters
        if is_globular:
            sasa, ss8 = calculate_sasa_ss8(filepath, struct)
        else:
            sasa, ss8 = None, None
        
        return filepath.name, is_globular, mean_plddt, sasa, ss8
        
    except Exception as e:
        print(f"Error processing {filepath}: {str(e)}")
        return filepath.name, False, 0, None, None

def process_batch(file_batch):
    """Process a batch of files"""
    return [check_globularity(f) for f in file_batch]

def main():
    # Directory with CIF files
    af_dir = Path('./tartest')  # Changed to your test directory
    all_files = list(af_dir.glob('*.cif'))
    print(f"Found {len(all_files)} CIF files to process")
     
    # For testing with small number of files, we can process them directly
    # without multiprocessing to better see any errors
    results = [check_globularity(f) for f in tqdm(all_files, desc="Processing structures")]
    
    # Save results
    with open('structure_analysis_results.txt', 'w') as f:
        f.write("filename\tis_globular\tplddt\tsasa\tss8\n")
        for name, is_glob, plddt, sasa, ss8 in results:
            sasa_str = ','.join(map(str, sasa)) if sasa is not None else ''
            ss8_str = ''.join(ss8) if ss8 is not None else ''
            f.write(f"{name}\t{is_glob}\t{plddt}\t{sasa_str}\t{ss8_str}\n")
    
    # Print summary
    total = len(results)
    passed = sum(1 for r in results if r[1])  # count globular structures
    print(f"\nProcessed {total} structures")
    print(f"Passed filters: {passed} ({passed/total*100:.1f}%)")

if __name__ == '__main__':
    main()