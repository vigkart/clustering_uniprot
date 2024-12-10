import plyvel
import os
import logging
from pathlib import Path
from tqdm import tqdm
import psutil
import json
from dataclasses import dataclass, asdict
from typing import Dict, List, Optional

@dataclass
class ActiveSite:
    position: int
    amino_acid: str
    description: str

@dataclass
class Interaction:
    partner: str
    interaction_type: str
    method: str
    evidence: str
    references: List[str]

@dataclass
class Cofactor:
    name: str
    description: str
    evidence: str
    chebi_id: Optional[str]

@dataclass
class UniprotEntry:
    # Primary fields
    id: str
    accessions: List[str]
    description: str
    gene_names: List[str]
    organism: str
    taxonomy: List[str]
    references: List[Dict]
    comments: List[Dict]
    features: List[Dict]
    keywords: List[str]
    sequence: str
    sequence_length: int
    sequence_mass: int
    
    # Function annotations
    ec_numbers: List[str]
    go_terms: List[Dict]
    database_refs: Dict[str, List[str]]
    catalytic_activity: List[str]
    pathway: List[str]
    
    # Structural features
    active_sites: List[ActiveSite]
    binding_sites: List[Dict]
    interactions: List[Interaction]
    cofactors: List[Cofactor]
    subunit_structure: str

class UniprotDatReader:
    def __init__(self, db_path: str, dat_file_path: str, batch_limit: int = 10000):
        # Setup logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler('uniprot_ingestion.log'),
                logging.StreamHandler()
            ]
        )
        
        self.db_path = Path(db_path)
        self.dat_file_path = Path(dat_file_path)
        self.batch_limit = batch_limit
        
        # Initialize LevelDB with optimized settings
        self.db = plyvel.DB(
            str(self.db_path),
            create_if_missing=True,
            write_buffer_size=4*1024*1024*1024,  # 4GB
            max_open_files=10000,
            lru_cache_size=64*1024*1024*1024,  # 64GB
            block_size=4*1024,
            compression=None  # Disable compression for speed
        )
        
        # State variables
        self.current_entry = None
        self.sequence_lines = []
        self.current_interaction = None
        self.current_cofactor = None
        self.parsing_interaction = False
        self.parsing_cofactor = False

    def _initialize_entry(self):
        """Initialize a new UniprotEntry with empty fields"""
        return UniprotEntry(
            id="",
            accessions=[],
            description="",
            gene_names=[],
            organism="",
            taxonomy=[],
            references=[],
            comments=[],
            features=[],
            keywords=[],
            sequence="",
            sequence_length=0,
            sequence_mass=0,
            ec_numbers=[],
            go_terms=[],
            database_refs={},
            catalytic_activity=[],
            pathway=[],
            active_sites=[],
            binding_sites=[],
            interactions=[],
            cofactors=[],
            subunit_structure=""
        )

    def _parse_id_line(self, line: str) -> None:
        """Parse ID line for entry name"""
        parts = line.split()
        self.current_entry.id = parts[1]

    def _parse_ac_line(self, line: str) -> None:
        """Parse AC line for accession numbers"""
        accessions = line[5:].replace(" ", "").rstrip().rstrip(";").split(";")
        self.current_entry.accessions.extend([acc for acc in accessions if acc])

    def _parse_de_line(self, line: str) -> None:
        """Parse DE line for description and EC numbers"""
        if "EC=" in line:
            ec_start = line.find("EC=")
            ec_end = line.find(";", ec_start)
            if ec_end == -1:
                ec_end = len(line)
            ec_numbers = line[ec_start+3:ec_end].split(" ")
            self.current_entry.ec_numbers.extend([ec for ec in ec_numbers if ec])
        else:
            desc_line = line[5:].rstrip()
            if self.current_entry.description:
                self.current_entry.description += " " + desc_line
            else:
                self.current_entry.description = desc_line

    def _parse_ft_line(self, line: str) -> None:
        """Parse feature table lines"""
        parts = line[5:].split()
        if not parts:
            return
            
        feature_type = parts[0]
        if feature_type == "ACT_SITE":
            # Handle possible position ranges and special positions
            position_str = parts[1]
            try:
                # Handle isoform-specific positions (e.g., "P12821-3:414")
                if ":" in position_str:
                    position_str = position_str.split(":")[1]
                    
                if ".." in position_str:
                    start_str, end_str = position_str.split("..")
                    # Remove special characters and convert to int
                    start = int(start_str.strip("<>"))
                    end = int(end_str.strip("<>"))
                    position = start  # Use start position for active sites
                else:
                    # Remove special characters for single positions
                    position = int(position_str.strip("<>"))
                    
                aa = parts[2] if len(parts) > 2 else ""
                description = " ".join(parts[3:]) if len(parts) > 3 else ""
                
                active_site = ActiveSite(
                    position=position,
                    amino_acid=aa,
                    description=description
                )
                self.current_entry.active_sites.append(active_site)
                
            except ValueError:
                logging.warning(f"Could not parse position: {position_str}")
                return
            
        elif feature_type in ["BINDING", "METAL"]:
            # Handle possible position ranges and special positions
            position_str = parts[1]
            try:
                # Handle isoform-specific positions (e.g., "P12821-3:414")
                if ":" in position_str:
                    position_str = position_str.split(":")[1]
                    
                if ".." in position_str:
                    start_str, end_str = position_str.split("..")
                    # Remove special characters and convert to int
                    start = int(start_str.strip("<>"))
                    end = int(end_str.strip("<>"))
                    position = start  # Use start position for binding sites
                else:
                    # Remove special characters for single positions
                    position = int(position_str.strip("<>"))
                    
                binding_site = {
                    "type": feature_type,
                    "position": position,
                    "amino_acid": parts[2] if len(parts) > 2 else "",
                    "ligand": " ".join(parts[3:]) if len(parts) > 3 else ""
                }
                self.current_entry.binding_sites.append(binding_site)
                
            except ValueError:
                logging.warning(f"Could not parse position: {position_str}")
                return

    def _parse_cc_line(self, line: str) -> None:
        """Parse comment lines"""
        comment = line[5:].rstrip()
        
        if "-!- CATALYTIC ACTIVITY:" in comment:
            activity = comment.replace("-!- CATALYTIC ACTIVITY:", "").strip()
            self.current_entry.catalytic_activity.append(activity)
        elif "-!- PATHWAY:" in comment:
            pathway = comment.replace("-!- PATHWAY:", "").strip()
            self.current_entry.pathway.append(pathway)
        elif "-!- COFACTOR:" in comment:
            self._parse_cc_cofactor(line)
        elif "-!- INTERACTION" in comment:
            self._parse_cc_interaction(line)
        elif "-!- SUBUNIT:" in comment:
            self.current_entry.subunit_structure = comment.replace("-!- SUBUNIT:", "").strip()
        else:
            self.current_entry.comments.append(comment)

    def _parse_cc_interaction(self, line: str) -> None:
        """Parse interaction data"""
        line = line[5:].strip()
        
        if line.startswith("-!- INTERACTION"):
            self.parsing_interaction = True
            self.current_interaction = {
                "partners": [],
                "type": "",
                "method": "",
                "evidence": "",
                "references": []
            }
        elif self.parsing_interaction:
            if line.startswith("-!-"):
                self.parsing_interaction = False
                if self.current_interaction:
                    interaction = Interaction(
                        partner="; ".join(self.current_interaction["partners"]),
                        interaction_type=self.current_interaction["type"],
                        method=self.current_interaction["method"],
                        evidence=self.current_interaction["evidence"],
                        references=self.current_interaction["references"]
                    )
                    self.current_entry.interactions.append(interaction)
            else:
                if "Xref=" in line:
                    self.current_interaction["references"].append(
                        line.split("Xref=")[1].strip()
                    )
                elif ";" in line:
                    parts = line.split(";")
                    for part in parts:
                        if ":" in part:
                            key, value = part.split(":", 1)
                            key = key.strip().lower()
                            value = value.strip()
                            if key == "with":
                                self.current_interaction["partners"].append(value)
                            elif key == "method":
                                self.current_interaction["method"] = value
                            elif key == "evidence":
                                self.current_interaction["evidence"] = value

    def _parse_cc_cofactor(self, line: str) -> None:
        """Parse cofactor information"""
        line = line[5:].strip()
        
        if line.startswith("-!- COFACTOR"):
            self.parsing_cofactor = True
            self.current_cofactor = {
                "name": "",
                "description": "",
                "evidence": "",
                "chebi_id": None
            }
        elif self.parsing_cofactor:
            if line.startswith("-!-"):
                self.parsing_cofactor = False
                if self.current_cofactor:
                    cofactor = Cofactor(
                        name=self.current_cofactor["name"],
                        description=self.current_cofactor["description"],
                        evidence=self.current_cofactor["evidence"],
                        chebi_id=self.current_cofactor["chebi_id"]
                    )
                    self.current_entry.cofactors.append(cofactor)
            else:
                if "Name=" in line:
                    self.current_cofactor["name"] = line.split("Name=")[1].split(";")[0].strip()
                elif "Xref=ChEBI:" in line:
                    self.current_cofactor["chebi_id"] = line.split("ChEBI:")[1].split(";")[0].strip()
                else:
                    self.current_cofactor["description"] += " " + line.strip()

    def _parse_dr_line(self, line: str) -> None:
        """Parse database reference lines"""
        parts = line[5:].split(';')
        db_name = parts[0].strip()
        refs = [p.strip().rstrip('.') for p in parts[1:]]
        
        if db_name == "GO":
            go_term = {
                'id': refs[0],
                'term': refs[1],
                'evidence': refs[2]
            }
            self.current_entry.go_terms.append(go_term)
        else:
            if db_name not in self.current_entry.database_refs:
                self.current_entry.database_refs[db_name] = []
            self.current_entry.database_refs[db_name].extend(refs)

    def _parse_sq_line(self, line: str) -> None:
        """Parse sequence lines"""
        if line.startswith("SQ"):
            parts = line.split()
            self.current_entry.sequence_length = int(parts[2])
            self.current_entry.sequence_mass = int(parts[4])
        else:
            # Sequence lines start with spaces and contain the actual sequence
            # Remove all spaces and numbers
            sequence_part = ''.join(c for c in line if c.isalpha())
            if sequence_part:
                self.sequence_lines.append(sequence_part)
            
    def _parse_gn_line(self, line: str) -> None:
        """Parse gene name lines"""
        parts = line[5:].split(';')
        for part in parts:
            if "Name=" in part:
                gene_name = part.split('=')[1].strip()
                self.current_entry.gene_names.append(gene_name)

    def _parse_ox_line(self, line: str) -> None:
        """Parse organism taxonomy line"""
        if "NCBI_TaxID=" in line:
            taxid = line.split('=')[1].strip().rstrip(';')
            self.current_entry.taxonomy.append(taxid)

    def _parse_os_line(self, line: str) -> None:
        """Parse organism species line"""
        organism = line[5:].strip().rstrip('.')
        if self.current_entry.organism:
            self.current_entry.organism += ' ' + organism
        else:
            self.current_entry.organism = organism

    def _validate_entry(self) -> bool:
        """Validate required fields before saving"""
        if not self.current_entry.id:
            logging.warning("Entry missing ID, skipping")
            return False
        if not self.current_entry.accessions:
            logging.warning(f"Entry {self.current_entry.id} missing accessions")
        if not self.current_entry.sequence:
            logging.warning(f"Entry {self.current_entry.id} missing sequence")
            return False
        return True

    def _log_memory_usage(self):
        """Log current memory usage"""
        process = psutil.Process(os.getpid())
        mem_usage = process.memory_info().rss / 1024 / 1024  # in MB
        logging.info(f"Current memory usage: {mem_usage:.2f} MB")

    def process_dat_file(self):
        """Process the UniProt DAT file"""
        batch = self.db.write_batch()
        entries_processed = 0
        
        logging.info(f"Starting to process: {self.dat_file_path}")
        
        with open(self.dat_file_path, 'r') as dat_file:
            for line in tqdm(dat_file):
                line = line.rstrip()
                
                if line.startswith("ID"):
                    if self.current_entry is not None:
                        # Process previous entry
                        if self.sequence_lines:
                            self.current_entry.sequence = "".join(self.sequence_lines)
                            logging.debug(f"Processed sequence for {self.current_entry.id}: length={len(self.current_entry.sequence)}")
                        else:
                            logging.debug(f"No sequence lines found for {self.current_entry.id}")
                    
                    # Initialize new entry
                    self.current_entry = self._initialize_entry()
                    self.sequence_lines = []
                    self._parse_id_line(line)
                
                # Parse different line types
                elif line.startswith("AC"):
                    self._parse_ac_line(line)
                elif line.startswith("DE"):
                    self._parse_de_line(line)
                elif line.startswith("DR"):
                    self._parse_dr_line(line)
                elif line.startswith("CC"):
                    self._parse_cc_line(line)
                elif line.startswith("FT"):
                    self._parse_ft_line(line)
                elif line.startswith("SQ") or (line and line[0] == " "):
                    self._parse_sq_line(line)
                elif line.startswith("//"):
                    # End of entry
                    if self.current_entry:
                        if self.sequence_lines:
                            self.current_entry.sequence = "".join(self.sequence_lines)
                            logging.debug(f"Final sequence for {self.current_entry.id}: length={len(self.current_entry.sequence)}")
                        if self._validate_entry():
                            entry_key = self.current_entry.id.encode()
                            entry_value = json.dumps(asdict(self.current_entry)).encode()
                            batch.put(entry_key, entry_value)
                            entries_processed += 1
                elif line.startswith("GN"):
                    self._parse_gn_line(line)
                elif line.startswith("OS"):
                    self._parse_os_line(line)
                elif line.startswith("OX"):
                    self._parse_ox_line(line)

        # Write final batch
        batch.write()
        logging.info(f"Successfully processed {entries_processed} entries")
        
    def close(self):
        """Close the database connection"""
        self.db.close()

def main():
    # Configuration
    DB_PATH = "./final_db"  # Local directory for LevelDB storage
    DAT_FILE = "./data/uniprot_trembl.dat"
    # DAT_FILE = "./data/uniprot_sprot.dat"  # Path to your SPROT DAT file 
    BATCH_SIZE = 50000  # Smaller batch size for testing
    
    # Initialize and run ingestion
    reader = UniprotDatReader(
        db_path=DB_PATH,
        dat_file_path=DAT_FILE,
        batch_limit=BATCH_SIZE
    )
    
    try:
        reader.process_dat_file()
    finally:
        reader.close()

if __name__ == "__main__":
    main()