import plyvel
import os
import logging
from pathlib import Path
from tqdm import tqdm
import psutil
import json
from dataclasses import dataclass, asdict, field
from typing import Dict, List, Optional, Union
import time

@dataclass
class ActiveSite:
    position: int
    amino_acid: str
    description: str
    evidence: str = ""
    annotations: List[str] = None

    def __post_init__(self):
        if self.annotations is None:
            self.annotations = []

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
    # Primary fields (no defaults)
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
    ec_numbers: List[str]
    go_terms: List[Dict]
    database_refs: Dict[str, List[str]]
    pathway: List[str]
    active_sites: List[ActiveSite]
    binding_sites: List[Dict]
    interactions: List[Interaction]
    cofactors: List[Cofactor]
    subunit_structure: str
    
    # Fields with defaults must come last
    catalytic_activity: List[Dict[str, Union[str, List[str]]]] = field(default_factory=list)

class UniprotDatReader:
    def __init__(self, db_path: str, dat_file_path: str, batch_limit: int = 10000, 
                 max_memory_percent: float = 80.0, resume: bool = True):
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
        self.max_memory_percent = max_memory_percent
        self.resume = resume
        
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
        
        # Resume state
        self.progress_file = self.db_path / "progress.json"
        self.last_processed_position = 0
        self._load_progress()

    def _load_progress(self):
        """Load progress from previous run if it exists"""
        if self.resume and self.progress_file.exists():
            try:
                with open(self.progress_file, 'r') as f:
                    progress = json.load(f)
                self.last_processed_position = progress.get('position', 0)
                logging.info(f"Resuming from position {self.last_processed_position}")
            except Exception as e:
                logging.warning(f"Could not load progress file: {e}")
                self.last_processed_position = 0

    def _save_progress(self, position: int):
        """Save current progress"""
        try:
            with open(self.progress_file, 'w') as f:
                json.dump({'position': position}, f)
        except Exception as e:
            logging.warning(f"Could not save progress: {e}")

    def _check_memory_usage(self) -> bool:
        """
        Check if memory usage is too high
        Returns True if we should pause processing
        """
        memory_percent = psutil.virtual_memory().percent
        if memory_percent > self.max_memory_percent:
            logging.warning(f"Memory usage too high ({memory_percent}%). Pausing processing.")
            return True
        return False

    def _wait_for_memory(self, initial_wait: int = 60, max_wait: int = 300):
        """Wait for memory usage to decrease"""
        wait_time = initial_wait
        while self._check_memory_usage():
            logging.info(f"Waiting {wait_time} seconds for memory to free up...")
            time.sleep(wait_time)
            wait_time = min(wait_time * 2, max_wait)  # Exponential backoff with cap

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
            # Strip ECO evidence codes
            desc_line = line[5:].rstrip()
            if "{ECO:" in desc_line:
                desc_line = desc_line.split("{ECO:")[0].strip()
            
            if self.current_entry.description:
                self.current_entry.description += " " + desc_line
            else:
                self.current_entry.description = desc_line

    def _parse_ft_line(self, line: str) -> None:
        """Parse feature table lines"""
        try:
            parts = line[5:].split(maxsplit=3)
            if len(parts) < 2:  # Need at least type and position
                return
            
            feature_type = parts[0]
            position_str = parts[1]
            description = parts[3] if len(parts) > 3 else ""
            aa = parts[2] if len(parts) > 2 else ""

            # Create general feature entry
            feature = {
                "type": feature_type,
                "position": position_str,  # Store original position string
                "amino_acid": aa,
                "description": description,
                "evidence": ""
            }

            # Only try to parse numeric positions for specific feature types
            if feature_type in ["ACT_SITE", "BINDING", "METAL"]:
                try:
                    # Handle isoform-specific positions
                    if ":" in position_str:
                        position_str = position_str.split(":")[1]
                    
                    # Parse position
                    if ".." in position_str:
                        start_str, end_str = position_str.split("..")
                        start = int(start_str.strip('<>'))
                        end = int(end_str.strip('<>'))
                        feature["position"] = f"{start}..{end}"
                    else:
                        feature["position"] = int(position_str.strip("<>"))
                except ValueError:
                    # Keep original position string if not numeric
                    pass

            # Parse evidence if present
            if "/evidence=" in description:
                evidence_part = description[description.find("/evidence="):].split()[0]
                feature["evidence"] = evidence_part.replace("/evidence=", "").strip('"')
                description = description.split("/evidence=")[0].strip()

            # Add to general features list
            self.current_entry.features.append(feature)

            # Handle specific feature types
            if feature_type == "ACT_SITE":
                try:
                    pos = int(position_str.strip("<>"))
                    active_site = ActiveSite(
                        position=pos,
                        amino_acid=aa,
                        description=description,
                        evidence=feature["evidence"]
                    )
                    self.current_entry.active_sites.append(active_site)
                except ValueError:
                    logging.warning(f"Skipping non-numeric active site position: {position_str}")
                    
            elif feature_type in ["BINDING", "METAL"]:
                binding_site = {
                    "type": feature_type,
                    "position": feature["position"],
                    "amino_acid": aa,
                    "ligand": description.split("/")[0].strip(),
                    "evidence": feature["evidence"],
                    "note": ""
                }
                
                if "/note=" in description:
                    note_part = description[description.find("/note="):].split(";")[0]
                    binding_site["note"] = note_part.replace("/note=", "").strip('"')
                    
                self.current_entry.binding_sites.append(binding_site)

        except Exception as e:
            logging.warning(f"Could not parse feature line '{line[:50]}...': {str(e)}")

    def _parse_cc_line(self, line: str) -> None:
        """Parse comment lines"""
        comment = line[5:].rstrip()
        
        if "CATALYTIC ACTIVITY:" in comment:
            if "-!- CATALYTIC ACTIVITY:" in comment:
                activity = comment.replace("-!- CATALYTIC ACTIVITY:", "").strip()
                
                # Check format
                if "Reaction=" in activity:
                    # Complex format
                    self.current_catalytic = {
                        "reaction": "",
                        "ec_number": "",
                        "description": activity,
                        "xrefs": []
                    }
                    self.parsing_catalytic = True
                else:
                    # Simple format
                    if "{ECO:" in activity:
                        activity = activity.split("{ECO:")[0].strip()
                    self.current_entry.catalytic_activity.append({
                        "description": activity,
                        "reaction": "",
                        "ec_number": "",
                        "xrefs": []
                    })
            elif self.parsing_catalytic:
                if line.startswith("CC       "):
                    if "Reaction=" in line:
                        reaction = line[line.find("Reaction="):].split(";")[0]
                        self.current_catalytic["reaction"] = reaction.replace("Reaction=", "").strip()
                    if "EC=" in line:
                        ec = line[line.find("EC="):].split(";")[0]
                        self.current_catalytic["ec_number"] = ec.replace("EC=", "").strip()
                    if "Xref=" in line:
                        xref_parts = line[line.find("Xref="):].split(";")[0].split(",")
                        for xref in xref_parts:
                            if "ChEBI:" in xref:
                                self.current_catalytic["xrefs"].append(xref.strip())
                else:
                    # End of section
                    self.parsing_catalytic = False
                    if self.current_catalytic:
                        self.current_entry.catalytic_activity.append(self.current_catalytic)
        elif "-!- PATHWAY:" in comment:
            pathway = comment.replace("-!- PATHWAY:", "").strip()
            self.current_entry.pathway.append(pathway)
        elif "-!- COFACTOR:" in comment:
            if "-!- COFACTOR:" in comment:
                self.parsing_cofactor = True
                self.current_cofactor = {
                    "name": "",
                    "description": "",
                    "evidence": "",
                    "chebi_id": None,
                    "note": ""
                }
            elif self.parsing_cofactor:
                if line.startswith("CC       "):
                    if "Name=" in line:
                        name_part = line[line.find("Name="):].split(";")[0]
                        self.current_cofactor["name"] = name_part.replace("Name=", "").strip()
                    if "Xref=" in line:
                        xref_parts = line[line.find("Xref="):].split(";")[0].split(",")
                        for xref in xref_parts:
                            if "ChEBI:" in xref:
                                self.current_cofactor["chebi_id"] = xref.replace("ChEBI:", "").strip()
                            self.current_cofactor["xrefs"].append(xref.strip())
                    if "Note=" in line:
                        note_part = line[line.find("Note="):].split(";")[0]
                        self.current_cofactor["note"] = note_part.replace("Note=", "").strip()
                        if self.current_cofactor["description"]:
                            self.current_cofactor["description"] += " " + self.current_cofactor["note"]
                        else:
                            self.current_cofactor["description"] = self.current_cofactor["note"]
                else:
                    # End of cofactor section
                    self.parsing_cofactor = False
                    if self.current_cofactor:
                        cofactor = Cofactor(
                            name=self.current_cofactor["name"],
                            description=self.current_cofactor["description"].strip(),
                            evidence=self.current_cofactor["evidence"],
                            chebi_id=self.current_cofactor["chebi_id"]
                        )
                        self.current_entry.cofactors.append(cofactor)
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

    def _parse_dr_line(self, line: str) -> None:
        """Parse database reference lines"""
        parts = line[5:].split(';')
        db_name = parts[0].strip()
        refs = [p.strip().rstrip('.') for p in parts[1:]]
        
        if db_name == "GO":
            go_id = refs[0]
            term = refs[1]
            evidence = refs[2] if len(refs) > 2 else ""
            qualifier = ""
            
            # Handle GO term qualifiers
            if "|" in term:
                qualifier, term = term.split("|", 1)
            
            go_term = {
                'id': go_id,
                'term': term,
                'evidence': evidence,
                'qualifier': qualifier
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
            # Initialize sequence for this entry
            if not hasattr(self.current_entry, '_sequence_parts'):
                self.current_entry._sequence_parts = []
        else:
            # Sequence lines start with spaces and contain the actual sequence
            # Remove all spaces and numbers
            sequence_part = ''.join(c for c in line if c.isalpha())
            if sequence_part:
                if not hasattr(self.current_entry, '_sequence_parts'):
                    self.current_entry._sequence_parts = []
                self.current_entry._sequence_parts.append(sequence_part)
                # Join all parts when we have them
                self.current_entry.sequence = ''.join(self.current_entry._sequence_parts)

    def _parse_gn_line(self, line: str) -> None:
        """Parse gene name lines"""
        parts = line[5:].split(';')
        for part in parts:
            if "Name=" in part:
                gene_name = part.split('=')[1].strip()
                # Strip evidence codes
                if "{ECO:" in gene_name:
                    gene_name = gene_name.split("{ECO:")[0].strip()
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

    def _parse_kw_line(self, line: str) -> None:
        """Parse keyword lines"""
        keywords = line[5:].replace(" ", "").rstrip().rstrip(";").split(";")
        self.current_entry.keywords.extend([kw for kw in keywords if kw])

    def _parse_rx_line(self, line: str) -> None:
        """Parse reference line"""
        ref = {}
        parts = line[5:].split(';')
        for part in parts:
            if "=" in part:
                key, value = part.split('=', 1)
                ref[key.strip()] = value.strip()
        if ref:
            self.current_entry.references.append(ref)

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

    def _process_line(self, line):
        """Process a line based on its type"""
        try:
            if not line or len(line) < 5:  # Ensure line has minimum length
                return
            
            line_type = line[:2].strip()  # Get first two chars for type
            
            # Initialize new entry when we see an ID line
            if line_type == "ID":
                self.current_entry = self._initialize_entry()
                self._parse_id_line(line)
                return
                
            # Skip processing if we don't have a current entry
            if not self.current_entry:
                return
                
            if line_type == "AC":
                self._parse_ac_line(line)
            elif line_type == "DE":
                self._parse_de_line(line)
            elif line_type == "DR":
                self._parse_dr_line(line)
            elif line_type == "CC":
                self._parse_cc_line(line)
            elif line_type == "FT":
                self._parse_ft_line(line)
            elif line_type == "KW":
                self._parse_kw_line(line)
            elif line_type == "GN":
                self._parse_gn_line(line)
            elif line_type == "OS":
                self._parse_os_line(line)
            elif line_type == "OX":
                self._parse_ox_line(line)
            elif line_type == "SQ" or line[0] == " ":
                self._parse_sq_line(line)
            elif line_type == "//":  # End of entry
                if self._validate_entry():
                    self._save_entry()
                self.current_entry = None
        except Exception as e:
            logging.warning(f"Error processing line '{line[:50]}...': {str(e)}")

    def process_dat_file(self):
        """Process the UniProt DAT file"""
        try:
            batch = self.db.write_batch()
            entries_processed = 0
            current_position = 0
            
            logging.info(f"Starting to process: {self.dat_file_path}")
            
            with open(self.dat_file_path, 'r') as dat_file:
                # Skip to last processed position if resuming
                if self.last_processed_position > 0:
                    dat_file.seek(self.last_processed_position)
                    # Skip partial entry if any
                    for line in dat_file:
                        if line.startswith("ID"):
                            break

                # Get file size for progress bar
                dat_file.seek(0, 2)  # Seek to end
                file_size = dat_file.tell()
                dat_file.seek(self.last_processed_position)  # Return to last position

                # Create progress bar
                pbar = tqdm(total=file_size, unit='B', unit_scale=True)
                pbar.update(self.last_processed_position)

                while True:
                    try:
                        current_position = dat_file.tell()  # Get position before reading
                        line = dat_file.readline()
                        
                        if not line:  # End of file
                            break
                            
                        line = line.rstrip()
                        pbar.update(len(line) + 1)  # +1 for newline character
                        
                        if not line:  # Skip empty lines
                            continue
                            
                        # Process line based on type
                        if line.startswith("ID"):
                            if self.current_entry is not None:
                                # Process previous entry
                                if self.sequence_lines:
                                    self.current_entry.sequence = "".join(self.sequence_lines)
                                
                            # Initialize new entry
                            self.current_entry = self._initialize_entry()
                            self.sequence_lines = []
                            self._parse_id_line(line)
                        elif line.startswith("//"):
                            # End of entry processing
                            if self.current_entry:
                                if self._validate_entry():
                                    entry_key = self.current_entry.id.encode()
                                    entry_value = json.dumps(asdict(self.current_entry)).encode()
                                    batch.put(entry_key, entry_value)
                                    entries_processed += 1
                                    
                                    if entries_processed % 10000 == 0:
                                        self._save_progress(current_position)
                                        self._log_memory_usage()
                                        
                                    if entries_processed % self.batch_limit == 0:
                                        batch.write()
                                        batch = self.db.write_batch()
                        else:
                            # Process other line types
                            self._process_line(line)
                            
                    except Exception as e:
                        logging.error(f"Error processing line: {e}")
                        continue
                        
        except Exception as e:
            logging.error(f"Error during processing: {e}")
            raise
        finally:
            pbar.close()
            # Write final batch and save progress
            batch.write()
            self._save_progress(current_position)
            logging.info(f"Successfully processed {entries_processed} entries")

    def close(self):
        """Close the database connection and cleanup"""
        self.db.close()
        if not self.resume:
            # Clean up progress file if not needed for resuming
            try:
                self.progress_file.unlink()
            except FileNotFoundError:
                pass

def main():
    # Configuration
    DB_PATH = "./uniprot_protein_db"
    DAT_FILE = "./data/uniprot_sprot.dat"
    BATCH_SIZE = 50000
    MAX_MEMORY_PERCENT = 80.0
    
    # Initialize and run ingestion
    reader = UniprotDatReader(
        db_path=DB_PATH,
        dat_file_path=DAT_FILE,
        batch_limit=BATCH_SIZE,
        max_memory_percent=MAX_MEMORY_PERCENT,
        resume=True
    )
    
    try:
        reader.process_dat_file()
    except KeyboardInterrupt:
        logging.info("Processing interrupted by user")
    except Exception as e:
        logging.error(f"Error during processing: {e}")
    finally:
        reader.close()

if __name__ == "__main__":
    main()