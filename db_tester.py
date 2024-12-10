import plyvel
import random
import json
from typing import List
from pathlib import Path

def has_interesting_features(entry: dict) -> bool:
    """Check if entry has any interesting structural or functional features."""
    features = [
        entry.get('catalytic_activity', []),
        entry.get('active_sites', []),
        entry.get('binding_sites', []),
        entry.get('interactions', []),
        entry.get('cofactors', [])
    ]
    # Return True if any of these lists are non-empty
    return any(feature for feature in features if feature)

def get_entries_with_features(db_path: str, num_samples: int = 5, max_attempts: int = 1000) -> List[dict]:
    """
    Retrieve entries that have interesting structural or functional features.
    
    Args:
        db_path: Path to the LevelDB database
        num_samples: Number of entries to retrieve
        max_attempts: Maximum number of random samples to try
    
    Returns:
        List of protein entries as dictionaries
    """
    db = plyvel.DB(db_path, create_if_missing=False)
    
    try:
        # Get all keys
        keys = [key for key in db.iterator(include_value=False)]
        total_entries = len(keys)
        print(f"Total entries in database: {total_entries}")
        
        # Keep sampling until we find enough interesting entries
        interesting_entries = []
        attempts = 0
        
        while len(interesting_entries) < num_samples and attempts < max_attempts:
            key = random.choice(keys)
            value = db.get(key)
            if value:
                entry = json.loads(value.decode())
                entry['id'] = key.decode()
                
                if has_interesting_features(entry):
                    interesting_entries.append(entry)
                    print(f"Found entry with features: {entry['id']}")
            
            attempts += 1
        
        print(f"\nFound {len(interesting_entries)} entries with features after {attempts} attempts")
        return interesting_entries
    
    finally:
        db.close()

def main():
    DB_PATH = "./uniprot_protein_db"
    OUTPUT_FILE = "interesting_entries.json"
    
    # Get entries with interesting features
    entries = get_entries_with_features(DB_PATH, num_samples=5)
    
    # Write to JSON file with nice formatting
    with open(OUTPUT_FILE, 'w') as f:
        json.dump(entries, f, indent=2)
    
    # Print a summary of what we found
    for entry in entries:
        print(f"\nEntry: {entry['id']}")
        print("Features found:")
        if entry.get('catalytic_activity'): print(f"- Catalytic activities: {len(entry['catalytic_activity'])}")
        if entry.get('active_sites'): print(f"- Active sites: {len(entry['active_sites'])}")
        if entry.get('binding_sites'): print(f"- Binding sites: {len(entry['binding_sites'])}")
        if entry.get('interactions'): print(f"- Interactions: {len(entry['interactions'])}")
        if entry.get('cofactors'): print(f"- Cofactors: {len(entry['cofactors'])}")
    
    print(f"\nWrote full entries to {OUTPUT_FILE}")

if __name__ == "__main__":
    main()