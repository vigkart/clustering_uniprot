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
                position = {"start": start, "end": end}
            else:
                # Single position
                pos = int(position_str.strip("<>"))
                position = {"start": pos, "end": pos}
                
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
        position_str = parts[1]
        try:
            # Handle isoform-specific positions
            if ":" in position_str:
                position_str = position_str.split(":")[1]
                
            if ".." in position_str:
                start_str, end_str = position_str.split("..")
                # Remove special characters and convert to int
                start = int(start_str.strip("<>"))
                end = int(end_str.strip("<>"))
                position = {"start": start, "end": end}
            else:
                # Single position
                pos = int(position_str.strip("<>"))
                position = {"start": pos, "end": pos}
                
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