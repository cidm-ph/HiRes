import logging
from scripts import assists

def pbp3_checking(results):
    pbp3_resistance = {
        ("R517H"): "Group I strains",
        ("N526K"): "Group IIa strains",
        ("N350D", "G490Q", "N526K", "A530S"): "Group IIa strains",
        ("N526K", "A502V"): "Group IIb strains",
        ("N350D", "M377I", "A502V", "N526K"): "Group IIb strains",
        ("N526K", "A502T"): "Group IIc strains",
        ("N350D", "A502T", "N526K"): "Group IIc strains",
        ("N526K", "I449V"): "Group IId strains",
        ("M377I", "S385T", "N526K"): "Group III strains",
        ("M377I", "S385T", "R517H"): "Group III-like strains",
        ("M377I", "S385T", "L389F", "N526K"): "Group III+ Strains",
        ("M377I", "S385T", "L389F", "R517H"): "Group III-like+ Strains",
    }

    pbp3_info = results[results['Locus_Tag'].str.contains('HI_1132')].reset_index()
    pbp3_info['PROT-POS-1'] = pbp3_info['PROT-POS'].apply(assists.convert)
    pbp3_result = set(tuple(pbp3_info['PROT-POS-1'].values))

    # Variables to track the best match
    best_match_group = None
    best_match_count = 0
    best_match_key = None

    # Loop through the pbp3_resistance dictionary and check for matches
    for mutations, group in pbp3_resistance.items():
        # Count how many mutations from the dictionary are in php3_result
        matched_count = len(set(mutations).intersection(pbp3_result))
        
        # Update the best match if this group has more matched mutations
        if matched_count > best_match_count:
            best_match_count = matched_count
            best_match_group = group
            best_match_key = '+'.join(sorted(mutations))

    # If a best match was found, print it, otherwise print "No match found"
    if best_match_group:
        logging.info(f"Best match: {best_match_group} with {best_match_count} mutations matched.")
        pbp3_res = "BLNAR"
    else:
        logging.info(f"No match found")
        pbp3_res = "BLNAS"

    return best_match_key, pbp3_res, best_match_group
