import logging
from scripts import assists

def cipro_checking(results):
    cipro_resistance = {
        ("gyrA:D88Y",): "Intermediate resistance",
        ("gyrA:D88N",): "Intermediate resistance",
        ("gyrA:S84L",): "Intermediate resistance",
        ("gyrA:S84F", "parC:S84R"): "Intermediate resistance (Reduced Susceptibility - 'Group B')",
        ("gyrA:S84F", "gyrA:D420N"): "Intermediate resistance (Reduced Susceptibility - 'Group B')",
        ("gyrA:S84L", "parC:S84R"): "Intermediate resistance (Reduced Susceptibility - 'Group B')",
        ("gyrA:S84L", "parC:S84I"): "Resistant; Predicted MIC: >0.5 ug/mL",
        ("gyrA:S84L", "gyrA:E88K"): "Resistant; Predicted MIC: 2 ug/mL",
        ("gyrA:S84L", "gyrA:D420N"): "Intermediate resistance (Reduced Susceptibility - 'Group B')",
        ("gyrA:D88N", "parC:S84I"): "Resistant; Predicted MIC: >0.5 ug/mL",
        ("gyrA:D88N", "gyrA:E88K"): "Resistant; Predicted MIC: >0.5 ug/mL",
        ("gyrA:S84Y", "gyrA:E88K"): "Resistant; Predicted MIC: >0.5 ug/mL",
        ("gyrA:S84Y", "gyrA:D88N", "gyrA:E88K"): "Resistant; Predicted MIC: >0.5 ug/mL",
        ("gyrA:S84Y", "gyrA:D88Y", "gyrA:D83G", "gyrA:S84A"): "Resistant; Predicted MIC: 2 ug/mL",
        ("gyrA:S84L", "gyrA:D88Y", "gyrA:S84N"): "Resistant; Predicted MIC: 8 ug/mL",
        ("gyrA:S84L", "gyrA:D88Y", "gyrA:E88K"): "Resistant; Predicted MIC: 8-16 ug/mL",
        ("gyrA:S84L", "gyrA:D88N", "gyrA:E82D"): "Resistant; Predicted MIC: 16 ug/mL",
        ("gyrA:S84L", "gyrA:D88N", "parC:S84R"): "Resistant; Predicted MIC: 16 ug/mL",
        ("gyrA:S84Y", "gyrA:D88G", "parC:S84I", "gyrA:D420N"): "Resistant; Predicted MIC: 16 ug/mL",
        ("gyrA:S84L", "gyrA:D88G", "parC:S84I", "gyrA:D420N"): "Resistant; Predicted MIC: 16 ug/mL",
        ("gyrA:S84L", "gyrA:D88N", "parC:S84I", "gyrA:D420N"): "Resistant; Predicted MIC: 16 ug/mL",
        ("gyrA:S84L", "gyrA:D88Y", "gyrA:G82C", "gyrA:E88K"): "Resistant; Predicted MIC: 16 ug/mL",
        ("gyrA:S84Y", "gyrA:D88Y", "parC:S84I", "gyrA:E88K"): "Resistant; Predicted MIC: 16 ug/mL",
        ("gyrA:D88Y", "parC:S84I"): "Resistant; Predicted MIC: >0.5 ug/mL",
        ("gyrA:D88Y", "gyrA:E88K"): "Resistant; Predicted MIC: 2 ug/mL",
        ("gyrA:E83C", "parC:S84R"): "Resistant; Predicted MIC: >0.5 ug/mL",
        ("gyrA:S84L", "gyrA:D88N", "parC:S84I"): "Resistant; Predicted MIC: >0.5 ug/mL",
        ("gyrA:A426V",): "Resistant",
    }
    
    gyra_info = results[results['Locus_Tag'].str.contains('HI_1264')].reset_index()
    parc_info = results[results['Locus_Tag'].str.contains('HI_1529')].reset_index()
    pare_info = results[results['Locus_Tag'].str.contains('HI_1528')].reset_index()

    info_list = [
        (gyra_info, "gyrA"),
        (parc_info, "parC"),
        (pare_info, "parE"),
    ]
    results_with_prefix = []  # Collect sets of prefixed mutations
    for info, prefix in info_list:
        # Convert PROT-POS values and add prefix
        info['PROT-POS-1'] = info['PROT-POS'].apply(assists.convert)
        info['PROT-POS-1'] = info['PROT-POS-1'].apply(lambda x: f"{prefix}:{x}")
        # Collect results as a set
        results_with_prefix.append(set(info['PROT-POS-1'].values))

    # Combine all prefixed results
    combined_results = set.union(*results_with_prefix)

    best_match_group = None
    best_match_count = 0
    best_match_key = None

    # Loop through the resistance dictionary
    for mutations, group in cipro_resistance.items():
        # Count how many mutations from the dictionary match combined_results
        matched_count = len(set(mutations).intersection(combined_results))
        
        # Update the best match if this group has more matched mutations
        if matched_count > best_match_count:
            best_match_count = matched_count
            best_match_group = group
            best_match_key = '+'.join(sorted(mutations))  # Concatenate mutations for key

    # Output the best match
    if best_match_group:
        logging.info(f"Best match: {best_match_group} with {best_match_count} mutations matched.")
        logging.info(f"Mutations (key): {best_match_key}")
    else:
        logging.info("No match found.")

    return best_match_group, best_match_key