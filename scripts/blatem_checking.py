import logging
import pandas as pd
from scripts import assists

def blatem_checking(assembly, analysis_outdir):
    blatem_result, blatem_res = None, None
    abricate_cmd = f"abricate --db card --quiet {assembly} > {analysis_outdir}/card.txt"
    assists.run_cmd(abricate_cmd)
    try:
        card_info = pd.read_csv(f"{analysis_outdir}/card.txt", sep="\t", header=0)
    except Exception as e:
        logging.warning(f"Failed to load one or more files: {e}")

    if not card_info.empty:
        blatem_info = card_info[card_info['GENE'].str.contains('TEM-1')].reset_index()
        if not blatem_info.empty:
            blatem_result = "blaTEM-1"
            blatem_res = "BLPAR"
    return blatem_result, blatem_res