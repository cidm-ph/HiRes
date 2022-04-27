import os
import argparse
import pandas as pd
import vcf2ciprores as v2cr


parser = argparse.ArgumentParser(description='HiRes')
parser.add_argument('input', help='Input file')
args = vars(parser.parse_args())
output_csvs= []

for file in os.listdir(args['input']):
    filename = os.path.join(args['input'], os.fsdecode(file))
    outfile = os.path.join(
        args['input'], os.path.splitext(os.path.basename(file))[0] + '.ciprores.vcf'
    )
    if filename.endswith((".vcf")) and os.stat(filename).st_size != 0:
        results = v2cr.reference_split(filename)
        if results is not None and results.empty is False:
            results.to_csv(outfile, sep='\t', index = False)

