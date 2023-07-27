import os
import argparse
import vcf2ciprores as v2cr


parser = argparse.ArgumentParser(description='HiRes')
parser.add_argument('--input', '-i', help='Input file/folder')
parser.add_argument('--output', '-o', help='Output folder')
args = vars(parser.parse_args())
output_csvs= []

def main():
    file_list = os.listdir(args['input'])
    if len(file_list) < 1:
        file_run(file)
    else:
        for file in os.listdir(args['input']):
            file_run(file)

def file_run(file):
    filename = os.path.join(args['input'], os.fsdecode(file))
    outfile = os.path.join(
        args['input'], os.path.splitext(os.path.basename(file))[0] + '.ciprores.vcf'
    )
    if filename.endswith((".vcf")) and os.stat(filename).st_size != 0:
        results = v2cr.reference_split(filename)
        if results is not None and results.empty is False:
            results.to_csv(outfile, sep='\t', index = False)

if __name__ == "__main__":
    main()
