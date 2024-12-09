import datetime
import logging
import os
import sys
import warnings
import glob
from scripts import arguments
from scripts import assists
from scripts import vcf2ciprores as v2cr
from scripts import pbp3_checking as pbp3
from scripts import cipro_checking as cipro
from scripts import blatem_checking as blatem

__version__ = "0.0.1"
warnings.simplefilter(action="ignore", category=FutureWarning)
logging.getLogger().setLevel(logging.INFO)
formatter = logging.Formatter(
    "hires:%(levelname)s:%(asctime)s: %(message)s", datefmt="%y/%m/%d %I:%M:%S %p"
)

dependency_list = ["abricate", "spades.py", "snippy"]
output_csvs= []

def main(args):
    now = datetime.datetime.now()
    date = now.strftime("%Y%m%d")
    is_assembly = bool(args.fasta is not None)
    is_reads = bool(args.R1 is not None)
    is_force = bool(args.force is not None)
    forced = ""
    if is_force is True:
        forced = "--force"

    

    # set outdir defaults - if no outdir is set, it will default to either the fasta or R1 location
    if args.outdir is None and args.fasta is not None:
        default = os.path.dirname(args.fasta)
        outdir = default
    elif args.outdir is None and args.R1 is not None:
        default = os.path.dirname(args.R1)
        outdir = default
    else:
        outdir = args.outdir

    # force creation of new folder within set outdir
    maindir = outdir 
    final_dict = {
        "Folder": maindir
    }
    
    #newdir = maindir + "/bams"
    folder_exists = os.path.exists(maindir)
    if not folder_exists:
        os.makedirs(maindir)
        logging.info("Making output folder")
    else:
        logging.info(f"Folder exists")

    # error log
    errorlog = os.path.join(outdir, "hires_" + date + ".log")

        # Clear existing handlers
    logger = logging.getLogger()
    if logger.hasHandlers():
        logger.handlers.clear()

    stdout_handler = logging.StreamHandler(sys.stdout)
    file_handler = logging.FileHandler(errorlog, mode="w+")
    for handler in [stdout_handler, file_handler]:
        handler.setLevel(logging.INFO)
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    # cmd checks
    if is_reads is True:
        if args.R2 is None:
            logging.error("R2 was not provided, please provide the paired reads")
            sys.exit(1)

    # launch line
    logging.info(
        "Launching hires v%s and writing output files to directory %s",
        __version__,
        outdir,
    )

    # checking file integrity and existence of output directory
    if is_reads is True and is_assembly is True:
        logging.critical(f"Can't handle both fasta and reads, pick one or the other.")

    # checking all the versions and installations of dependencies.
    logging.info("Checking installs of dependencies")
    for dependency in dependency_list:
        assists.check_dependencies(dependency)

    snippy_outdir = maindir + "/snippy"
    folder_exists = os.path.exists(snippy_outdir)
    name = os.path.basename(maindir)
    snippy_result = False

    # Testing for PBP3 Mutations
    ref = os.path.join(os.path.dirname(__file__), "databases/RdKW20.gbk")
    if folder_exists:
        snippy_result = assists.check_snippy_finished(snippy_outdir, name)
    elif not folder_exists:
        if is_reads and is_assembly is False and snippy_result is False:
            snippy = f"snippy --prefix {name} --outdir {snippy_outdir} --ref {ref} --R1 {args.R1} --R2 {args.R1} --report {forced}"
            assists.run_cmd(snippy)
        elif is_assembly and is_reads is False and snippy_result is False:
            snippy = f"snippy --prefix {name} --outdir {snippy_outdir} --ref {ref} --ctgs {args.fasta} --report {forced}"
            assists.run_cmd(snippy)
        
    filename = f"{snippy_outdir}/{name}.vcf"
    outfile = f"{snippy_outdir}/{name}.ciprores.vcf"
    if os.stat(filename).st_size != 0:
       results = v2cr.reference_split(filename)
       if results is not None and results.empty is False:
           results.to_csv(outfile, sep='\t', index = False)
    
    pbp3_mutations, pbp3_designation, pbp3_group = pbp3.pbp3_checking(results)

    # checking for blaTEM 
    if is_assembly:
        assembly = args.fasta
    else:
        spades_outdir = maindir + "/spades"
        folder_exists = os.path.exists(spades_outdir)
        if not folder_exists:
            os.makedirs(spades_outdir)
            logging.info("Making spades output folder")
        else:
            logging.info(f"Spades folder exists")

        spades_result = assists.check_spades_finished(spades_outdir)
        if spades_result is False:
            spades = f"spades.py --careful --only-assembler --pe1-1 {args.R1} --pe1-2 {args.R2} -o {maindir}/spades"
            assists.run_cmd(spades)
        else:
            logging.info("Spades has already finished for this sample. Skipping.")
        assembly = spades_outdir + "/contigs.fasta"
    
    analysis_outdir = maindir + "/analysis"
    folder_exists = os.path.exists(analysis_outdir)
    if not folder_exists:
        os.makedirs(analysis_outdir)
        logging.info("Making analysis output folder")
    else:
        logging.info(f"Analysis folder exists")
    blatem_gene, blatem_designation = blatem.blatem_checking(assembly, analysis_outdir)
    cipro_designation, cipro_mutations = cipro.cipro_checking(results)

    outfile = os.path.join(maindir, f"{name}.hires.tsv")
    
    result_dict = {
        'PHP3 Mutations': pbp3_mutations,
        'PHP3 Group': pbp3_group,
        'PHP3 Designation': pbp3_designation,
        'BlaTEM Gene': blatem_gene,
        'BlaTEM Designation': blatem_designation,
        'Fluoroquinolone Mutations': cipro_mutations,
        'Fluoroquinolone Designation': cipro_designation
    }
    final_dict.update(result_dict)

    headers = list(final_dict.keys())
    values = [[str(item) if item is not None else 'None' for item in sublist] 
          if isinstance(sublist, list) else [str(sublist) if sublist is not None else 'None'] 
          for sublist in final_dict.values()]
    rows = ["\t".join(row) for row in zip(*values)] 

    tsv_lines = ["\t".join(headers)] + rows
    tsv_string = "\n".join(tsv_lines)

    with open(outfile, 'w') as output_file:
        output_file.write(tsv_string)
        logging.info(f"Writing information to {outfile}")
    logging.info(f"Complete!")

if __name__ == "__main__":
    parser = arguments.create_parser()  # pylint: disable=E1101
    args = parser.parse_args()
    main(args)
