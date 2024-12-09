import os
import subprocess
import sys
import logging
import shutil
import os.path
import pkg_resources
from Bio.Data.IUPACData import protein_letters_3to1

bor_vfdb_db = os.path.join(os.path.dirname(os.path.dirname(__file__)), "databases")

def run_cmd(command):
    """
    Run commands with error outputs.
    """
    logging.info("Running command: %s", command)
    result = subprocess.run(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True
    )
    result.stdout = result.stdout.decode()
    result.stderr = result.stderr.decode()
    if result.returncode != 0:
        logging.critical("Failed to run command: %s", result.args)
        logging.critical("stdout: %s", result.stdout)
        logging.critical("stderr: %s", result.stderr)
        sys.exit(1)
    return result


def check_files(file):
    """
    Check input files if they exist and have contents
    """

    if os.path.isfile(file) is True and os.stat(file).st_size != 0:
        truemsg = file + " exists and not empty, continuing..."
        logging.info(truemsg)
    else:
        msg = (
            file
            + " either file is does not exist or is empty, please check files. Exiting."
        )
        logging.critical(msg)
        sys.exit(1)


def check_folders(folder):
    """
    Check the output folder if it exists, if not make new directory.
    """
    if os.path.exists(folder) is True:
        truemsg = folder + " output folder exists"
        logging.info(truemsg)
    else:
        os.makedirs(folder)
        msg = folder + " does not exist, making output folder"
        logging.info(msg)


def check_dependencies(cmd_exec):
    cmd_path = shutil.which(cmd_exec)
    vcmd = subprocess.run([cmd_exec, "--version"], capture_output=True, text=True)
    if vcmd.returncode == 0 and vcmd.stdout != '':
        result = vcmd.stdout.splitlines()
    else:
        result = vcmd.stderr.splitlines()
    if cmd_exec == "abricate":
        version = " ".join(result).replace("abricate ", "v")
        if pkg_resources.parse_version(version) < pkg_resources.parse_version(
            "1.0.1"
        ):
            logging.critical("Abricate version too old, please upgrade to v1.0.0+")
            sys.exit(1)
    elif cmd_exec == "spades.py":
        version = result[0].replace("SPAdes genome assembler ", "")
    elif cmd_exec == 'minimap2':
        version = "v" + result[0]
    elif cmd_exec == 'snippy':
        version = result[0].replace("snippy ", "v")
    if cmd_exec == "samtools":
            version = result[0].replace("samtools ", "")
            if pkg_resources.parse_version(version) < pkg_resources.parse_version(
                "1.10"
            ):
                logging.critical("Samtools version too old, please upgrade to v1.10.0+")
                sys.exit(1)
    if cmd_exec == "bcftools":
        version = result[0].replace("bcftools ", "")
    if cmd_path is not None:
        msg = "Located " + cmd_exec + " " + version + " in " + cmd_path
        logging.info(msg)
    else:
        msg = cmd_exec + " was not found, please check installation on your device"
        logging.critical(msg)
        sys.exit(1)

        
def check_snippy_finished(snippy_outdir, name):
    result = "### rm -f "
    prokka_log = name + ".log"
    prokka_log = os.path.join(snippy_outdir, prokka_log)
    
    if os.path.isfile(prokka_log) and os.stat(prokka_log).st_size != 0:
        with open(prokka_log, 'r') as log:
            for line in log:
                if result in line:
                    return True
        return False
    else:
        return False
    
def convert(mutation):
    ref, pos, alt = mutation[:3], mutation[3:-3], mutation[-3:]
    ref_one = protein_letters_3to1.get(ref.capitalize(), ref)
    alt_one = protein_letters_3to1.get(alt.capitalize(), alt)
    return f"{ref_one}{pos}{alt_one}"

def check_spades_finished(spades_outdir):
    result = "======= SPAdes pipeline finished."
    result2 = "======= SPAdes pipeline finished WITH WARNINGS!"
    spades_log = os.path.join(spades_outdir, "spades.log")
    
    if os.path.isfile(spades_log) and os.stat(spades_log).st_size != 0:
        with open(spades_log, 'r') as log:
            for line in log:
                if result in line or result2 in line:
                    return True
        return False
    else:
        return False