import logging
import os
import sys
import click
from Bio import Entrez
from Bio import SeqIO
import urllib

from src.config.configure import EMAIL


logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')


# Some Functions
def parse_list(f_name: str):
    """
    Parse the target file list

    :param f_name: The input target name
    """
    res = set()
    with open(f_name, 'r') as IN:
        for line in IN:
            res.add(line.strip())
    return res


# Main
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version='1.0.0')
@click.option("-t", "--target",
              required=True,
              type=click.Path(),
              help="The file contain the target name list")
@click.option("-d", "--db",
              default="nucleotide",
              show_default=True,
              type=click.Choice(["nucleotide", "protein", "gene", "genome", "assembly"]),
              help="The database used to search")
@click.option("--format",
              default="fasta",
              show_default=True,
              type=click.Choice(["fasta", "gb"]),
              help="The format to download")
@click.option("-p", "--prefix",
              default="res",
              show_default=True,
              help="The out put prefix")
def main(target, db, format, prefix):
    """
    使用NCBI的entrez进行数据的下载
    """
    d_out = os.path.dirname(os.path.abspath(prefix))
    os.makedirs(d_out, exist_ok=True)
    f_out = f"{prefix}.fasta" if format == "fasta" else f"{prefix}.gbk"
    logging.info(f"Read the target list file {target}")
    ids = parse_list(target)
    logging.info(f"Start to fetch the seqinfo")
    Entrez.email = EMAIL
    try:
        handle = Entrez.efetch(db=db, id=ids, rettype="gb", retmode="text")
    except urllib.error.HTTPError:
        logging.error("Can't get info through Entrez, please check your input and database")
        sys.exit(1)
    f_format = "genbank" if format == "gb" else "fasta"
    record = SeqIO.parse(handle, "genbank")
    SeqIO.write(record, f_out, f_format)
    handle.close()


if __name__ == "__main__":
    main()
