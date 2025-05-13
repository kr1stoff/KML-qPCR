import re
import click
import logging
from pathlib import Path
from multiprocessing import Pool

logging.basicConfig(level=logging.DEBUG,
                    format="%(asctime)s - %(levelname)s - %(filename)s - %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")


@click.command()
@click.option("--workdir", default="qpcr_analysis", show_default=True, help="引物分析工作目录.")
@click.option("--threads", default=8, type=int, show_default=True, help="最大线程数.")
@click.help_option(help="显示帮助信息.")
def main(workdir, threads):
    """
    输入输出同目录，接上 qpcr_primer3_design.py 的工作目录，解析 primer3 输出引物表格
    """
    logging.info("开始解析 primer3 输出引物表格")
    primer3_outdir = Path(workdir).resolve().joinpath("primer3_out")
    primer3_pardir = Path(workdir).resolve().joinpath("primer3_parse")
    primer3_pardir.mkdir(parents=True, exist_ok=True)
    out_par_tuple = []
    for p3out in Path(primer3_outdir).glob("*.out"):
        p3par = primer3_pardir.joinpath(p3out.name.replace(".out", ".tsv"))
        out_par_tuple.append((p3out, p3par))
    with Pool(threads) as pool:
        pool.map(parse_primer3_out, out_par_tuple)
    logging.info("解析 primer3 输出表格完成")


def parse_primer3_out(input_tuple):
    """
    解析 primer3 输出表格

    :param primer3_out: primer3 输出表格
    :param parse_tab: 解析后的表格
    """
    primer3_out, parse_tab = input_tuple
    content = open(primer3_out, "r").read()
    chrom = re.findall(r"SEQUENCE_ID=(\S+)_sliding", content)[0]
    forward_pos = re.findall(r"PRIMER_LEFT_\d+=(\d+),(\d+)", content)
    forward_seq = re.findall(r"PRIMER_LEFT_\d+_SEQUENCE=(\S+)", content)
    forward_tm = re.findall(r"PRIMER_LEFT_\d+_TM=(\S+)", content)
    forward_gc = re.findall(r"PRIMER_LEFT_\d+_GC_PERCENT=(\S+)", content)
    reverse_pos = re.findall(r"PRIMER_RIGHT_\d+=(\d+),(\d+)", content)
    reverse_seq = re.findall(r"PRIMER_RIGHT_\d+_SEQUENCE=(\S+)", content)
    reverse_tm = re.findall(r"PRIMER_RIGHT_\d+_TM=(\S+)", content)
    reverse_gc = re.findall(r"PRIMER_RIGHT_\d+_GC_PERCENT=(\S+)", content)
    amplicon_len = re.findall(r"PRIMER_PAIR_\d+_PRODUCT_SIZE=(\S+)", content)
    amplicon_tm = re.findall(r"PRIMER_PAIR_\d+_PRODUCT_TM=(\S+)", content)
    probe_pos = re.findall(r"PRIMER_INTERNAL_\d+=(\d+),(\d+)", content)
    probe_seq = re.findall(r"PRIMER_INTERNAL_\d+_SEQUENCE=(\S+)", content)
    probe_tm = re.findall(r"PRIMER_INTERNAL_\d+_TM=(\S+)", content)
    probe_gc = re.findall(r"PRIMER_INTERNAL_\d+_GC_PERCENT=(\S+)", content)

    # ! 没有找到引物就不生成结果了
    if not forward_pos:
        return

    with open(parse_tab, "w") as f:
        header = "forward_sequence\treverse_sequence\tprobe_sequence\t"\
        "chromosome\tforward_length\tforward_tm\tforward_gc\treverse_length\treverse_tm\treverse_gc\t"\
        "probe_length\tprobe_tm\tprobe_gc\tamplicon_tm\tamplicon_length\n"
        f.write(header)
        for i in range(len(forward_pos)):
            # * 输出结果使用 samtools faidx 1-based 坐标
            outlist = [forward_seq[i], reverse_seq[i], probe_seq[i],
                       chrom, forward_pos[i][1], forward_tm[i], forward_gc[i], reverse_pos[i][1], reverse_tm[i],
                       reverse_gc[i], probe_pos[i][1], probe_tm[i], probe_gc[i], amplicon_tm[i], amplicon_len[i]]
            f.write("\t".join(outlist) + "\n")


if __name__ == "__main__":
    main()
