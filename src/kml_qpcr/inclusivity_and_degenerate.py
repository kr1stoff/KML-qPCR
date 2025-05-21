import click
from pathlib import Path
from functools import reduce
import pandas as pd
from subprocess import run, PIPE
from collections import Counter
from multiprocessing import Pool
from Bio.Seq import Seq

from src.config.cnfg_software import AMPLICON_EXTRACTOR, MAMBA


@click.command()
@click.option("--ref-seqs", required=True, help="输入参考序列集, 比如多个基因组的同源基因,CDS,转录本,保守区域等.")
@click.option("--workdir", default="qpcr_analysis", show_default=True, help="引物分析工作目录, 在该目录下读 primer3_parse 作为输入.")
@click.option("--threads", default=8, type=int, show_default=True, help="最大线程数.")
@click.help_option(help="显示帮助信息.")
def main(ref_seqs, workdir, threads):
    workdir = Path(workdir)
    # MAIN
    prsdir = workdir.joinpath("primer3_parse")
    shelldir = workdir.joinpath("shell")
    shelldir.mkdir(parents=True, exist_ok=True)
    incldir = workdir.joinpath("inclusivity_and_degenerate")
    incldir.mkdir(parents=True, exist_ok=True)
    # 合并所有引物探针，然后去重
    prs_files = list(prsdir.glob("*.tsv"))
    # ! 每个区域只要前 5 个, 如果敏感性不够再加
    dfs = [pd.read_csv(pf, sep='\t', usecols=['forward_sequence', 'reverse_sequence',
                                              'probe_sequence']).head(5) for pf in prs_files]
    rdc_df = reduce(lambda up, down: pd.concat([up, down]), dfs)
    rdc_df = rdc_df.drop_duplicates().reset_index(drop=True)

    # ! 耗时久调试过程中已完成就跳过, 并行批量运行扩增子/探针提取步骤
    # run_amplicon_command_batch(rdc_df, ref_seqs, incldir, shelldir, threads)

    # 批量计算包含性和简并引物探针
    pargs = [(item.forward_sequence, item.reverse_sequence, item.probe_sequence, ref_seqs, incldir)
             for item in rdc_df.itertuples()]
    with Pool(threads) as pool:
        pool.map(analysis_inclusivity, pargs)


def get_degenerate_primer(inseqs) -> tuple[str, int]:
    """
    根据输入序列集生成简并引物
    :param inseqs: 输入序列集
    :return: 简并引物序列和简并位点数量
    """
    # 简并碱基字典
    dgnrt_dict = {
        "A": "A", "C": "C", "G": "G", "T": "T",
        "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
        "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
    }
    # 反向简并字典
    rvs_dgnrt_dict = {v: k for k, v in dgnrt_dict.items()}
    # 生成简并引物, 并记录简并位点数量
    dgnrt_primer, dgnrt_num = "", 0
    # 遍历序列的每个位置
    for i in range(len(inseqs[0])):
        # 获取所有序列在当前位置的碱基
        ibs = [inseqs[j][i] for j in range(len(inseqs))]
        # 参考序列中存在简并碱基. 将简并碱基转换为普通碱基. 如果碱基在ATGC中则保持不变,否则根据简并碱基字典转换
        newibs = [nuc for nuc in (ib if ib in "ATGC" else dgnrt_dict[ib] for ib in ibs)]
        # 统计每个碱基的出现次数
        bscnt_dict = Counter(newibs)
        # 筛选出频率大于 0.2 的碱基
        fnlbss = [bs for bs in bscnt_dict if bscnt_dict[bs] / len(ibs) > 0.2]
        # 将筛选出的碱基排序并合并,然后查找对应的简并碱基
        fnlbs = rvs_dgnrt_dict["".join(sorted(fnlbss))]
        # 将简并碱基添加到结果中
        dgnrt_primer += fnlbs
        # 如果不是ATGC中的碱基,则简并位点数量加1
        if fnlbs not in "ATGC":
            dgnrt_num += 1
    return dgnrt_primer, dgnrt_num


def calc_inclusivity(sqs, curdir) -> float:
    """
    计算包含性
    :Param sqs: 参考序列文件
    :Param curdir: 当前输出目录
    :Return: 包含性
    """
    def count_sequence(file_path: str) -> str:
        return run(f"grep -c '>' {file_path}", shell=True, stdout=PIPE, encoding="utf-8").stdout.strip()
    # 包含引物探针的序列数 / 参考序列数
    nsqs = count_sequence(sqs)
    nincl = count_sequence(f"{curdir}/incl.probe.fa")
    # ! 如果没有抓取到扩增子和探针序列, 则返回 0
    if nsqs and nincl:
        return int(nincl) / int(nsqs)
    return 0


def analyze_degenerate_primer_probe(fwd, rvs, curdir) -> tuple[str, int, str, int, str, int]:
    """
    分析简并引物和探针
    :param fwd: 正向引物序列
    :param rvs: 反向引物序列
    :param curdir: 输出目录
    :return: 简并引物序列和简并位点数量
    """
    def read_sequences(file_path: str) -> list[str]:
        # ! seqkit amplicon 会输出很多空行, 是假阳性, 要删掉
        raw_sqs = run(f"grep -v '>' {curdir}/{file_path}", shell=True,
                      check=True, stdout=PIPE, encoding="utf-8").stdout.strip().split("\n")
        return [seq for seq in raw_sqs if seq]
    # 获取目标扩增子集序列
    amplc_sqs = read_sequences("dgnrt.amplicon.fa")
    # 截取上游和下游引物
    fwd_sqs = [seq[:len(fwd)] for seq in amplc_sqs]
    # * 下游引物需要反向互补
    rvs_sqs = [Seq(seq).reverse_complement()[:len(rvs)] for seq in amplc_sqs]
    prb_sqs = read_sequences("dgnrt.probe.fa")
    # 获取引物和探针的简并序列和简并位点数量
    fwd_dgnrt_seq, fwd_dgnrt_num = get_degenerate_primer(fwd_sqs)
    rvs_dgnrt_seq, rvs_dgnrt_num = get_degenerate_primer(rvs_sqs)
    prb_dgnrt_seq, prb_dgnrt_num = get_degenerate_primer(prb_sqs)
    return (fwd_dgnrt_seq, fwd_dgnrt_num,
            rvs_dgnrt_seq, rvs_dgnrt_num,
            prb_dgnrt_seq, prb_dgnrt_num)


def analysis_inclusivity(pargs) -> None:
    """
    包容性分析
    :param pargs: 并行分析参数元组
        - fwd: 正向引物序列
        - rvs: 反向引物序列
        - prb: 探针序列
        - sqs: 参考序列集
        - incldir: 包含性分析结果输出目录
    """
    # 解析参数
    fwd, rvs, prb, sqs, incldir = pargs
    # 创建当前引物探针组工作目录
    curdir = incldir.joinpath(f"{fwd}-{rvs}-{prb}")
    curdir.mkdir(parents=True, exist_ok=True)
    # 计算包容性
    inclusivity = calc_inclusivity(sqs, curdir)
    # ! 包容性小于 50% 就不往下做了
    if inclusivity < 0.5:
        return
    # 分析简并引物和探针
    rtntpl = analyze_degenerate_primer_probe(fwd, rvs, curdir)
    fwd_dgnrt_seq, fwd_dgnrt_num, rvs_dgnrt_seq, rvs_dgnrt_num, prb_dgnrt_seq, prb_dgnrt_num = rtntpl
    # 写入当前文件夹下的包含性结果
    with open(curdir.joinpath("inclusivity_and_degenerate.txt"), "w") as f:
        f.write("\t".join([fwd_dgnrt_seq, str(fwd_dgnrt_num),
                           rvs_dgnrt_seq, str(rvs_dgnrt_num),
                           prb_dgnrt_seq, str(prb_dgnrt_num),
                           str(inclusivity)]) + "\n")


def run_amplicon_command_batch(rdc_df, refsqs, incldir: Path, shelldir: Path, threads: int) -> None:
    """
    扩增子提取命令行拆分出来跑, 避免 "error: libmamba Could not set lock" 和 "反复运行"
    :param rdc_df: 去重后的上下游引物+探针序列3列数据框
    :param refsqs: 参考序列集
    :param incldir: 包含性分析结果输出目录
    :param shelldir: 拆分出来的 shell 脚本目录
    :param threads: 线程数
    :return: None
    """
    amplc_extct_shelldir = shelldir.joinpath("amplicon_extractor")
    amplc_extct_shelldir.mkdir(parents=True, exist_ok=True)
    total_shell = []
    # * 如果要做简并引物的话 amplicon-extractor 没有 seqkit amplicon 对齐的好; 如果只做包容性分析用 amplicon-extractor 更准确
    for item in rdc_df.itertuples():
        fwd, rvs, prb = item.forward_sequence, item.reverse_sequence, item.probe_sequence
        # 创建当前引物探针组工作目录
        curdir = incldir.joinpath(f"{fwd}-{rvs}-{prb}")
        cmd = f"""
        mkdir -p {curdir}
        {AMPLICON_EXTRACTOR} -F {refsqs} -m 6 -f {fwd} -r {rvs} -o {curdir}/incl.amplicon.fa
        {MAMBA} run -n basic seqkit amplicon --only-positive-strand --output-mismatches --line-width 0 \
            --max-mismatch 6 --forward {prb} {curdir}/incl.amplicon.fa > {curdir}/incl.probe.fa
        {MAMBA} run -n basic seqkit amplicon --only-positive-strand --output-mismatches --line-width 0 \
            --max-mismatch 3 --forward {fwd} --reverse {rvs} {refsqs} > {curdir}/dgnrt.amplicon.fa
        {MAMBA} run -n basic seqkit amplicon --only-positive-strand --output-mismatches --line-width 0 \
            --max-mismatch 6 --forward {prb} {curdir}/dgnrt.amplicon.fa > {curdir}/dgnrt.probe.fa
        """
        # 写入当前文件夹下的包含性结果
        with open(amplc_extct_shelldir.joinpath(f"{fwd}-{rvs}-{prb}.sh"), "w") as f:
            f.write(cmd)
        total_shell.append(f"bash {amplc_extct_shelldir}/{fwd}-{rvs}-{prb}.sh")
    # 写入总 shell 脚本, 然后运行
    with open(shelldir.joinpath("amplicon_extractor.sh"), "w") as f:
        f.write("\n".join(total_shell))
    run(f"cat {shelldir}/amplicon_extractor.sh | {MAMBA} run -n basic parallel -j {threads}",
        shell=True, check=True)


if __name__ == "__main__":
    main()
