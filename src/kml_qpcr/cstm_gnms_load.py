import logging
from pathlib import Path
from subprocess import run

from src.utils.util_command import multi_run_command


def load_customer_genomes(cstm_gnm_dir: str, sci_name: str, gnm_set_dir: str, threads: int) -> None:
    """
    导入客户基因组集, 格式化成符合项目结构的目录结构.

    :param cstm_gnm_dir: 客户基因组目录
    :param sci_name: 物种学名
    :param gnm_set_dir: 项目基因组集目录
    """
    logging.info(f"开始导入客户基因组集: {cstm_gnm_dir}, 物种学名: {sci_name}, 基因组集目录: {gnm_set_dir}")
    all_dir = Path(gnm_set_dir) / sci_name.replace(" ", "_") / "all"
    cmds = []
    for fna in Path(cstm_gnm_dir).glob("*.fna"):
        # 获取基因组名称
        gnm_name = fna.stem
        # 目标目录
        target_dir = all_dir / gnm_name
        target_dir.mkdir(parents=True, exist_ok=True)
        # 复制文件到目标目录
        copy_cmd = f"cp --force {fna} {target_dir}/{gnm_name}.fna"
        cmds.append(copy_cmd)
        # 执行所有复制命令
    if cmds:
        multi_run_command(cmds, threads=threads)
