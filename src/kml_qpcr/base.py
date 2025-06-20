from pathlib import Path


class BaseQPCR():
    def __init__(self, sci_name: str, genome_set_dir: str, threads: int, force: bool):
        """
        qPCR 项目中公用的参数
        :sci_name: 目标微生物科学名
        :genome_set_dir: 基因组集合目录
        :threads: 线程数
        :force: 是否强制执行
        """
        self.sci_name = sci_name
        self.genome_set_dir = genome_set_dir
        self.threads = threads
        self.force = force
        # 当前微生物基因组目录
        self.gnm_dir = Path(genome_set_dir).joinpath(sci_name.replace(" ", "_"))
