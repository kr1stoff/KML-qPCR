import logging
from src.kml_qpcr.cli import mycli

# 日志规则
logging.basicConfig(level=logging.DEBUG,
                    format="%(asctime)s - %(levelname)s - %(filename)s - %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")

# 命令行参数和主流程在 mycli 中
if __name__ == "__main__":
    mycli(obj={})
