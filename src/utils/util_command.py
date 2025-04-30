import logging
from subprocess import run, CalledProcessError


def execute_command(cmd: str) -> str:
    """
    执行 shell 命令并返回标准输出。

    :param cmd: 要执行的命令
    :return: 命令的标准输出
    :raises RuntimeError: 如果命令执行失败
    """
    try:
        logging.debug(f"执行命令: {cmd}")
        res = run(cmd, shell=True, capture_output=True, text=True, check=True)
        return res.stdout.strip()
    except CalledProcessError as e:
        raise RuntimeError(f"命令执行失败: {cmd}\n错误信息: {e.stderr}")
