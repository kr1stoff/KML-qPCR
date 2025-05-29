from subprocess import run
from multiprocessing import Pool
from functools import partial
from subprocess import run, CalledProcessError


def execute_cmd_and_get_stdout(cmd: str) -> str:
    """
    执行 shell 命令并返回标准输出。
    :param cmd: 要执行的命令
    :return: 命令的标准输出
    :raises RuntimeError: 如果命令执行失败
    """
    try:
        res = run(cmd, shell=True, capture_output=True,
                  text=True, check=True, executable="/bin/bash")
        return res.stdout.strip()
    except CalledProcessError as e:
        raise RuntimeError(f"命令执行失败: {cmd}\n错误信息: {e.stderr}")


def multi_run_command(cmds: list[str], thread: int) -> None:
    """
    使用多线程执行一组 shell 命令.

    :param cmds: 要执行的命令列表.
    :param thread: 线程数.
    """
    _run = partial(run, shell=True, check=True)
    with Pool(processes=thread) as pool:
        pool.map(_run, cmds)
