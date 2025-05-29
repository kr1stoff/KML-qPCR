def list2txt(list, filename) -> None:
    """
    将列表写入文本文件，每行一个条目
    """
    with open(filename, "w") as f:
        for item in list:
            f.write(f"{item}\n")



def count_file_lines(file_path) -> int:
    """
    获取文件的行数
    :param file_path: 文件路径
    :return: 行数
    """
    with open(file_path, 'r') as file:
        return sum(1 for _ in file)
