def list2txt(list, filename) -> None:
    """
    将列表写入文本文件，每行一个条目
    """
    with open(filename, "w") as f:
        for item in list:
            f.write(f"{item}\n")
