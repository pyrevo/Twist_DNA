def extract_chr(file):
    chr = None
    with open(file) as lines:
        chr = [line.split("\t")[0] for line in lines]
    return chr


def extract_stub(path, sep):
    parts = path.split(sep)
    if len(parts) <= 1:
        raise Exception("Unable to split path: " + path + " using separator " + sep)
    return parts[0]
