def extract_chr(file):
    chr = None
    with open(file) as lines:
        chr = [line.split("\t")[0] for line in lines]
    return chr


def create_chr_entries_for_vff_header(file_fai, assembly):
    chr = []
    with open(file_fai) as lines:
        for line in lines:
            row = line.split("\t")
            chr.append([("ID", row[0]), ("length", row[1]), ("assembly", assembly)])
    return chr


def extract_stub(path, sep):
    parts = path.split(sep)
    if len(parts) <= 1:
        raise Exception("Unable to split path: " + path + " using separator " + sep)
    return parts[0]
