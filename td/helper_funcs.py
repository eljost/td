import os

IRREPS_REPL = {
    "'": "_",
    "\"": "__",
}
FLOAT_RE = "([-\.\deEdD]+)"
CONV_DICT = {
        "s": str,
        "i": int,
        "f": float,
}
THIS_DIR = os.path.dirname(os.path.realpath(__file__))


def conv(to_convert, fmt_str):
    return [CONV_DICT[t](item) for item, t in zip(to_convert, fmt_str)]


def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i+n]
