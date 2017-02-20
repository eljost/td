HARTREE2EV = 27.211386
HARTREE2NM = 45.56335
EV2NM = 1.239841e3
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


def conv(to_convert, fmt_str):
    return [CONV_DICT[t](item) for item, t in zip(to_convert, fmt_str)]
