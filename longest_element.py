import re
from chv1r6180 import *

def parse_type(s: str):
    s = s.strip().upper()
    m = re.fullmatch(r'([A-G])\s*([0-9]+)', s)
    if not m:
        raise ValueError("Type format must look like A6, B4, C5, D7, E6, E7, E8, F4, G2")

    ctype = m.group(1)
    k = int(m.group(2))

    if ctype == "A":
        n_param = k + 1
        rank = k
        if rank < 1:
            raise ValueError("A1 or above required.")
    elif ctype in {"B", "C", "D"}:
        n_param = k
        rank = k
        if ctype == "D" and rank < 3:
            raise ValueError("D_n requires n >= 3.")
        if rank < 2 and ctype in {"B", "C"}:
            raise ValueError(f"{ctype}_n requires n >= 2.")
    elif ctype == "E":
        if k not in {6, 7, 8}:
            raise ValueError("E must be E6, E7, or E8.")
        n_param = k
        rank = k
    elif ctype == "F":
        if k != 4:
            raise ValueError("F must be F4.")
        n_param = 4
        rank = 4
    elif ctype == "G":
        if k != 2:
            raise ValueError("G must be G2.")
        n_param = 2
        rank = 2
    else:
        raise ValueError("Unsupported type.")
    return ctype, rank, n_param

def get_longest_w(ctype: str, rank: int):
    W = coxeter(ctype, rank)
    w0 = longestperm(W)
    w_0 = W.permtoword(w0)

    if ctype in {"B", "C", "D"}:
        w = [rank - x for x in w_0]
    else:
        w = [x + 1 for x in w_0]
    return w, w_0


def main():
    coxeter_type = input(
        "Coxeter type is (An, Bn, Cn, Dn, E6, E7, E8, F4, G2): "
    )
    ctype, rank, n = parse_type(coxeter_type)

    W = coxeter(ctype, rank)
    w0 = longestperm(W)
    w_0 = W.permtoword(w0)

    if ctype in {"B", "C", "D"}:
        w = [rank - x for x in w_0]
        print(f"Longest element is (Geck's notation form) : {w_0}")
        print(f"Longest element is (simple reflection form) : {w}")
    else:
        w = [x + 1 for x in w_0]
        if ctype in {"A", "E", "F", "G"}:
            print(f"Longest element is (Geck's notation form) : {w_0}")
        print(f"Longest element is (simple reflection form) : {w}")

if __name__ == "__main__":
    main()



