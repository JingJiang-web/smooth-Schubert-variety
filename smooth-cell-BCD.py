from chv1r6180 import *
from package_simple_ref_to_oneline import calculate_oneline_notation
import itertools
from typing import Tuple, Dict, Set, Iterable, Optional
import colorama
colorama.init(autoreset=True)

RED = "\033[91m"
GREEN = '\033[92m'
RESET = "\033[0m"

print('  ###############################################################################################')
print('  ## Program name : smooth-cell-BCD                                                            ##')
print('  ##                                                                                           ##')
print(f"   ## Output there are two colors : {RED}red{RESET}, {GREEN}green{RESET}                                                 ##")
print('  ## Red one is "w = s_iw_0" we input                                                          ##')
print('  ## Green one is smooth element                                                               ##')
print('  ##                                                                                           ##')
print('  ## Remark : In type C_n, if we input w = s_1w_0, it is both red and green since smoothness   ##')
print('  ###############################################################################################')
print('                                                                                                 ')



DEFAULT_FORBIDDEN: Dict[str, Set[Tuple[int, ...]]] = {
    # 'B': {
    #     (-1,2,-3),(1,-2,-3),(1,2,-3),(1,-3,-2),(-2,-1,-3),(-2,1,-3),(2,-1,-3),
    #     (2,-3,-1),(-3,1,-2),(-3,-2,-1),(-3,-2,1),(-3,2,-1),(3,-2,-1),(3,-2,1),
    #     (-2,-4,3,1),(2,-4,3,1),(-3,-4,-1,-2),(-3,4,-1,2),(-3,4,1,2),(3,4,-1,2),(3,4,1,2),
    #     (4,-1,3,-2),(4,1,3,-2),(-4,2,3,1),(4,2,3,-1),(4,2,3,1),(-2,-1)
    # },
    # 'C': {
    #     (-1,2,-3),(1,-2,-3),(1,2,-3),(1,-3,-2),(-2,-1,-3),(-2,1,-3),(2,-1,-3),
    #     (2,-3,-1),(-3,1,-2),(-3,-2,-1),(-3,-2,1),(-3,2,-1),(3,-2,-1),(3,-2,1),
    #     (-2,-4,3,1),(2,-4,3,1),(-3,-4,-1,-2),(-3,4,-1,2),(-3,4,1,2),(3,4,-1,2),(3,4,1,2),
    #     (4,-1,3,-2),(4,1,3,-2),(-4,2,3,1),(4,2,3,-1),(4,2,3,1),(1,-2)
    # },

    'D': {
        (1,2,-3),(-1,2,-3),(-1,-3,-2),(1,-3,-2),(-2,-1,-3),(-3,-2,-1),(-1,4,-3,2),(-2,1,-3,-4),(2,-1,-3,-4),(2,1,-3,-4),
        (-2,-3,1,-4),(2,-3,1,-4),(2,-4,3,1),(-2,-4,3,-1),(-2,4,-3,-1),(2,4,-3,-1),(2,-4,3,-1),(-2,4,-3,1),(-2,-4,3,1),
        (3,-1,-2,-4),(3,1,-2,-4),(3,-2,1,-4),(3,-2,-4,1),(-3,-4,1,-2),(3,-4,-1,-2),(-3,4,1,2),(3,4,-1,2),(-3,4,-1,2),
        (3,4,1,2),(3,-4,1,-2),(-3,-4,-1,-2),(-3,-4,-2,1),(3,4,-2,-1),(-3,4,-2,1),(3,-4,-2,1),(-4,-1,-3,2),(4,1,3,-2),
        (-4,-1,3,-2),(4,-1,3,-2),(-4,1,3,-2),(4,-1,-3,2),(-4,1,-3,2),(4,-2,1,-3),(4,-2,-3,-1),(-4,-2,-3,-1),(-4,2,3,1),
        (4,2,3,-1),(-4,2,3,-1),(4,2,3,1),(4,-2,-3,1),(4,-3,-1,-2),(4,-3,-1,2),(-4,-3,1,2),(4,-3,1,-2),(4,-3,-2,1)
    },

    'B': {
        (3,4,1,2),(4,2,3,1),(-1,4,-3,2),(3,-2,1),(-3,-2,1),(-3,-2,-1),(3,-2,-1),(-3,2,-1),(2,-3,-1),(1,-3,-2),
        (-3,1,-2),(1,2,-3),(-1,2,-3),(2,-1,-3),(-2,1,-3),(1,-2,-3),(-2,-1,-3),(-2,-1)
    },

    'C': {
        (3,4,1,2),(4,2,3,1),(-1,4,-3,2),(3,-2,1),(-3,-2,1),(-3,-2,-1),(3,-2,-1),(-3,2,-1),(2,-3,-1),(1,-3,-2),
        (-3,1,-2),(1,2,-3),(-1,2,-3),(2,-1,-3),(-2,1,-3),(1,-2,-3),(-2,-1,-3),(1,-2)
    },

    # 'D': {
    #     (3,4,1,2),(4,2,3,1),(-1,4,-3,2),(3,-2,1),(-3,-2,1),(-3,-2,-1),(3,-2,-1),(-3,2,-1),(2,-3,-1),(1,-3,-2),
    #     (-3,1,-2),(1,2,-3),(-1,2,-3),(2,-1,-3),(-2,1,-3),(1,-2,-3),(-2,-1,-3)
    # },
}

def _group_by_length(patterns: Iterable[Iterable[int]]):
    grouped = {}
    for p in patterns:
        pt = tuple(p)
        grouped.setdefault(len(pt), set()).add(pt)
    return grouped
    #print(grouped)

def _normalize_by_abs_with_sign_tiebreak(subperm: Tuple[int, ...], indices: Tuple[int, ...]) -> Tuple[int, ...]:
    k = len(subperm)
    if k != len(indices):
        raise ValueError("The length of subperm should keep pace with indices")

    tmp = []
    for pos, (x, idx) in enumerate(zip(subperm, indices)):
        tmp.append((abs(x), idx, pos, 1 if x > 0 else -1))

    sorted_tmp = sorted(tmp, key=lambda t: (t[0], t[1]))

    rank_map = {}
    for rank, (absval, idx, pos, sign) in enumerate(sorted_tmp, start=1):
        rank_map[(absval, idx, pos)] = rank

    normalized = []
    for pos, (x, idx) in enumerate(zip(subperm, indices)):
        r = rank_map[(abs(x), idx, pos)]
        normalized.append(r if x > 0 else -r)
    return tuple(normalized)

def is_oneline_smooth(oneline: Tuple[int, ...], Coxeter_type: str,
                      forbidden_patterns: Optional[Dict[str, Set[Tuple[int, ...]]]] = None
                      ) -> Tuple[bool, Optional[dict]]:

    if Coxeter_type not in ('B','C','D'):
        raise ValueError("Coxeter_type must be 'B','C' or 'D'")

    if forbidden_patterns is None:
        forbidden = DEFAULT_FORBIDDEN
    else:
        forbidden = {k: set(v) for k, v in DEFAULT_FORBIDDEN.items()}
        for k, v in forbidden_patterns.items():
            forbidden[k] = set(tuple(x) for x in v)

    patterns = forbidden.get(Coxeter_type, set())
    if not patterns:
        return True, None

    patterns_by_len = _group_by_length(patterns)
    n = len(oneline)

    for k in (2,3,4):
        if k > n:
            continue
        patk = patterns_by_len.get(k)
        if not patk:
            continue

        for indices in itertools.combinations(range(n), k):
            sub = tuple(oneline[i] for i in indices)
            norm = _normalize_by_abs_with_sign_tiebreak(sub, indices)
            if norm in patk:
                match_info = {
                    'length': k,
                    'subsequence': sub,
                    'pattern': norm
                }
                return False, match_info
    return True, None

def replace_zeros(seq, n):
    return [n if x == 0 else x for x in seq]

def is_given_by_square(W, word1, word2):
    w1 = word1
    w2 = word2
    try:
        concat = w1[::-1] + w2
        square = W.mattoword(W.wordtomat(concat))
        if square == []:
            return True, f"{RED}{w2}{RESET}"
        else:
            return False, str(w2)
    except Exception as e:

        print(f"Warning: failed square check for words {w1} and {w2}: {e}")
        return False, str(w2)

coxeter_type = input("Coxeter type is (e.g.,B, C, D): ").strip().upper()
n = int(input("n : "))

if coxeter_type == 'A':
    W = coxeter(coxeter_type, n - 1)
elif coxeter_type in ['B', 'C', 'D']:
    W = coxeter(coxeter_type, n)
else:
    print(f"It is an error!")
    exit(1)


w_0 = [0, 1, 0, 1, 2, 1, 0, 1, 2, 3, 2, 1, 0, 1, 2, 3]
w_input = [1]
w_input = [n - x for x in w_input]
w1 = w_0 + w_input
elms = leftcellelm(W, w1)['elms']

for idx, elm in enumerate(elms, start=1):
    w = W.coxelmtoword(elm)
    w_processed = replace_zeros(w, n)
    #w_processed = w_processed[::-1]
    current_seq_3, _ = calculate_oneline_notation(coxeter_type, n, w_processed)
    tuple_seq = tuple(current_seq_3)
    is_smooth_flag, match = is_oneline_smooth(tuple_seq, coxeter_type)
    is_given, w_display = is_given_by_square(W, w1, w)


    if is_smooth_flag:
        if is_given:
            print(f"{RED}{GREEN}Smooth Weyl element #{idx} : {w_display}{RESET}")
        else:
            print(f"{GREEN}Smooth Weyl element #{idx} : {w_display}")
        print(f"\nOne-line notation : {tuple_seq}")

    else:
        if is_given:
            print(f"{RED}Non-smooth Weyl element #{idx} : {w_display}{RESET}")
        else:
            print(f"Non-smooth Weyl element #{idx} : {w_display}")
        print(f"\nOne-line notation : {tuple_seq}")
        print(f"\nMatched:", match)
    print("=============================================================================================")
