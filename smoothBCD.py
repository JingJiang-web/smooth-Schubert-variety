import itertools
from typing import Iterable, Tuple, Dict, Set, Optional

DEFAULT_FORBIDDEN: Dict[str, Set[Tuple[int, ...]]] = {
    'B': {
        (3,4,1,2),(4,2,3,1),(-1,4,-3,2),
        (3,-2,1),(3,-2,-1),(-3,-2,1),(-3,-2,-1),
        (-3,2,-1),(-2,-1,-3),(1,2,-3),(-1,2,-3),
        (2,-1,-3),(-2,1,-3),(1,-2,-3),(-2,-1)
    },
    'C': {
        (3,4,1,2),(4,2,3,1),(-1,4,-3,2),
        (3,-2,1),(3,-2,-1),(-3,-2,1),(-3,-2,-1),
        (-3,2,-1),(-2,-1,-3),(1,2,-3),(-1,2,-3),
        (2,-1,-3),(-2,1,-3),(1,-2,-3),(-2,-1)
    },
    'D': {
        (3,4,1,2),(4,2,3,1),(-1,4,-3,2)
    },
}


def _group_by_length(patterns: Iterable[Iterable[int]]) -> Dict[int, Set[Tuple[int, ...]]]:
    grouped: Dict[int, Set[Tuple[int, ...]]] = {}
    for p in patterns:
        pt = tuple(p)
        grouped.setdefault(len(pt), set()).add(pt)
    return grouped

def _normalize_by_abs_with_sign_tiebreak(subperm: Tuple[int, ...], indices: Tuple[int, ...]) -> Tuple[int, ...]:
    k = len(subperm)
    if k != len(indices):
        raise ValueError("The length of subperm and indices should keep the same")

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

def is_smooth(
    perm: Iterable[int],
    Coxeter_type: str,
    forbidden_patterns: Optional[Dict[str, Set[Tuple[int, ...]]]] = None,
    verbose: bool = False
) -> Tuple[bool, Optional[Dict[str, object]]]:

    if Coxeter_type not in ('B','C','D'):
        raise ValueError("Coxeter_type must be one of 'B','C','D'")

    perm_list = list(perm)
    n = len(perm_list)
    if n < 2:
        return True, None

    if forbidden_patterns is None:
        forbidden = {k: set(v) for k, v in DEFAULT_FORBIDDEN.items()}
    else:
        forbidden = {k: set(v) for k, v in DEFAULT_FORBIDDEN.items()}
        if not isinstance(forbidden_patterns, dict):
            raise ValueError("forbidden_patterns must be a dict mapping type->set_of_patterns")
        for k, v in forbidden_patterns.items():
            forbidden[k] = set(tuple(x) for x in v)

    forbidden_set = forbidden.get(Coxeter_type, set())
    if not forbidden_set:

        return True, None

    forbidden_by_len = _group_by_length(forbidden_set)

    for k in (2,3,4):
        if k > n-1:
            continue
        patterns_k = forbidden_by_len.get(k)
        if not patterns_k:
            continue
        for indices in itertools.combinations(range(n), k):
            sub = tuple(perm_list[i] for i in indices)
            norm = _normalize_by_abs_with_sign_tiebreak(sub, indices)
            if norm in patterns_k:
                match_info = {
                    'Coxeter_type': Coxeter_type,
                    'indices': indices,
                    'original_subsequence': sub,
                    'normalized': norm,
                    'matched_pattern': norm
                }
                if verbose:
                    print("This element is non-smooth since the subsequence", sub, "corresponding to", norm)
                return False, match_info

    return True, None

if __name__ == "__main__":

    perm = (2, -3, -1)

    smooth_flag, info = is_smooth(perm, 'B', verbose=True)
    print("This element is smooth :", smooth_flag)


