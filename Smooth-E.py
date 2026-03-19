from chv1r6180 import *
import itertools
from fractions import Fraction
from typing import List, Sequence
from collections import deque


RED = "\033[91m"
GREEN = '\033[92m'
RESET = "\033[0m"

coxeter_type = input("Coxeter type is (E6, E7, E8) : ").strip().upper()
n = {"E6":6, "E7":7, "E8":8}[coxeter_type]

if coxeter_type == 'E6':
    W = coxeter('E', 6)
if coxeter_type == 'E7':
    W = coxeter('E', 7)
if coxeter_type == 'E8':
    W = coxeter('E', 8)

if coxeter_type == "E6":
    w_0 = [0, 1, 2, 0, 3, 1, 2, 0, 3, 2, 4, 3, 1, 2, 0, 3, 2, 4, 3, 1, 5, 4, 3, 1, 2, 0, 3,
           2, 4, 3, 1, 5, 4, 3, 2, 0]
elif coxeter_type == "E7":
    w_0 = [0, 1, 2, 0, 3, 1, 2, 0, 3, 2, 4, 3, 1, 2, 0, 3, 2, 4, 3, 1, 5, 4, 3, 1, 2, 0, 3,
           2, 4, 3, 1, 5, 4, 3, 2, 0, 6, 5, 4, 3, 1, 2, 0, 3, 2, 4, 3, 1, 5, 4, 3, 2, 0, 6,
           5, 4, 3, 1, 2, 3, 4, 5, 6]
elif coxeter_type == "E8":
    w_0 = [0, 1, 2, 0, 3, 1, 2, 0, 3, 2, 4, 3, 1, 2, 0, 3, 2, 4, 3, 1, 5, 4, 3, 1, 2, 0, 3,
           2, 4, 3, 1, 5, 4, 3, 2, 0, 6, 5, 4, 3, 1, 2, 0, 3, 2, 4, 3, 1, 5, 4, 3, 2, 0, 6,
           5, 4, 3, 1, 2, 3, 4, 5, 6, 7, 6, 5, 4, 3, 1, 2, 0, 3, 2, 4, 3, 1, 5, 4, 3, 2, 0,
           6, 5, 4, 3, 1, 2, 3, 4, 5, 6, 7, 6, 5, 4, 3, 1, 2, 0, 3, 2, 4, 3, 1, 5, 4, 3, 2,
           0, 6, 5, 4, 3, 1, 2, 3, 4, 5, 6, 7]


w_input = [2,4,5,6] #Normal
w_input = [x - 1 for x in w_input]
w1_0 = w_0 + w_input

all_roots = [tuple(int(x) for x in r) for r in W.roots]
all_root_set = set(all_roots)

def is_positive(r):
    return all(x >= 0 for x in r) and any(x != 0 for x in r)

positive_roots = [r for r in all_roots if is_positive(r)]
print(f"Total roots: {len(all_roots)}, positive roots: {len(positive_roots)}")

C = cartanmat('E', n)
Cmat = [[int(x) for x in row] for row in C]
n = len(Cmat)
B = [[Fraction(Cmat[i][j], 1) for j in range(n)] for i in range(n)]


def print_matrix(mat):
    rows = [[str(x) for x in row] for row in mat]
    ncols = len(rows[0])
    colw = [max(len(row[j]) for row in rows) for j in range(ncols)]

    def fmt(row):
        return "  ".join(row[j].rjust(colw[j]) for j in range(ncols))

    if len(rows) == 1:
        print("⎡ " + fmt(rows[0]) + " ⎤")
        return

    print("⎡ " + fmt(rows[0]) + " ⎤")
    for row in rows[1:-1]:
        print("⎢ " + fmt(row) + " ⎥")
    print("⎣ " + fmt(rows[-1]) + " ⎦")


# print("Cartan matrix :")
# print_matrix(C)

def inner(u, v):
    s = Fraction(0, 1)
    for i in range(n):
        for j in range(n):
            s += Fraction(u[i], 1) * B[i][j] * Fraction(v[j], 1)
    return s

def cartan_entry(alpha, beta):
    return (2 * inner(beta, alpha)) / inner(alpha, alpha)

def reflect(alpha, beta):
    k = cartan_entry(alpha, beta)
    if k.denominator != 1:
        raise ValueError(f"Non-integral reflection coefficient {k}")
    k = int(k)
    return tuple(beta[i] - k * alpha[i] for i in range(n))

def generate_subsystem_from_simple(S):
    X = set(S)
    X |= {tuple(-x for x in r) for r in list(X)}
    changed = True
    while changed:
        changed = False
        current = list(X)
        for a in S:
            for b in current:
                g = reflect(a, b)
                if g in all_root_set and g not in X:
                    X.add(g)
                    X.add(tuple(-x for x in g))
                    changed = True
    return frozenset(X)

def rank_Q(vectors):
    if not vectors:
        return 0
    M = [[Fraction(x, 1) for x in v] for v in vectors]
    r = 0
    m = len(M)
    k = len(M[0])
    for col in range(k):
        pivot = None
        for i in range(r, m):
            if M[i][col] != 0:
                pivot = i
                break
        if pivot is None:
            continue
        M[r], M[pivot] = M[pivot], M[r]
        piv = M[r][col]
        M[r] = [x / piv for x in M[r]]
        for i in range(m):
            if i != r and M[i][col] != 0:
                f = M[i][col]
                M[i] = [M[i][j] - f * M[r][j] for j in range(k)]
        r += 1
        if r == m:
            break
    return r

def simple_root_1(j: int, n: int) -> List[int]:
    v = [0] * n
    v[j - 1] = 1
    return v


def reflect_on_root(i: int, root: List[int], A: Sequence[Sequence[int]]) -> List[int]:
    n = len(root)
    ii = i - 1
    dot = 0
    for j in range(n):
        dot += root[j] * A[j][ii]
    new_root = root[:]
    new_root[ii] -= dot
    return new_root


def gammas_from_word(w, A, *, use_reversed_prefix=True):
    inv_word = list(reversed(w))
    gammas = []
    prefix_reflections = []

    for idx in inv_word:
        r = simple_root_1(idx, n)

        if use_reversed_prefix:
            it = reversed(prefix_reflections)
        else:
            it = prefix_reflections

        for s in it:
            r = reflect_on_root(s, r, A)

        gammas.append(r)
        prefix_reflections.append(idx)

    return gammas











def is_edge(a, b):
    x = cartan_entry(a, b)
    y = cartan_entry(b, a)
    return x.denominator == 1 and y.denominator == 1 and int(x) == -1 and int(y) == -1


def dynkin_neighbors(S):
    S = list(S)
    N = {a: [] for a in S}
    for i, a in enumerate(S):
        for j in range(i + 1, len(S)):
            b = S[j]
            if is_edge(a, b):
                N[a].append(b)
                N[b].append(a)
    return N

def stable_root_key(r):

    r = tuple(int(x) for x in r)
    nnz = sum(1 for x in r if x != 0)
    height = sum(r)
    if nnz == 1 and 1 in r:
        i = r.index(1) + 1
    else:
        i = 10**9
    return (nnz, i, height, tuple(r))

def relabel_A3(S):

    S = [tuple(int(x) for x in r) for r in S]
    N = dynkin_neighbors(S)
    deg = {a: len(N[a]) for a in S}

    mid = [a for a in S if deg[a] == 2]
    ends = [a for a in S if deg[a] == 1]
    if len(mid) != 1 or len(ends) != 2:
        raise ValueError(f"Not an A3 chain by edges: degrees={deg}")

    beta2 = mid[0]
    e1, e2 = sorted(ends, key=stable_root_key)
    return (e1, beta2, e2)

def relabel_D4(S):
    S = [tuple(int(x) for x in r) for r in S]
    N = dynkin_neighbors(S)
    deg = {a: len(N[a]) for a in S}

    center = [a for a in S if deg[a] == 3]
    leaves = [a for a in S if deg[a] == 1]
    if len(center) != 1 or len(leaves) != 3:
        raise ValueError(f"Not a D4 star by edges: degrees={deg}")

    beta2 = center[0]
    l1, l2, l3 = sorted(leaves, key=stable_root_key)
    return (l1, beta2, l2, l3)

def relabel_by_type(S, t):
    if t == "A3":
        return relabel_A3(S)
    if t == "D4":
        return relabel_D4(S)
    return tuple(tuple(int(x) for x in r) for r in S)

def is_orthogonal(a, b):
    x = cartan_entry(a, b)
    y = cartan_entry(b, a)
    return x == 0 and y == 0

def canonical_generators(S):
    return tuple(sorted(S))

def classify_subsystem(X):
    r = rank_Q(list(X))
    size = len(X)
    if r == 3 and size == 12:
        return "A3"
    if r == 4 and size == 24:
        return "D4"
    return None


def find_subsystems_E(return_roots=True, return_pos_roots=True):
    results = {t: {} for t in ["A3", "D4"]}

    for S in itertools.combinations(positive_roots, 3):
        if rank_Q(S) != 3:
            continue
        X = generate_subsystem_from_simple(list(S))
        t = classify_subsystem(X)
        if t == "A3":
            if frozenset(X) not in results["A3"]:
                item = {"gens": tuple(S)}
                if return_roots:
                    item["roots"] = frozenset(X)
                if return_pos_roots:
                    item["pos_roots"] = tuple(sorted([r for r in X if is_positive(r)]))
                results["A3"][frozenset(X)] = item

    pos = list(positive_roots)

    def is_edge(a, b):
        x = cartan_entry(a, b)
        y = cartan_entry(b, a)
        return x.denominator == 1 and y.denominator == 1 and int(x) == -1 and int(y) == -1

    for c in pos:
        neigh = [r for r in pos if r != c and is_edge(c, r)]
        for r1, r2, r3 in itertools.combinations(neigh, 3):
            if inner(r1, r2) != 0 or inner(r1, r3) != 0 or inner(r2, r3) != 0:
                continue
            S = [r1, c, r2, r3]
            if rank_Q(S) != 4:
                continue
            X = generate_subsystem_from_simple(S)
            t = classify_subsystem(X)
            if t == "D4":
                if frozenset(X) not in results["D4"]:
                    item = {"gens": tuple(S)}
                    if return_roots:
                        item["roots"] = frozenset(X)
                    if return_pos_roots:
                        item["pos_roots"] = tuple(sorted([r for r in X if is_positive(r)]))
                    results["D4"][frozenset(X)] = item

    out = {t: list(results[t].values()) for t in results}
    return out

subs = find_subsystems_E()

def solve_in_beta_basis(betas, r):
    betas = [tuple(int(x) for x in b) for b in betas]
    r = tuple(int(x) for x in r)
    k = len(betas)
    n_dim = len(r)

    M = [[Fraction(betas[j][i], 1) for j in range(k)] + [Fraction(r[i], 1)]
         for i in range(n_dim)]

    row = 0
    pivcols = []
    for col in range(k):
        pivot = None
        for i in range(row, n_dim):
            if M[i][col] != 0:
                pivot = i
                break
        if pivot is None:
            continue
        M[row], M[pivot] = M[pivot], M[row]
        piv = M[row][col]
        M[row] = [x / piv for x in M[row]]

        for i in range(n_dim):
            if i != row and M[i][col] != 0:
                f = M[i][col]
                M[i] = [M[i][j] - f * M[row][j] for j in range(k + 1)]
        pivcols.append(col)
        row += 1
        if row == n_dim:
            break

    for i in range(n_dim):
        if all(M[i][j] == 0 for j in range(k)) and M[i][k] != 0:
            raise ValueError(f"No solution in span(betas) for r={r}")

    x = [Fraction(0, 1) for _ in range(k)]
    for i, col in enumerate(pivcols):
        x[col] = M[i][k]

    if any(xx.denominator != 1 for xx in x):
        raise ValueError(f"Non-integral beta-coeffs for r={r}: {x}")

    return [int(xx) for xx in x]


def beta_expr_from_coeffs(coeffs, sym="β"):
    terms = []
    for i, c in enumerate(coeffs, start=1):
        if c == 0:
            continue
        if c == 1:
            terms.append(f"{sym}{i}")
        elif c == -1:
            terms.append(f"-{sym}{i}")
        else:
            terms.append(f"{c}{sym}{i}")
    if not terms:
        return "0"
    s = " + ".join(terms)
    return s.replace("+ -", "- ")


def root_to_expr(beta, sym="α"):
    terms = []
    for i, c in enumerate(beta, start=1):
        if c == 0:
            continue
        if c == 1:
            terms.append(f"{sym}{i}")
        else:
            terms.append(f"{c}*{sym}{i}")
    return " + ".join(terms) if terms else "0"

def print_root_subsystem(X, title=None):

    if title:
        print(title)

    def sort_key(r):
        height = sum(r)
        nnz = sum(1 for x in r if x != 0)
        return (height, nnz, r)

    X_sorted = sorted(X, key=sort_key)

    for r in X_sorted:
        print("   ", root_to_expr(r, sym="α"))



def pos_root_sort_key(r):
    height = sum(r)
    nnz = sum(1 for x in r if x != 0)
    return (height, nnz, r)

print("A3 subsystems :", len(subs["A3"]))
print("D4 subsystems :", len(subs["D4"]))



w1 = [x + 1 for x in w1_0]

any_red = False

for t in ["A3", "D4"]:
    print(f"\nType {t} :")
    print("=" * 60)
    for i, item in enumerate(subs[t][:2000], start=1):
        gens = item["gens"]
        pos_roots = item["pos_roots"]
        betas = list(relabel_by_type(gens, t))
        print(f"\n{i}: simple roots :")
        for j, b in enumerate(betas, start=1):
            print(f"    β{j} = {root_to_expr(b, sym='α')}")

        print("-" * 40)
        # print(f"    positive roots : {len(pos_roots)}")
        # for r in sorted(pos_roots, key=pos_root_sort_key):
        #     coeffs = solve_in_beta_basis(betas, r)
        #     print(f"    {root_to_expr(r, sym='α')} = {beta_expr_from_coeffs(coeffs)}")


        forbidden_patterns = {
            "A3": [
                [2, 1, 3, 2],
                [1, 2, 3, 2, 1],
            ],
            "D4": [
                [2, 1, 3, 4, 2],
            ],
        }


        def root_len_sq(v):
            return inner(v, v)


        if w1_0 == w_input:
            gs = gammas_from_word(w1, Cmat, use_reversed_prefix=True)
        else:
            gs = gammas_from_word(w1, Cmat, use_reversed_prefix=False)

        inv_w = list(reversed(w1))
        gamma_set = {tuple(g) for g in gs}
        pos_roots = item["pos_roots"]
        sub_pos_set = {tuple(r) for r in pos_roots}
        intersection = sub_pos_set & gamma_set

        def cartan_matrix_from_betas(betas):
            k = len(betas)
            A = [[0] * k for _ in range(k)]
            for i in range(k):
                for j in range(k):
                    x = cartan_entry(betas[i], betas[j])
                    if x.denominator != 1:
                        raise ValueError("Non-integral Cartan entry from betas!")
                    A[i][j] = int(x)
            return A


        def apply_simple_reflection(i, root_vec, A_sub):
            return reflect_on_root(i, root_vec, A_sub)


        def state_key_from_word(word, A_sub, k):
            def e(i):
                v = [0] * k
                v[i - 1] = 1
                return v

            imgs = []
            for i in range(1, k + 1):
                v = e(i)
                for s in word:
                    v = apply_simple_reflection(s, v, A_sub)
                imgs.append(tuple(v))
            return tuple(imgs)


        _WSUB_ENUM_CACHE = {}


        def enumerate_Wsub_words_from_Asub(A_sub):
            key = tuple(tuple(row) for row in A_sub)
            if key in _WSUB_ENUM_CACHE:
                return _WSUB_ENUM_CACHE[key]

            k = len(A_sub)

            def e(i):
                v = [0] * k
                v[i - 1] = 1
                return v

            start = tuple(tuple(e(i)) for i in range(1, k + 1))  # identity images
            q = deque([start])
            word_of = {start: []}

            while q:
                imgs = q.popleft()
                wword = word_of[imgs]
                for s in range(1, k + 1):
                    new_imgs = []
                    for img in imgs:
                        new_imgs.append(tuple(apply_simple_reflection(s, list(img), A_sub)))
                    new_imgs = tuple(new_imgs)
                    if new_imgs not in word_of:
                        word_of[new_imgs] = wword + [s]
                        q.append(new_imgs)

            _WSUB_ENUM_CACHE[key] = word_of
            return word_of


        _REFLECT_WORD_CACHE = {}


        def reduce_word_in_group(word, A_sub):
            k = len(A_sub)
            enum = enumerate_Wsub_words_from_Asub(A_sub)
            st = state_key_from_word(word, A_sub, k)
            return enum[st]


        def reflection_word_from_root_in_subsystem(coeffs, A_sub):
            Akey = tuple(tuple(row) for row in A_sub)
            key = (Akey, tuple(coeffs))
            if key in _REFLECT_WORD_CACHE:
                return _REFLECT_WORD_CACHE[key]

            k = len(A_sub)
            enum = enumerate_Wsub_words_from_Asub(A_sub)

            target = list(coeffs)
            best = None  # (len(w_word), w_word, i)

            for imgs, w_word in enum.items():
                for i in range(1, k + 1):
                    if list(imgs[i - 1]) == target:
                        cand = (len(w_word), w_word, i)
                        if best is None or cand[0] < best[0]:
                            best = cand
                        break

            if best is None:
                raise ValueError(f"Cannot realize root {coeffs} as w(simple_i) in this subsystem basis")

            _, w_word, i = best
            refl = w_word + [i] + list(reversed(w_word))

            red = reduce_word_in_group(refl, A_sub)
            _REFLECT_WORD_CACHE[key] = red
            return red


        def build_I_sigma_and_reflection_product(
                betas, intersection, t, *,
                print_I_sigma=True, print_sigma=True, check_forbidden=True
        ):
            I_sigma_roots = list(intersection)

            def sort_key(r):
                return (sum(r), sum(1 for x in r if x != 0), r)

            I_sigma_roots.sort(key=sort_key)

            if print_I_sigma:
                if I_sigma_roots:
                    beta_forms = []
                    for r in I_sigma_roots:
                        coeffs = solve_in_beta_basis(betas, r)
                        beta_forms.append(beta_expr_from_coeffs(coeffs, sym="β"))
                    print(f"    {GREEN}I(sigma){RESET} = {{{', '.join(beta_forms)}}}")
                else:
                    print(f"    {GREEN}I(sigma){RESET} = ∅")

            A_sub = cartan_matrix_from_betas(betas)

            sigma_concat = []
            for r in I_sigma_roots:
                coeffs = solve_in_beta_basis(betas, r)
                rw = reflection_word_from_root_in_subsystem(coeffs, A_sub)
                sigma_concat.extend(rw)

            sigma = reduce_word_in_group(sigma_concat, A_sub)

            is_bad = False
            hit_pat = None

            if check_forbidden and t in forbidden_patterns:
                sigma_inv = list(reversed(sigma))
                for pat in forbidden_patterns[t]:
                    test = sigma_inv + pat
                    red = reduce_word_in_group(test, A_sub)
                    if red == []:
                        is_bad = True
                        hit_pat = pat
                        break

            if print_sigma:
                if is_bad:
                    print(f"    sigma is {RED}{sigma}{RESET}")
                else:
                    print(f"    sigma is {sigma}")

            return I_sigma_roots, sigma, is_bad


        _, sigma, is_bad = build_I_sigma_and_reflection_product(betas, intersection, t)
        if is_bad:
                any_red = True

print("=" * 90)
print("X(w) is non-smooth" if any_red else "X(w) is smooth")