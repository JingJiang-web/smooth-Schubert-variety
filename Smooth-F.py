from chv1r6180 import *
import itertools
from fractions import Fraction
from typing import List, Sequence
from collections import deque
import traceback

RED = "\033[91m"
GREEN = "\033[92m"
RESET = "\033[0m"

coxeter_type = "F4"

if coxeter_type == 'F4':
    W = coxeter('F', 4)
    w_0 = [0, 1, 0, 2, 1, 0, 2, 1, 2, 3, 2, 1, 0, 2, 1, 2, 3, 2, 1, 0, 2, 1, 2, 3]

w_input = [1]
w_input = [x - 1 for x in w_input]
w1_0 = w_0 + w_input

try:
    try:
        C = cartanmat(coxeter_type)
    except TypeError:
        rank0 = {"F4": 4}[coxeter_type]
        C = cartanmat(coxeter_type, rank0)

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

    Cmat = [[int(x) for x in row] for row in C]
    n = len(Cmat)

    # ---- roots ----
    all_roots = [tuple(int(x) for x in r) for r in W.roots]
    all_root_set = set(all_roots)

    def positive_root(root):
        # lexicographic positive convention (same as original Smooth-F.py)
        for comp in root:
            if comp > 0:
                return True
            if comp < 0:
                return False
        return False

    positive_roots = [r for r in all_roots if positive_root(r)]
    print(f"Total roots: {len(all_roots)}, positive roots: {len(positive_roots)}")

    # ---- symmetrized bilinear form (for non-simply-laced) ----
    d = [None] * n
    d[0] = Fraction(1, 1)
    changed = True
    while changed:
        changed = False
        for i in range(n):
            if d[i] is None:
                continue
            for j in range(n):
                aij, aji = Cmat[i][j], Cmat[j][i]
                if i != j and aij != 0 and aji != 0:
                    if d[j] is None:
                        d[j] = d[i] * Fraction(aij, aji)
                        changed = True

    for i in range(n):
        if d[i] is None:
            d[i] = Fraction(1, 1)

    den_lcm = 1
    for x in d:
        den_lcm = den_lcm * x.denominator // __import__("math").gcd(den_lcm, x.denominator)
    dint = [int(x * den_lcm) for x in d]
    g = dint[0]
    for x in dint[1:]:
        g = __import__("math").gcd(g, x)
    dint = [x // g for x in dint]

    B = [[Fraction(dint[i] * Cmat[i][j], 1) for j in range(n)] for i in range(n)]

    def inner(u, v):
        s = Fraction(0, 1)
        for i in range(n):
            for j in range(n):
                s += Fraction(u[i], 1) * B[i][j] * Fraction(v[j], 1)
        return s

    def cartan_entry(alpha, beta):
        return (2 * inner(beta, alpha)) / inner(alpha, alpha)

    def reflect(alpha, beta):
        # s_alpha(beta) = beta - <alpha^vee, beta> alpha
        k = cartan_entry(alpha, beta)
        if k.denominator != 1:
            raise ValueError(f"Non-integral reflection coefficient {k} for alpha={alpha}, beta={beta}")
        k = int(k)
        return tuple(beta[i] - k * alpha[i] for i in range(n))

    def generate_subsystem_from_simple(S):
        X = set(S)
        X |= {tuple(-x for x in r) for r in list(X)}
        changed2 = True
        while changed2:
            changed2 = False
            current = list(X)
            for a in S:
                for b in current:
                    g2 = reflect(a, b)
                    if g2 in all_root_set and g2 not in X:
                        X.add(g2)
                        X.add(tuple(-x for x in g2))
                        changed2 = True
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

    def classify_subsystem(X):
        vecs = list(X)
        r = rank_Q(vecs)
        size = len(X)

        # F4 subsystems by rank/size + long/short counts
        if r == 2 and size == 8:
            return "B2"
        if r == 3 and size == 12:
            return "A3"
        if r == 3 and size == 18:
            lens = [inner(v, v) for v in vecs]
            mn, mx = min(lens), max(lens)
            num_short = sum(1 for L in lens if L == mn)
            num_long = sum(1 for L in lens if L == mx)
            if num_long == 12 and num_short == 6:
                return "B3"
            if num_long == 6 and num_short == 12:
                return "C3"
            return "B3/C3(?)"
        return None

    def find_subsystems_F4(return_roots=True, return_pos_roots=True):
        results = {t: {} for t in ["A3", "B2", "B3", "C3"]}

        # rank 2 (B2)
        for S in itertools.combinations(positive_roots, 2):
            if rank_Q(S) != 2:
                continue
            X = generate_subsystem_from_simple(list(S))
            t = classify_subsystem(X)
            if t in results and frozenset(X) not in results[t]:
                item = {"gens": tuple(S)}
                if return_roots:
                    item["roots"] = frozenset(X)
                if return_pos_roots:
                    item["pos_roots"] = tuple(sorted([r for r in X if positive_root(r)]))
                results[t][frozenset(X)] = item

        # rank 3 (A3, B3, C3)
        for S in itertools.combinations(positive_roots, 3):
            if rank_Q(S) != 3:
                continue
            X = generate_subsystem_from_simple(list(S))
            t = classify_subsystem(X)
            if t in results and frozenset(X) not in results[t]:
                item = {"gens": tuple(S)}
                if return_roots:
                    item["roots"] = frozenset(X)
                if return_pos_roots:
                    item["pos_roots"] = tuple(sorted([r for r in X if positive_root(r)]))
                results[t][frozenset(X)] = item

        return {t: list(results[t].values()) for t in results}

    # ---- pretty printing ----
    def root_to_expr(beta, sym="α"):
        terms = []
        for i, c in enumerate(beta, start=1):
            if c == 0:
                continue
            if c == 1:
                terms.append(f"{sym}{i}")
            elif c == -1:
                terms.append(f"-{sym}{i}")
            else:
                terms.append(f"{c}{sym}{i}")
        return " + ".join(terms).replace("+ -", "- ") if terms else "0"

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
        return " + ".join(terms).replace("+ -", "- ") if terms else "0"

    # ---- ambient simple-reflection action on coefficients ----
    def simple_root_1(j: int, n0: int) -> List[int]:
        v = [0] * n0
        v[j - 1] = 1
        return v

    def reflect_on_root(i: int, root: List[int], A: Sequence[Sequence[int]]) -> List[int]:
        n0 = len(root)
        ii = i - 1
        dot = 0
        for j in range(n0):
            dot += root[j] * A[j][ii]
        new_root = root[:]
        new_root[ii] -= dot
        return new_root

    def is_simple_root_coeffs(coeffs):
        """Return i (1-based) if coeffs is a simple root basis vector, else None."""
        ones = [k for k, c in enumerate(coeffs) if c == 1]
        if len(ones) == 1 and all(c == 0 for k, c in enumerate(coeffs) if k != ones[0]):
            return ones[0] + 1
        return None


    def reflection_word_from_beta_coeffs(t, coeffs):
        """
        coeffs: tuple/list in beta basis (Bourbaki order beta1,beta2,(beta3))
        Return a simple-reflection word for s_{beta}.
        """
        coeffs = tuple(int(x) for x in coeffs)

        # simple root: s_{beta_i} = s_i
        i = is_simple_root_coeffs(coeffs)
        if i is not None:
            return [i]

        # non-simple: use your explicit table
        try:
            return list(REFLECTION_WORDS[t][coeffs])
        except KeyError:
            raise KeyError(f"No reflection word for type {t}, beta-coeffs={coeffs}. "
                           f"Either add it to REFLECTION_WORDS or check beta numbering.")


    def sigma_from_ordered_I_sigma(t, I_sigma_coeffs_ordered, A_sub):
        """
        I_sigma_coeffs_ordered: e.g. [(1,0,0),(1,1,0),(1,1,1)] in C3, LEFT->RIGHT order.
        Returns reduced word for sigma = Π_{beta in I(sigma)} s_beta (left to right),
        reduced in the subgroup defined by A_sub.
        """
        concat = []
        for coeffs in I_sigma_coeffs_ordered:
            concat.extend(reflection_word_from_beta_coeffs(t, coeffs))

        # reduce in the subgroup
        return reduce_word_in_group(concat, A_sub)


    def I_sigma_to_beta_coeffs_ordered(I_sigma_roots_ordered, betas):
        I_list = []
        for r in I_sigma_roots_ordered:
            coeffs = solve_in_beta_basis(betas, r)
            I_list.append(tuple(coeffs))
        return I_list


    def gammas_from_word(w, A, *, use_reversed_prefix=True):
        inv_word = list(reversed(w))
        gammas = []
        prefix_reflections = []

        for idx in inv_word:
            r = simple_root_1(idx, len(A))
            it = reversed(prefix_reflections) if use_reversed_prefix else prefix_reflections
            for s in it:
                r = reflect_on_root(s, r, A)
            gammas.append(r)
            prefix_reflections.append(idx)

        return gammas

    # ---- solve coefficients of an ambient root in the beta-basis ----
    def solve_in_beta_basis(betas, r):
        betas = [tuple(int(x) for x in b) for b in betas]
        r = tuple(int(x) for x in r)
        k = len(betas)
        n_dim = len(r)

        M = [[Fraction(betas[j][i], 1) for j in range(k)] + [Fraction(r[i], 1)] for i in range(n_dim)]

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

    # ---- build subsystem Coxeter group purely from its Cartan matrix ----
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


    def beta_positive(coeffs):
        """Positivity test in beta-coordinates (lexicographic, consistent with positive_root)."""
        for c in coeffs:
            if c > 0:
                return True
            if c < 0:
                return False
        return False

    def sigma_from_inversion_set_beta(I_coeffs, A_sub):
        """Recover sigma from its inversion set I(sigma) given in beta-coordinates.

        Uses the standard 'peel simple roots' recursion:
          if beta_i ∈ I(w), then w = w' s_i and I(w') = s_i(I(w) \ {beta_i}).
        Returns a reduced word for w (generators numbered 1..k).
        """
        k = len(A_sub)

        def e(i):
            v = [0] * k
            v[i - 1] = 1
            return tuple(v)

        I = {tuple(int(x) for x in a) for a in I_coeffs}
        peeled = []  # rightmost letters, in peeling order

        while I:
            i_found = None
            for i in range(1, k + 1):
                if e(i) in I:
                    i_found = i
                    break
            if i_found is None:
                raise ValueError(f"Inversion set contains no simple beta root; cannot peel. I={sorted(I)}")

            peeled.append(i_found)
            I.remove(e(i_found))

            newI = set()
            for alpha in I:
                img = tuple(apply_simple_reflection(i_found, list(alpha), A_sub))
                # The recursion guarantees positivity, but normalize defensively.
                if not beta_positive(img):
                    img = tuple(-x for x in img)
                newI.add(img)
            I = newI

        # peeled letters are from right to left, reverse to get word
        return list(reversed(peeled))
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

        start = tuple(tuple(e(i)) for i in range(1, k + 1))
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


    def reduce_word_in_group(word, A_sub):
        k = len(A_sub)
        enum = enumerate_Wsub_words_from_Asub(A_sub)
        st = state_key_from_word(word, A_sub, k)
        return enum[st]

    STANDARD_CARTAN = {
        "A3": [
            [2, -1, 0],
            [-1, 2, -1],
            [0, -1, 2],
        ],
        "B2": [
            [2, -2],
            [-1, 2],
        ],
        "C2": [
            [2, -1],
            [-2, 2],
        ],
        "B3": [
            [2, -1, 0],
            [-1, 2, -2],
            [0, -1, 2],
        ],
        "C3": [
            [2, -1, 0],
            [-1, 2, -1],
            [0, -2, 2],
        ],
    }

    REFLECTION_WORDS = {

        "B2": {
            (1, 1): [1, 2, 1],

            (1, 2): [2, 1, 2],
        },

        "A3": {
            (1, 1, 0): [1, 2, 1],

            (0, 1, 1): [2, 3, 2],

            (1, 1, 1): [1, 2, 3, 2, 1],
        },

        "B3": {
            (1, 1, 0): [1, 2, 1],

            (0, 1, 1): [2, 3, 2],

            (1, 1, 1): [1, 2, 3, 2, 1],

            (0, 1, 2): [3, 2, 3],

            (1, 1, 2): [1, 3, 2, 3, 1],

            (1, 2, 2): [2, 1, 3, 2, 3, 1, 2],
        },

        "C3": {
            (1, 1, 0): [1, 2, 1],

            (0, 1, 1): [3, 2, 3],

            (1, 1, 1): [1, 3, 2, 3, 1],

            (1, 2, 1): [2, 1, 3, 2, 3, 1, 2],

            (0, 2, 1): [2, 3, 2],

            (2, 2, 1): [1, 2, 3, 2, 1],
        },
    }


    def permute_matrix(A, perm):
        return [[A[perm[i]][perm[j]] for j in range(len(perm))] for i in range(len(perm))]


    def canonicalize_betas(betas, t):
        """
        Reorder subsystem simple roots (betas) so that beta indices match the Dynkin diagram
        conventions in your picture:

          A3 :  β1—β2—β3

          Bn :  β1—...—β(n-1) => βn   (arrow points to SHORT root, so βn is short)
                In particular B2: β1 => β2  (β2 short)

          Cn :  β1—...—β(n-1) <= βn   (arrow points to SHORT root, so β(n-1) is short, βn is long)
                In particular C3: β1—β2<=β3  (β2 short, β3 long)

        Strategy:
          (1) Try permuting the given generators to match the target Cartan matrix exactly.
          (2) If not matched, fall back to chain-graph detection + oriented double edge from Cartan asymmetry.
        """
        betas = list(betas)
        k = len(betas)

        # build Cartan matrix in the given beta order
        A = cartan_matrix_from_betas(betas)

        # ---- Target Cartan matrices under YOUR desired numbering ----
        # Key point:
        #   - B2 with β2 short has Cartan [[2,-1],[-2,2]]  (this is often called "C2" in some conventions)
        #   - C3 with β2 short, β3 long has double edge oriented with -2 in row 2 col 3.
        TARGET_CARTAN = {
            "A3": [
                [2, -1, 0],
                [-1, 2, -1],
                [0, -1, 2],
            ],
            "B2": [  # enforce β2 short
                [2, -1],
                [-2, 2],
            ],
            "B3": [  # enforce β3 short
                [2, -1, 0],
                [-1, 2, -2],
                [0, -1, 2],
            ],
            "C3": [  # enforce β2 short, β3 long (arrow points to β2)
                [2, -1, 0],
                [-1, 2, -2],
                [0, -1, 2],
            ],
        }


        def permute_matrix(M, perm):
            return [[M[perm[i]][perm[j]] for j in range(len(perm))] for i in range(len(perm))]

        # ---------- (1) Try exact target match by permutation ----------
        T = TARGET_CARTAN.get(t)
        if T is not None and len(T) == k:
            for perm in itertools.permutations(range(k)):
                Aperm = permute_matrix(A, perm)
                if Aperm == T:
                    return [betas[i] for i in perm]

        # ---------- Helpers for fallback ----------
        def neighbors(i):
            return [j for j in range(k) if i != j and A[i][j] != 0]

        deg = [len(neighbors(i)) for i in range(k)]

        def is_chain_rank3():
            return k == 3 and sorted(deg) == [1, 1, 2]

        def chain_mid_and_ends():
            # assumes rank-3 chain
            mid = deg.index(2)
            ends = [i for i in range(3) if i != mid]
            return mid, ends

        def find_double_edge_pair():
            # returns (i,j) with i<j if unique pair has abs entry 2
            pairs = []
            for i in range(k):
                for j in range(i + 1, k):
                    if abs(A[i][j]) == 2 or abs(A[j][i]) == 2:
                        pairs.append((i, j))
            if len(pairs) == 1:
                return pairs[0]
            return None

        def stable_end_order(ends):
            # deterministic choice for A3 endpoints if needed
            return sorted(ends, key=lambda i: tuple(int(x) for x in betas[i]))

        # ---------- (2) Fallback by type ----------
        if t == "B2" and k == 2:
            # Force β2 short: need A[1][0] = -2 and A[0][1] = -1 (in 0-based indices)
            if A[1][0] == -2 and A[0][1] == -1:
                return [betas[0], betas[1]]
            if A[0][1] == -2 and A[1][0] == -1:
                return [betas[1], betas[0]]
            # if something degenerate happens, keep deterministic:
            return sorted(betas, key=lambda r: tuple(int(x) for x in r))

        if t == "A3" and is_chain_rank3():
            mid, ends = chain_mid_and_ends()
            e1, e2 = stable_end_order(ends)
            return [betas[e1], betas[mid], betas[e2]]

        if t in ("B3", "C3") and is_chain_rank3():
            mid, ends = chain_mid_and_ends()
            dbl = find_double_edge_pair()

            # identify which end is connected to mid by the double edge
            double_end = None
            if dbl is not None:
                i, j = dbl
                if mid in (i, j):
                    other = j if i == mid else i
                    # other should be one of the ends
                    if other in ends:
                        double_end = other

            # if can't detect, just do deterministic chain ordering
            if double_end is None:
                e1, e2 = stable_end_order(ends)
                return [betas[e1], betas[mid], betas[e2]]

            other_end = ends[0] if ends[1] == double_end else ends[1]

            # Candidate order: [other_end, mid, double_end]
            ordered = [betas[other_end], betas[mid], betas[double_end]]
            A2 = cartan_matrix_from_betas(ordered)

            if t == "B3":
                # enforce β3 short: A2[1][2] = -2, A2[2][1] = -1 (0-based: row 1 col 2 is -2)
                if A2[1][2] == -2 and A2[2][1] == -1:
                    return ordered
                # otherwise swap the two ends
                return [betas[double_end], betas[mid], betas[other_end]]

            if t == "C3":
                # enforce β2 short, β3 long: A2[1][2] = -2, A2[2][1] = -1
                # (same oriented double edge as we enforce for your C3 diagram)
                if A2[1][2] == -2 and A2[2][1] == -1:
                    return ordered
                return [betas[double_end], betas[mid], betas[other_end]]

        # ---------- final deterministic fallback ----------
        return sorted(betas, key=lambda r: tuple(int(x) for x in r))


    def sigma_from_ordered_I_sigma_product(t, I_sigma_beta_coeffs_ordered, A_sub):
        """
        I_sigma_beta_coeffs_ordered: 已按 LEFT->RIGHT 排序的 beta坐标正根列表
          例如 C3: [(1,0,0),(1,1,0),(1,1,1)]
        返回：sigma = Π_{γ in I(sigma) (left->right)} s_γ 的（在子群中化简后的）简约词
        """
        concat = []
        for coeffs in I_sigma_beta_coeffs_ordered:
            concat.extend(reflection_word_from_beta_coeffs(t, coeffs))

        # 在子Coxeter群里化简
        return reduce_word_in_group(concat, A_sub)


    # ---- main enumeration ----
    subs = find_subsystems_F4()
    print("\n==== Root subsystems in F4 ====")
    print("A3:", len(subs["A3"]))
    print("B2:", len(subs["B2"]))
    print("B3:", len(subs["B3"]))
    print("C3:", len(subs["C3"]))


    w = [x + 1 for x in w1_0]

    forbidden_patterns = {
        "B2": [
            [2, 1, 2],
        ],
        "A3": [
            [2, 1, 3, 2],
            [1, 2, 3, 2, 1],
        ],
        "B3": [
            [2, 1, 3, 2],
            [1, 2, 3, 2, 1],
            [1, 2, 3, 2, 1, 3],
            [1, 2, 3, 2, 1, 3, 2],
            [1, 2, 3, 2, 1, 2, 3],
            [1, 2, 3, 2, 1, 2, 3, 2],
        ],
        "C3": [
            [2, 1, 3, 2],
            [3, 2, 1, 3, 2],
            [2, 1, 3, 2, 3],
            [3, 2, 1, 3, 2, 3],
            [3, 2, 1, 2, 3],
            [1, 2, 3, 2, 1, 3, 2, 3],
        ],
    }



    def build_I_sigma_and_reflection_product(betas, intersection, t, *,
                                             print_I_sigma=True, print_sigma=True, check_forbidden=True):
        # intersection is a set of (ambient) positive roots in the subsystem
        # intersection is an ORDERED list of ambient roots (LEFT->RIGHT)
        I_sigma_roots = list(intersection)

        if t == "B2":
            Wsub = coxeter('F', 4)
        elif t == "A3":
            Wsub = coxeter('A', 3)
        elif t == "B3":
            Wsub = coxeter('B', 3)
        elif t == "C3":
            Wsub = coxeter('C', 3)

        if print_I_sigma:
            if I_sigma_roots:
                beta_forms = []
                for r in I_sigma_roots:
                    coeffs = solve_in_beta_basis(betas, r)
                    beta_forms.append(beta_expr_from_coeffs(coeffs, sym="β"))
                print(f"     {GREEN}I(sigma){RESET} = {{{', '.join(beta_forms)}}}")
            else:
                print(f"     {GREEN}I(sigma){RESET} = ∅")

        A_sub = cartan_matrix_from_betas(betas)

        I_list = I_sigma_to_beta_coeffs_ordered(I_sigma_roots, betas)

        sigma = sigma_from_ordered_I_sigma_product(t, I_list, A_sub)
        sigma = [x - 1 for x in sigma]
        sigma = Wsub.mattoword(Wsub.wordtomat(sigma))
        sigma = [x + 1 for x in sigma]




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
                print(f"     sigma is {RED}{sigma}{RESET}")
            else:
                print(f"     sigma is {sigma}")

        return I_sigma_roots, sigma, is_bad, hit_pat

    any_red = False

    for t in ["B2", "A3", "B3", "C3"]:
        if not subs[t]:
            continue
        print(f"\nType {t} :")
        print("=" * 60)

        for i, item in enumerate(subs[t][:2000], start=1):
            gens = item["gens"]
            pos_roots_sub = item["pos_roots"]

            betas = canonicalize_betas(gens, t)

            print(f"\n{i}: simple roots :")
            for j, b in enumerate(betas, start=1):
                print(f"     β{j} = {root_to_expr(b, sym='α')}")

            gs = gammas_from_word(w, Cmat, use_reversed_prefix=True)
            # print(gs)
            sub_pos_set = {tuple(r) for r in pos_roots_sub}

            intersection_ordered = [tuple(g) for g in gs if tuple(g) in sub_pos_set]

            _, _, is_bad, _ = build_I_sigma_and_reflection_product(betas, intersection_ordered, t)

            if is_bad:
                any_red = True

    print("=" * 90)
    print("X(w_0C) is non-smooth" if any_red else "X(w_0C) is smooth")

except Exception:
    print("=== Exception occurred ===")
    traceback.print_exc()