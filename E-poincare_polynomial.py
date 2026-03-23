from chv1r6180 import *

coxeter_type = input("Coxeter type is (E6, E7, E8) : ").strip().upper()
n = {"E6": 6, "E7": 7, "E8": 8}[coxeter_type]

RED = "\033[91m"
GREEN = "\033[92m"
RESET = "\033[0m"

W = coxeter(coxeter_type, n)

if coxeter_type == "E6":
    w_0 = [0, 1, 2, 0, 3, 1, 2, 0, 3, 2, 4, 3, 1, 2, 0, 3, 2, 4, 3, 1, 5, 4, 3, 1, 2, 0, 3, 2, 4, 3, 1, 5, 4, 3, 2, 0]
elif coxeter_type == "E7":
    w_0 = [0, 1, 2, 0, 3, 1, 2, 0, 3, 2, 4, 3, 1, 2, 0, 3, 2, 4, 3, 1, 5, 4, 3, 1, 2, 0, 3, 2, 4, 3, 1, 5, 4, 3, 2, 0, 6, 5, 4, 3, 1, 2, 0, 3, 2, 4, 3, 1, 5, 4, 3, 2, 0, 6, 5, 4, 3, 1, 2, 3, 4, 5, 6]
elif coxeter_type == "E8":
    w_0 = [0, 1, 2, 0, 3, 1, 2, 0, 3, 2, 4, 3, 1, 2, 0, 3, 2, 4, 3, 1, 5, 4, 3, 1, 2, 0, 3, 2, 4, 3, 1, 5, 4, 3, 2, 0, 6, 5, 4, 3, 1, 2, 0, 3, 2, 4, 3, 1, 5, 4, 3, 2, 0, 6, 5, 4, 3, 1, 2, 3, 4, 5, 6, 7, 6, 5, 4, 3, 1, 2, 0, 3, 2, 4, 3, 1, 5, 4, 3, 2, 0, 6, 5, 4, 3, 1, 2, 3, 4, 5, 6, 7, 6, 5, 4, 3, 1, 2, 0, 3, 2, 4, 3, 1, 5, 4, 3, 2, 0, 6, 5, 4, 3, 1, 2, 3, 4, 5, 6, 7]

def word_length(word):
    return len(word)

def poincare_coeff_list(W, target):

    target_word = W.coxelmtoword(target)
    max_len = len(target_word)

    coeff = [0] * (max_len + 1)

    for y in redleftcosetreps(W):
        if bruhat(W, y, target):
            ell = len(W.coxelmtoword(y))
            coeff[ell] += 1

    return coeff

def is_smooth_from_coeff_list(coeff_list):
    return coeff_list == coeff_list[::-1]

w_input = [1]
w_input0 = [x - 1 for x in w_input]
target = W.wordtocoxelm(w_0 + w_input0)

print(f"w_input = {w_0 +w_input0}")
print("=" * 120)

coeff_list = poincare_coeff_list(W, target)
smooth = is_smooth_from_coeff_list(coeff_list)

line = (
    f"P_w coeffs = {coeff_list}\n"
    f"{'X(w) is smooth' if smooth else 'X(w) is non-smooth'}"
)

if smooth:
    print(GREEN + line + RESET)
else:
    print(RED + line + RESET)

print("-" * 120)