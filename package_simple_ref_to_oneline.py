RED = "\033[91m"
GREEN = '\033[92m'
RESET = "\033[0m"

print('                                                                                       ')
print('  #####################################################################################')
print('  ## Program name : package_simple_ref_to_oneline                                    ##')
print('  ##                                                                                 ##')
print('  ## Input one of following four types A(n-1), B(n), C(n), D(n)                      ##')
print('  ## For w, if we want to convert simple reflection form to one-line notation form   ##')
print('  ## We need input simple reflection form                                            ##')
print('  ##                                                                                 ##')
print(f"  ## {RED}For Example:{RESET}                                                                    ##")
print(f"  ## {RED}In type B(3), let s_i be the simple reflection{RESET}                                  ##")
print(f"  ## {RED}Then s_1s_2s_3s_1s_2s_1s_3 = [1,2,3,1,2,1,3]  (simple reflection form){RESET}          ##")
print(f"  ## {RED}We need input '1,2,3,1,2,1,3'{RESET}                                                   ##")
print('  #####################################################################################')
print('                                                                                       ')

def replace_zeros(seq, n):
    return [n if x == 0 else x for x in seq]

def apply_reflection(seq, reflection, coxeter_type, n):
    seq = seq.copy()
    if coxeter_type == 'A':
        idx = reflection - 1
        if 0 <= idx < len(seq) - 1:
            seq[idx], seq[idx + 1] = seq[idx + 1], seq[idx]
    elif coxeter_type in ['B', 'C']:
        if reflection < n:
            idx = reflection - 1
            if 0 <= idx < len(seq) - 1:
                seq[idx], seq[idx + 1] = seq[idx + 1], seq[idx]
        else:
            seq[0] = -seq[0]
    elif coxeter_type == 'D':
        if reflection < n:
            idx = reflection - 1
            if 0 <= idx < len(seq) - 1:
                seq[idx], seq[idx + 1] = seq[idx + 1], seq[idx]
        else:
            seq[0], seq[1] = -seq[1], -seq[0]
    return seq


def calculate_oneline_notation(coxeter_type, n, w):
    if coxeter_type == 'A':
        seq_length = n
    else:
        seq_length = n

    initial_seq = list(range(1, seq_length + 1))

    reversed_w = w[::-1]
    current_seq = initial_seq.copy()

    steps = []
    for r in reversed_w:
        current_seq = apply_reflection(current_seq, r, coxeter_type, n)
        steps.append(current_seq.copy())

    return current_seq, steps


def main():
    """交互式命令行入口"""
    valid_types = {'A', 'B', 'C', 'D'}
    coxeter_type = input("Coxeter type is (e.g., A, B, C, D): ").strip().upper()

    min_n = 2 if coxeter_type in ['A', 'B', 'C', 'D'] else 3
    while True:
        try:
            n = int(input(f"n ( type {coxeter_type} needs n ≥ {min_n}): ").strip())
            if n >= min_n:
                break
            print(f"错误：{coxeter_type}型要求n≥{min_n}")
        except ValueError:
            print("input valid n")

    if coxeter_type == 'A':
        seq_length = n
        max_reflection = seq_length - 1
    else:
        seq_length = n
        max_reflection = n

    while True:
        w_input = input(f"w is (e.g., 1,2,1,3,2 (simple reflection)): ").strip()
        try:
            w = [int(i) for i in w_input.split(',') if i.strip()]
            w = w[::-1]
            if coxeter_type == 'A':
                pass
            else:
                w = [n - x for x in w]
                w = replace_zeros(w, n)
            valid = all(1 <= r <= max_reflection for r in w)
            if valid:
                break
            print(f"number should between 1 and {max_reflection}")

        except ValueError:
            print("input valid sequence")

            valid = all(1 <= r <= max_reflection for r in w)
            if valid:
                break
            print(f"number should between 1 and {max_reflection}")
        except ValueError:
            print("input valid sequence")

    current_seq, _ = calculate_oneline_notation(coxeter_type, n, w)
    print(f"One-line notation corresponding to w is : {tuple(current_seq)}")


if __name__ == "__main__":
    main()
