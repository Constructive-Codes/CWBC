
perm = [0, 1, 2, 21, 11, 13, 5, 8, 17, 10, 20, 7, 14, 6, 12, 22, 15, 19, 23, 9, 4, 16, 3, 18]

def read_code(path):
    with open(path) as f:
        return [list(map(int, line.split())) for line in f if line.strip()]

A = read_code("code960-from-1666.txt")
B = read_code("code960-24-8-10.txt")

A_perm = [[row[j] for j in perm] for row in A]

assert {tuple(r) for r in A_perm} == {tuple(r) for r in B}
print("Success")

