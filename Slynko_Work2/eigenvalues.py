import numpy as np

dim = 100
L_plus_D = np.zeros((dim, dim))
for i in range(0, dim):
    L_plus_D[i][i] = 20.0
for i in range(0, 4):
    for j in range(0, i):
        L_plus_D[i][j] = 1.0
for i in range(4, dim):
    for j in range(i - 4, i):
        L_plus_D[i][j] = 1.0

U = np.zeros((dim, dim))
for i in range(0, dim - 4):
    for j in range(i + 1, i + 5):
        U[i][j] = 1.0
for i in range(dim - 4, dim):
    for j in range(i + 1, dim):
        U[i][j] = 1.0

try:
    tmp = np.linalg.inv(L_plus_D)
except np.linalg.LinAlgError:
    print("Something went wrong with matrix initialization")
    exit(1)
else:
    B = np.multiply(np.matmul(tmp, U), -1.0)
    w, v = np.linalg.eig(B)
    np.sort(w)
    print("The condition number is {}".format(np.linalg.norm(L_plus_D + U, np.inf)
                                              * np.linalg.norm(np.linalg.inv(L_plus_D + U), np.inf)))
    print("The minimum eigenvalue is {}".format(w[0]))
    print("The maximum eigenvalue is {}".format(w[len(w) - 1]))
