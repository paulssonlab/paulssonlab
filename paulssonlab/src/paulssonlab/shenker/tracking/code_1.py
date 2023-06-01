import numpy as np


def Make_d_nodes(array):
    node_arr = []
    for node in array:
        node = str(node)
        node1 = node
        node = [node + "a"]
        node1 = [node1 + "b"]
        node_arr += node + node1

    return node_arr


def divide_arr(array, m):
    m = Mmax

    l = len(array)
    arr_live = ["0"] * Mmax

    if l > Mmax:
        for i in range(0, m):
            arr_live[i] = array[i]

    else:
        arr_live = array
    return arr_live


M = input("Enter number of cells at t0 : ")
Mmax = input("Enter maximum possible number of cells : ")
T = input("Enter the last time point : ")
T = int(T)
M = int(M)
Mmax = int(Mmax)

out = []
for j in range(0, M):
    if j > (Mmax - 1):
        break
    out += [str(j)]

out1 = Make_d_nodes(out)
final_out = [[]] * T
final_out[0] = out

for i in range(1, T):
    final_out[i] = Make_d_nodes(final_out[i - 1])


for i in range(0, len(final_out)):
    final_out[i] = divide_arr(final_out[i], Mmax)


print(final_out)
