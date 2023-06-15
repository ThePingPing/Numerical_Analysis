import numpy as np
import math


def s_i(row):
    new_row = [abs(item) for item in row[:-1]]
    return max(new_row)


#  detect when there is no solution
def detect_false_line(A):
    for index, row in enumerate(A):
        if A[index][-1] != 0 and row == [0] * len(row):
            return True
    return False


def delete_zero_lines(A):
    # return A after change and the non-zero-index
    amount_of_delete = 0
    index = 0
    while index < len(A):
        if A[index] == [0] * len(A[index]):
            print(f'delete one zero line now')
            del A[index]
            amount_of_delete += 1
        else:
            index += 1
    return A, amount_of_delete


def subtract_rows(A, i, j, multiply):
    # make the operation R(j) <- R(j) - multiply * R(i)
    A[j] = [item_j - multiply * item_i for item_i, item_j in zip(A[i], A[j])]
    return A


def swap_rows(A, i, j):
    print(f'i is: {i} and j is: {j}')
    A[i], A[j] = A[j][:], A[i][:]

    print(f"after swap A is: {A}")
    return A


def calc_scale(A):
    scale = []
    for i in range(len(A)):
        current_scale = max([abs(item) for item in A[i][:-1]])
        scale.append(current_scale)
    return scale


def one_elimination(A, k, scale):
    a_to_s_ratio = [0 for i in range(len(A)-k)]
    for index, row in enumerate(A[k:]):
        current_scale = abs(A[index + k][k]) / scale[index+k]
        a_to_s_ratio[index] = current_scale

    max_ratio = max(a_to_s_ratio)
    max_pivot = a_to_s_ratio.index(max_ratio)

    if max_pivot != 0:
        A = swap_rows(A, k, k+max_pivot)

    for index, row in enumerate(A[k+1:]):
        multiply = A[index+k+1][k] / A[k][k]
        A = subtract_rows(A, k, index+k+1, multiply)


def solve_equation(A):
    k = 0
    scale = calc_scale(A)
    print(f'scale is: {scale}')

    while k < len(A):
        if detect_false_line(A):
            print("error there is no solution")
            return False

        delete_zero_lines(A)

        if len(A) < len(A[0]) - 1:
            print(A)
            print("there is no single solution")
            return "infinite"

        one_elimination(A, k, scale)
        A = np.array(A)
        A = np.around(A, 10)
        A = A.tolist()
        print(f"with k at: {k} the matrix is: {A}")
        k += 1

    if detect_false_line(A):
        print("error there is no solution")
        return False

    delete_zero_lines(A)

    if len(A) < len(A[0]) - 1:
        print("there is infinite amount of solutions")
        return "infinite"

    print(f"before back_substitution the matrix is: {A}")
    return back_substitution(A)


def back_substitution(A):
    A = np.array(A)

    n = A[:][-1].size - 1
    x = np.zeros_like(A[0][:-1])

    if A[n-1, n-1] == 0:
        print("error there is no solution")
        return False

    for i in range(n-1, -1, -1):
        x[i] = A[i][-1]
        for j in range(i+1, n):
            x[i] -= A[i][j] * x[j]
        x[i] /= A[i][i]

    return x


def get_column(A, i):
    col = []
    for row in A:
        col.append(row[i])
    return col


def remove_last_column(A):
    for i in range(len(A)):
        del A[i][-1]
    return A


def append_column(A, column):
    for i in range(len(A)):
        A[i].append(column[i])
    return A


def iterative_refinement(A, e):
    n = len(A)
    estimate = [0 for i in range(n)]
    error = math.inf
    initial_b = np.array(get_column(A, len(A[0])-1))
    current_b = get_column(A, len(A[0])-1)
    A = remove_last_column(A)
    i = 0

    while error > e or i > 10:
        A = append_column(A, current_b)

        x_roof = solve_equation(A)
        x_roof.astype('float')

        A = remove_last_column(A)

        A, x_roof = np.array(A), np.array(x_roof)
        b_roof = A.dot(x_roof)

        remainder = np.subtract(b_roof, current_b)
        current_b = remainder

        estimate = np.array(estimate)
        estimate = np.add(estimate, x_roof)

        estimate_b = A.dot(estimate)
        error = np.subtract(initial_b, estimate_b)
        error = np.linalg.norm(error, 2)
        A = A.tolist()
        i += 1

    return estimate


def main():
    A = [[2.11, -4.21, 0.921, 2.01], [4.01, 10.2, -1.12, -3.09], [1.09, 0.987, 0.832, 4.21]]

    x = solve_equation(A)
    print(f"solution x is: {x}")

    iterative_refinement(A, 2)


if __name__ == '__main__':
    main()