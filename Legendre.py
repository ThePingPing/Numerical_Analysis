from scipy.special import roots_legendre
import numpy as np
from math import pi as PI


def simpson(a, b):
    h = (b-a)/2
    total = h / 3 * (np.sin(a) + 4 * np.sin(a + h) + np.sin(b))
    return total


def second_script_adaptive():
    error_tol = 0.00001
    all_simpsons = simpson(0, PI)
    points = [0, PI]
    index = 0

    while abs(all_simpsons - 2) >= error_tol and index <= 10:
        print(f'i is: {index}')
        all_simpsons = 0
        for i in range(len(points)-1):
            all_simpsons += simpson(points[i], points[i+1])

        new_points = []
        for i in range(len(points)-1):
            first, second = points[i], points[i+1]
            new_points.append(first)
            new_points.append((first+second)/2)
            new_points.append(second)
        points = new_points
        index += 1

    print(f'total all_simpsons is: {all_simpsons}')


def first_script_legendre():
    N = int(input("please enter N sample parameter: "))
    sample_roots = roots_legendre(N)[0]
    print(f'sample roots are: {sample_roots}')

    b = [0 if i % 2 == 1 else 2 / (i + 1) for i in range(2*N)]
    b = np.array(b)

    A = [[sample_roots[j] ** i for j in range(N)] for i in range(2*N)]
    A = np.array(A)

    C = np.linalg.lstsq(A, b)[0]
    print(f"the C-i elements are: {C}")

    estimate = 0
    for i in range(N):
        estimate += C[i] * np.sin(PI/2*(sample_roots[i]+1))

    estimate *= PI/2

    print(type(estimate))
    print(f'final estimate is legender: {estimate}')


if __name__ == '__main__':
    first_script_legendre()
    second_script_adaptive()
