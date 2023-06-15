from sympy.abc import t, y
from sympy import Function
from matplotlib import pyplot as plt
import math


def euler(f, y0, h, a, b, x_vals, y_vals): # f(t,y(t)), y0 of edge point - y(a), h - partition parameter, [a,b] - segment
    N = int((b - a) / h)  # num iterations
    w = y0
    for i in range(N + 1):
        ti = a + i * h
        x_vals.append(ti)
        y_vals.append(w)
        w = w + h * f.subs([(t, ti), (y, w)])
    return w


def check_sol(y1, f):
    f1 = y1.diff(t)
    f2 = f.subs([(y,y1)])
    return f1.compare(f2) == 0


def plot_graph(y1, t_vals, w_vals):
    y_vals = [y1.subs(t,ti) for ti in t_vals]
    plt.plot(t_vals, w_vals)
    plt.plot(t_vals, y_vals)
    plt.legend(["Estimated", "Real"])
    plt.show()


def calculate_error(t_vals, w_vals, h, a, y1):
    real_d_vals = [((h * 2) / 2 * 1) * (math.e ** ((ti - a) * 1) - 1) for ti in t_vals]
    d_vals = [abs(y1.subs(t, t_vals[i]) - w_vals[i]) for i in range(len(t_vals))]
    plt.plot(t_vals, d_vals)
    plt.plot(t_vals, real_d_vals)
    plt.legend(["Differences", "Differences Bound"])
    plt.show()


def linear_interpolation(x0, w0, x1, w1, val):
    l = ((w1 - w0)/(x1 - x0)) * (t - x0) + w1
    return l.subs([(t, val)])


def hermit_interpolation(x0, w0, x1, w1, val, f):
    v1 = f.subs([(t, x0), (y, w0)])
    m = (w1 - w0) / (x1 - x0)
    v2 = f.subs([(t, x1), (y, w1)])

    v1t = (m - v1) / (x1 - x0)
    v2t = (v2 - m) / (x1 - x0)

    v1tt = (v2t - v1t) / (x1 - x0)

    yt = w0 + v1 * (t - x0) + v1t * ((t - x0) ** 2) + v1tt * ((t - x0) ** 2) * (t - x1)
    return yt.subs([(t, val)])


def calculate_h():
    return 5 / (100 * (math.e - 1))


def taylor(N, f, yf, y0, h, a, b, tx_vals, ty_vals): # f(t,y(t)), y0 of edge point - y(a), h - partition parameter, [a,b] - segment
    n = int((b - a) / h)  # num iterations
    w = y0
    for i in range(n + 1):
        ti = a + i * h
        tx_vals.append(ti)
        ty_vals.append(w)
        T = sum([n_der(f, yf, j).subs([(t, ti), (y, w)]) * (h ** j) / (math.factorial(j+1)) for j in range(N)])
        w = w + h * T
    return w

def n_der(f, yf, n):
    dif = f
    for i in range(n):
        dif = dif.diff(t)
        dif = dif.subs([(yf.diff(t), f)])
    return dif.subs([(yf, y)])

def rk(f, y0, h, a, b, x_vals, y_vals):
    N = int((b - a) / h)  # num iterations
    w = y0
    for i in range(N + 1):
        ti = a + i * h
        x_vals.append(ti)
        y_vals.append(w)
        k1 = h * f.subs([(t, ti), (y, w)])
        k2 = h * f.subs([(t, ti + h / 2), (y, w + k1 / 2)])
        k3 = h * f.subs([(t, ti + h / 2), (y, w + k2 / 2)])
        k4 = h * f.subs([(t, ti + h), (y, w + k3)])
        w = w + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return w



def main():
    t_vals1 = []
    w_vals1 = []
    t_vals2 = []
    w_vals2 = []
    print("Q1")
    f = (1/t) ** 2 - (y/t) - y ** 2
    yf = Function('yf')(t)
    ft = (1/t) ** 2 - (yf/t) - yf ** 2
    print("euler with h = 0.05:")
    print(euler(f, -1, 0.05, 1, 2, t_vals1, w_vals1))
    print("euler with h = 0.5:")
    print(euler(f, -1, 0.5, 1, 2, t_vals2, w_vals2))
    print("Q2")
    y1 = - 1 / t
    print(check_sol(y1, f))
    print("Q3")
    print("plot euler with h = 0.05")
    plot_graph(y1, t_vals1, w_vals1)
    print("plot euler with h = 0.5")
    plot_graph(y1, t_vals2, w_vals2)
    print("Q4")
    print("plot differences with h = 0.05")
    calculate_error(t_vals1, w_vals1, 0.05, 1, y1)
    print("plot differences with h = 0.5")
    calculate_error(t_vals2, w_vals2, 0.5, 1, y1)
    print("Q5")
    print("linear interpolation:")
    print(linear_interpolation(t_vals1[11], w_vals1[11], t_vals1[12], w_vals1[12], 1.555))
    print("Q6")
    print("hermit interpolation:")
    print(hermit_interpolation(t_vals1[11], w_vals1[11], t_vals1[12], w_vals1[12], 1.555, f))
    print("Q7")
    print("h value:")
    print(calculate_h())
    print("Q8")
    taylor_x1_vals = []
    taylor_y1_vals = []
    print("taylor with h = 0.05:")
    print(taylor(2, ft, yf, -1, 0.05, 1, 2, taylor_x1_vals, taylor_y1_vals))
    taylor_x2_vals = []
    taylor_y2_vals = []
    print("taylor with h = 0.5:")
    print(taylor(2, ft, yf, -1, 0.5, 1, 2, taylor_x2_vals, taylor_y2_vals))
    print("plot taylor with h = 0.05")
    plot_graph(y1, taylor_x1_vals, taylor_y1_vals)
    print("plot taylor with h = 0.5")
    plot_graph(y1, taylor_x2_vals, taylor_y2_vals)
    rk_x1_vals = []
    rk_y1_vals = []
    rk_x2_vals = []
    rk_y2_vals = []
    print("rk with h = 0.05:")
    print(rk(f, -1, 0.05, 1, 2, rk_x1_vals, rk_y1_vals))
    print("rk with h = 0.5:")
    print(rk(f, -1, 0.5, 1, 2, rk_x2_vals, rk_y2_vals))
    print("plot rk with h = 0.05")
    plot_graph(y1, rk_x1_vals, rk_y1_vals)
    print("plot rk with h = 0.5")
    plot_graph(y1, rk_x2_vals, rk_y2_vals)





if __name__ == "__main__":
    main()

