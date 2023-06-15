from sympy.abc import y, t
from sympy import log
from matplotlib import pyplot as plt


def rkf(f, a, b, alpha, tol, hmin, hmax, x_vals, y_vals):
    h = hmax
    ti = a
    w = alpha
    x_vals.append(ti)
    y_vals.append(w)

    while ti < b:
        k1 = h * f.subs([(t, ti), (y, w)])
        k2 = h * f.subs([(t, ti + h/4), (y, w + k1/4)])
        k3 = h * f.subs([(t, ti + 3*h/8), (y, w + (3*k1+9*k2)/32)])
        k4 = h * f.subs([(t, ti + 12*h/13), (y, w + (1932*k1 - 7200*k2 + 7296*k3)/2197)])
        k5 = h * f.subs([(t, ti + h), (y, w + 439*k1/216 - 8*k2 + 3680*k3/513 - 845*k4/4104)])
        k6 = h * f.subs([(t, ti + h/2), (y, w - 8*k1/27 + 2*k2 - 3544*k3/2565 + 1859*k4/4104 - 11*k5/40)])
        r = (1/h)*abs(k1/360 - 128*k3/4275 - 2197*k4/75240 + k5/50 + 2*k6/55)
        if r <= tol:
            ti += h
            w += 25*k1/216 + 1408*k3/2565 + 2197*k4/4104 - k5/5
            x_vals.append(ti)
            y_vals.append(w)
        delta = 0.84 * ((tol / r) ** (1 / 4))
        if delta <= 0.1:
            h = 0.1*h
        elif delta >= 4:
            h = 4.0*h
        else:
            h = delta*h
        if h > hmax:
            h = hmax
        if ti >= b:
            break
        elif ti + h > b:
            h = b - ti
        elif h < hmin:
            break

def plot_graph(y_sol, t_vals, w_vals):
    y_vals = [y_sol.subs(t, ti) for ti in t_vals]
    plt.plot(t_vals, w_vals)
    plt.plot(t_vals, y_vals)
    plt.legend(["Estimated", "Real"])
    plt.show()

def plot_error(t_vals, w_vals, h, a, y1):
    d_vals = [abs(y1.subs(t, t_vals[i]) - w_vals[i]) for i in range(len(t_vals))]
    plt.plot(t_vals, d_vals)
    plt.legend("Differences")
    plt.show()

f = -((y/t)**2)+(y/t)
# a = float(input("enter left end point, a:\n"))
# b = float(input("enter right end point, b:\n"))
# alpha = float(input("enter the initial condition, alpha:\n"))
# tol = float(eval(input("enter tolerance, tol:\n")))
# hmin = float(input("enter minimum mesh, hmin:\n"))
# hmax = float(input("enter maximum mesh, hmax:\n"))
x_vals = []
y_vals = []
rkf(f, 1, 4, 1, 10**(-6), 0.05, 0.5, x_vals, y_vals)
y_sol = t/(1+log(t))
plot_graph(y_sol, x_vals, y_vals)
plot_error(x_vals, y_vals, 0.05, 1, y_sol)
