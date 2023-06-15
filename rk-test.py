import math
from sympy.abc import y, t, z, x, c
from sympy import log
from matplotlib import pyplot as plt

# Sorry Dr" Altar I did not manage to solve the problem for c = 18, I would do better for the next code
# for c = 9 behavior and chaotic


def rk(x_tag, y_tag, z_tag, h, a, b, x0, y0, z0, x_vals, y_vals, z_vals):
    vals = []
    ax = plt.axes(projection='3d')
    m = 0
    for c0 in [4, 6, 8.5, 8.7, 9, 12]:
        x_vals = []
        y_vals = []
        z_vals = []
        z_tag_i = z_tag.subs([(c, c0)])
        w1 = x0
        w2 = y0
        w3 = z0
        t0 = a
        print("the c0 is:"+str(c0))
        while t0 <= b:
            print(t0)
            x_vals.append(w1)
            y_vals.append(w2)
            z_vals.append(w3)
            k11 = h * x_tag.subs([(t, t0), (x, w1), (y, w2), (z, w3)])
            k12 = h * y_tag.subs([(t, t0), (x, w1), (y, w2), (z, w3)])
            k13 = h * z_tag_i.subs([(t, t0), (x, w1), (y, w2), (z, w3)])
            k21 = h * x_tag.subs([(t, t0+h/2), (x, w1 + 0.5 * k11), (y, w2 + 0.5 * k12), (z, w3 + 0.5 * k13)])
            k22 = h * y_tag.subs([(t, t0+h/2), (x, w1 + 0.5 * k11), (y, w2 + 0.5 * k12), (z, w3 + 0.5 * k13)])
            k23 = h * z_tag_i.subs([(t, t0+h/2), (x, w1 + 0.5 * k11), (y, w2 + 0.5 * k12), (z, w3 + 0.5 * k13)])
            k31 = h * x_tag.subs([(t, t0+h/2), (x, w1 + 0.5 * k21), (y, w2 + 0.5 * k22), (z, w3 + 0.5 * k23)])
            k32 = h * y_tag.subs([(t, t0+h/2), (x, w1 + 0.5 * k21), (y, w2 + 0.5 * k22), (z, w3 + 0.5 * k23)])
            k33 = h * z_tag_i.subs([(t, t0+h/2), (x, w1 + 0.5 * k21), (y, w2 + 0.5 * k22), (z, w3 + 0.5 * k23)])
            k41 = h * x_tag.subs([(t, t0+h), (x, w1 + k31), (y, w2 + k32), (z, w3 + k33)])
            k42 = h * y_tag.subs([(t, t0+h), (x, w1 + k31), (y, w2 + k32), (z, w3 + k33)])
            k43 = h * z_tag_i.subs([(t, t0+h), (x, w1 + k31), (y, w2 + k32), (z, w3 + k33)])
            t0 += h
            w1 += ((1/6) * (k11 + 2 * k21 + 2 * k31 + k41))
            w2 += ((1/6) * (k12 + 2 * k22 + 2 * k32 + k42))
            w3 += ((1/6) * (k13 + 2 * k23 + 2 * k33 + k43))
        vals.append([x_vals[:], y_vals[:], z_vals[:]])
        if c0 == 9:
            ax.plot3D(vals[m][0], vals[m][1], vals[m][2], 'green')
        else:
            ax.plot3D(vals[m][0], vals[m][1], vals[m][2], 'gray')
        m += 1

    plt.show()


if __name__ == "__main__":
    x_tag = - y - z
    y_tag = x + 0.1 * y
    z_tag = 0.1 + z * (x - c)
    rk(x_tag, y_tag, z_tag, 0.1, 0, 200, 2, 5, 0.1, [], [], [])
