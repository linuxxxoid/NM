import matplotlib.pyplot as plt
import numpy as np

def fet(a, b, c, d, x):
    #print("dd")
    tmp = a + b * x + c * (x ** 2) + d * (x ** 3)
    #print(tmp)
    #print("TTTTTTTT")
    return tmp

def converter(f):
    return [float(f[i]) for i in range(0, len(f))]


def draw_plot(points, values, a, b, c, d):
    x, y = [], []
    n = len(points) - 1
    for i in range(n):
        print(i)
        x1 = np.linspace(points[i], points[i + 1], 10, endpoint = True)
        for j in x1:
            y1 = (fet(a[i], b[i], c[i], d[i], j - points[i]))
        #y1 = f(a[i], b[i], c[i], d[i], j - points[i]) for j in x1 
        #print("QQQQQQQ")
        #print(y1)

        x.append(x1)
        y.append(y1)

    print(x)
    print(y)
    plt.scatter(points, values, color='r')
    for i in range(n):
        plt.plot(x[i], y[i], color='b')

    plt.show()


if __name__ == '__main__':
    f = open('log.txt', 'r')
    points = f.readline().split()
    values = f.readline().split()

    a = f.readline().split()

    b = f.readline().split()
    c = f.readline().split()
    d = f.readline().split()

    points = converter(points)
    values = converter(values)
    a = converter(a)
    b = converter(b)
    c = converter(c)
    d = converter(d)
    
    draw_plot(points, values, a, b, c, d)
