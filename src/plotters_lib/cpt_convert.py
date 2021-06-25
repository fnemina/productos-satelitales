import numpy as np
import colorsys


def load_cpt(path):
    try:
        f = open(path)
    except Exception:
        print("File ", path, "not found")
        return None

    lines = f.readlines()

    f.close()

    x = np.array([])
    r = np.array([])
    g = np.array([])
    b = np.array([])

    color_model = 'RGB'

    for l in lines:
        ls = l.split()
        if l[0] == '#':
            if ls[-1] == 'HSV':
                color_model = 'HSV'
                continue
            else:
                continue
        if ls[0] == 'B' or ls[0] == 'F' or ls[0] == 'N':
            pass
        else:
            x = np.append(x, float(ls[0]))
            r = np.append(r, float(ls[1]))
            g = np.append(g, float(ls[2]))
            b = np.append(b, float(ls[3]))
            xtemp = float(ls[4])
            rtemp = float(ls[5])
            gtemp = float(ls[6])
            btemp = float(ls[7])

        x = np.append(x, xtemp)
        r = np.append(r, rtemp)
        g = np.append(g, gtemp)
        b = np.append(b, btemp)

    if color_model == 'HSV':
        for i in range(r.shape[0]):
            rr, gg, bb = colorsys.hsv_to_rgb(r[i] / 360., g[i], b[i])
        r[i] = rr
        g[i] = gg
        b[i] = bb

    if color_model == 'RGB':
        r = r / 255.0
        g = g / 255.0
        b = b / 255.0

    x_norm = (x - x[0]) / (x[-1] - x[0])

    red = []
    blue = []
    green = []

    for i in range(len(x)):
        red.append([x_norm[i], r[i], r[i]])
        green.append([x_norm[i], g[i], g[i]])
        blue.append([x_norm[i], b[i], b[i]])

    color_dict = {'red': red, 'green': green, 'blue': blue}

    return color_dict
