def poly2latex(poly, variable="x", width=2):
    t = ["{0:0.{width}f}"]
    t.append(t[-1] + " {variable}")
    t.append(t[-1] + "^{1}")

    def f():
        for i, v in enumerate(reversed(poly)):
            idx = i if i < 2 else 2
            yield t[idx].format(v, i, variable=variable, width=width)

    return "${}$".format("+".join(f()))
