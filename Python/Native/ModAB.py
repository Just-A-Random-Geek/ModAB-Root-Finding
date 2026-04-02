def modAB_root(f, x1, x2, y, xtol=1e-14, ytol=1e-16, maxiter=200):
    """
    Finds the root of f(x) = target within [x1, x2] using
    modified Anderson-Björk method (Ganchovski, Traykov).
    f(x) must be continuous and sign(f(x1) - y) ≠ sign(f(x2) - y).
    """
    if x2 < x1:
        x1, x2 = x2, x1 
        
    epsy = ytol * max(abs(y), 1)
    y1 = f(x1) - y
    if abs(y1) <= epsy:
        return x1

    y2 = f(x2) - y
    if abs(y2) <= epsy:
        return x2

    epsx = xtol * (x2 - x1)
    side = 0
    bisection = True
    threshold = x2 - x1  # Threshold to fall back to bisection if AB fails to shrink the interval enough
    for _ in range(maxiter):
        if bisection:
            x3 = (x1 + x2) * 0.5
            y3 = f(x3) - y  # Function value at midpoint
            ym = (y1 + y2) * 0.5  # Ordinate of chord at midpoint
            dy = y2 - y1
            r = 1 - abs(ym / dy)  # Symmetry factor
            # Check if the function is close enough to linear
            if abs(ym - y3) < r * r * (abs(y3) + abs(ym)):
                bisection = False
                threshold = (x2 - x1) * 16  # Safety factor: 4 bisection iterations = 2^4
        else:
            x3 = (x1 * y2 - y1 * x2) / (y2 - y1)
            if x3 <= x1:
                x3, y3 = x1, y1
            elif x3 >= x2:
                x3, y3 = x2, y2
            else:
                y3 = f(x3) - y
                
            threshold *= 0.5

        if abs(y3) <= epsy or x2 - x1 <= epsx:
            return x3

        if (y1 > 0) == (y3 > 0):  # Same sign check
            if side == 1:
                m = 1 - y3 / y1
                y2 = y2 * m if m > 0 else y2 * 0.5
            elif not bisection:
                side = 1
            x1, y1 = x3, y3
        else:
            if side == -1:
                m = 1 - y3 / y2
                y1 = y1 * m if m > 0 else y1 * 0.5
            elif not bisection:
                side = -1
            x2, y2 = x3, y3

        if x2 - x1 > threshold:  # AB failed to shrink the interval enough
            bisection = True
            side = 0

    return float('nan')
