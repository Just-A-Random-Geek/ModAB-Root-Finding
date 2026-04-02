A fast and robust root-finding library using the Modified Anderson-Bjork method (Ganchovski, Traykov), written in C for Python.


## Installation

```bash
pip install pymodab
```

## Usage

```python
import math
from pymodab import find_root, find_default, get_evaluation_count

# Find the root of cos(x) - x = 0 in [0, 1]
root = find_root(lambda x: math.cos(x) - x, 0, 1, 100)
print(f"Root: {root}")  # 0.7390851332151607

# Using default tolerances
root = find_default(lambda x: x**2 - 2, 1, 2)
print(f"sqrt(2) = {root}")  # 1.4142135623730951

# Get the number of function evaluations
print(f"Evaluations: {get_evaluation_count()}")
```

## API

### `find_root(f, x1, x2, atol=1e-14, rtol=1e-14, max_iter=200)`

Find the root of `f(x) = 0` within the interval `[x1, x2]`.

**Parameters:**
- `f`: A continuous function of one variable
- `x1`, `x2`: Bracket interval endpoints (must satisfy `f(x1) * f(x2) < 0`)
- `atol`: Absolute tolerance (default: 1e-14)
- `rtol`: Relative tolerance (default: 1e-14)
- `max_iter`: Maximum iterations (default: 200)

**Returns:** The root, or `NaN` if not found.

### `find_default(f, x1, x2)`

Convenience wrapper for `find_root` with default tolerances.

### `get_evaluation_count()`

Returns the number of function evaluations from the last root-finding call.

## Algorithm

Modified Anderson-Björck's method is a new robust and efficient bracketing root-finding algorithm. It combines bisection with Anderson-Björk's method to achieve both fast performance and worst-case optimality.

References:

Ganchovski N.; Traykov A. Modified Anderson-Björck's method for solving non-linear equations in structural mechanics. IOP Conference Series: Materials Science and Engineering 2023, 1276 (1) 012010, IOP Publishing. 
https://iopscience.iop.org/article/10.1088/1757-899X/1276/1/012010/pdf

Ganchovski, N.; Smith, O.; Rackauckas, C.; Tomov, L.; Traykov, A. Improvements of the Modified Anderson-Björck (modAB) Root-Finding Algorithm. Preprints 2026, 2026032190.
https://www.preprints.org/manuscript/202603.2190

## License

MIT License
