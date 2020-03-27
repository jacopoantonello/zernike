# Zernike

Python code for Zernike polynomials. This module was part of
[enzpy](https://github.com/jacopoantonello/enzpy) but will be further developed
here instead.

![](./media/table.png)

# Example

```python
import numpy as np
import matplotlib.pyplot as plt
from zernike import RZern

cart = RZern(6)
L, K = 200, 250
ddx = np.linspace(-1.0, 1.0, K)
ddy = np.linspace(-1.0, 1.0, L)
xv, yv = np.meshgrid(ddx, ddy)
cart.make_cart_grid(xv, yv)

c = np.zeros(cart.nk)
plt.figure(1)
for i in range(1, 10):
    plt.subplot(3, 3, i)
    c *= 0.0
    c[i] = 1.0
    Phi = cart.eval_grid(c, matrix=True)
    plt.imshow(Phi, origin='lower', extent=(-1, 1, -1, 1))
    plt.axis('off')

plt.show()
```

# Installation

## Linux

```bash
python setup.py develop --user
```

## Windows

- You should first install the following requirements:
    - [Anaconda for Python 3](https://www.anaconda.com/download). This includes
      Python as well as some necessary scientific libraries.
    - [Git](https://git-scm.com/download/win). This is necessary for the
      automatic version numbering of this package. Also, make sure you choose
      *Git from the command line and also 3rd-party software* in *Adjusting
      your PATH environment*.
- *Clone* this repository using Git. From any folder in File Explorer,
  right-click and hit *Git Bash here*. Paste `git clone
  https://github.com/jacopoantonello/zernike` and hit enter. Do not use
  GitHub's *Download ZIP* button above, as the installation script will not
  work in that case.
- Finally, double-click on `install.bat`.

## Run tests

```bash
cd tests
nosetests -v -x --pdb *.py
```
