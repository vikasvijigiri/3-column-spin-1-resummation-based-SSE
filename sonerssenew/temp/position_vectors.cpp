#include "variables.hpp"

void POSvectors::pbc() {
  int k, x1, y1, i, d1, d2;

  for (k = 1; k <= 2 * lx * ly; ++k) {
    // Block 1(Bottom left)
    if (k <= lx * ly) {
      if (k % lx <= lx / 2 && k <= (int((k - 1) / lx) + 1) * lx && k <= (Np / 2) + lx) {
        if (k % lx != 0) {
          x1 = (k % lx) - 1;
          y1 = int((k - 1) / lx);
          d1 = x1 + y1;
          d2 = y1 - x1;
        } else {
          x1 = (-k % lx) - 1;
          y1 = int((k - 1) / lx);
          d1 = x1 + y1;
          d2 = y1 - x1;
        }
      }
      // Block 2(Botton right)
      else if (k % lx > lx / 2 && k <= (int((k - 1) / lx) + 1) * lx && k <= (Np / 2) + lx) {
        x1 = (k % lx) - 1 - lx;
        y1 = int((k - 1) / lx);
        d1 = x1 + y1;
        d2 = y1 - x1;
      }
      // Block 3(top left)
      else if (k % lx <= lx / 2 && k <= (int((k - 1) / lx) + 1) * lx && k >= (Np / 2) + lx) {
        x1 = (k % lx) - 1;
        y1 = int((k - 1) / lx) - lx;
        d1 = x1 + y1;
        d2 = y1 - x1;
      }
      // Block 4(top right)
      else if (k % lx > lx / 2 && k <= (int((k - 1) / lx) + 1) * lx && k >= (Np / 2) + lx) {
        x1 = (k % lx) - 1 - lx;
        y1 = int((k - 1) / lx) - lx;
        d1 = x1 + y1;
        d2 = y1 - x1;
      }
      psite[k - 1][0] = d1;
      psite[k - 1][1] = d2;

    } else {
      if ((k - 1) % lx <= lx / 2 && k <= (int((k - 1) / lx) + 1) * lx && k <= 2 * Np - Np / 2) {
        x1 = k % lx - 1 + int(k / lx) - lx;
        y1 = int(((k - Np) - 1) / lx) + 1 - (k - 1) % lx;
      } else if ((k - 1) % lx > lx / 2 && k <= (int((k - 1) / lx) + 1) * lx && k <= 2 * Np - Np / 2) {
        x1 = (k - 1) % lx - lx + int((k - Np - 1) / lx);
        y1 = lx - (k - 1) % lx + int((k - Np - 1) / lx) + 1;
      } else if ((k - 1) % lx <= lx / 2 && k <= (int((k - 1) / lx) + 1) * lx && k >= 2 * Np - Np / 2) {
        x1 = (k - Np) % lx - 1 + int((k - Np) / lx) - lx;
        y1 = int(((k - 2 * Np)) / lx) - (k - Np - 1) % lx;
      } else if ((k - 1) % lx > lx / 2 && k <= (int((k - 1) / lx) + 1) * lx && k >= 2 * Np - Np / 2) {
        x1 = (k - 1) % lx - lx - int((2 * Np - k + 1) / lx) - 1;
        y1 = -(k - 1) % lx + lx - int((2 * Np - k + 1) / lx);
      }

      psite[k - 1][0] = x1;
      psite[k - 1][1] = y1;
    }
  } // for loop
} // end func pbc
