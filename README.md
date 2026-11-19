# caustic_finder

Tools for computing gravitational microlensing caustics and critical curves,
based on the `causticfinder` program by Robert W. Schmidt (Astronomisches
Rechen-Institut, Heidelberg), which implements Hans Witt's method
(PhD thesis, Hamburg, 1991).

## Directory layout

| Directory | Description |
|---|---|
| `mod_code_v1/` | Original single-file version of `causticfinder` (C++). Computes critical and caustic curves for a set of point-mass lenses. |
| `mod_code_v2/` | Modified two-part version: `causticfinder_init` (initialization) and `causticfinder_cmap` (OpenMP-parallelized map computation). |
| `performance/` | Jupyter notebooks and timing data (`cmaps_time.dat`) used to benchmark the v1/v2 code. |

## Requirements

- C++ compiler (`g++`)
- [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/) with development headers
- OpenMP support (v2 only)

On Arch Linux: `sudo pacman -S gsl`

## Building

Each code version has its own `Makefile`:

```bash
cd mod_code_v1 && make        # builds ./causticfinder
cd mod_code_v2 && make        # builds ./causticfinder_init and ./causticfinder_cmap
```

Note: the Makefiles link against `-L/usr/lib64`; adjust `LIBS` if your GSL
libraries live elsewhere (e.g. `/usr/lib`).

## Usage

Star positions are given in lens-plane Einstein radii,
$x = \sqrt{4GM/c^2 D_d D_{ds}/D_s}$; all stars currently have the same mass.

Example run (v1):

```bash
./causticfinder stars.txt critcaust.txt 2 0.0 -0.5
```

- `stars.txt` — input file with star positions (one per line, e.g. `350, 500`)
- `critcaust.txt` — output file
- `2` — number of stars
- `0.0` — $\sigma_c$ (if zero, gamma is reduced shear)
- `-0.5` — gamma

The output contains 4 columns: x and y coordinates of the critical curves
(columns 1–2) and caustic curves (columns 3–4). By default 200 points are
computed per curve segment; change `M` at the top of the source to adjust.

See `mod_code_v1/README` and `mod_code_v2/README` for the original
documentation.

## License

The code is distributed under the GNU General Public License v3 or later
(see `COPYING` in each code directory).
