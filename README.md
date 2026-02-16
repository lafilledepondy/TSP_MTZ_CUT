# TSP_MTZ_CUT

Solve the Asymmetric Traveling Salesman Problem (ATSP) using a MILP formulation with Miller-Tucker-Zemlin (MTZ) subtour elimination and the Gurobi optimizer.

## Requirements

- C++14 compiler
- CMake 3.14+
- Gurobi (with a valid license)

The default Gurobi path is set in CMake. Update it if your install differs:

- Linux: `/opt/gurobi1301/linux64`
- Windows: `C:/gurobi952/win64`
- macOS: `/Library/gurobi952/mac64`

## Build

From the project root:

```bash
mkdir -p build
cd build
cmake ..
make TSP_Gurobi
```

## Run

From the build directory:

```bash
./TSP_Gurobi data/att48.tsp
```

The `data/` symlink is created automatically, so instances can be referenced with a relative path.

### Modes

You can choose the solve mode as the second argument:

```bash
# MTZ (default)
./TSP_Gurobi data/br17.atsp MTZ

# CUT (integer MIP with lazy subtour cuts)
./TSP_Gurobi data/br17.atsp CUT

# CUT_LP (fractional LP with iterative cut generation)
./TSP_Gurobi data/br17.atsp CUT_LP
```

### Summary Output (for scripts)

Use `--summary` to print a single-line, machine-readable result:

```bash
./TSP_Gurobi data/att48.tsp MTZ --summary
```

This prints a line like:

```
RESULT instance=att48.tsp mode=MTZ obj=... bound=... nodes=... cuts=... status=... time=...
```

### Generate LaTeX Results Table

The script [scripts/generate_results.py](scripts/generate_results.py) runs all instances in `data/` and writes a LaTeX table to `results.tex`:

```bash
python3 scripts/generate_results.py
```

## Outputs

During execution, the solver writes these files in the build directory:

- `model.lp`: exported MILP model
- `solution.sol`: solver solution (if found)
- `atsp_mtz.log`: Gurobi log

## Notes

- The time limit is set to 180 seconds in the solver.
- Single-threaded solve is enforced.
