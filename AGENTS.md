# AGENTS.md — phyex1d

## Dev setup
```bash
python3 -m venv venv
source venv/bin/activate
pip install -e .
```
Requires Python >= 3.10.

## Runtime PHYEX dependency
The PHYEX Fortran package must be compiled before running physics:
```bash
source PHYEX/tools/env.sh
phyex1d ...
```
The compiled `pyphyex` module must be importable (via `PYTHONPATH` or `--exp pyphyex=<path>`).

## CLI
Entry point: `phyex1d` → `phyex1d.cli:main`
```bash
source PHYEX/tools/env.sh
phyex1d ./DEPHY-SCM/ARMCU/REF/ARMCU_REF_SCM_driver.nc --exp dt=60 --plot rc.png rc y_var=P
```
Experiment classes (via `--exp class=..`): `PhysicsAromeTQ` (default), `PhysicsAromeThetaR`, `SimplePhysicsTQ`, `SimplePhysicsThetaR`.

## Programmatic API
Main orchestrator:
```python
from phyex1d import ExperimentRunner, CaseXarray, CaseCommonFormat
runner = ExperimentRunner(case, experiments, output_dir)
runner.run()
```
`CaseCommonFormat` takes a netCDF path; `CaseXarray` takes an `xarray.Dataset` (see `examples/case_xarray_demo.py`).

## Defaults
- **Grid**: `L90arome` (looked up in `~/.phyex1d/grids/` then `src/phyex1d/grids/`; `.grid` files are JSON)
- **Namelist**: `default` (looked up in `~/.phyex1d/namelists/` then `src/phyex1d/namelists/`; f90nml format)

## Tests
- No tests yet. Add new tests in `tests/` using **pytest**.

## Lint
```bash
flake8 src/ --max-line-length=100
pylint src/phyex1d/
```

## Style
- **Docstrings**: NumPy style (Doxygen also used for API docs)
- **Line length**: max 100 (flake8 enforced)
- **No type hints** used
- **Naming**: PEP8
- **Formatting**: no auto-formatter
- **Pre-commit**: none

## CI
Only Doxygen docs deploy to GitHub Pages on push to `main` (`.github/workflows/doxygen.yaml`).
