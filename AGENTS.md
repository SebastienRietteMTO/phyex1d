# AGENTS.md — phyex1d

## Dev setup
```bash
python3 -m venv venv
source venv/bin/activate
pip install -e .
```

## Tests
- No tests yet. Add new tests in `tests/` using **pytest**.

## Lint
```bash
flake8 src/
pylint src/phyex1d/
```

## Coding conventions
- **Docstrings**: NumPy style
- **Line length**: max 100 (enforced by flake8)
- **No type hints** currently used
- **Naming**: PEP8 (snake_case for functions/vars, PascalCase for classes)
- **Formatting**: no auto-formatter configured

See `pyproject.toml` for dependencies and entry point.
