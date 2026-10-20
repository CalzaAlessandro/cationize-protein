# Working with the private repository

Repository: https://github.com/CalzaAlessandro/cationize-protein

Only the owner and explicitly granted collaborators can access this private
repository. Making it public is a separate future decision.

The **Code** tab contains the current files. `README.md` is the starting page;
`docs/` explains methods and outputs; `src/` contains the program; `tests/`
checks behavior; `parameters/` holds the accepted fitted libraries.

The **Issues** tab can track requested improvements or reproducible problems.
**Commits** records saved development steps. A commit is a named snapshot;
a push copies local commits to GitHub without making a private repository
public or installing the software on another computer.

## Local use

```sh
git clone https://github.com/CalzaAlessandro/cationize-protein.git
cd cationize-protein
python3 install.py
.venv/bin/cationize-protein --help
```

Cloning requires an authenticated account with repository access. Use
`git status` and `git diff` to inspect local changes. Commit only intended
source, documentation, tests and parameters. Environments, external binaries,
`.env` files and `.development/` are ignored; production trajectories and
private quantum calculations stay outside Git.

The optional GitHub Actions workflow is manual: open **Actions**, select
**Python tests**, then **Run workflow**. It installs the Python package and
runs unit tests. It does not install licensed software, run HSA MD, or
establish force-field accuracy on arbitrary proteins.

Before a public release, choose the licence and citation/authorship details,
review redistribution terms for reused parameters, and qualify supported
installation platforms. This initial private repository is a development
workspace.
