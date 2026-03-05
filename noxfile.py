import os
import pathlib
import shutil

import nox

ROOT = pathlib.Path(__file__).parent


nox.options.sessions = ["lint", "test", "test-notebooks"]


@nox.session
def test(session: nox.Session) -> None:
    """Run the tests."""
    session.install("-r", "requirements-testing.txt")
    session.install(".")

    args = [
        "-n",
        "auto",
        "--cov",
        "cascade",
        "-vvv",
    ] + session.posargs

    if "CI" in os.environ:
        args.append(f"--cov-report=xml:{ROOT.absolute()!s}/coverage.xml")
    session.run("pytest", *args)

    if "CI" not in os.environ:
        session.run("coverage", "report", "--ignore-errors", "--show-missing")


@nox.session(name="test-notebooks")
def test_notebooks(session: nox.Session) -> None:
    """Run the notebooks."""
    args = [
        "pytest",
        "notebooks",
        "--nbmake",
        "--nbmake-kernel=python3",
        "--nbmake-timeout=3000",
        "-n",
        "auto",
        "-vvv",
    ] + session.posargs

    session.install("-r", "requirements-testing.txt")
    session.install("nbmake")
    session.install("-r", "notebooks/requirements.in")
    session.install(".")

    session.run(*args)


@nox.session
def lint(session: nox.Session) -> None:
    """Look for lint."""
    session.install("pre-commit")
    session.run("pre-commit", "run", "--all-files")


@nox.session
def build(session: nox.Session) -> None:
    """Build sdist and wheel dists."""
    session.install("pip")
    session.install("build")
    session.run("python", "--version")
    session.run("pip", "--version")
    session.run("python", "-m", "build", "--outdir", "./build/wheelhouse")


@nox.session
def release(session):
    """Tag, build and publish a new release to PyPI."""
    session.install("zest.releaser[recommended]")
    session.run("fullrelease")


@nox.session(name="publish-testpypi")
def publish_testpypi(session):
    """Publish wheelhouse/* to TestPyPI."""
    session.install("twine")
    session.run("twine", "check", "build/wheelhouse/*")
    session.run(
        "twine",
        "upload",
        "--skip-existing",
        "--repository-url",
        "https://test.pypi.org/legacy/",
        "build/wheelhouse/*.tar.gz",
    )


@nox.session(name="publish-pypi")
def publish_pypi(session):
    """Publish wheelhouse/* to PyPI."""
    session.install("twine")
    session.run("twine", "check", "build/wheelhouse/*")
    session.run(
        "twine",
        "upload",
        "--skip-existing",
        "build/wheelhouse/*.tar.gz",
    )


@nox.session(python=False)
def clean(session):
    """Remove all .venv's, build files and caches in the directory."""
    folders_to_remove = (
        "coastal_cascade.egg-info",
        ".pytest_cache",
        ".venv",
        "build",
        "dist",
    )
    folders_to_clean = ("cascade", "notebooks", "tests")
    patterns_to_clean = (
        "*.py[co]",
        "__pycache__",
        "*.c",
        "*.so",
        "*-checkpoint.ipynb",
        ".ipynb_checkpoints",
    )

    with session.chdir(ROOT):
        for folder in folders_to_remove:
            session.log(f"rm -r {folder}")
            shutil.rmtree(folder, ignore_errors=True)

        for folder in folders_to_clean:
            with session.chdir(folder):
                for pattern in patterns_to_clean:
                    session.log(f"rm -r {folder}/**/{pattern}")
                    _clean_rglob(pattern)


def _clean_rglob(pattern):
    nox_dir = pathlib.Path(".nox")

    for p in pathlib.Path(".").rglob(pattern):
        if nox_dir in p.parents:
            continue
        if p.is_dir():
            p.rmdir()
        else:
            p.unlink()
