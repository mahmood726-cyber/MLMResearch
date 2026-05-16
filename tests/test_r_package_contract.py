from pathlib import Path
import shutil
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[1]


def test_package_metadata_contract():
    description = (ROOT / "DESCRIPTION").read_text(encoding="utf-8")

    assert "Package: MLMResearch" in description
    assert "Title: Multilevel Meta-Analysis Method Evaluation and Improvement" in description
    assert "License: MIT + file LICENSE" in description
    assert "Config/testthat/edition:" not in description or "testthat" in description
    assert (ROOT / "NAMESPACE").exists()
    assert (ROOT / "README.md").exists()


def test_r_surface_and_testthat_files_present():
    r_files = sorted((ROOT / "R").glob("*.R"))
    testthat_files = sorted((ROOT / "tests" / "testthat").glob("test-*.R"))

    assert len(r_files) >= 5
    assert (ROOT / "R" / "adaptive_rve.R").exists()
    assert (ROOT / "tests" / "testthat.R").exists()
    assert testthat_files


def test_native_testthat_optional():
    rscript = shutil.which("Rscript")
    if not rscript:
        pytest.skip("Rscript is not installed in this environment")

    completed = subprocess.run(
        [rscript, "-e", "testthat::test_dir('tests/testthat')"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
