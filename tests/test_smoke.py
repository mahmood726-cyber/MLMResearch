import shutil
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


def _find_rscript():
    candidates = [
        shutil.which("Rscript"),
        "/mnt/c/Program Files/R/R-4.5.2/bin/Rscript.exe",
        "/mnt/c/Program Files/R/R-4.5.2/bin/x64/Rscript.exe",
        "/mnt/c/Program Files/R/R-4.5.1/bin/Rscript.exe",
        "/mnt/c/Program Files/R/R-4.5.1/bin/x64/Rscript.exe",
    ]
    for candidate in candidates:
        if candidate and Path(candidate).exists():
            return candidate
    return None


def test_core_package_files_exist():
    required = [
        ROOT / "DESCRIPTION",
        ROOT / "NAMESPACE",
        ROOT / "README.md",
        ROOT / "R" / "adaptive_rve.R",
        ROOT / "R" / "improved_tau2.R",
        ROOT / "R" / "methods_registry.R",
        ROOT / "tests" / "testthat.R",
        ROOT / "tests" / "testthat" / "test-adaptive_rve.R",
        ROOT / "e156-submission" / "index.html",
    ]
    missing = [str(path.relative_to(ROOT)) for path in required if not path.exists()]
    assert not missing, f"Missing required MLMResearch files: {missing}"


def test_description_and_namespace_contract():
    description = (ROOT / "DESCRIPTION").read_text(encoding="utf-8")
    namespace = (ROOT / "NAMESPACE").read_text(encoding="utf-8")

    assert "Package: MLMResearch" in description
    assert "License: MIT + file LICENSE" in description
    for export in [
        "export(adaptive_rve)",
        "export(get_rve_recommendation)",
        "export(improved_tau2_estimator)",
        "export(validate_method_comparison)",
    ]:
        assert export in namespace


def test_r_test_suite_is_present_and_named():
    testthat_entry = (ROOT / "tests" / "testthat.R").read_text(encoding="utf-8")
    test_file = (ROOT / "tests" / "testthat" / "test-adaptive_rve.R").read_text(
        encoding="utf-8"
    )

    assert 'test_check("MLMResearch")' in testthat_entry
    assert "adaptive_rve selects correct method" in test_file
    assert "get_rve_recommendation" in test_file


def test_r_test_suite_runs_when_rscript_is_available():
    rscript = _find_rscript()
    if not rscript:
        pytest.skip("Rscript is not installed in this environment")

    result = subprocess.run(
        [rscript, "--vanilla", "-e", "testthat::test_local('.')"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )
    assert result.returncode == 0, result.stderr or result.stdout

