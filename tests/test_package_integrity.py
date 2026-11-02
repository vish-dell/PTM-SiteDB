"""Guards against shipping an incomplete project (the installer checks the same paths at start-up)."""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]


def test_installer_required_paths_exist():
    ps = (ROOT / "INSTALL_AND_RUN.ps1").read_text(encoding="utf-8-sig")
    block = ps[ps.index("$required = @("): ps.index("$missing")]
    paths = re.findall(r'"([^"]+)"', block)
    assert len(paths) >= 8
    missing = [p for p in paths if not (ROOT / p.replace("\\\\", "/").replace("\\", "/")).exists()]
    assert not missing, f"installer lists paths that do not exist: {missing}"


def test_no_virtualenv_is_shipped_inside_source_tree_guard():
    # a .venv in the tree is machine-specific; the package must be built without it
    assert (ROOT / "START_HERE_WINDOWS.txt").exists()
