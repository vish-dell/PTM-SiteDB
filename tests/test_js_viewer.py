"""Runs the jsdom viewer test (blank-viewer regression) when node + jsdom are available; otherwise skips."""
import os
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _jsdom():
    if not shutil.which("node"):
        return None
    for cand in (os.environ.get("JSDOM_PATH"), "jsdom"):
        if not cand:
            continue
        r = subprocess.run(["node", "-e", f"require({cand!r})"], capture_output=True, cwd=ROOT)
        if r.returncode == 0:
            return cand
    return None


@pytest.mark.skipif(_jsdom() is None, reason="needs node and `npm install jsdom` (see README_FIXES)")
def test_spatial_viewer_rebinds_to_live_container_after_rerun():
    env = dict(os.environ, JSDOM_PATH=_jsdom())
    r = subprocess.run(["node", "tests/js/spatial_viewer_rebind.test.js", "static/js/features.js",
                        "static/js/ptmviewer.js"], capture_output=True, text=True, cwd=ROOT, env=env, timeout=60)
    assert r.returncode == 0, r.stdout + r.stderr
    assert '"liveContainerHasCanvas":true' in r.stdout
