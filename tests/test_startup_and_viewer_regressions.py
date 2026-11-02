from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_main_app_is_the_complete_startup_entrypoint():
    app = (ROOT / "app.py").read_text(encoding="utf-8")
    assert 'app.register_blueprint(features_bp)' in app
    assert 'DispatcherMiddleware' in app
    assert '"/research": research_app' in app
    assert 'app.run(debug=False, host="127.0.0.1", port=5000, use_reloader=False)' in app


def test_catalog_protein_viewer_does_not_require_ptm_rows():
    app = (ROOT / "app.py").read_text(encoding="utf-8")
    start = app.index('@app.route("/viewer/<uniprot>")')
    end = app.index('@app.route("/api/search")', start)
    route = app[start:end]
    assert 'if not rows and not cat' in route
    assert 'ptms=[dict(r) for r in rows]' in route


def test_every_protein_detail_and_feature_page_offers_the_structure_viewer():
    protein = (ROOT / "templates" / "protein.html").read_text(encoding="utf-8")
    feature = (ROOT / "templates" / "feature_base.html").read_text(encoding="utf-8")
    assert 'Open 3D Structure / PTM Viewer' in protein
    assert 'feature-viewer-link' in feature
    assert 'v=20261010-evidence-fix' in feature


def test_busy_state_does_not_reference_uninitialized_hypothesis_payload():
    js = (ROOT / "static" / "js" / "features.js").read_text(encoding="utf-8")
    start = js.index('const busy =')
    end = js.index('async function api(path', start)  # signature gained a timeout parameter
    assert 'j.' not in js[start:end]
    assert 'Array.isArray(j.hypotheses)' in js


def test_genomic_preparation_script_uses_project_root_virtualenv_and_resolves_inputs_first():
    script = (ROOT / "research" / "PREPARE_GENOMIC_INPUTS.ps1").read_text(encoding="utf-8-sig")
    assert '$projectRoot = Split-Path -Parent $PSScriptRoot' in script
    assert r'Join-Path $projectRoot ".venv312\Scripts\python.exe"' in script
    assert script.index('$callerDirectory = (Get-Location).Path') < script.index('Set-Location -LiteralPath $PSScriptRoot')


def test_viewer_checks_residue_identity_and_labels_predicted_coordinates():
    js = (ROOT / "static" / "js" / "viewer.js").read_text(encoding="utf-8")
    assert "function expectedResidue(site)" in js
    assert "function verifiedResidue(mapping, ptm)" in js
    assert "Predicted AlphaFold coordinates—not an experimentally solved structure" in js
    assert "No coordinates or PTM mappings were fabricated" in js
    assert "v6.pdb" in js and "v4.pdb" in js
