import json
from scripts import cohort_guard as G


def _summary(d, cohort, sha):
    (d / f"{cohort}_analysis_summary.json").write_text(json.dumps({"cohort": cohort, "input_provenance": {"maf_sha256": sha}}))


def test_identical_maf_under_two_labels_is_detected(tmp_path):
    maf = tmp_path / "a.maf"; maf.write_text("Hugo_Symbol\nTP53\n")
    _summary(tmp_path, "COHORT_A", G.file_sha256(maf))
    assert G.find_duplicate_input(maf, "COHORT_B", tmp_path) == "COHORT_A"
    assert G.find_duplicate_input(maf, "COHORT_A", tmp_path) is None        # same label re-run is fine


def test_different_maf_is_not_flagged(tmp_path):
    a = tmp_path / "a.maf"; a.write_text("x\n"); b = tmp_path / "b.maf"; b.write_text("y\n")
    _summary(tmp_path, "A", G.file_sha256(a))
    assert G.find_duplicate_input(b, "B", tmp_path) is None


def test_duplicate_groups():
    assert G.duplicate_groups({"A": "h", "B": "h", "C": "z"}) == [["A", "B"]]
    assert G.duplicate_groups({"A": "h", "C": "z"}) == []
