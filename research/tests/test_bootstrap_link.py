from types import SimpleNamespace

import scripts.bootstrap_sources as bootstrap


def test_intogen_zip_link_in_query_parameter(monkeypatch):
    fake = SimpleNamespace(
        url="https://www.intogen.org/download",
        text='<a href="/download?file=IntOGen-Drivers-20240920.zip">download</a>',
        raise_for_status=lambda: None,
    )
    monkeypatch.setattr(bootstrap.requests, "get", lambda *a, **k: fake)
    found = bootstrap.find_intogen_zip_from_page()
    assert found == "https://www.intogen.org/download?file=IntOGen-Drivers-20240920.zip"
