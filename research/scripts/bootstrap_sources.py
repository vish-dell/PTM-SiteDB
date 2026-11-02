from __future__ import annotations

import gzip
import hashlib
import json
import re
import shutil
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urljoin, urlparse

import requests
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
RAW.mkdir(parents=True, exist_ok=True)
PROCESSED.mkdir(parents=True, exist_ok=True)

DBPTM_URL = "https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Phosphorylation.zip"
INTOGEN_PAGE = "https://www.intogen.org/download"
INTOGEN_FALLBACK_NAME = "IntOGen-Drivers-20240920.zip"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, out: Path, timeout: int = 180) -> None:
    print(f"Downloading {url}")
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(out.suffix + ".partial")
    try:
        with requests.get(url, stream=True, timeout=timeout, headers={"User-Agent": "PTM-SiteDB-Integrated/3.0"}) as r:
            r.raise_for_status()
            total = int(r.headers.get("content-length") or 0)
            done = 0
            with tmp.open("wb") as f:
                for chunk in r.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        f.write(chunk)
                        done += len(chunk)
                        if total:
                            print(f"  {done/total:.0%}", end="\r")
        if tmp.stat().st_size < 1000:
            raise RuntimeError("Downloaded file is unexpectedly small; refusing to keep it as a source archive.")
        with zipfile.ZipFile(tmp) as z:
            if z.testzip() is not None:
                raise RuntimeError("ZIP integrity check failed.")
        tmp.replace(out)
    except Exception:
        if tmp.exists():
            tmp.unlink()
        raise
    print()


def find_intogen_zip_from_page() -> str | None:
    """Find the official IntOGen archive even when the ZIP filename is a query parameter."""
    try:
        r = requests.get(INTOGEN_PAGE, timeout=45, headers={"User-Agent": "PTM-SiteDB-Integrated/1.0"})
        r.raise_for_status()
    except Exception as exc:
        print(f"IntOGen page download failed: {exc}")
        return None
    hrefs = re.findall(r'href=["\']([^"\']+)["\']', r.text, flags=re.I)
    candidates = []
    for href in hrefs:
        absolute = urljoin(r.url, href)
        # IntOGen commonly serves files as /download?file=...zip, so checking
        # only urlparse(...).path incorrectly misses the actual archive link.
        url_text = absolute.lower()
        if "driver" in url_text and ".zip" in url_text:
            candidates.append(absolute)
    if candidates:
        # Prefer the specifically named driver compendium if multiple ZIPs are linked.
        candidates.sort(key=lambda value: ("intogen-drivers-20240920.zip" not in value.lower(), value))
        return candidates[0]
    # Current official-style parameter URL fallback. Download validation still checks ZIP integrity.
    return urljoin(r.url, "/download?file=IntOGen-Drivers-20240920.zip")


def _safe_extract(z: zipfile.ZipFile, destination: Path) -> None:
    """Extract official archives without allowing member paths to escape the target directory."""
    base = destination.resolve()
    for member in z.infolist():
        target = (base / member.filename).resolve()
        if target != base and base not in target.parents:
            continue
        z.extract(member, base)


def _expand_nested_zips(root: Path, max_depth: int = 2) -> None:
    """Some source downloads wrap the actual table in a second ZIP."""
    for _ in range(max_depth):
        nested = [p for p in root.rglob("*.zip") if p.is_file() and not any(part.endswith("_nested") for part in p.parts)]
        changed = False
        for zpath in nested:
            target = zpath.parent / f"{zpath.stem}_nested"
            if target.exists():
                continue
            try:
                target.mkdir(parents=True, exist_ok=True)
                with zipfile.ZipFile(zpath) as z:
                    if z.testzip() is not None:
                        continue
                    _safe_extract(z, target)
                changed = True
            except (OSError, zipfile.BadZipFile):
                continue
        if not changed:
            break


def _decompress_tabular_files(root: Path) -> None:
    """Unpack gzip-compressed text tables shipped inside ZIP files."""
    for source in list(root.rglob("*")):
        if not source.is_file() or source.suffix.lower() not in {".gz", ".bgz"}:
            continue
        target = source.with_suffix("")
        if target.exists():
            continue
        try:
            with gzip.open(source, "rb") as src, target.open("wb") as dst:
                shutil.copyfileobj(src, dst)
        except (OSError, EOFError):
            try:
                target.unlink(missing_ok=True)
            except OSError:
                pass


def _candidate_header(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        try:
            frame = pd.read_excel(path, nrows=8, dtype=str)
            header = "\t".join(map(str, frame.columns))
            preview = frame.head(3).to_csv(index=False, sep="\t")
            return (header + "\n" + preview)[:12000].lower()
        except Exception as exc:
            return f"excel_read_error:{exc}".lower()
    try:
        with path.open("rb") as stream:
            raw = stream.read(12000)
        return raw.decode("utf-8-sig", errors="ignore").lower()
    except OSError:
        return ""


def _convert_spreadsheet(path: Path) -> Path:
    """Convert a supported spreadsheet table to UTF-8 TSV for the normal parser."""
    frame = pd.read_excel(path, dtype=str)
    if frame.empty or len(frame.columns) < 2:
        raise RuntimeError(f"Spreadsheet contains no usable tabular data: {path.name}")
    out = path.with_name(path.stem + "_normalized.tsv")
    frame.to_csv(out, sep="\t", index=False)
    return out


def discover_table_from_zip(zip_path: Path, target_name: str, preferred_keywords: list[str],
                            exclude_keywords: tuple[str, ...] = ()) -> Path:
    """Discover a real source table across text, compressed-text and spreadsheet archive layouts.

    The selection is based on source filenames and table headers/content; documentation files
    are excluded. On failure, report the archive inventory rather than only saying no table exists.
    """
    extract_dir = PROCESSED / "_extract" / target_name
    if extract_dir.exists():
        shutil.rmtree(extract_dir)
    extract_dir.mkdir(parents=True)
    with zipfile.ZipFile(zip_path) as z:
        _safe_extract(z, extract_dir)
    _expand_nested_zips(extract_dir)
    _decompress_tabular_files(extract_dir)

    table_suffixes = {".txt", ".tsv", ".csv", ".tab", ".dat", ".psv", ".xlsx", ".xls"}
    all_files = [p for p in extract_dir.rglob("*") if p.is_file()]
    ignored_words = ("readme", "license", "citation", "metadata", "changelog", "manual", "help")
    scored = []
    for path in all_files:
        if path.suffix.lower() not in table_suffixes and path.suffix:
            continue
        if any(word in path.name.lower() for word in ignored_words):
            continue
        if any(word in path.name.lower() for word in exclude_keywords):
            continue
        if path.stat().st_size < 100:
            continue
        name = path.name.lower()
        head = _candidate_header(path)
        score = sum(10 for keyword in preferred_keywords if keyword.lower() in name)
        score += sum(3 for keyword in preferred_keywords if keyword.lower() in head)
        # A named data table with a readable multi-column header can be a valid headerless
        # source, but documentation or an unrelated table must not pass without source clues.
        if score:
            scored.append((score, path))

    if not scored:
        inventory = [f"{path.relative_to(extract_dir)} ({path.stat().st_size} bytes)" for path in all_files[:80]]
        detail = "; ".join(inventory) or "archive contained no extractable files"
        raise RuntimeError(
            f"No candidate table matched source keywords in {zip_path.name}. "
            f"Archive inventory (first 80 files): {detail}. "
            "Check the official archive version and inspect its table format; no records were fabricated."
        )
    scored.sort(key=lambda item: (-item[0], -item[1].stat().st_size, str(item[1]).lower()))
    winner = scored[0][1]
    if winner.suffix.lower() in {".xlsx", ".xls"}:
        try:
            winner = _convert_spreadsheet(winner)
        except Exception as exc:
            raise RuntimeError(f"Candidate spreadsheet {winner.name} could not be converted: {exc}") from exc
    return winner


def declared_intogen_release(extract_root: Path) -> str | None:
    """Read the release declared inside the downloaded archive; do not trust filename dates."""
    for p in extract_root.rglob("*"):
        if p.is_file() and p.name.lower().startswith("readme"):
            try:
                text = p.read_text(encoding="utf-8", errors="ignore")
            except Exception:
                continue
            m = re.search(r"IntOGen\s+RELEASE\s+([0-9]{8}|[0-9]{4}\.[0-9]{2}\.[0-9]{2})", text, re.I)
            if m:
                return m.group(1)
    return None


def main() -> int:
    manifest = {"generated_utc": datetime.now(timezone.utc).isoformat(), "sources": {}}
    dbptm_zip = RAW / "dbptm_phosphorylation.zip"
    intogen_zip = RAW / INTOGEN_FALLBACK_NAME

    # A prior verified ZIP is reused. Failed remote retrieval is recorded and the website may still launch.
    if not dbptm_zip.exists():
        try:
            download(DBPTM_URL, dbptm_zip)
        except Exception as exc:
            print(f"dbPTM download failed: {exc}")
            manifest["sources"]["dbPTM_phosphorylation"] = {"status": "download_failed", "url": DBPTM_URL, "error": str(exc)}
    if dbptm_zip.exists():
        try:
            with zipfile.ZipFile(dbptm_zip) as z:
                if z.testzip() is not None:
                    raise RuntimeError("ZIP integrity check failed")
            manifest["sources"]["dbPTM_phosphorylation"] = {"url": DBPTM_URL, "path": str(dbptm_zip.relative_to(ROOT)), "sha256": sha256(dbptm_zip)}
        except Exception as exc:
            manifest["sources"]["dbPTM_phosphorylation"] = {"status": "invalid_archive", "error": str(exc)}
            try: dbptm_zip.unlink()
            except OSError: pass

    if not intogen_zip.exists():
        url = find_intogen_zip_from_page()
        if not url:
            print(f"Could not discover the official IntOGen ZIP. Place {INTOGEN_FALLBACK_NAME} in {RAW} and rerun.")
        else:
            try:
                download(url, intogen_zip)
            except Exception as exc:
                print(f"Automatic IntOGen download failed: {exc}")
                print(f"Place {INTOGEN_FALLBACK_NAME} in {RAW} and rerun.")
    if intogen_zip.exists():
        manifest["sources"]["IntOGen"] = {"download_page": INTOGEN_PAGE, "path": str(intogen_zip.relative_to(ROOT)), "sha256": sha256(intogen_zip)}
    else:
        manifest["sources"]["IntOGen"] = {"status": "pending_manual_download", "download_page": INTOGEN_PAGE}

    if dbptm_zip.exists() and manifest["sources"].get("dbPTM_phosphorylation", {}).get("sha256"):
        try:
            ptm_source = discover_table_from_zip(dbptm_zip, "dbptm_phosphorylation", ["phosphorylation", "uniprot", "position", "modified"])
            ptm_dest = PROCESSED / "dbptm_phosphorylation_source.txt"
            shutil.copy2(ptm_source, ptm_dest)
            manifest["sources"]["dbPTM_phosphorylation"]["extracted_table"] = str(ptm_dest.relative_to(ROOT))
        except Exception as exc:
            manifest["sources"]["dbPTM_phosphorylation"].update({"status": "table_discovery_failed", "error": str(exc)})

    if intogen_zip.exists() and manifest["sources"].get("IntOGen", {}).get("sha256"):
        try:
            intogen_extract = PROCESSED / "_extract" / "intogen_drivers"
            # The archive holds BOTH the curated catalogue (Compendium_Cancer_Genes.tsv, ~630 drivers) and
            # Unfiltered_drivers.tsv (every gene tested, ~7,700). Matching on the word "drivers" picked the unfiltered
            # table, which made nearly every gene a "driver". Name the curated catalogue and exclude the other.
            table = discover_table_from_zip(intogen_zip, "intogen_drivers",
                                            ["compendium_cancer_genes", "compendium", "symbol"],
                                            exclude_keywords=("unfiltered",))
            declared = declared_intogen_release(intogen_extract)
            drivers_dest = PROCESSED / "intogen_drivers_source.tsv"
            shutil.copy2(table, drivers_dest)
            src = manifest["sources"]["IntOGen"]
            src["extracted_table"] = str(drivers_dest.relative_to(ROOT))
            src["release_declared_inside_archive"] = declared
            src["archive_name_date_is_not_authoritative"] = True
            if declared and re.fullmatch(r"[0-9]{8}", declared):
                try:
                    archive_dt = datetime.strptime(declared, "%Y%m%d").date()
                    src["release_date_interpretation"] = str(archive_dt)
                except ValueError:
                    pass
        except Exception as exc:
            manifest["sources"]["IntOGen"].update({"status": "table_discovery_failed", "error": str(exc)})

    (PROCESSED / "source_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    if not (PROCESSED / "dbptm_phosphorylation_source.txt").exists():
        print("WARNING: dbPTM table unavailable; gene panel cannot be built until source ingestion succeeds.")
    if not (PROCESSED / "intogen_drivers_source.tsv").exists():
        print(f"WARNING: IntOGen table unavailable; add the official ZIP to {RAW} and rerun this script.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
