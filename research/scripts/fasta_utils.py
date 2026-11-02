"""Standard-library FASTA indexing and random access for Windows/Linux/macOS.

Supports ordinary uncompressed, line-wrapped FASTA and the standard samtools .fai
five-column index. Coordinates passed to fetch() are 0-based, half-open. The code
creates an index when missing, avoiding the native HTSlib build required by pysam.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class FaiRecord:
    length: int
    offset: int
    line_bases: int
    line_width: int


class IndexedFasta:
    def __init__(self, path: str | Path):
        self.path = Path(path)
        if not self.path.is_file():
            raise FileNotFoundError(f"FASTA file not found: {self.path}")
        self.index_path = Path(str(self.path) + ".fai")
        if not self.index_path.exists():
            self.build_index()
        self.index = self._read_index()
        self.references = list(self.index.keys())
        if not self.references:
            raise ValueError(f"No FASTA entries indexed in {self.path}")
        self._handle = self.path.open("rb")

    def build_index(self, progress: bool = True) -> None:
        """Create a standard .fai index from a consistently line-wrapped FASTA.

        Streams the file once with constant memory and prints progress, because a whole genome
        (~3 GB) takes a few minutes in pure Python and must not look frozen.
        """
        import sys
        import time

        rows: list[tuple[str, int, int, int, int]] = []
        names_seen: set[str] = set()
        current_name: str | None = None
        total_len = 0
        first_offset = 0
        line_bases = 0
        line_width = 0
        prev_len = 0
        prev_width = 0
        n_seq_lines = 0

        def finish() -> None:
            nonlocal current_name, total_len, first_offset, line_bases, line_width, prev_len, prev_width, n_seq_lines
            if current_name is None:
                return
            if n_seq_lines == 0 or total_len <= 0:
                raise ValueError(f"FASTA record {current_name!r} has no sequence")
            if prev_len > line_bases:
                raise ValueError(f"FASTA record {current_name!r} has an invalid final sequence line")
            rows.append((current_name, total_len, first_offset, line_bases, line_width))
            current_name = None
            total_len = first_offset = line_bases = line_width = prev_len = prev_width = n_seq_lines = 0

        size = self.path.stat().st_size
        step = 256 * 1024 * 1024
        next_report = step
        started = time.time()
        offset = 0
        if progress and size > step:
            print(f"Indexing FASTA ({size / 1e9:.1f} GB). A whole genome can take a few minutes - please do not press Ctrl+C.",
                  flush=True)
        with self.path.open("rb", buffering=1024 * 1024) as handle:
            for raw in handle:
                start = offset
                offset += len(raw)
                if progress and offset >= next_report:
                    print(f"  ...{offset / 1e6:,.0f} / {size / 1e6:,.0f} MB "
                          f"({100 * offset / size:.0f}%), {time.time() - started:.0f}s elapsed", flush=True)
                    next_report += step
                if raw[:1] == b">":
                    finish()
                    header = raw[1:].strip().split()
                    if not header:
                        raise ValueError("FASTA contains an empty header")
                    current_name = header[0].decode("utf-8", errors="replace")
                    if current_name in names_seen:
                        raise ValueError(f"Duplicate FASTA sequence identifier: {current_name}")
                    names_seen.add(current_name)
                    continue
                if current_name is None:
                    if raw.strip():
                        raise ValueError("FASTA sequence content occurs before the first '>' header")
                    continue
                seq_len = len(raw.rstrip(b"\r\n"))
                if seq_len == 0:
                    continue
                if line_bases == 0:
                    first_offset = start
                    line_bases = seq_len
                    line_width = len(raw)
                else:
                    # every sequence line except the last of a record must have the full width
                    if prev_len != line_bases:
                        raise ValueError(
                            f"FASTA record {current_name!r} has inconsistent sequence line lengths; "
                            "reformat it to a consistent wrapped FASTA before indexing."
                        )
                    if prev_width != line_width:
                        raise ValueError(
                            f"FASTA record {current_name!r} has inconsistent line endings/widths; "
                            "reformat it to a consistent wrapped FASTA before indexing."
                        )
                prev_len = seq_len
                prev_width = len(raw)
                n_seq_lines += 1
                total_len += seq_len
        finish()

        tmp = Path(str(self.index_path) + ".tmp")
        with tmp.open("w", encoding="utf-8", newline="\n") as out:
            for row in rows:
                out.write("\t".join(map(str, row)) + "\n")
        tmp.replace(self.index_path)
        if progress and size > step:
            print(f"FASTA index written: {len(rows)} sequences in {time.time() - started:.0f}s.", flush=True)

    def _read_index(self) -> dict[str, FaiRecord]:
        out: dict[str, FaiRecord] = {}
        with self.index_path.open("r", encoding="utf-8") as handle:
            for line in handle:
                cols = line.rstrip("\r\n").split("\t")
                if len(cols) < 5:
                    raise ValueError(f"Invalid FASTA index row in {self.index_path}: {line!r}")
                name = cols[0]
                try:
                    rec = FaiRecord(int(cols[1]), int(cols[2]), int(cols[3]), int(cols[4]))
                except ValueError as exc:
                    raise ValueError(f"Invalid numeric fields in FASTA index for {name}") from exc
                if rec.length < 0 or rec.offset < 0 or rec.line_bases <= 0 or rec.line_width < rec.line_bases:
                    raise ValueError(f"Invalid FASTA index geometry for {name}")
                out[name] = rec
        return out

    def fetch(self, contig: str, start: int, end: int) -> str:
        if contig not in self.index:
            raise KeyError(f"FASTA contig not found: {contig}")
        if start < 0 or end < start:
            raise ValueError(f"Invalid FASTA slice: {contig}:{start}-{end}")
        rec = self.index[contig]
        if end > rec.length:
            raise ValueError(f"FASTA slice exceeds contig length: {contig}:{start}-{end} (length {rec.length})")
        chunks: list[bytes] = []
        pos = start
        while pos < end:
            col = pos % rec.line_bases
            count = min(end - pos, rec.line_bases - col)
            byte_pos = rec.offset + (pos // rec.line_bases) * rec.line_width + col
            self._handle.seek(byte_pos)
            chunk = self._handle.read(count)
            if len(chunk) != count:
                raise IOError(f"Unexpected end of FASTA while reading {contig}:{start}-{end}")
            chunks.append(chunk)
            pos += count
        return b"".join(chunks).decode("ascii", errors="strict").upper()

    def close(self) -> None:
        handle = getattr(self, "_handle", None)
        if handle is not None and not handle.closed:
            handle.close()

    def __enter__(self) -> "IndexedFasta":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()
