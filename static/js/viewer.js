(() => {
  "use strict";
  const cfg = window.PTM_VIEWER || {};
  const statusEl = document.getElementById("viewer-status");
  const listEl = document.getElementById("ptm-list");
  const noteEl = document.getElementById("mapping-note");
  const accession = String(cfg.uniprot || "").trim();
  const baseAcc = accession.split("-")[0];            // isoform suffix is never used by RCSB/SIFTS
  const ptms = Array.isArray(cfg.ptms) ? cfg.ptms : [];

  const MAX_ENTRIES = 50;      // PDB entries pulled from RCSB search
  const MAX_STRUCTURES = 6;    // best-covering entries whose coordinates are downloaded
  const MOTIF_OK = 0.8;        // fraction of flanking residues that must agree with the structure
  const BFAC_HIGH = 80;        // experimental B-factor above this = poorly ordered residue

  let viewer = null, expRes = null, afRes = null, active = null, searchInfo = { ids: [], coveredBy: {} };

  // ------------------------------------------------------------------ residue tables
  const ONE_TO_THREE = {
    A: "ALA", R: "ARG", N: "ASN", D: "ASP", C: "CYS", Q: "GLN", E: "GLU", G: "GLY",
    H: "HIS", I: "ILE", L: "LEU", K: "LYS", M: "MET", F: "PHE", P: "PRO", S: "SER",
    T: "THR", W: "TRP", Y: "TYR", V: "VAL"
  };
  const THREE_TO_ONE = Object.fromEntries(Object.entries(ONE_TO_THREE).map(([o, t]) => [t, o]));
  Object.assign(THREE_TO_ONE, { SEP: "S", TPO: "T", PTR: "Y", MSE: "M", ALY: "K", M3L: "K", MLY: "K",
    MLZ: "K", M2L: "K", KCX: "K", CSO: "C", HYP: "P", SNN: "N" });
  // residue names that mean the PTM itself is present in the deposited structure
  const PTM_RESIDUES = new Set(["SEP", "TPO", "PTR", "ALY", "M3L", "MLY", "MLZ", "M2L", "KCX", "CSO", "HYP"]);
  const THREE_RE = /^(Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val)\s*[-_:]?\s*(\d+)/i;

  function setStatus(message) { if (statusEl) statusEl.textContent = message; }
  function escapeHtml(value) {
    return String(value == null ? "" : value).replace(/[&<>"']/g, c => ({
      "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
  }

  // ------------------------------------------------------------------ site label parsing
  // Accepts S473, pS473, Ser473, Ser-473, T 308 ... Returns null residue when unreadable.
  function expectedResidue(site) {
    const raw = String(site || "").trim();
    let m = raw.match(THREE_RE);
    if (m) return THREE_TO_ONE[m[1].toUpperCase()] || null;
    m = raw.match(/^(?:ph|p)?([ACDEFGHIKLMNPQRSTVWY])\s*[-_:]?\s*\d+/);
    return m ? m[1] : null;
  }
  function parseSite(ptm) {
    const raw = String(ptm.site || "").trim();
    const aa = expectedResidue(raw);
    const m = raw.match(/(\d+)\s*$/) || raw.match(/(\d+)/);
    const pos = m ? Number(m[1]) : NaN;
    if (!aa || !Number.isFinite(pos)) return { aa: null, pos: null, issue: `site label "${raw}" is not letter+number` };
    const dbPos = Number(ptm.position);
    const issue = Number.isFinite(dbPos) && dbPos !== pos ? `position column (${dbPos}) differs from label (${pos})` : null;
    return { aa, pos, issue };
  }
  const sites = ptms.map((ptm, idx) =>
    Object.assign({ idx, ptm, flank: String(ptm.flanking || "").replace(/\s+/g, "") }, parseSite(ptm)));

  // ------------------------------------------------------------------ PDB text parser
  // Reads coordinates AND the experimental evidence records that the old viewer ignored:
  // EXPDTA/RESOLUTION, REMARK 465 (missing residues), SEQADV (engineered mutations), MODRES.
  function parsePdb(text) {
    const out = { res: new Map(), chains: new Set(), missing: new Set(), seqadv: new Map(),
                  modres: new Map(), method: "", resolution: null };
    for (let line of String(text).split("\n")) {
      line = line.replace(/\r$/, "");
      const rec = line.slice(0, 6);
      if (rec === "ENDMDL") break;                       // NMR: first model only
      if (rec === "ATOM  " || rec === "HETATM") {
        const resn = line.slice(17, 20).trim();
        if (rec === "HETATM" && !THREE_TO_ONE[resn]) continue;
        const alt = line[16];
        if (alt !== " " && alt !== "A" && alt !== "1") continue;
        if (line[26] && line[26].trim()) continue;       // insertion code: numbering ambiguous
        const chain = (line[21] || "").trim();
        const resi = parseInt(line.slice(22, 26), 10);
        if (Number.isNaN(resi)) continue;
        const key = chain + "|" + resi;
        let r = out.res.get(key);
        if (!r) { r = { resn, aa: THREE_TO_ONE[resn] || null, b: null }; out.res.set(key, r); }
        if (line.slice(12, 16).trim() === "CA") r.b = parseFloat(line.slice(60, 66));
        out.chains.add(chain);
      } else if (rec === "SEQADV") {
        out.seqadv.set((line[16] || "").trim() + "|" + parseInt(line.slice(18, 22), 10), line.slice(49).trim());
      } else if (rec === "MODRES") {
        out.modres.set((line[16] || "").trim() + "|" + parseInt(line.slice(18, 22), 10),
          `${line.slice(12, 15).trim()} (parent ${line.slice(24, 27).trim()}) ${line.slice(29).trim()}`);
      } else if (rec === "EXPDTA") {
        out.method = line.slice(10).trim();
      } else if (rec === "REMARK") {
        let m = line.match(/^REMARK\s+2 RESOLUTION\.\s+([\d.]+)/);
        if (m) out.resolution = parseFloat(m[1]);
        m = line.match(/^REMARK 465\s+([A-Z0-9]{3})\s+(\S)\s+(-?\d+)\s*$/);
        if (m) out.missing.add(m[2] + "|" + parseInt(m[3], 10));
      }
    }
    return out;
  }

  // ------------------------------------------------------------------ SIFTS (UniProt -> PDB numbering)
  // The old code used only the FIRST block that covered a site, so a residue that is disordered in
  // chain A but solved in chain B was reported as "unmapped". Every covering block is now returned.
  function siftsCandidates(blocks, pos) {
    const out = [];
    for (const b of blocks || []) {
      const us = Number(b.unp_start), ue = Number(b.unp_end);
      const start = b.start || {}, end = b.end || {};
      if (!(us <= pos && pos <= ue) || start.author_residue_number == null) continue;
      const resi = Number(start.author_residue_number) + (pos - us);
      const contiguous = end.author_residue_number != null &&
        Number(end.author_residue_number) - Number(start.author_residue_number) === ue - us &&
        !start.author_insertion_code && !end.author_insertion_code;
      out.push({ chain: String(b.chain_id), resi, contiguous });
    }
    return out;
  }

  // ------------------------------------------------------------------ the correctness checks
  // Compare the record's flanking sequence with the residues the structure really has around the
  // mapped position. One matching letter proves little (about 1 residue in 14 is Ser); a matching
  // 15-mer proves the numbering is right.
  function motifScore(parsed, chain, resi, s) {
    const f = s.flank;
    if (!f || f.length < 5) return null;
    const mid = (f.length - 1) / 2;
    let best = null;
    for (let c = 0; c < f.length; c++) {                 // the site residue should sit at (or within 2 of) the centre
      if (f[c].toUpperCase() !== s.aa || Math.abs(c - mid) > 2) continue;
      let cmp = 0, hit = 0;
      for (let i = 0; i < f.length; i++) {
        if (i === c) continue;
        const q = f[i].toUpperCase();
        if (!ONE_TO_THREE[q]) continue;
        const r = parsed.res.get(chain + "|" + (resi + i - c));
        if (!r || !r.aa) continue;                       // unresolved neighbours are skipped, not counted against
        cmp++; if (r.aa === q) hit++;
      }
      if (!best || hit / Math.max(cmp, 1) > best.hit / Math.max(best.cmp, 1)) best = { cmp, hit };
    }
    if (!best) return { state: "unrecognised" };
    if (best.cmp < 3) return { state: "too_few", cmp: best.cmp, hit: best.hit };
    return { state: best.hit / best.cmp >= MOTIF_OK ? "match" : "mismatch", cmp: best.cmp, hit: best.hit };
  }

  function evaluateResidue(parsed, chain, resi, s, kind, contiguous) {
    const key = chain + "|" + resi;
    const r = parsed.res.get(key);
    const rec = { chain, resi, kind, status: "", notes: [] };
    if (!r) {
      if (kind === "alphafold") {
        rec.status = "unresolved";
        rec.detail = `position ${resi} is outside the model (wrong isoform/species numbering, or label error)`;
      } else if (parsed.missing.has(key)) {
        rec.status = "disordered";
        rec.detail = "listed as missing in REMARK 465 (no electron density: disordered region)";
      } else {
        rec.status = "unresolved";
        rec.detail = "no coordinates at this residue (outside the crystallised construct or not modelled)";
      }
      return rec;
    }
    rec.resn = r.resn; rec.b = r.b;
    if (r.aa !== s.aa) {
      rec.status = "mismatch";
      rec.detail = `structure has ${r.resn} where the record says ${s.aa}: wrong position, isoform, species or a mutant construct`;
      return rec;
    }
    rec.modified = PTM_RESIDUES.has(r.resn);
    rec.motif = motifScore(parsed, chain, resi, s);
    rec.status = "ok";
    const adv = parsed.seqadv.get(key);
    if (adv) {
      rec.notes.push(`SEQADV: ${adv}`);
      if (/MUTATION|CONFLICT/i.test(adv)) { rec.status = "suspect"; rec.detail = "residue is an engineered mutation/conflict, not the natural residue"; }
    }
    const mr = parsed.modres.get(key);
    if (mr) rec.notes.push(`MODRES: ${mr}`);
    if (rec.status === "ok" && rec.motif && rec.motif.state === "mismatch") {
      rec.status = "suspect";
      rec.detail = `letter matches but only ${rec.motif.hit}/${rec.motif.cmp} flanking residues agree` +
        (contiguous === false ? " (numbering gap/insertion code in this SIFTS block, offset may be off)" : ": mapping is probably wrong");
    }
    return rec;
  }

  function analyseExperimental(pdb, text, blocks) {
    const parsed = parsePdb(text);
    const perChain = new Map();
    for (const s of sites) {
      if (!s.aa) continue;
      for (const c of siftsCandidates(blocks, s.pos)) {
        const rec = evaluateResidue(parsed, c.chain, c.resi, s, "experimental", c.contiguous);
        if (!perChain.has(c.chain)) perChain.set(c.chain, new Map());
        const m = perChain.get(c.chain), old = m.get(s.idx);
        if (!old || (rec.status === "ok" && old.status !== "ok")) m.set(s.idx, rec);
      }
    }
    let bestChain = null, bestN = -1, bestSize = -1;
    for (const [chain, m] of perChain) {
      const n = [...m.values()].filter(r => r.status === "ok").length;
      if (n > bestN || (n === bestN && m.size > bestSize)) { bestChain = chain; bestN = n; bestSize = m.size; }
    }
    return { kind: "experimental", id: pdb.toUpperCase(), chain: bestChain, text, parsed,
             recs: perChain.get(bestChain) || new Map(), okCount: Math.max(bestN, 0), perChain };
  }

  function analyseAlphaFold(id, text) {
    const parsed = parsePdb(text);
    let best = null;
    for (const chain of parsed.chains) {
      const recs = new Map();
      for (const s of sites) if (s.aa) recs.set(s.idx, evaluateResidue(parsed, chain, s.pos, s, "alphafold", true));
      const n = [...recs.values()].filter(r => r.status === "ok").length;
      if (!best || n > best.okCount) best = { kind: "alphafold", id, chain, text, parsed, recs, okCount: n };
    }
    return best;
  }

  // ------------------------------------------------------------------ network
  async function jget(url, opts) {
    const r = await fetch(url, opts);
    if (!r.ok) throw new Error(`HTTP ${r.status}`);
    return r.json();
  }
  async function findExperimentalIds() {
    const term = (attribute, value) => ({ type: "terminal", service: "text",
      parameters: { attribute, operator: "exact_match", value } });
    const query = {
      query: { type: "group", logical_operator: "and", nodes: [
        term("rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name", "UniProt"),
        term("rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession", baseAcc) ] },
      return_type: "entry",
      request_options: { paginate: { start: 0, rows: MAX_ENTRIES }, results_content_type: ["experimental"] }
    };
    const r = await fetch("https://search.rcsb.org/rcsbsearch/v2/query", {
      method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(query) });
    if (r.status === 204) return [];                     // RCSB: no hits
    if (!r.ok) throw new Error(`RCSB search returned HTTP ${r.status}`);
    const data = await r.json();
    return (data.result_set || []).map(row => row.identifier).filter(Boolean);
  }
  async function getSiftsBlocks(pdb) {
    try {
      const d = await jget(`https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/${encodeURIComponent(pdb.toLowerCase())}`);
      const e = d[pdb.toLowerCase()];
      const set = e && (e.UniProt || e.uniprot);
      const blocks = set && set[baseAcc] && set[baseAcc].mappings;
      return Array.isArray(blocks) ? blocks : null;
    } catch (_) { return null; }
  }
  async function getPdbText(pdb) {
    const r = await fetch(`https://files.rcsb.org/download/${encodeURIComponent(pdb.toUpperCase())}.pdb`);
    if (!r.ok) throw new Error(`PDB ${pdb} returned HTTP ${r.status}`);
    const t = await r.text();
    if (!/^ATOM/m.test(t)) throw new Error(`PDB ${pdb} did not contain coordinate records`);
    return t;
  }
  async function pool(items, n, fn) {
    const out = new Array(items.length); let i = 0;
    await Promise.all(Array.from({ length: Math.min(n, items.length) }, async () => {
      while (i < items.length) { const k = i++; out[k] = await fn(items[k]); }
    }));
    return out;
  }
  async function fetchAlphaFoldModel() {
    const tried = [];
    try {                                                // ask AlphaFold DB for the CURRENT file instead of guessing a version
      const d = await jget(`https://alphafold.ebi.ac.uk/api/prediction/${encodeURIComponent(baseAcc)}`);
      if (Array.isArray(d) && d[0] && d[0].pdbUrl) tried.push({ url: d[0].pdbUrl, id: d[0].entryId || `AF-${baseAcc}-F1` });
    } catch (_) { /* fall through to fixed versions */ }
    const base = `https://alphafold.ebi.ac.uk/files/AF-${encodeURIComponent(baseAcc)}-F1-model_`;
    for (const v of ["v6.pdb", "v4.pdb", "v3.pdb"]) tried.push({ url: base + v, id: `AF-${baseAcc}-F1` });
    for (const c of tried) {
      try {
        const r = await fetch(c.url);
        if (!r.ok) continue;
        const text = await r.text();
        if (/^ATOM\s/m.test(text) && text.length > 1000) return { id: c.id, text };
      } catch (_) { /* next */ }
    }
    return null;
  }

  // ------------------------------------------------------------------ viewer helpers
  const plddtLabel = b => b >= 90 ? "very high" : b >= 70 ? "confident" : b >= 50 ? "low" : "very low";
  const plddtColor = a => a.b >= 90 ? "#0053D6" : a.b >= 70 ? "#65CBF3" : a.b >= 50 ? "#FFDB13" : "#FF7D45";
  function baseStyle(res) {
    return res.kind === "alphafold" ? { cartoon: { colorfunc: plddtColor } } : { cartoon: { color: "spectrum" } };
  }
  function sel(rec) { const s = { resi: rec.resi }; if (rec.chain) s.chain = rec.chain; return s; }
  function atomsFor(selection) {
    try { return viewer.selectedAtoms(selection) || []; } catch (_) { return []; }
  }
  // Final live cross-check inside the 3Dmol model itself (independent of the text parser above).
  function verifiedResidue(mapping, ptm) {
    const expected = expectedResidue(ptm.site);
    if (!expected || !viewer) return false;
    const names = new Set(atomsFor(sel(mapping)).map(a => String(a.resn || "").toUpperCase()).filter(Boolean));
    return [...names].some(n => THREE_TO_ONE[n] === expected);
  }
  function usable(rec) { return rec && (rec.status === "ok" || rec.status === "suspect"); }
  function showModel(res) {
    active = res;
    const el = document.getElementById("viewer3d");
    viewer = window.PTMViewer.mount(el);                  // always the viewer of the element that is on the page
    const n = window.PTMViewer.load(viewer, res.text);    // throws if the coordinate file parsed to 0 atoms
    highlightAll();
    const fin = window.PTMViewer.finish(viewer, el, chainSel());     // fit, resize, draw and confirm a canvas really exists
    if (!fin.ok) throw new Error(fin.reason);
    return n;
  }
   function chainSel() { return active && active.chain ? { chain: active.chain } : {}; }
  function highlightAll() {
    viewer.setStyle({}, {});
    viewer.setStyle(chainSel(), baseStyle(active));
    for (const [idx, rec] of active.recs) {
      if (!usable(rec)) continue;
      const colour = rec.status === "ok" ? "yellow" : "orange";
      viewer.setStyle(sel(rec), { stick: { radius: 0.24, color: colour }, cartoon: { color: colour } });
    }
  }
  function headline(res) {
    if (res.kind === "alphafold")
      return `AlphaFold PREDICTED model ${res.id} · chain ${res.chain} · ${res.okCount}/${ptms.length} PTM sites residue-verified`;
    const p = res.parsed;
    const meth = p.method ? ` · ${p.method}${p.resolution ? " " + p.resolution + " Å" : ""}` : "";
    return `Experimental structure ${res.id}${meth} · chain ${res.chain} · ${res.okCount}/${ptms.length} PTM sites verified`;
  }
  function noteFor(res) {
    return res.kind === "alphafold"
      ? "Predicted AlphaFold coordinates—not an experimentally solved structure. Cartoon is coloured by pLDDT (blue = reliable, orange = unreliable). PTM markers appear only when the position exists and its residue identity matches the site label. Missing or mismatched residues remain unmarked."
      : "Experimental coordinates. A site is highlighted yellow only if SIFTS maps it, the residue letter agrees AND the flanking sequence agrees. Orange = letter agrees but the mapping looks wrong (see report). Unresolved, disordered or mismatched residues stay unhighlighted.";
  }

  // ------------------------------------------------------------------ report (trace table)
  function describe(rec, kind) {
    if (!rec) return { icon: "–", cls: "na", text: "not mapped" };
    const loc = kind === "alphafold" ? `${rec.resn || ""} ${rec.resi}` : `${rec.resn || ""} ${rec.chain}:${rec.resi}`;
    if (rec.status === "ok") {
      const bits = [loc.trim()];
      let cls = "ok";
      if (kind === "alphafold") {
        bits.push(`pLDDT ${Math.round(rec.b)} (${plddtLabel(rec.b)})`);
        if (rec.b < 70) cls = "warn";
      } else {
        if (rec.modified) bits.push(`${rec.resn} present = PTM seen in the structure`);
        if (rec.b != null && rec.b > BFAC_HIGH) { bits.push(`B=${Math.round(rec.b)} poorly ordered`); cls = "warn"; }
      }
      if (rec.motif && rec.motif.state === "match") bits.push(`motif ${rec.motif.hit}/${rec.motif.cmp}`);
      else if (rec.motif && rec.motif.state === "too_few") bits.push("motif not checkable (flank mostly unresolved)");
      else if (rec.motif && rec.motif.state === "unrecognised") bits.push("motif not checked (site not centred in flanking string)");
      else if (!rec.motif) bits.push("no flanking sequence to check");
      return { icon: "✔", cls, text: bits.join(" · ") };
    }
    if (rec.status === "suspect") return { icon: "⚠", cls: "warn", text: `${loc.trim()} · ${rec.detail}` };
    if (rec.status === "mismatch") return { icon: "✘", cls: "bad", text: rec.detail };
    return { icon: "○", cls: "na", text: rec.detail || rec.status };
  }
  function expCell(s) {
    if (!s.aa) return { icon: "✘", cls: "bad", text: s.issue };
    const rec = expRes && expRes.recs.get(s.idx);
    if (rec) {
      const d = describe(rec, "experimental");
      const alt = (expRes.alsoOk[s.idx] || []).filter(x => x !== expRes.id);
      if (rec.status !== "ok" && alt.length) d.text += ` · but verified in ${alt.join(", ")}`;
      return d;
    }
    if (!searchInfo.ids.length) return { icon: "–", cls: "na", text: "no experimental PDB entry exists for this protein" };
    const cov = searchInfo.coveredBy[s.idx] || [];
    if (!cov.length) return { icon: "–", cls: "na", text: `outside every solved segment (${searchInfo.ids.length} entries checked): no experimental structure for this region` };
    const alt = expRes ? (expRes.alsoOk[s.idx] || []) : [];
    if (alt.length) return { icon: "○", cls: "na", text: `not in ${expRes.id}, but verified in ${alt.join(", ")}` };
    return { icon: "○", cls: "na", text: `covered by ${cov.slice(0, 4).join(", ")} but not verified (see their chains)` };
  }
  function afCell(s) {
    if (!s.aa) return { icon: "✘", cls: "bad", text: s.issue };
    if (!afRes) return { icon: "–", cls: "na", text: "no AlphaFold model available" };
    return describe(afRes.recs.get(s.idx), "alphafold");
  }
  const COLORS = { ok: "#0b7a4b", warn: "#a15c00", bad: "#b3261e", na: "#5f7377" };
  function renderReport() {
    let host = document.getElementById("verify-report");
    if (!host) {
      host = document.createElement("div");
      host.id = "verify-report";
      host.style.cssText = "margin:0 0 10px;font-size:11px";
      noteEl.insertAdjacentElement("afterend", host);
    }
    if (!sites.length) { host.innerHTML = ""; return; }
    const rows = sites.map(s => {
      const e = expCell(s), a = afCell(s);
      const td = c => `<td style="padding:5px 6px;border-top:1px solid #e3ecec;color:${COLORS[c.cls]}">${c.icon} ${escapeHtml(c.text)}</td>`;
      const warn = s.issue && s.aa ? `<br><span style="color:#a15c00">⚠ ${escapeHtml(s.issue)}</span>` : "";
      return `<tr><td style="padding:5px 6px;border-top:1px solid #e3ecec;font-weight:800;white-space:nowrap">${escapeHtml(s.ptm.site || "—")}${warn}</td>${td(e)}${td(a)}</tr>`;
    }).join("");
    host.innerHTML = `<details open><summary style="cursor:pointer;font-weight:800;margin-bottom:6px">Mapping &amp; correctness report</summary>
      <div style="overflow-x:auto"><table style="min-width:0;width:100%;border-collapse:collapse"><thead><tr>
      <th style="position:static">Site</th><th style="position:static">Experimental PDB</th><th style="position:static">AlphaFold (predicted)</th></tr></thead>
      <tbody>${rows}</tbody></table></div>
      <p style="color:#5f7377;margin:6px 0 0">✔ residue letter and flanking motif agree · ⚠ check · ✘ contradicts the record · ○/– no usable coordinates. Motif = the record's flanking residues compared with the residues actually in the structure.</p></details>`;
    renderSwitch(host);
  }
  function renderSwitch(host) {
    if (!(expRes && afRes)) return;
    const old = document.getElementById("source-switch");
    if (old) old.remove();
    const box = document.createElement("div");
    box.id = "source-switch";
    box.style.cssText = "display:flex;gap:6px;margin-bottom:8px";
    for (const [label, res] of [[`Experimental ${expRes.id}`, expRes], [`AlphaFold`, afRes]]) {
      const b = document.createElement("button");
      b.type = "button"; b.className = "btn btn-secondary btn-sm"; b.textContent = label;
      if (active === res) b.style.borderColor = "#0b7878";
      b.addEventListener("click", () => { switchTo(res); });
      box.appendChild(b);
    }
    host.insertAdjacentElement("beforebegin", box);
  }
  function switchTo(res) {
    let nAtoms;
    try { nAtoms = showModel(res); }
    catch (e) {
      setStatus(`The ${res.kind === "alphafold" ? "AlphaFold" : "experimental"} structure ${res.id} was retrieved but could not be displayed: ${e.message}. No highlight was drawn.`);
      noteEl.textContent = "Nothing is shown rather than a misleading blank scene.";
      renderCards(); renderReport();
      return;
    }
    noteEl.textContent = noteFor(res);
    setStatus(headline(res) + ` · ${nAtoms} atoms drawn`);
    if (res.kind === "alphafold" && statusEl) statusEl.insertAdjacentHTML("beforeend", " · " + window.PTMViewer.PLDDT_LEGEND);
    renderCards(); renderReport();
  }
  function addResetButton() {
    const card = document.querySelector(".viewer-card");
    if (!card || document.getElementById("viewer-reset")) return;
    const b = document.createElement("button");
    b.type = "button"; b.id = "viewer-reset"; b.className = "btn btn-secondary btn-sm";
    b.style.cssText = "position:absolute;right:10px;top:10px;z-index:2";
    b.textContent = "Reset view";
    b.addEventListener("click", () => { if (viewer) { viewer.zoomTo(chainSel()); viewer.render(); } });
    card.appendChild(b);
  }

  // ------------------------------------------------------------------ cards + selection
  function renderCards() {
    listEl.innerHTML = "";
    if (!ptms.length) {
      const empty = document.createElement("p");
      empty.className = "empty";
      empty.textContent = "No real PTM records are loaded for this protein. The structure viewer can still show an available protein model; no PTM markers are invented.";
      listEl.appendChild(empty);
      return;
    }
    sites.forEach(s => {
      const ptm = s.ptm, rec = active && active.recs.get(s.idx);
      const mapped = usable(rec);
      const el = document.createElement("button");
      el.type = "button";
      el.className = "ptm-card" + (mapped ? " mapped" : "");
      const state = mapped ? describe(rec, active.kind) : { icon: "○", text: rec && rec.detail ? rec.detail : "not mapped in this structure", cls: "na" };
      el.innerHTML = `<span class="ptm-site">${escapeHtml(ptm.site || "—")}</span>
        <span><b>PTM / evidence</b><br>${escapeHtml(ptm.ptm_type || "—")} · ${escapeHtml(ptm.evidence || "No evidence field")}</span>
        <span><b>Context / confidence</b><br>${escapeHtml(ptm.disease_context || "—")} · ${escapeHtml(ptm.confidence || "—")}</span>
        <span style="color:${COLORS[state.cls]}">${state.icon} ${escapeHtml(state.text)}</span>`;
      el.addEventListener("click", e => { e.preventDefault(); select(s); });
      listEl.appendChild(el);
    });
  }
  function select(s) {
    const rec = active && active.recs.get(s.idx);
    if (!usable(rec) || !verifiedResidue(rec, s.ptm)) {
      setStatus(`${s.ptm.site || "This PTM site"} is not reliably mapped in the loaded structure: ${rec && rec.detail ? rec.detail : "no coordinates"}`);
      return;
    }
    const selection = sel(rec);
    highlightAll();
    viewer.setStyle(selection, { stick: { radius: 0.35 }, sphere: { radius: 0.55, color: "red" }, cartoon: { color: "red" } });
    viewer.zoomTo(selection, 500);
    viewer.render();
    setStatus(`${s.ptm.site} · ${s.ptm.ptm_type || "PTM"} · ${active.kind === "alphafold" ? "AlphaFold prediction" : "Experimental PDB " + active.id} · residue ${rec.resi}`);
  }

  // ------------------------------------------------------------------ main flow
  async function loadExperimental() {
    let ids = [];
    try { ids = await findExperimentalIds(); } catch (_) { return null; }
    searchInfo.ids = ids;
    if (!ids.length) return null;
    const blocksList = await pool(ids, 8, getSiftsBlocks);           // all SIFTS maps, in parallel
    const cand = [];
    ids.forEach((pdb, i) => {
      const blocks = blocksList[i];
      if (!blocks) return;
      let n = 0;
      for (const s of sites) {
        if (s.aa && siftsCandidates(blocks, s.pos).length) {
          n++; (searchInfo.coveredBy[s.idx] = searchInfo.coveredBy[s.idx] || []).push(pdb.toUpperCase());
        }
      }
      if (n) cand.push({ pdb, blocks, n });
    });
    cand.sort((a, b) => b.n - a.n);                                  // rank by SIFTS coverage, not RCSB order
    let best = null; const alsoOk = {};
    for (const c of cand.slice(0, MAX_STRUCTURES)) {
      let text;
      try { text = await getPdbText(c.pdb); } catch (_) { continue; }
      const res = analyseExperimental(c.pdb, text, c.blocks);
      for (const [idx, rec] of res.recs) if (rec.status === "ok") (alsoOk[idx] = alsoOk[idx] || []).push(res.id);
      if (!best || res.okCount > best.okCount) best = res;
      if (best.okCount === sites.filter(s => s.aa).length) break;
    }
    if (best) best.alsoOk = alsoOk;
    return best;
  }

  async function load() {
    if (!accession || !window.$3Dmol || !window.PTMViewer || !statusEl || !listEl || !noteEl) {
      setStatus("3D viewer could not initialize. Check the browser connection and the 3Dmol.js resource.");
      if (noteEl) noteEl.textContent = "No structural result is being substituted.";
      return;
    }
    try { viewer = window.PTMViewer.mount(document.getElementById("viewer3d")); }
    catch (e) { setStatus(`3D viewer could not initialize: ${e.message}`); return; }
    addResetButton();
    renderCards();
    setStatus("Checking available experimental structure mappings…");
    try { expRes = await loadExperimental(); } catch (_) { expRes = null; }
    if (expRes && !expRes.alsoOk) expRes.alsoOk = {};
    setStatus("Checking AlphaFold…");
    try {
      const model = await fetchAlphaFoldModel();
      afRes = model ? analyseAlphaFold(model.id, model.text) : null;
    } catch (_) { afRes = null; }

    // Prefer a structure that verifies at least one PTM site. If none does, still show whichever real structure
    // exists (experimental first) with no highlights, rather than claiming there is no structure.
    const chosen = expRes && expRes.okCount > 0 ? expRes : (afRes && afRes.okCount > 0 ? afRes : (expRes || afRes));
    if (!chosen) {
      active = null;
      renderCards(); renderReport();
      setStatus("Structure unavailable from the configured public structure sources.");
      noteEl.innerHTML = `No coordinate model could be loaded. Check internet access or open <a href="https://alphafold.ebi.ac.uk/entry/AF-${escapeHtml(baseAcc)}-F1" target="_blank" rel="noopener">AlphaFold DB</a> and <a href="https://www.rcsb.org/" target="_blank" rel="noopener">RCSB PDB</a>. No coordinates or PTM mappings were fabricated.`;
      return;
    }
    switchTo(chosen);
    if (chosen === afRes && expRes && statusEl) {
      statusEl.insertAdjacentText("beforeend", " · (no experimental structure verified any site)");
    }
    if (ptms.length && chosen.okCount === 0 && statusEl) {
      statusEl.insertAdjacentText("beforeend", " · NO PTM site could be residue-verified in this structure, so none is highlighted");
    }
  }
  if (cfg.__test) window.__ptmViewerTest = { parseSite, parsePdb, siftsCandidates, evaluateResidue, motifScore,
    analyseExperimental, analyseAlphaFold, sites, describe };
  load();
})();