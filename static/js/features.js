(() => {
const F = window.FEATURE || {};
const $ = (s, r = document) => r.querySelector(s);
const esc = s => String(s == null ? "" : s).replace(/[&<>"']/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
const TYPE_COLOR = {Phosphorylation:"#18a6a6",Acetylation:"#e8871e",Ubiquitination:"#8a5cf6",Methylation:"#3aa8e6",SUMOylation:"#e0507a",Glycosylation:"#5cb85c"};
const tcol = t => TYPE_COLOR[t] || "#657b80";
const out = $("#result");
const up = () => ($("#protein-select") || {}).value || F.selected;
const val = (id, d) => { const e = $("#" + id); return e && e.value !== "" ? e.value : d; };
const pm = p => `<a href="https://pubmed.ncbi.nlm.nih.gov/${esc(p)}/" target="_blank" rel="noopener">${esc(p)}</a>`;
const chips = (a, cls = "") => (a && a.length ? a.map(x => `<span class="chip ${cls}">${esc(x)}</span>`).join("") : '<span class="muted">—</span>');
const cards = items => `<div class="kpi-grid">${items.map(([v, l, n]) => `<div class="kpi"><div class="kpi-v">${esc(v)}</div><div class="kpi-l">${esc(l)}</div>${n ? `<div class="kpi-n">${esc(n)}</div>` : ""}</div>`).join("")}</div>`;
const method = m => m ? `<details class="method"><summary>Method &amp; assumptions</summary><p>${esc(m)}</p></details>` : "";
const table = (heads, rows, empty = "Nothing to show.") => rows.length ?
  `<div class="table-wrap"><table><thead><tr>${heads.map(h => `<th>${h}</th>`).join("")}</tr></thead><tbody>${rows.map(r => `<tr>${r.map(c => `<td>${c}</td>`).join("")}</tr>`).join("")}</tbody></table></div>` :
  `<div class="empty">${esc(empty)}</div>`;
const pfmt = p => p == null ? "—" : (p < 0.001 ? p.toExponential(1) : p.toFixed(3));
const busy = m => out.innerHTML = `<div class="busy"><span class="spin"></span>${esc(m || "Working…")}</div>`;
async function api(path, timeoutMs = 120000) {
  const ctrl = new AbortController();
  const timer = setTimeout(() => ctrl.abort(), timeoutMs);
  try {
    let r;
    try {
      r = await fetch(path, { signal: ctrl.signal });
    } catch (e) {
      if (e.name === "AbortError")
        throw new Error(`The server did not answer within ${timeoutMs / 1000} s. Ensembl is probably slow; check the Flask terminal and try again.`);
      throw new Error("Could not reach the Flask server. Is it still running?");
    }
    if (!r.ok) throw new Error("Server error " + r.status);
    const j = await r.json();
    if (j.empty) { out.innerHTML = `<div class="notice warn">${esc(j.message)}</div>`; return null; }
    if (j.available === false) { out.innerHTML = `<div class="notice warn"><b>Not available.</b> ${esc(j.reason || "")}${(j.notes||[]).map(n=>"<br>"+esc(n)).join("")}<br><span class="muted">Nothing was estimated or substituted.</span></div>` + retrievalHtml(j.retrieval); return null; }
    if (j.error) { out.innerHTML = `<div class="notice warn">${esc(j.error)}</div>`; return null; }
    return j;
  } finally {
    clearTimeout(timer);
  }
}

// cBioPortal retrieval audit: which study/profile/sample list was really used and whether it is complete.
function retrievalHtml(r) {
  if (!r) return "";
  const bits = [
    r.study_id ? `study <b>${esc(r.study_id)}</b>` : "",
    r.profile_id ? `mutation profile <code>${esc(r.profile_id)}</code>` : "mutation profile: <b>not found</b>",
    r.sample_list_id ? `sample list <code>${esc(r.sample_list_id)}</code> (${r.n_sequenced ?? "?"} sequenced samples)` : "",
    r.pages ? `${r.pages} page(s), ${r.complete ? "complete" : "<b>INCOMPLETE</b>"}` : "",
    r.n_malformed_records ? `${r.n_malformed_records} malformed record(s) ignored` : ""
  ].filter(Boolean).join(" · ");
  return `<div class="notice"><b>cBioPortal retrieval</b> (status: ${esc(r.status)}): ${bits}` +
    (r.warnings || []).map(w => `<br>⚠ ${esc(w)}`).join("") + `</div>`;
}

// ------------------------------------------------------- shared SVG track
function track(L, sites, extra = "", h = 90) {
  const W = 900, x = p => 20 + (p / L) * (W - 40);
  let s = `<svg viewBox="0 0 ${W} ${h}" width="100%" class="track"><line x1="20" x2="${W-20}" y1="${h-22}" y2="${h-22}" stroke="#9bb" stroke-width="2"/>`;
  const step = [10,20,25,50,100,200,250,500,1000].find(v => L / v <= 10) || 1000;
  for (let p = 0; p <= L; p += step) s += `<text x="${x(p)}" y="${h-4}" font-size="9" fill="#789" text-anchor="middle">${p}</text>`;
  s += extra;
  for (const u of sites) s += `<line x1="${x(u.pos)}" x2="${x(u.pos)}" y1="${h-30}" y2="${h-22}" stroke="${tcol(u.ptm_type)}" stroke-width="2"><title>${esc(u.site)} ${esc(u.ptm_type)}</title></line>`;
  return s + "</svg>";
}
const legend = () => `<div class="legend">${Object.entries(TYPE_COLOR).map(([k, c]) => `<span><i style="background:${c}"></i>${k}</span>`).join("")}</div>`;

// ---------------------------------------------------------------- HOTSPOTS
async function hotspots() {
  busy("Scanning windows…");
  const j = await api(`/api/f/hotspots/${up()}?window=${val("p-window",15)}&min_sites=${val("p-min",3)}`); if (!j) return;
  const L = j.length, sig = j.regions.filter(r => r.adj_p < 0.05);
  const maxc = Math.max(1, ...j.profile.map(p => p[1])), W = 900, H = 90;
  const x = p => 20 + (p / L) * (W - 40);
  let area = "", regs = "";
  if (j.profile.length) area = `<polyline fill="rgba(24,166,166,.18)" stroke="#18a6a6" stroke-width="1" points="${x(1)},${H-22} ${j.profile.map(p => `${x(p[0])},${H-22-(p[1]/maxc)*52}`).join(" ")} ${x(L)},${H-22}"/>`;
  for (const r of j.regions) regs += `<rect x="${x(r.start)}" y="6" width="${Math.max(3, x(r.end) - x(r.start))}" height="8" rx="3" fill="${r.adj_p < 0.05 ? "#d64545" : "#e0a020"}"><title>${r.start}–${r.end}: ${r.n_residues} residues, adj.p=${pfmt(r.adj_p)}</title></rect>`;
  const sites = []; const seen = new Set();
  j.regions.forEach(r => r.sites.forEach(s => { const k = s.ptm_type + s.pos; if (!seen.has(k)) { seen.add(k); sites.push(s); } }));
  out.innerHTML = cards([[j.n_positions, "Modified residues", `${j.n_units} residue×PTM-type units`], [j.regions.length, "Candidate regions", `window ${j.params.window}, ≥${j.params.min_sites} residues`], [sig.length, "Significant (adj. p<0.05)"], [L, "Protein length (aa)"]]) +
    `<div class="panel"><h3>Density profile — modified residues per ${j.params.window}-residue window</h3>${track(L, [], area + regs)}<div class="legend"><span><i style="background:#d64545"></i>significant region</span><span><i style="background:#e0a020"></i>candidate (not significant)</span></div></div>` +
    `<h3>Regions</h3>` + table(["Region", "Residues", "Best window", "Expected", "Fold", "p", "adj. p", "PTM types", "Kinases", "Disease context", "PMIDs"],
      j.regions.map(r => [`<b>${r.start}–${r.end}</b>`, r.n_residues, `${r.best_window[0]}–${r.best_window[1]} (${r.best_window_count})`, r.expected_in_window, r.fold_enrichment ?? "—", pfmt(r.p_value), `<b>${pfmt(r.adj_p)}</b>`,
        Object.entries(r.ptm_types).map(([t, n]) => `<span class="chip" style="border-color:${tcol(t)}">${n}× ${esc(t)}</span>`).join(""), chips(r.kinases), chips(r.diseases), r.pmids.map(pm).join(", ") || "—"]),
      "No region reaches the minimum number of modified residues in one window for this protein.") +
    `<p class="caveat">A hotspot is a <b>candidate PTM-dense region</b>. Curated databases over-represent well-studied regions, so enrichment is a screening signal, not proof of regulatory importance.</p>` + method(j.method);
}
// --------------------------------------------------------------- CROSSTALK
// --------------------------------------------------------------- CROSSTALK
const CLOSE_DIST = 10;   // pairs this close (in residues) are flagged as "close"

async function crosstalk() {
  busy("Analysing PTM pairs…");

  const distance = Number(val("p-dist", 10));
  if (!Number.isFinite(distance) || distance < 0) {
    out.innerHTML = '<div class="notice warn">Enter a valid non-negative residue distance.</div>';
    return;
  }

  const sameType = $("#p-same").checked ? 1 : 0;
  const j = await api(`/api/f/crosstalk/${encodeURIComponent(up())}?max_dist=${distance}&same_type=${sameType}`);
  if (!j) return;

  const pairs = Array.isArray(j.pairs) ? j.pairs : [];
  const nCross = pairs.filter(p => p.a.ptm_type !== p.b.ptm_type).length;
  const nClose = pairs.filter(p => Number(p.distance) <= CLOSE_DIST).length;
  const L = Math.max(1, Number(j.protein?.length) || 1);
  const W = 900, H = 240;
  const x = pos => 20 + (Number(pos) / L) * (W - 40);
  const baseline = H - 22;
  const skey = s => `${s.pos}|${s.ptm_type}`;

  // ---- arcs: height grows with the span
  let arcs = "";
  pairs.forEach((p, i) => {
    const a = x(p.a.pos), b = x(p.b.pos);
    const color = p.class === "same-residue competition" ? "#d64545" : tcol(p.a.ptm_type);
    if (Number(p.a.pos) === Number(p.b.pos)) {
      arcs += `<circle cx="${a}" cy="${baseline - 9 - (i % 3) * 5}" r="5" fill="none" stroke="${color}" stroke-width="2">` +
        `<title>${esc(p.a.site)} ${esc(p.a.ptm_type)} ↔ ${esc(p.b.site)} ${esc(p.b.ptm_type)} (same residue)</title></circle>`;
      return;
    }
    const w = Number(p.distance) <= CLOSE_DIST ? 2.6 : 1.6;
    const offset = Math.min(baseline - 14, 14 + Math.abs(b - a) * 1.1);
    arcs += `<path class="xt-arc" data-i="${i}" data-w="${w}" d="M ${a} ${baseline} Q ${(a + b) / 2} ${baseline - offset} ${b} ${baseline}" ` +
      `fill="none" stroke="${color}" stroke-opacity="0.55" stroke-width="${w}" style="cursor:pointer">` +
      `<title>${esc(p.a.site)} ${esc(p.a.ptm_type)} ↔ ${esc(p.b.site)} ${esc(p.b.ptm_type)} (${esc(p.distance)} residues)</title></path>`;
  });

  // ---- one dot + one hidden label per modified residue
  const sites = [], seen = new Set();
  pairs.forEach(p => [p.a, p.b].forEach(s => {
    if (!seen.has(skey(s) + s.site)) { seen.add(skey(s) + s.site); sites.push({ pos: Number(s.pos), ptm_type: s.ptm_type, site: s.site }); }
  }));
  const dots = sites.map(s => `<circle class="xt-dot" data-k="${esc(skey(s))}" cx="${x(s.pos)}" cy="${baseline}" r="3.4" fill="${tcol(s.ptm_type)}" stroke="#fff" stroke-width="1"><title>${esc(s.site)} ${esc(s.ptm_type)} (position ${s.pos})</title></circle>`).join("");
  const lbls = sites.map(s => `<text class="xt-lbl" data-k="${esc(skey(s))}" x="${x(s.pos)}" y="${baseline - 12}" font-size="10" font-weight="700" text-anchor="middle" fill="#17363a" opacity="0" style="pointer-events:none">${esc(s.site)}</text>`).join("");

  // ---- legend: types that do not occur in this result are greyed out
  const present = new Set(sites.map(s => s.ptm_type));
  const xtLegend = `<div class="legend">${Object.entries(TYPE_COLOR).map(([k, c]) =>
    `<span style="${present.has(k) ? "" : "opacity:.35"}" title="${present.has(k) ? "" : "No pair of this type in this result"}"><i style="background:${c}"></i>${k}</span>`).join("")}</div>`;

  // ---- what the protein actually has (explains every zero)
  const ubt = j.units_by_type || null;
  const typeKeys = ubt ? Object.keys(ubt) : [];
  const typeLine = ubt
    ? `<p class="muted" style="margin:0 0 8px">This protein has <b>${esc(j.n_units)}</b> modified residue(s) in the database: ` +
      `${typeKeys.map(k => `${esc(ubt[k])} ${esc(k)}`).join(", ") || "none"}.</p>`
    : "";

  // ---- mapping self-check: arcs must match the recorded residue positions
  const siteNum = s => { const m = String(s).match(/(\d+)/); return m ? Number(m[1]) : null; };
  const problems = [];
  pairs.forEach(p => {
    [p.a, p.b].forEach(s => {
      const n = siteNum(s.site);
      if (n !== null && n !== Number(s.pos)) problems.push(`${s.site} is drawn at position ${s.pos}`);
      if (Number(s.pos) > L) problems.push(`${s.site} (position ${s.pos}) lies beyond the protein length ${L}`);
    });
    if (Math.abs(Number(p.b.pos) - Number(p.a.pos)) !== Number(p.distance))
      problems.push(`${p.a.site}–${p.b.site}: Δ ${p.distance} does not equal the position difference`);
  });
  const check = !pairs.length ? "" : (problems.length
    ? `<div class="notice warn"><b>Mapping check found ${problems.length} problem(s):</b><br>${[...new Set(problems)].slice(0, 8).map(esc).join("<br>")}</div>`
    : `<div class="notice ok"><b>Mapping check passed.</b> All ${pairs.length} pair(s): each site number equals the position where its arc end is drawn, Δ equals the position difference, and every site lies within the ${L}-residue protein.</div>`);

  // ---- honest framing: crosstalk means DIFFERENT PTM types
  let notice = "";
  if (pairs.length && nCross === 0) {
    notice = `<div class="notice warn"><b>0 pairs of different PTM types.</b> All ${pairs.length} returned pair(s) link the same modification type, ` +
      `so they are <b>proximal same-type clusters</b>, not crosstalk between different modifications.</div>`;
  }
  const heading = pairs.length && nCross === 0 ? "Proximal same-type PTM clusters" : "Crosstalk arcs";

  let emptyMsg;
  if (sameType) {
    emptyMsg = `No two modified residues of this protein lie within ${distance} residues of each other.`;
  } else if (ubt && typeKeys.length < 2) {
    emptyMsg = `<b>0 pairs is correct here.</b> ${typeKeys.length ? `Only <b>${esc(typeKeys[0])}</b> is recorded for this protein` : "No PTM is recorded for this protein"}, ` +
      `so no pair of <b>different</b> PTM types can exist. Tick “Include same-type pairs” to see neighbouring ${esc(typeKeys[0] || "")} sites.`;
  } else {
    emptyMsg = `No pair of <b>different</b> PTM types lies within ${distance} residues. Try a larger distance, or tick “Include same-type pairs”.`;
  }

  // ---- table helpers
  const hint = (s, t) => `${s} <span class="muted" title="${esc(t)}" style="cursor:help">ⓘ</span>`;
  const bothPmids = p => { const s = new Set(p.b.pmids || []); return (p.a.pmids || []).filter(v => s.has(v)); };
  const pmidCell = p => {
    const both = bothPmids(p), either = p.pmids || [];
    if (both.length) return `<b>both sites:</b> ${both.map(pm).join(", ")}`;
    return either.length ? `<span class="muted">one site each:</span> ${either.map(pm).join(", ")}` : "—";
  };
  const kinaseCell = p => {
    const shared = p.shared_kinases || [], ka = p.a.kinases || [], kb = p.b.kinases || [];
    if (shared.length) return chips(shared);
    if (!ka.length && !kb.length) return '<span class="muted" title="No kinase is recorded for either site">not recorded</span>';
    return `<span class="muted" title="Kinases recorded per site; none is shared">${esc(p.a.site)}: ${esc(ka.join(", ") || "—")}<br>${esc(p.b.site)}: ${esc(kb.join(", ") || "—")}</span>`;
  };

  out.innerHTML =
    cards([
      [pairs.length, "Candidate pairs", `≤ ${distance} residues`],
      [nCross, "Different-type pairs", "true crosstalk candidates"],
      [pairs.filter(p => Number(p.distance) === 0).length, "Same-residue competition", "two PTM types on one residue"],
      [nClose, `Close pairs (≤ ${CLOSE_DIST} residues)`],
      [pairs.filter(p => (p.shared_diseases || []).length > 0).length, "Share a disease context"]
    ]) +
    notice +
    `<div class="panel">
       <h3>${heading}</h3>
       ${typeLine}
       ${pairs.length ? check + track(L, [], arcs + dots + lbls, H) +
         `<div class="muted" style="text-align:center;font-size:12px">Residue position</div>` + xtLegend
         : `<div class="empty">${emptyMsg}</div>`}
       <p class="caveat">Each arc is one returned pair; thicker arcs are ≤ ${CLOSE_DIST} residues apart. Arc height grows with the
       distance between the sites and does not represent interaction strength. Hover an arc or a table row to highlight it and label its two sites.
       Proximity suggests a candidate relationship, not proven biochemical crosstalk.</p>
     </div>` +
    `<h3>PTM-type combinations</h3><div>${
      Object.entries(j.type_pair_counts || {}).map(([k, v]) => `<span class="chip big">${esc(k)} <b>${esc(v)}</b></span>`).join("")
      || '<span class="muted">None detected</span>'}</div>` +
    `<h3>Pairs</h3>` +
    table(
      ["Site A", "PTM A", "Site B", "PTM B", "Δ residues", "Class",
       hint("Shared kinase", "“not recorded” means no kinase is recorded in this database for the site(s), not that none exists."),
       hint("Shared disease", "“—” means no disease context is recorded for both sites, not that none exists."),
       hint("PMIDs", "“both sites” = a paper recorded for both residues. “one site each” = papers recorded separately, which does not show the sites were studied together.")],
      pairs.map(p => [
        `<b>${esc(p.a.site)}</b>`, esc(p.a.ptm_type), `<b>${esc(p.b.site)}</b>`, esc(p.b.ptm_type), esc(p.distance),
        `<span class="badge">${esc(p.class)}</span>` +
          (Number(p.distance) > 0 && Number(p.distance) <= CLOSE_DIST ? ` <span class="badge red" title="Within ${CLOSE_DIST} residues">close</span>` : ""),
        kinaseCell(p), chips(p.shared_diseases || []), pmidCell(p)
      ]),
      "No eligible PTM pairs were found for the selected settings."
    ) +
    method(j.method);

  // ---- hover linking: arc <-> table row <-> the two endpoint dots (enlarged and labelled)
  const arcEls = [...out.querySelectorAll(".xt-arc")];
  const dotEls = [...out.querySelectorAll(".xt-dot")];
  const lblEls = [...out.querySelectorAll(".xt-lbl")];
  const rows = [...out.querySelectorAll(".table-wrap tbody tr")];
  const on = i => {
    const p = pairs[i]; if (!p) return;
    const ka = skey(p.a), kb = skey(p.b);
    arcEls.forEach(el => {
      const hit = Number(el.dataset.i) === i;
      el.style.opacity = hit ? 1 : 0.12;
      el.setAttribute("stroke-width", hit ? 4 : el.dataset.w);
    });
    dotEls.forEach(d => d.setAttribute("r", d.dataset.k === ka || d.dataset.k === kb ? 6 : 3.4));
    lblEls.forEach(l => {
      const isA = l.dataset.k === ka, isB = l.dataset.k === kb;
      l.setAttribute("opacity", isA || isB ? 1 : 0);
      l.setAttribute("y", isB && !isA ? baseline - 28 : baseline - 12);
    });
  };
  const off = () => {
    arcEls.forEach(el => { el.style.opacity = ""; el.setAttribute("stroke-width", el.dataset.w); });
    dotEls.forEach(d => d.setAttribute("r", 3.4));
    lblEls.forEach(l => l.setAttribute("opacity", 0));
  };
  arcEls.forEach(el => {
    const i = Number(el.dataset.i);
    el.addEventListener("mouseenter", () => { on(i); if (rows[i]) rows[i].style.background = "#eaf6f6"; });
    el.addEventListener("mouseleave", () => { off(); if (rows[i]) rows[i].style.background = ""; });
  });
  rows.forEach((r, i) => {
    r.addEventListener("mouseenter", () => on(i));
    r.addEventListener("mouseleave", off);
  });
}
// ----------------------------------------------------------------- NETWORK
async function network() {
  busy("Building network…");
  const scope = val("p-scope", "protein");
  const j = await api(`/api/f/network/${scope === "all" ? "ALL" : up()}`); if (!j) return;
  const dis = $("#p-disease").checked;
  const els = [], seen = new Set();
  const node = (id, label, type, extra = {}) => { if (!seen.has(id)) { seen.add(id); els.push({data: {id, label, type, ...extra}}); } };
  j.paths.forEach((p, i) => {
    node("k:" + p.kinase, p.kinase, "kinase"); node("p:" + p.protein, p.protein, "protein");
    node("s:" + p.site_id, p.site.split(" ")[1] || p.site, "site", {ptm: p.ptm_type, full: p.site, pmids: p.pmids.join(", "), conf: p.confidence});
    els.push({data: {id: "e1" + i, source: "k:" + p.kinase, target: "s:" + p.site_id, label: "modifies"}});
    if (!seen.has("es" + p.site_id)) { seen.add("es" + p.site_id); els.push({data: {id: "es" + p.site_id, source: "s:" + p.site_id, target: "p:" + p.protein, label: "on"}}); }
    if (dis) p.diseases.forEach(d => { node("d:" + d, d, "disease"); if (!seen.has("ed" + p.site_id + d)) { seen.add("ed" + p.site_id + d); els.push({data: {id: "ed" + p.site_id + d, source: "s:" + p.site_id, target: "d:" + d, label: "context"}}); } });
  });
  out.innerHTML = cards([[j.kinases.length, "Kinases named", "in records"], [j.n_sites_with_kinase, "Sites with a kinase"], [j.n_sites_without_kinase, "Sites with no kinase recorded", "not drawn"], [j.cascades.length, "Cascade links", "kinase is itself a database protein"]]) +
    `<div class="panel"><div id="cy"></div><div class="legend"><span><i style="background:#f2c14e"></i>kinase</span><span><i style="background:#0d2b2e"></i>protein</span><span><i style="background:#d64545"></i>PTM site</span><span><i style="background:#8a5cf6"></i>cancer/disease context</span></div><div id="cy-info" class="muted">Click a node for details.</div></div>` +
    `<h3>Kinases</h3>` + table(["Kinase", "Sites", "Proteins", "Publications", "In database as protein", "Disease context"], j.kinases.map(k => [`<b>${esc(k.kinase)}</b>`, k.n_sites, k.n_proteins, k.n_pmids, k.also_in_database ? "✓" : "—", chips(k.diseases)]), "No record names a kinase.") +
    (j.cascades.length ? `<h3>Cascade candidates</h3>` + table(["Upstream kinase", "Modifies", "Which is kinase", "…that modifies", "PMIDs"], j.cascades.map(c => [esc(c.upstream_kinase), esc(c.via_site), esc(c.kinase_protein), c.downstream_sites.map(esc).join(", "), c.pmids.map(pm).join(", ")])) : "") + method(j.method);
  if (!window.cytoscape) { $("#cy").innerHTML = '<div class="notice warn">Cytoscape.js could not be loaded (offline?). Tables below are still complete.</div>'; return; }
  const col = {kinase: "#f2c14e", protein: "#0d2b2e", site: "#d64545", disease: "#8a5cf6"};
  const cy = cytoscape({container: $("#cy"), elements: els, wheelSensitivity: .3,
    style: [{selector: "node", style: {label: "data(label)", "background-color": e => col[e.data("type")], color: "#17363a", "font-size": 9, "text-valign": "bottom", "text-margin-y": 3, width: e => e.data("type") === "protein" ? 34 : 20, height: e => e.data("type") === "protein" ? 34 : 20}},
      {selector: "edge", style: {width: 1.2, "line-color": "#b6c7c9", "target-arrow-color": "#8aa0a4", "target-arrow-shape": "triangle", "curve-style": "bezier"}},
      {selector: ".hl", style: {"line-color": "#18a6a6", "target-arrow-color": "#18a6a6", "border-width": 3, "border-color": "#18a6a6", width: 2.5}}],
    layout: {name: "cose", animate: false, nodeRepulsion: 9000, idealEdgeLength: 70}});
  cy.on("tap", "node", e => { const n = e.target; cy.elements().removeClass("hl"); n.closedNeighborhood().addClass("hl"); const d = n.data();
    $("#cy-info").innerHTML = `<b>${esc(d.full || d.label)}</b> · ${esc(d.type)}${d.ptm ? " · " + esc(d.ptm) : ""}${d.conf ? " · " + esc(d.conf) + " confidence" : ""}${d.pmids ? " · PMID " + esc(d.pmids) : ""}<br>Connected: ${n.neighborhood("node").map(m => esc(m.data("label"))).join(", ")}`; });
}

// ------------------------------------------------------------- CONSERVATION
async function conservation() {
  busy("Fetching Ensembl orthologue alignments…");
  const j = await api(`/api/f/conservation/${up()}`); if (!j) return;
  const ok = j.sites.filter(s => s.status === "ok");
  const full = ok.filter(s => s.n_aligned && s.n_identical === s.n_aligned).length;
  out.innerHTML = cards([[j.n_orthologs, "Orthologues aligned", j.species.join(", ")], [ok.length, "Sites checked"], [full, "Identical in every aligned species"], [j.sites.length - ok.length, "Sites not checkable", "residue mismatch / out of range"]]) +
    `<div class="table-wrap"><table><thead><tr><th>Site</th><th>PTM</th><th>Human</th><th>Conserved</th>${j.species.map(s => `<th class="rot">${esc(s)}</th>`).join("")}</tr></thead><tbody>` +
    j.sites.map(s => `<tr><td><b>${esc(s.site)}</b></td><td>${esc(s.ptm_type)}</td><td class="mono">${esc(s.human_residue || "—")}</td><td>${s.status === "ok" ? `<b>${s.n_identical}/${s.n_aligned}</b>` : `<span class="muted" title="${esc(s.status)}">n/a</span>`}</td>` +
      j.species.map(sp => { const o = s.orthologs.find(x => x.species === sp); if (!o) return `<td class="cons none">·</td>`; const c = o.residue === null ? "none" : o.gap ? "gap" : o.identical ? "same" : "diff"; return `<td class="cons ${c}" title="${esc(sp)} · ${esc(o.orthology || "")} · flank identity ${o.flank_identity ?? "n/a"}">${esc(o.residue ?? "·")}</td>`; }).join("") + "</tr>").join("") +
    `</tbody></table></div><div class="legend"><span><i style="background:#1f9d55"></i>identical</span><span><i style="background:#e0a020"></i>different residue</span><span><i style="background:#8a97a0"></i>gap in alignment</span></div>` +
    `<p class="caveat">Conservation of a residue is consistent with functional constraint but does not prove the modification itself is conserved. Species with no returned alignment are left blank, not imputed.</p>` + method(j.method);
}

// --------------------------------------------------------------- MUTATIONS
let studiesLoaded = false;
async function loadStudies() {
  if (studiesLoaded) return; studiesLoaded = true;
  try { const r = await fetch("/api/cbio/studies"); const j = await r.json();
    if (j.available === false) { studiesLoaded = false; fallbackStudies("cBioPortal unavailable: " + j.reason); return; }
    $("#study-list").innerHTML = j.studies.map(s => `<option value="${esc(s.id)}">${esc(s.name)} (${s.samples ?? "?"} samples)</option>`).join("");
    $("#study-status").textContent = `${j.studies.length} studies loaded from cBioPortal.`;
  } catch (e) { studiesLoaded = false; fallbackStudies("Could not load the study list."); }
}
function fallbackStudies(why) {
  const ids = ["brca_tcga", "luad_tcga", "lusc_tcga", "coadread_tcga", "prad_tcga", "skcm_tcga", "stad_tcga", "ucec_tcga", "gbm_tcga", "blca_tcga"];
  $("#study-list").innerHTML = ids.map(i => `<option value="${i}">`).join("");
  $("#study-status").textContent = why + " Common study IDs are suggested (availability not checked); you can also type any ID. Reload to retry.";
}
async function mutations() {
  loadStudies();
  const w = val("p-window", 0), study = val("p-study", "");
  busy("Fetching real variants…");
  const j = await api(`/api/f/variants/${up()}?window=${w}`);
  let html = "";
  if (j) { const s = j.summary;
    html += `<h3>A · Curated & large-scale variants (UniProt / EBI Proteins API)</h3>` + cards([[s.n_variants_total, "Single-residue variants known"], [s.n_at_modified_residue, "Alter a modified residue"], [s.n_cancer_associated_hits, "Cancer-associated hits"], [`${s.sites_with_hits}/${s.sites_total}`, "PTM sites with a variant"]]) +
      table(["PTM site", "PTM", "Variant", "Offset", "Consequence", "Cancer-associated", "Disease annotation", "COSMIC", "ClinVar", "Site kinase", "Site disease"],
        j.hits.slice(0, 400).map(h => [`<b>${esc(h.site)}</b>`, esc(h.ptm_type), `<b>${esc(h.variant)}</b>${h.at_modified_residue ? ' <span class="badge red">at residue</span>' : ""}`, h.offset, esc(h.consequence || "—"), h.cancer_associated ? "✓" : "—", chips(h.diseases), chips(h.cosmic.slice(0, 3)), chips(h.clinvar.slice(0, 2)), chips(h.site_kinases), chips(h.site_diseases)]),
        "No known variant lies at or within the chosen window of a recorded PTM residue.") + method(j.method);
  }
  if (study) {
    const c = await fetchJson(`/api/f/cohort/${up()}?study=${encodeURIComponent(study)}&window=${w}`);
    html += `<h3>B · Cohort mutations in <code>${esc(study)}</code> (cBioPortal)</h3>` + (c.available === false ? `<div class="notice warn">${esc(c.reason)}</div>` :
      cards([[c.n_sequenced ?? "?", "Sequenced samples"], [c.n_mutations, "Mutations in this gene"], [c.n_mut_on_ptm_residues, "On PTM residues", `expected ${c.expected_on_ptm_uniform} if uniform`], [pfmt(c.p_enrichment), "Enrichment p (binomial)"]]) +
      table(["PTM site", "PTM", "Mutations at residue", "Samples at residue", `Samples ±${c.window}`, "Frequency", "Most frequent changes", "Types"],
        c.sites.filter(s => s.n_mut_window).map(s => [`<b>${esc(s.site)}</b>`, esc(s.ptm_type), s.n_mut_at_residue, s.n_samples_at_residue, s.n_samples_window, s.freq_window == null ? "—" : (s.freq_window * 100).toFixed(2) + "%", s.top_changes.map(t => `${esc(t.change)}×${t.n}`).join(", "), Object.entries(s.types).map(([k, v]) => `${esc(k)}:${v}`).join(", ")]), "No mutations at recorded PTM residues in this study.") + method(c.method));
  } else html += `<div class="notice">Select a cBioPortal cancer study above to add cohort-level mutation counts. No cohort is assumed.</div>`;
  out.innerHTML = html || out.innerHTML;
  out.insertAdjacentHTML("beforeend", `<p class="caveat">Positional overlap shows where a mutation and a modified residue coincide; whether the mutation changes the modification must be tested experimentally.</p>`);
}
async function fetchJson(u) { const r = await fetch(u); return r.json(); }

// ----------------------------------------------------------------- SPATIAL

 // ----------------------------------------------------------------- SPATIAL
let viewer = null, viewerElement = null, lastSpatial = null, spatialModelText = null, spatialModelKey = null;

async function spatial() {
  busy("Locating real structures and verifying residues (first run can take a minute)…");

  const j = await api(
    `/api/f/spatial/${encodeURIComponent(up())}?threshold=${val("p-thr", 10)}&alphafold=${$("#p-af").checked ? 1 : 0}`
  );
  if (!j) return;

  lastSpatial = j;

  const cov = j.coverage || {};
  const exp = cov.experimental || {
    mapped: 0, total: j.n_units || 0, sites: [], unmapped: []
  };
  const af = cov.alphafold || {
    mapped: 0, total: j.n_units || 0, sites: [], unmapped: []
  };

  const sourceText = j.source === "alphafold"
    ? "AlphaFold predicted model"
    : "Experimental PDB";

  const mappingRows = [];
  const allSites = new Set([
    ...(exp.sites || []),
    ...(af.sites || []),
    ...(j.unmapped || [])
  ]);

  for (const site of [...allSites].sort()) {
    const inExp = (exp.sites || []).includes(site);
    const inAf = (af.sites || []).includes(site);
    const unmapped = (j.unmapped || []).includes(site);

    mappingRows.push([
      `<b>${esc(site)}</b>`,
      inExp
        ? `<span class="badge">Experimental${exp.structure_id ? ` · ${esc(exp.structure_id)}:${esc(exp.chain || "")}` : ""}</span>`
        : "—",
      inAf
        ? `<span class="badge warn">AlphaFold${af.structure_id ? ` · ${esc(af.structure_id)}` : ""}</span>`
        : "—",
      unmapped
        ? '<span class="badge bad">Not mapped</span>'
        : '<span class="badge ok">Mapped</span>'
    ]);
  }

  const noteHtml = (j.notes || [])
    .map(n => `<div class="notice">${esc(n)}</div>`)
    .join("");

  out.innerHTML =
    `<div class="notice ${j.source === "alphafold" ? "warn" : "ok"}">
      <b>${esc(sourceText)}: ${esc(j.structure_id)}</b>
      chain ${esc(j.chain)} — ${esc(j.basis)}
    </div>` +
    noteHtml +
    cards([
      [`${j.n_mapped}/${j.n_units}`, "PTM sites mapped in selected model",
        (j.unmapped || []).length
          ? "Unmapped: " + j.unmapped.join(", ")
          : "All selected-model sites mapped"],
      [`${exp.mapped}/${exp.total}`, "Best experimental PDB coverage",
        exp.structure_id ? `${exp.structure_id}:${exp.chain || ""}` : "No verified experimental chain"],
      [`${af.mapped}/${af.total}`, "AlphaFold coverage",
        af.structure_id ? `${af.structure_id}:${af.chain || ""}` : "No usable AlphaFold model"],
      [(j.pairs || []).length, "Pairs measured"],
      [j.n_close, `Within ${j.threshold} Å (Cα–Cα)`],
      [j.n_long_range_close, "Far in sequence, close in 3D", "Δ ≥ 20 aa"]
    ]) +
    `<h3>PTM structural mapping</h3>` +
    table(
      ["PTM site", "Experimental evidence", "AlphaFold evidence", "Selected analysis"],
      mappingRows,
      "No PTM residue mapping could be verified from available real structures."
    ) +
    `<div class="split spatial-split">
      <div class="viewer-card small">
        <div id="v3d"></div>
        <button type="button" id="v3d-reset" class="btn btn-secondary btn-sm" style="position:absolute;right:10px;top:10px;z-index:2">Reset view</button>
        <div id="v3d-status" class="viewer-status">Loading structure coordinates…</div>
      </div>
      <div>` +
    table(
      ["A", "B", "Δ aa", "Cα–Cα Å", "Closest atoms Å", "Flags"],
      (j.pairs || []).map((p, i) => [
        `<a href="#" data-i="${i}" class="pairlink">${esc(p.site_a)}</a>`,
        `<a href="#" data-i="${i}" class="pairlink">${esc(p.site_b)}</a>`,
        p.seq_sep,
        `<b>${p.ca_dist}</b>`,
        p.min_atom_dist,
        [
          p.cross_type ? '<span class="badge">cross-type</span>' : "",
          p.long_range_close ? '<span class="badge red">long-range close</span>' : "",
          !p.both_confident ? '<span class="badge warn">pLDDT&lt;70</span>' : ""
        ].join(" ")
      ]),
      "No verified residue pairs are available for this protein."
    ) +
    `</div></div>
    <h3>Selected-model provenance</h3>` +
    table(
      ["Site", "Mapping status", "Source", "Structure", "Chain", "Residue"],
      (j.mapping || []).map(m => [
        esc(m.site),
        esc(m.mapping_status),
        esc(m.structure_source || "—"),
        esc(m.structure_id || "—"),
        esc(m.chain || "—"),
        m.structure_residue_number ?? "—"
      ]),
      "No mapping records returned."
    ) +
    `<p class="caveat">${esc(j.note || "")} ${esc((j.notes || []).join(" "))}</p>`;

  out.querySelectorAll(".pairlink").forEach(a => {
    a.addEventListener("click", e => {
      e.preventDefault();
      show3d(Number(a.dataset.i));
    });
  });

  const rb = $("#v3d-reset");
  if (rb) rb.addEventListener("click", () => {
    if (viewer && viewerElement === $("#v3d")) { viewer.zoomTo(); viewer.render(); }
  });

  // Load the selected structure immediately, without waiting for a pair click.
  await loadSpatialOverview();
}

function spatialModelCacheKey(j) {
  return `${up()}|${j.source || ""}|${j.structure_id || ""}`;
}

async function loadSpatialOverview() {
  const j = lastSpatial;
  const st = $("#v3d-status");

  if (!j || !st) return false;

  if (!window.$3Dmol) {
    st.textContent =
      "3Dmol.js did not load. Check internet/CDN access and refresh.";
    return false;
  }

  if (!j.structure_id || !j.source) {
    st.textContent = "No verified structure is available for this protein.";
    return false;
  }

  st.textContent = "Loading selected structure coordinates…";

  try {
    const key = spatialModelCacheKey(j);
    let text = spatialModelKey === key ? spatialModelText : null;

    if (!text) {
      const url =
        `/api/f/spatial/${encodeURIComponent(up())}/model` +
        `?source=${encodeURIComponent(j.source)}` +
        `&id=${encodeURIComponent(j.structure_id)}`;

      const r = await fetch(url);

      if (!r.ok) {
        throw new Error(`Coordinate endpoint returned HTTP ${r.status}`);
      }

      text = await r.text();

      if (
        !/^\s*(?:HEADER|TITLE|COMPND|ATOM|MODEL|REMARK|CRYST1)/m.test(text) ||
        !/^ATOM  /m.test(text)
      ) {
        throw new Error("Response does not contain a valid PDB coordinate model");
      }

      spatialModelText = text;
      spatialModelKey = key;
    }

    const element = $("#v3d");
    // The result area is rebuilt on every run, so the viewer MUST be bound to the element that is on the
    // page now. (A viewer cached from an earlier run drew into a detached <div> -> blank area.)
    viewer = PTMViewer.mount(element);
    viewerElement = element;
    const nAtoms = PTMViewer.load(viewer, text);
    viewer.setStyle({}, j.source === "alphafold"
      ? { cartoon: { colorfunc: PTMViewer.plddtColor } }
      : { cartoon: { color: "#cfd8dc" } });

    let shown = 0;

    for (const m of (j.mapping || [])) {
      if (
        m.mapping_status !== "mapped" ||
        m.structure_residue_number == null
      ) continue;

      const selection = {
        chain: String(m.chain || j.chain || ""),
        resi: Number(m.structure_residue_number)
      };

      if (!viewer.selectedAtoms({ ...selection, atom: "CA" }).length) {
        continue;
      }

      viewer.setStyle(selection, {
        cartoon: { color: "#18a6a6" },
        stick: { radius: 0.18, colorscheme: "cyanCarbon" }
      });

      shown += 1;
    }

    const fin = PTMViewer.finish(viewer, element);
    if (!fin.ok) throw new Error(fin.reason);

    const kind = j.source === "alphafold"
      ? "AlphaFold PREDICTED model (cartoon coloured by pLDDT)"
      : "Experimental PDB";

    st.innerHTML =
      `${esc(kind)} ${esc(j.structure_id)}, chain ${esc(j.chain)}: ${nAtoms} atoms drawn. ` +
      `${shown} of ${(j.mapping || []).length} PTM site(s) are residue-verified and marked. ` +
      `Select a pair to focus its distance. ` +
      (j.source === "alphafold" ? PTMViewer.PLDDT_LEGEND : "");

    return true;
  } catch (error) {
    st.textContent =
      `3D structure could not be displayed: ${error.message}. ` +
      "Check connection diagnostics and refresh.";

    return false;
  }
}

async function show3d(i) {
  const j = lastSpatial;
  const p = j && Array.isArray(j.pairs) ? j.pairs[i] : null;
  const st = $("#v3d-status");

  if (!j || !p || !st) return;

  if (!window.$3Dmol) {
    st.textContent = "3Dmol.js did not load. Check internet access.";
    return;
  }

  if (!viewer || viewerElement !== $("#v3d") || spatialModelKey !== spatialModelCacheKey(j)) {
    const loaded = await loadSpatialOverview();
    if (!loaded || !viewer) return;
  }

  try {
    viewer.removeAllShapes();
    viewer.removeAllLabels();
    viewer.setStyle({}, { cartoon: { color: "#cfd8dc" } });

    for (const m of (j.mapping || [])) {
      if (
        m.mapping_status !== "mapped" ||
        m.structure_residue_number == null
      ) continue;

      const selection = {
        chain: String(m.chain || j.chain || ""),
        resi: Number(m.structure_residue_number)
      };

      if (viewer.selectedAtoms({ ...selection, atom: "CA" }).length) {
        viewer.setStyle(selection, {
          cartoon: { color: "#18a6a6" },
          stick: { radius: 0.18, colorscheme: "cyanCarbon" }
        });
      }
    }

    const A = { chain: j.chain, resi: p.resi_a };
    const B = { chain: j.chain, resi: p.resi_b };

    viewer.setStyle(A, {
      cartoon: { color: "#d64545" },
      stick: { colorscheme: "redCarbon" }
    });
    viewer.setStyle(B, {
      cartoon: { color: "#1f6fd6" },
      stick: { colorscheme: "blueCarbon" }
    });

    const a = viewer.selectedAtoms({ ...A, atom: "CA" })[0];
    const b = viewer.selectedAtoms({ ...B, atom: "CA" })[0];

    if (!a || !b) {
      st.textContent =
        "Selected pair is not present in the loaded structure chain. " +
        "No distance marker was drawn.";
      return;
    }

    viewer.addCylinder({
      start: { x: a.x, y: a.y, z: a.z },
      end: { x: b.x, y: b.y, z: b.z },
      radius: 0.12,
      color: "black",
      dashed: true
    });

    viewer.addLabel(`${p.ca_dist} Å`, {
      position: {
        x: (a.x + b.x) / 2,
        y: (a.y + b.y) / 2,
        z: (a.z + b.z) / 2
      },
      backgroundColor: "black",
      fontColor: "white",
      fontSize: 12
    });

    viewer.zoomTo({ or: [A, B] });
    if (typeof viewer.resize === "function") viewer.resize();
    viewer.render();

    st.textContent =
      `${p.site_a} (red) ↔ ${p.site_b} (blue): ` +
      `Cα–Cα ${p.ca_dist} Å, closest atoms ${p.min_atom_dist} Å`;
  } catch (error) {
    st.textContent = `Could not render this pair: ${error.message}`;
  }
}


// -------------------------------------------------------------- HYPOTHESES
async function hypotheses() {
  busy("Collecting evidence layers…");
  const q = ["spatial", "variants", "conservation"].filter(k => $("#l-" + k).checked).map(k => k + "=1").join("&");
  const j = await api(`/api/f/hypotheses/${up()}?${q}`); if (!j) return;
  if (!Array.isArray(j.hypotheses)) {
    out.innerHTML = '<div class="notice warn">Analysis waiting for data. The server did not return a valid hypothesis list; check the source records and server log.</div>';
    return;
  }
  const kinds = [...new Set(j.hypotheses.map(h => h.kind))];
  const layerTxt = Object.entries(j.layers).map(([k, v]) => `<span class="chip ${v.used ? "ok" : "bad"}" title="${esc(v.reason || "")}">${esc(k)} ${v.used ? "✓" : "✗ " + esc(v.reason || "")}</span>`).join("");
  const draw = () => { const k = $("#h-kind").value; const list = j.hypotheses.filter(h => !k || h.kind === k);
    $("#h-list").innerHTML = list.map(h => `<article class="hyp"><header><span class="badge">${esc(h.kind)}</span><span class="hid">${esc(h.id)}</span><span class="hscore" title="Fraction of applicable evidence criteria met">${h.criteria_met}/${h.criteria_total}</span></header>
      <h4>${esc(h.statement)}</h4><div class="bar"><i style="width:${h.score * 100}%"></i></div>
      <p class="muted">${h.basis.map(esc).join(" ")}</p>
      <ul class="crit">${h.criteria.map(c => `<li class="${!c.applicable ? "na" : c.met ? "met" : "unmet"}">${!c.applicable ? "–" : c.met ? "✓" : "✗"} ${esc(c.name)}</li>`).join("")}</ul>
      <details><summary>Suggested experiments & evidence</summary><ul>${h.suggested_assays.map(a => `<li>${esc(a)}</li>`).join("")}</ul>
      <p>Sites: ${chips(h.evidence.sites)} Kinases: ${chips(h.evidence.kinases)} Context: ${chips(h.evidence.diseases)}<br>PMIDs: ${h.evidence.pmids.map(pm).join(", ") || "—"}</p></details>
      <footer>${esc(h.disclaimer)}</footer></article>`).join("") || '<div class="empty">No hypotheses could be derived from the current evidence layers.</div>'; };
  out.innerHTML = cards([[j.n, "Research questions"], [kinds.length, "Question types"]]) + `<div class="panel"><b>Evidence layers:</b> ${layerTxt}</div>` +
    `<div class="filter-bar"><div class="form-group"><label>Type</label><select id="h-kind"><option value="">All</option>${kinds.map(k => `<option>${esc(k)}</option>`).join("")}</select></div><div class="form-actions"><button class="btn btn-secondary btn-sm" id="h-dl">↓ Export JSON</button></div></div><div id="h-list"></div><p class="caveat">${esc(j.policy)} The score is only the share of applicable criteria satisfied — every criterion is shown.</p>`;
  $("#h-kind").addEventListener("change", draw); draw();
  $("#h-dl").addEventListener("click", () => { const a = document.createElement("a"); a.href = URL.createObjectURL(new Blob([JSON.stringify(j, null, 2)], {type: "application/json"})); a.download = `${j.protein.gene}_hypotheses.json`; a.click(); });
}

// ---------------------------------------------------------------- PRIORITY
const CRIT = {confidence: "High confidence", in_vivo: "in vivo evidence", kinase: "Kinase known", disease: "Disease context", replicated: "≥2 publications", hotspot: "In enriched hotspot", crosstalk: "Crosstalk partner"};
async function priority() {
  busy();
  const w = Object.keys(CRIT).map(k => `${k}=${val("w-" + k, 1)}`).join("&");
  const j = await api(`/api/f/priority/${val("p-scope", "protein") === "all" ? "ALL" : up()}?${w}`); if (!j) return;
  out.innerHTML = cards([[j.n_sites, "Sites scored"], [j.sites.filter(s => s.score === 1).length, "Meet every criterion"]]) +
    table(["#", "Gene", "Site", "PTM", "Score", ...Object.values(CRIT), "Kinases", "Disease context"], j.sites.map((s, i) => [i + 1, `<a href="/protein/${esc(s.uniprot)}">${esc(s.gene)}</a>`, `<b>${esc(s.site)}</b>`, esc(s.ptm_type), `<b>${s.score.toFixed(2)}</b>`,
      ...Object.keys(CRIT).map(k => s.criteria[k] ? '<span class="dot on"></span>' : '<span class="dot"></span>'), chips(s.kinases), chips(s.diseases)])) +
    `<p class="caveat">PhosphoRank is a checklist for triage: each criterion is 0/1 from stored evidence; weights default to equal and are yours to change. It is not a probability of functional importance.</p>`;
}

// --------------------------------------------------------------- INTEGRITY
async function integrity() {
  busy("Comparing records with the UniProt sequence…");
  const j = await api(`/api/f/integrity/${up()}`); if (!j) return;
  const s = j.summary;
  out.innerHTML = cards([[j.n_records, "Records checked"], [s.residue_ok || 0, "Residue letter matches UniProt"], [s.residue_bad || 0, "Residue mismatch"], [s.flank_not_found || 0, "Flanking sequence not found in UniProt"], [j.length_matches_catalog === null ? "—" : j.length_matches_catalog ? "✓" : "✗", "Catalogue length matches UniProt", `${j.length_uniprot} aa`]]) +
    (s.residue_bad ? `<div class="notice warn"><b>${s.residue_bad} residue letter(s) disagree with UniProt.</b> Check the numbering/isoform in the original source before using them.</div>` : `<div class="notice ok">All checked residue letters match UniProt.</div>`) +
    (s.flank_not_found ? `<div class="notice"><b>${s.flank_not_found} stored flanking window(s) could not be matched</b> (needs ≥90% identity). The residue letters can still be fine — compare the stored and UniProt windows below; a different flank format in the catalogue is a common cause.</div>` : "") +
    `<h3>Records needing review</h3>` + table(["Site", "PTM", "PMID", "UniProt residue", "Residue OK", "Flank match", "Stored flank", "UniProt window"], j.bad.map(r => [`<b>${esc(r.site)}</b>`, esc(r.ptm_type), esc(r.pmid), esc(r.uniprot_residue || "—"), r.residue_matches_uniprot === null ? "—" : r.residue_matches_uniprot ? "✓" : "✗", r.flank_identity == null ? "—" : (r.flank_identity * 100).toFixed(0) + "%", `<code>${esc(r.stored_flank || "—")}</code>`, `<code>${esc(r.uniprot_window || "—")}</code>`]), "None.") +
    `<details><summary>All ${j.n_records} records</summary>` + table(["Site", "PTM", "UniProt residue", "Residue OK", "Flank found"], j.records.map(r => [esc(r.site), esc(r.ptm_type), esc(r.uniprot_residue || "—"), String(r.residue_matches_uniprot), String(r.flank_found_in_uniprot)])) + `</details>`;
}


// A cohort with no mutation in this gene is a property of the cohort, not an error. Offer real alternatives.
async function findMutatedStudies(thenRun) {
  const box = $("#study-suggest");
  box.innerHTML = '<span class="muted">Checking large cBioPortal studies for mutations in this gene (can take up to a minute)…</span>';
  try {
    const r = await fetch(`/api/cbio/suggest/${up()}`); const j = await r.json();
    if (j.available === false) { box.innerHTML = `<div class="notice warn">${esc(j.reason || "Unavailable")}</div>`; return; }
    if (!j.studies.length) { box.innerHTML = `<div class="notice">None of the ${j.checked} large studies checked has a missense mutation in ${esc(j.gene)}. Type another study ID manually.</div>`; return; }
    box.innerHTML = `<p class="muted">Studies with real ${esc(j.gene)} missense mutations (${j.checked} large studies checked):</p>` +
      j.studies.map(s => `<button type="button" class="btn btn-secondary btn-sm" data-study="${esc(s.id)}" style="margin:2px">${esc(s.name)} · ${s.n_missense} missense · ${s.samples} samples</button>`).join("");
    box.querySelectorAll("button[data-study]").forEach(b => b.addEventListener("click", () => { $("#p-study").value = b.dataset.study; thenRun(); }));
  } catch (e) { box.innerHTML = `<div class="notice warn">Could not search studies: ${esc(e.message)}</div>`; }
}

// "Not estimable" with the specific reason; a missing statistic is never shown as 0, 1 or a bare dash.
const neCell = (v, why, f = pfmt) => v == null
  ? `<span class="badge warn" title="${esc(why || "reason not reported")}">Not estimable ⓘ</span>` : f(v);
const kv = obj => Object.entries(obj || {}).map(([k, v]) => `${esc(k)}: <b>${v}</b>`).join("<br>") || "none";

let pxViewer = null, pxElement = null, pxModelKey = null, pxModelText = null;
async function drawProximityModel(j) {
  const st = $("#px3d-status"), el = $("#px3d");
  if (!st || !el) return;
  try {
    if (!window.PTMViewer || !window.$3Dmol) throw new Error("3Dmol.js did not load (local copy and CDN fallbacks both failed)");
    const key = `${up()}|${j.structure_id}`;
    if (pxModelKey !== key) {
      const r = await fetch(`/api/f/spatial/${encodeURIComponent(up())}/model?source=alphafold&id=${encodeURIComponent(j.structure_id)}`);
      if (!r.ok) throw new Error(`coordinate endpoint returned HTTP ${r.status}`);
      const t = await r.text();
      if (!/^ATOM  /m.test(t)) throw new Error("the response is not a PDB coordinate file");
      pxModelText = t; pxModelKey = key;
    }
    pxViewer = PTMViewer.mount(el); pxElement = el;               // bound to the element that is on the page now
    const nAtoms = PTMViewer.load(pxViewer, pxModelText);
    pxViewer.setStyle({}, { cartoon: { colorfunc: PTMViewer.plddtColor } });
    const chain = j.highlight.chain, rejected = [];
    let nP = 0, nM = 0;
    const mark = (m, colour, label) => {
      if (!PTMViewer.verified(pxViewer, chain, m.resi, m.aa)) { rejected.push(label); return false; }
      pxViewer.setStyle({ chain, resi: m.resi }, { cartoon: { color: colour }, stick: { radius: 0.25, color: colour } });
      const ca = pxViewer.selectedAtoms({ chain, resi: m.resi, atom: "CA" })[0];
      if (ca) pxViewer.addLabel(label, { position: { x: ca.x, y: ca.y, z: ca.z }, backgroundColor: colour, fontColor: "white", fontSize: 11 });
      return true;
    };
    for (const m of j.highlight.ptm) if (mark(m, "#18a6a6", m.site)) nP++;
    for (const m of j.highlight.mutations) if (mark(m, "#d64545", m.change || `${m.aa}${m.resi}`)) nM++;
    const fin = PTMViewer.finish(pxViewer, el);
    if (!fin.ok) throw new Error(fin.reason);
    st.innerHTML = `${esc(j.structure_label)} ${esc(j.structure_id)}, chain ${esc(chain)} · ${nAtoms} atoms · ` +
      `<span style="color:#7fe3df">${nP} verified PTM residue(s)</span>, <span style="color:#ff9a9a">${nM} qualifying mutated residue(s)</span>` +
      (rejected.length ? ` · not highlighted (identity not verified in the loaded model): ${esc(rejected.join(", "))}` : "") +
      ` · ${PTMViewer.PLDDT_LEGEND}`;
  } catch (e) {
    st.textContent = `3D structure could not be displayed: ${e.message}. No residue was highlighted.`;
  }
}
function focusResidue(pos) {
  if (!pxViewer || pxElement !== $("#px3d")) return;
  pxViewer.zoomTo({ resi: pos }, 600); pxViewer.render();
}

async function proximity() {
  loadStudies();
  const study = val("p-study", "");
  if (!study) { out.innerHTML = '<div class="notice">Select a cBioPortal study first. No cohort is assumed.</div>'; return; }
  busy("Retrieving mutations and mapping them onto the AlphaFold structure…");
  const qs = `study=${encodeURIComponent(study)}&threshold=${val("p-threshold", 10)}&min_sep=${val("p-minsep", 20)}` +
    `&dist_mode=${encodeURIComponent(val("p-distmode", "ca"))}&plddt_min=${val("p-plddt", 70)}&n_perm=${val("p-nperm", 5000)}&seed=${val("p-seed", 1)}`;
  const j = await api(`/api/f/proximity/${up()}?${qs}`);
  if (!j) return;
  const ck = j.cbio_check || {};
  const noMut = ck.n_records_raw === 0;
  const ex = j.exclusions || {};
  const exRows = [
    ...Object.entries(ex.missense_mutations || {}).map(([k, v]) => ["Missense mutation", esc(k), v]),
    ...(j.ptm_exclusions || []).map(e => ["PTM site", `${esc(e.site)}: ${esc(e.reason)}`, 1])
  ];
  out.innerHTML =
    (noMut ? `<div class="notice warn"><b>This study has no mutations in this gene</b> (0 records in ${ck.n_sequenced ?? "?"} sequenced samples). cBioPortal answered correctly; this is a property of the cohort, not an error and not evidence of no effect. The analysis is <b>not estimable</b>. <button type="button" id="find-studies" class="btn btn-primary btn-sm">Find studies where this gene is mutated</button><div id="study-suggest" style="margin-top:8px"></div></div>` : "") +
    retrievalHtml(j.retrieval) +
    `<div class="notice warn"><b>${esc(j.structure_label || "AlphaFold PREDICTED model")}</b>: ${esc(j.structure)}; chain ${esc(j.chain)}; mean pLDDT ${j.mean_plddt ?? "?"}. Distance = ${esc(j.distance_definition || "")} (≤ ${j.threshold} Å; sequence separation ≥ ${j.min_sep}). This is a prediction, not an experimental structure.</div>` +
    (j.note ? `<div class="notice warn">${esc(j.note)}</div>` : "") +
    (j.mapping_warning ? `<div class="notice warn">⚠ ${esc(j.mapping_warning)}</div>` : "") +
    (ck.n_records_raw != null ? `<div class="notice">cBioPortal returned <b>${ck.n_records_raw}</b> mutation records in <b>${ck.n_distinct_samples}</b> samples (of ${ck.n_sequenced ?? "?"} sequenced); types: ${esc(Object.entries(ck.by_type || {}).map(([k, v]) => k + " " + v).join(", ") || "none")}. <a href="${esc(ck.api_url || "#")}" target="_blank" rel="noopener">Open the raw API response</a> to compare.</div>` : "") +
    cards([
      [j.n_missense, "Missense records", `${j.n_missense_mapped} pass mapping checks`],
      [j.n_mutated_positions, "Mutated eligible positions", `${j.n_mutated_sequence_proximal ?? 0} sequence-proximal (not tested by design)`],
      [`${j.n_qualifying_residues}/${j.n_eligible_residues}`, "Qualifying / eligible residues", "null universe = eligible residues"],
      [j.n_observed, "3D-close, sequence-far", j.expected == null ? "expected: not estimable" : `expected ${j.expected}`],
      [j.fold_enrichment == null ? "Not estimable" : j.fold_enrichment.toFixed(2) + "×", "Fold enrichment",
        j.fold_enrichment == null ? (j.fold_reason || "") :
          (j.fold_ci95 ? `95% CI ${j.fold_ci95[0].toFixed(2)}–${j.fold_ci95[1].toFixed(2)} (approx.)` : `CI not estimable: ${j.ci_reason || ""}`)],
      [j.p_value == null ? "Not estimable" : pfmt(j.p_value), "Exact p (hypergeometric)", j.p_value == null ? (j.p_value_reason || "") : ""],
      [j.p_perm == null ? "Not estimable" : pfmt(j.p_perm), "Permutation p", j.p_perm == null ? (j.p_perm_note || "") : `${j.n_permutations} valid permutations, seed ${j.seed}`]
    ]) +
    (j.negative_control ? `<div class="notice">Synonymous/silent-mutation calibration control: <b>${j.n_synonymous_positions || 0}</b> distinct positions, <b>${j.negative_control.n_observed || 0}</b> 3D-close/sequence-far. ${j.negative_control_note ? esc(j.negative_control_note) : `Exact p: <b>${j.negative_control.p_value == null ? "not estimable" : pfmt(j.negative_control.p_value)}</b>.`} This is a positional control, not proof that synonymous mutations are neutral.</div>` : "") +
    `<h3>Mutations meeting both criteria</h3>` +
    table(["Position", "Residue", "Change(s)", "Mutations", "Nearest PTM", "3D dist (Å)", "Seq separation", "pLDDT"],
      (j.hits || []).map(h => [`<a href="#" class="hitlink" data-pos="${h.pos}"><b>${h.pos}</b></a>`, esc(h.residue || "—"), esc(h.protein_change.join(", ")), h.n_mutations, h.nearest_ptm_pos, h.dist, h.min_seq_sep, h.plddt]),
      j.n_positions_all ? "No verified missense mutation meets both criteria in this study (a valid negative result for this study, not an error)." : "No mutation could be tested (see the notes above).") +
    `<h3>3D structure</h3><div class="viewer-card small"><div id="px3d"></div>
      <button type="button" id="px3d-reset" class="btn btn-secondary btn-sm" style="position:absolute;right:10px;top:10px;z-index:2">Reset view</button>
      <div id="px3d-status" class="viewer-status">Loading predicted model…</div></div>` +
    `<h3>Excluded from the primary analysis</h3>` +
    table(["Kind", "Reason", "Count"], exRows, "Nothing was excluded.") +
    `<h3>Per-PTM-site 3D neighbourhoods</h3>` + table(["Site", "PTM", "Contacts ≤ radius", "Universe (≥ sep)", "Neighbourhood", "Mutated in universe", "Mutated in neighbourhood", "Expected", "Fold", "95% CI", "Exact p", "Permutation p", "q (exact)", "q (perm)", "Mutated residues"],
      (j.sites || []).map(x => [`<b>${esc(x.site)}</b>`, esc(x.ptm_type), x.n_contacts, x.n_universe, x.n_neighbourhood, x.n_mutated_in_universe, x.n_mutated_in_neighbourhood, x.expected == null ? "—" : x.expected,
        x.fold_enrichment == null ? "—" : Number(x.fold_enrichment).toFixed(2) + "×", x.fold_ci95 ? `${Number(x.fold_ci95[0]).toFixed(2)}–${Number(x.fold_ci95[1]).toFixed(2)}` : "—",
        neCell(x.p_value, x.p_value_reason), neCell(x.p_perm, x.p_perm_note), neCell(x.q_bh_protein, "needs a valid exact p"), neCell(x.q_perm_protein, "needs a valid permutation p"), esc(x.hit_positions.join(", ") || "—")]),
      "No verified PTM residue has confident coordinates.") +
    `<p class="caveat">“Not estimable” means the statistic cannot be validly computed (hover for the reason); it is never zero and never evidence of no effect. q-values are Benjamini–Hochberg over the valid p-values of this protein’s sites only. Permutation statistics come from a composition-matched null drawn from eligible residues; they are not additional biological observations.</p>` + method(j.method);
  const fb = $("#find-studies"); if (fb) fb.addEventListener("click", () => findMutatedStudies(run));
  out.querySelectorAll(".hitlink").forEach(a => a.addEventListener("click", e => { e.preventDefault(); focusResidue(Number(a.dataset.pos)); }));
  const rb = $("#px3d-reset"); if (rb) rb.addEventListener("click", () => { if (pxViewer && pxElement === $("#px3d")) { pxViewer.zoomTo(); pxViewer.render(); } });
  await drawProximityModel(j);
}
const RUN = {proximity, hotspots, crosstalk, network, conservation, mutations, spatial, hypotheses, priority, integrity};
async function run() { try { await RUN[F.page](); } catch (e) { out.innerHTML = `<div class="notice warn">Something went wrong: ${esc(e.message)}</div>`; } }
$("#run").addEventListener("click", run);
$("#protein-select") && $("#protein-select").addEventListener("change", () => {
  const selected = up();
  history.replaceState(null, "", "?uniprot=" + encodeURIComponent(selected));
  const viewerLink = $("#feature-viewer-link");
  if (viewerLink) viewerLink.href = "/viewer/" + encodeURIComponent(selected);
});
if (!F.dbEmpty) run();
})();
