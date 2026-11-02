// Runs the REAL static/js/features.js (+ ptmviewer.js when present) in jsdom with a fake 3Dmol, and checks that
// after TWO analysis runs the viewer that draws is bound to the element currently on the page.
// usage: node spatial_viewer_rebind.test.js <features.js> [<ptmviewer.js>]
const { JSDOM } = require(process.env.JSDOM_PATH || "jsdom");
const fs = require("fs");
const [featuresPath, helperPath] = process.argv.slice(2);

const PDB = ["ATOM      1  CA  SER A   1       0.000   0.000   0.000  1.00 90.00           C",
             "ATOM      2  CA  ALA A  50       5.000   0.000   0.000  1.00 90.00           C"].join("\n");
const SPATIAL = { protein: { gene: "G" }, available: true, source: "alphafold", structure_id: "AF-X-F1", chain: "A",
  threshold: 10, n_mapped: 2, n_units: 2, unmapped: [], notes: [], basis: "b", note: "n",
  coverage: { experimental: { mapped: 0, total: 2, sites: [] }, alphafold: { mapped: 2, total: 2, sites: ["S1", "S50"] } },
  mapping: [{ site: "S1", mapping_status: "mapped", structure_residue_number: 1, chain: "A" },
            { site: "A50", mapping_status: "mapped", structure_residue_number: 50, chain: "A" }],
  pairs: [{ site_a: "S1", site_b: "A50", seq_sep: 49, ca_dist: 5, min_atom_dist: 4, both_confident: true,
            resi_a: 1, resi_b: 50, cross_type: false, long_range_close: true }], n_close: 1, n_long_range_close: 1 };

(async () => {
  const dom = new JSDOM(`<body><select id="protein-select"><option value="P1" selected>x</option></select>
    <input id="p-thr" value="10"><input id="p-af" type="checkbox" checked><button id="run"></button><div id="result"></div></body>`, { runScripts: "outside-only", pretendToBeVisual: true });
  const w = dom.window;
  Object.defineProperty(w.HTMLElement.prototype, "clientWidth", { get() { return 600; } });
  Object.defineProperty(w.HTMLElement.prototype, "clientHeight", { get() { return 480; } });
  w.FEATURE = { page: "spatial", selected: "P1", dbEmpty: true };       // dbEmpty: do not auto-run; we click ourselves
  w.history.replaceState = () => {};
  w.fetch = async url => ({ ok: true, status: 200,
    json: async () => SPATIAL, text: async () => PDB });
  w.AbortController = AbortController;

  const created = [];                                                    // every viewer + the element it was bound to
  w.$3Dmol = { createViewer(el) {
    const v = { el, atoms: 0, model: false,
      removeAllModels() { this.model = false; }, removeAllShapes() {}, removeAllLabels() {},
      addModel() { this.model = true; this.atoms = 2; }, setStyle() {}, zoomTo() {}, resize() {},
      selectedAtoms(sel) { return this.model ? [{ resn: "SER", x: 0, y: 0, z: 0 }] : []; },
      addCylinder() {}, addLabel() {},
      render() { if (this.model && !this.el.querySelector("canvas")) this.el.appendChild(this.el.ownerDocument.createElement("canvas")); } };
    created.push(v); return v; } };
  if (helperPath) w.eval(fs.readFileSync(helperPath, "utf8"));
  w.eval(fs.readFileSync(featuresPath, "utf8"));

  const run = async () => { w.document.querySelector("#run").click(); await new Promise(r => setTimeout(r, 150)); };
  await run();                                                           // 1st run
  await run();                                                           // 2nd run: result area rebuilt via innerHTML
  const live = w.document.querySelector("#v3d");
  const drawn = !!(live && live.querySelector("canvas"));
  const status = (w.document.querySelector("#v3d-status") || {}).textContent;
  console.log(JSON.stringify({ viewersCreated: created.length, liveContainerHasCanvas: drawn, status: (status || "").slice(0, 90) }));
  process.exit(drawn ? 0 : 1);
})().catch(e => { console.error("HARNESS ERROR", e); process.exit(2); });
