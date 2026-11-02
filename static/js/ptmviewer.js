/* Shared 3Dmol helper for the Analysis Suite viewers.
 *
 * Why this exists: the pages rebuild their result area with innerHTML on every run. A 3Dmol viewer that was
 * created once and cached in a module variable stays attached to the OLD, now-detached <div>, so later runs
 * drew into nothing: the status line said "structure loaded" while the visible area stayed blank.
 * mount() therefore keys viewers by the live DOM element and never reuses one for a different element, and
 * finish() verifies that a canvas of non-zero size really exists before anything is reported as displayed.
 */
(function (root) {
  "use strict";
  const viewers = new WeakMap();      // DOM element -> 3Dmol viewer created for exactly that element

  function lib() {
    const L = root.$3Dmol;
    if (!L) throw new Error("3Dmol.js is not loaded (local copy and CDN fallbacks both failed)");
    return L;
  }

  /** Viewer bound to this exact, currently attached element. */
  function mount(element) {
    if (!element) throw new Error("3D viewer container is missing");
    if (element.isConnected === false) throw new Error("3D viewer container is not attached to the page");
    let v = viewers.get(element);
    if (!v) {
      element.textContent = "";
      v = lib().createViewer(element, { backgroundColor: "white" });
      viewers.set(element, v);
    }
    return v;
  }

  /** Replace the model; throws if the text does not parse to atoms (a blank scene is never "loaded"). */
  function load(viewer, pdbText) {
    viewer.removeAllModels();
    viewer.removeAllShapes();
    viewer.removeAllLabels();
    viewer.addModel(pdbText, "pdb");
    const n = (viewer.selectedAtoms({}) || []).length;
    if (!n) throw new Error("the coordinate file parsed to 0 atoms");
    return n;
  }

  /** Fit, resize, draw, then confirm something drawable exists. Returns {ok, reason}. */
function finish(viewer, element, selection) {
    viewer.zoomTo(selection || {});
    if (typeof viewer.resize === "function") viewer.resize();
    viewer.render();
    const w = element.clientWidth, h = element.clientHeight;
    if (w === 0 || h === 0) return { ok: false, reason: "the viewer container has zero size (hidden or collapsed panel)" };
    const canvas = element.querySelector ? element.querySelector("canvas") : true;
    if (!canvas) return { ok: false, reason: "WebGL canvas was not created (browser/GPU WebGL may be disabled)" };
    return { ok: true };
  }

  const T2O = { ALA: "A", ARG: "R", ASN: "N", ASP: "D", CYS: "C", GLN: "Q", GLU: "E", GLY: "G", HIS: "H", ILE: "I",
    LEU: "L", LYS: "K", MET: "M", PHE: "F", PRO: "P", SER: "S", THR: "T", TRP: "W", TYR: "Y", VAL: "V",
    SEP: "S", TPO: "T", PTR: "Y", MSE: "M", ALY: "K", M3L: "K", MLY: "K", MLZ: "K", KCX: "K", CSO: "C", HYP: "P" };

  /** One-letter identity of the residue actually loaded at chain:resi (null if absent or ambiguous). */
  function residueLetter(viewer, chain, resi) {
    const sel = chain ? { chain, resi } : { resi };
    const names = [...new Set((viewer.selectedAtoms(sel) || []).map(a => String(a.resn || "").toUpperCase()).filter(Boolean))];
    return names.length === 1 ? (T2O[names[0]] || null) : null;
  }
  /** True only if the loaded model really has the expected amino acid there. Unverified residues are never highlighted. */
  function verified(viewer, chain, resi, expected) {
    return !!expected && residueLetter(viewer, chain, resi) === expected;
  }

  const plddtBand = b => b >= 90 ? "very high (>=90)" : b >= 70 ? "confident (70-90)" : b >= 50 ? "low (50-70)" : "very low (<50)";
  const plddtColor = a => a.b >= 90 ? "#0053D6" : a.b >= 70 ? "#65CBF3" : a.b >= 50 ? "#FFDB13" : "#FF7D45";

  const PLDDT_LEGEND = '<span style="display:inline-flex;gap:8px;flex-wrap:wrap;font-size:11px">' +
    [["#0053D6", "pLDDT ≥ 90"], ["#65CBF3", "70–90"], ["#FFDB13", "50–70"], ["#FF7D45", "< 50"]]
      .map(([c, t]) => `<span><i style="display:inline-block;width:10px;height:10px;background:${c};margin-right:4px;border-radius:2px"></i>${t}</span>`).join("") +
    "</span>";

  const api = { mount, load, finish, residueLetter, verified, plddtBand, plddtColor, PLDDT_LEGEND };
  root.PTMViewer = api;
  if (typeof module !== "undefined" && module.exports) module.exports = api;
})(typeof window !== "undefined" ? window : globalThis);
