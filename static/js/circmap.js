
(() => {
const DATA = window.PROTEOMIC_MAP;
if (!DATA || !DATA.length) return;

const PTM_COLORS = {
  Phosphorylation: "#18a6a6", Acetylation: "#e8871e", Ubiquitination: "#8a5cf6",
  Methylation: "#3aa8e6", SUMOylation: "#e0507a", Glycosylation: "#5cb85c",
};
const PTM_DEFAULT = "#657b80";
const DOMAIN_PALETTE = ["#3aa8e6","#e8871e","#5cb85c","#e0c23a","#8a5cf6","#e0507a"];
const CONF_RANK = { High: 3, Medium: 2, Low: 1 };
const CONF_COLOR = { High: "#1f9d55", Medium: "#e0a020", Low: "#8a97a0" };

const svg = document.getElementById("pm-svg");
const linSvg = document.getElementById("pm-linear-svg");
const legendEl = document.getElementById("pm-domain-legend");
const selPanel = document.getElementById("pm-selected-panel");
const kinPanel = document.getElementById("pm-kinase-panel");
const diseasePanel = document.getElementById("pm-disease-panel");
const summaryEl = document.getElementById("pm-summary");

const CX = 240, CY = 240, R_BACKBONE = 160, R_DOMAIN_OUT = 205, R_DOMAIN_IN = 190, R_TICK = 225;
let domains = []; // real UniProt features only, or [] if none found
let selectedId = null;

function ptmColor(type){ return PTM_COLORS[type] || PTM_DEFAULT; }
function escapeHtml(s){ return String(s==null?"":s).replace(/[&<>"']/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c])); }

function angleFor(pos){ return (pos / DATA.length) * 2*Math.PI - Math.PI/2; }
function pointOnRing(pos, r){
  const a = angleFor(pos);
  return { x: CX + r*Math.cos(a), y: CY + r*Math.sin(a), a };
}
function donutSlicePath(startFrac, endFrac, rOuter, rInner){
  const a0 = startFrac*2*Math.PI - Math.PI/2, a1 = endFrac*2*Math.PI - Math.PI/2;
  const large = (endFrac-startFrac) > 0.5 ? 1 : 0;
  const p = (r,a)=>({x:CX+r*Math.cos(a), y:CY+r*Math.sin(a)});
  const A=p(rOuter,a0), B=p(rOuter,a1), C=p(rInner,a1), D=p(rInner,a0);
  return `M ${A.x} ${A.y} A ${rOuter} ${rOuter} 0 ${large} 1 ${B.x} ${B.y} L ${C.x} ${C.y} A ${rInner} ${rInner} 0 ${large} 0 ${D.x} ${D.y} Z`;
}
function niceTickStep(length){
  const steps=[5,10,20,25,50,100,200,250,500,1000,2000,5000];
  for(const s of steps){ if(length/s <= 9) return s; }
  return Math.ceil(length/9);
}
function svgEl(tag, attrs){
  const el=document.createElementNS("http://www.w3.org/2000/svg",tag);
  for(const k in attrs) el.setAttribute(k, attrs[k]);
  return el;
}

// --- Real domain annotations, fetched live from UniProt. Never fabricated:
// if the fetch fails or returns nothing usable, the ring just stays a plain
// backbone with no domain segments and no legend.
async function fetchDomains(){
  try{
    const r = await fetch(`https://rest.uniprot.org/uniprotkb/${DATA.uniprot}.json`);
    if(!r.ok) return [];
    const j = await r.json();
    const feats = (j.features||[]).filter(f => f.type === "Domain");
    const source = feats.length ? feats : (j.features||[]).filter(f => f.type === "Region");
    const parsed = source.map(f => ({
      start: f.location && f.location.start && f.location.start.value,
      end: f.location && f.location.end && f.location.end.value,
      description: f.description || f.type,
    })).filter(d => d.start && d.end && d.end > d.start);
    parsed.sort((a,b)=> (b.end-b.start) - (a.end-a.start));
    const kept=[];
    for(const d of parsed){
      const containedInExisting = kept.some(k => d.start>=k.start && d.end<=k.end);
      if(!containedInExisting) kept.push(d);
      if(kept.length>=6) break;
    }
    return kept.map((d,i)=>({...d, color: DOMAIN_PALETTE[i % DOMAIN_PALETTE.length]}));
  }catch(e){ return []; }
}

function renderCircle(){
  svg.innerHTML = "";
  // Backbone PTM ring
  svg.appendChild(svgEl("circle",{cx:CX,cy:CY,r:R_BACKBONE,fill:"none",stroke:"#dce7e8","stroke-width":10}));

  // Domain ring: gray "Other" base (a full stroked circle, not a 360-degree
  // donut path -- an arc whose start and end angles coincide degenerates to
  // nothing in SVG), then real domain segments drawn on top as donut slices
  const domRingMid = (R_DOMAIN_OUT + R_DOMAIN_IN) / 2, domRingWidth = R_DOMAIN_OUT - R_DOMAIN_IN;
  svg.appendChild(svgEl("circle",{cx:CX,cy:CY,r:domRingMid,fill:"none",stroke:"#e3ebec","stroke-width":domRingWidth}));
  domains.forEach(d=>{
    const path = svgEl("path",{d:donutSlicePath(d.start/DATA.length, d.end/DATA.length, R_DOMAIN_OUT, R_DOMAIN_IN), fill:d.color});
    const title = svgEl("title",{}); title.textContent = `${d.description} (${d.start}-${d.end})`;
    path.appendChild(title);
    svg.appendChild(path);
  });

  // Tick labels
  const step = niceTickStep(DATA.length);
  for(let pos=0; pos<DATA.length; pos+=step){
    const {x,y} = pointOnRing(Math.max(pos,1), R_TICK);
    const t = svgEl("text",{x,y,"text-anchor":"middle","dominant-baseline":"middle","font-size":10,fill:"#8aa0a4"});
    t.textContent = pos; svg.appendChild(t);
  }

  // PTM markers
  DATA.ptms.forEach(p=>{
    if(!p.position) return;
    const {x,y} = pointOnRing(p.position, R_BACKBONE);
    const isSel = p.id === selectedId;
    const c = svgEl("circle",{cx:x,cy:y,r:isSel?10:7,fill:ptmColor(p.ptm_type),stroke:"white","stroke-width":isSel?2.5:1.5,cursor:"pointer"});
    c.addEventListener("click", ()=>selectSite(p.id));
    const title = svgEl("title",{}); title.textContent = `${p.site} — ${p.ptm_type}`;
    c.appendChild(title);
    svg.appendChild(c);
    if(isSel){
      const lbl = svgEl("text",{x, y:y - 14, "text-anchor":"middle","font-size":11,"font-weight":800,fill:"#17363a"});
      lbl.textContent = p.site; svg.appendChild(lbl);
    }
  });

  // Center label
  const g1=svgEl("text",{x:CX,y:CY-6,"text-anchor":"middle","font-size":19,"font-weight":800,fill:"#17363a"}); g1.textContent=DATA.gene;
  const g2=svgEl("text",{x:CX,y:CY+14,"text-anchor":"middle","font-size":11,fill:"#657b80"}); g2.textContent=`${DATA.uniprot} \u00b7 ${DATA.length.toLocaleString()} aa`;
  svg.appendChild(g1); svg.appendChild(g2);

  // Legend
  legendEl.innerHTML = "";
  if(domains.length){
    domains.forEach(d=>{
      const span=document.createElement("span"); span.className="legend-item";
      span.innerHTML = `<i style="background:${d.color}"></i>${escapeHtml(d.description)}`;
      legendEl.appendChild(span);
    });
    const other=document.createElement("span"); other.className="legend-item";
    other.innerHTML = `<i style="background:#e3ebec"></i>Other`;
    legendEl.appendChild(other);
  }
}

function renderLinear(){
  linSvg.innerHTML="";
  const x0=30, usable=840;
  const xFor = pos => x0 + (pos/DATA.length)*usable;
  linSvg.appendChild(svgEl("rect",{x:x0,y:28,width:usable,height:14,rx:7,fill:"#e3ebec"}));
  domains.forEach(d=>{
    linSvg.appendChild(svgEl("rect",{x:xFor(d.start),y:28,width:Math.max(xFor(d.end)-xFor(d.start),1),height:14,fill:d.color}));
  });
  const step = niceTickStep(DATA.length);
  for(let pos=0; pos<=DATA.length; pos+=step){
    linSvg.appendChild(svgEl("text",{x:xFor(pos),y:58,"text-anchor":"middle","font-size":9,fill:"#8aa0a4"})).textContent=pos;
  }
  DATA.ptms.forEach(p=>{
    if(!p.position) return;
    const x=xFor(p.position); const isSel=p.id===selectedId;
    const c=svgEl("circle",{cx:x,cy:22,r:isSel?6:4,fill:ptmColor(p.ptm_type),stroke:"white","stroke-width":1.2,cursor:"pointer"});
    c.addEventListener("click",()=>selectSite(p.id));
    const title=svgEl("title",{}); title.textContent=`${p.site} — ${p.ptm_type}`; c.appendChild(title);
    linSvg.appendChild(c);
  });
}

function renderSelected(){
  const p = DATA.ptms.find(x=>x.id===selectedId);
  if(!p){ selPanel.innerHTML = `<h4>Selected PTM Site</h4><div class="pm-empty">Click a PTM marker to see details.</div>`; kinPanel.innerHTML = `<h4>Kinase Relationship</h4><div class="pm-empty">Select a PTM site first.</div>`; return; }
  const rows = [
    ["Site", `<b>${escapeHtml(p.site)}</b>`],
    ["Type", escapeHtml(p.ptm_type)],
    ["Kinase", escapeHtml(p.kinase || "—")],
    ["Evidence", escapeHtml(p.evidence || "—")],
    ["Confidence", escapeHtml(p.confidence || "—")],
    ["Cancer Association", escapeHtml(p.disease_context || "—")],
    ["PMID", p.pmid ? `<a href="https://pubmed.ncbi.nlm.nih.gov/${encodeURIComponent(p.pmid)}/" target="_blank">${escapeHtml(p.pmid)}</a>` : "—"],
    ["Position", p.position ?? "—"],
  ];
  selPanel.innerHTML = `<h4>Selected PTM Site</h4><table class="pm-kv">${rows.map(([k,v])=>`<tr><td>${k}</td><td>${v}</td></tr>`).join("")}</table>`;

  const hasKinase = p.kinase && p.kinase.trim() && p.kinase.trim() !== "—";
  if(hasKinase){
    kinPanel.innerHTML = `<h4>Kinase Relationship</h4>
      <div class="pm-flow">
        <div class="pm-flow-box pm-flow-kinase">${escapeHtml(p.kinase)}</div>
        <div class="pm-flow-arrow">\u2193</div>
        <div class="pm-flow-box pm-flow-protein">${escapeHtml(DATA.gene)} (${escapeHtml(DATA.uniprot)})</div>
        <div class="pm-flow-arrow">\u2193</div>
        <div class="pm-flow-box pm-flow-site">${escapeHtml(p.site)} \u2014 ${escapeHtml(p.ptm_type)}</div>
      </div>`;
  } else {
    kinPanel.innerHTML = `<h4>Kinase Relationship</h4><div class="pm-empty">No known upstream kinase recorded for this site.</div>`;
  }
}

function selectSite(id){
  selectedId = id;
  renderCircle(); renderLinear(); renderSelected();
}

function renderDiseasePanel(){
  const counts = new Map(); // disease -> {n, bestConf}
  DATA.ptms.forEach(p=>{
    if(!p.disease_context) return;
    p.disease_context.split(";").map(d=>d.trim()).filter(Boolean).forEach(d=>{
      const rank = CONF_RANK[p.confidence] || 0;
      const cur = counts.get(d) || {n:0, rank:0};
      counts.set(d, {n: cur.n+1, rank: Math.max(cur.rank, rank)});
    });
  });
  const list = [...counts.entries()].map(([d,v])=>({disease:d, n:v.n, rank:v.rank}))
    .sort((a,b)=> b.n-a.n);
  if(!list.length){ diseasePanel.innerHTML = `<h4>Cancer / Disease Associations</h4><div class="pm-empty">No disease associations recorded.</div>`; return; }
  const rankLabel = r => r===3?"High":r===2?"Medium":r===1?"Low":"—";
  const rankColor = r => CONF_COLOR[rankLabel(r)] || "#8a97a0";
  const shown = list.slice(0,4);
  const rest = list.slice(4);
  let html = `<h4>Cancer / Disease Associations</h4>` + shown.map(d=>
    `<div class="pm-disease-row"><span>${escapeHtml(d.disease)}</span><span class="pm-badge" style="background:${rankColor(d.rank)}22;color:${rankColor(d.rank)}">${rankLabel(d.rank)}</span></div>`
  ).join("");
  if(rest.length){
    html += `<button class="pm-viewall" id="pm-viewall-diseases">View all (${list.length})</button>
      <div class="pm-disease-more" id="pm-disease-more" style="display:none">${rest.map(d=>
        `<div class="pm-disease-row"><span>${escapeHtml(d.disease)}</span><span class="pm-badge" style="background:${rankColor(d.rank)}22;color:${rankColor(d.rank)}">${rankLabel(d.rank)}</span></div>`
      ).join("")}</div>`;
  }
  diseasePanel.innerHTML = html;
  const btn = document.getElementById("pm-viewall-diseases");
  if(btn) btn.addEventListener("click", ()=>{
    const more = document.getElementById("pm-disease-more");
    const open = more.style.display !== "none";
    more.style.display = open ? "none" : "block";
    btn.textContent = open ? `View all (${list.length})` : "Show less";
  });
}

function renderSummary(){
  const total = DATA.ptms.length;
  const byType = {};
  DATA.ptms.forEach(p=>{ byType[p.ptm_type] = (byType[p.ptm_type]||0)+1; });
  const kinases = new Set(DATA.ptms.filter(p=>p.kinase && p.kinase.trim() && p.kinase.trim()!=="—").map(p=>p.kinase.trim()));
  const diseases = new Set();
  DATA.ptms.forEach(p=>{ if(p.disease_context) p.disease_context.split(";").forEach(d=>{ d=d.trim(); if(d) diseases.add(d); }); });
  const withEvidence = DATA.ptms.filter(p=>p.evidence && p.evidence.trim()).length;
  const evidencePct = total ? Math.round((withEvidence/total)*100) : 0;
  const confCounts = {High:0, Medium:0, Low:0};
  DATA.ptms.forEach(p=>{ if(confCounts[p.confidence] !== undefined) confCounts[p.confidence]++; });

  const typeChips = Object.entries(byType).map(([t,n])=>
    `<span class="pm-type-chip" style="color:${ptmColor(t)}">${n}<small>${escapeHtml(t.slice(0,6))}</small></span>`).join("");

  const kinaseList = [...kinases]; const diseaseList = [...diseases];
  const shortList = (arr,max=3) => arr.length<=max ? arr.map(escapeHtml).join(", ") : arr.slice(0,max).map(escapeHtml).join(", ") + ` + ${arr.length-max} more`;

  summaryEl.innerHTML = `
    <div class="pm-sum-card"><div class="pm-sum-value">${total}</div><div class="pm-sum-label">Total PTM Sites</div><div class="pm-sum-chips">${typeChips}</div></div>
    <div class="pm-sum-card"><div class="pm-sum-value">${kinases.size}</div><div class="pm-sum-label">Unique Kinases</div><div class="pm-sum-note">${shortList(kinaseList) || "None recorded"}</div></div>
    <div class="pm-sum-card"><div class="pm-sum-value">${diseases.size}</div><div class="pm-sum-label">Diseases</div><div class="pm-sum-note">${shortList(diseaseList) || "None recorded"}</div></div>
    <div class="pm-sum-card"><div class="pm-sum-value" style="color:#1f9d55">${evidencePct}%</div><div class="pm-sum-label">Sites with Evidence</div><div class="pm-sum-note">${withEvidence}/${total} recorded</div></div>
    <div class="pm-sum-card"><div class="pm-sum-label" style="margin-bottom:8px">Confidence</div>${["High","Medium","Low"].map(l=>
      `<div class="pm-conf-row"><span>${l}</span><span>${confCounts[l]}</span></div>`).join("")}</div>`;
}

function exportPNG(){
  const xml = new XMLSerializer().serializeToString(svg);
  const svg64 = btoa(unescape(encodeURIComponent(xml)));
  const img = new Image();
  img.onload = ()=>{
    const canvas=document.createElement("canvas");
    canvas.width=960; canvas.height=960;
    const ctx=canvas.getContext("2d");
    ctx.fillStyle="#ffffff"; ctx.fillRect(0,0,canvas.width,canvas.height);
    ctx.drawImage(img,0,0,canvas.width,canvas.height);
    const a=document.createElement("a");
    a.download = `${DATA.gene || "protein"}_circular_map.png`;
    a.href = canvas.toDataURL("image/png");
    a.click();
  };
  img.src = "data:image/svg+xml;base64," + svg64;
}

async function init(){
  renderCircle(); renderLinear(); renderSelected(); renderDiseasePanel(); renderSummary();
  const exportBtn = document.getElementById("pm-export");
  if(exportBtn) exportBtn.addEventListener("click", exportPNG);

  domains = await fetchDomains();
  renderCircle(); renderLinear(); // re-render once real domain data (or lack of it) is known
}
init();
})();
