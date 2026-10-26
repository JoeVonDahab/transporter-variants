/* Gene page: header, track switch, protein map, variant table. */
"use strict";

(() => {
  const { esc, fmt, fmtLR, evChip, EV_LABEL, CONS_LABEL, getJSON, proteinKey, el } = TV;
  const P = new URLSearchParams(location.search);
  const GENE = (P.get("g") || "").trim();
  const app = document.getElementById("app");
  const PAGE = 50;

  const S = {
    track: P.get("t") === "routed" ? "routed" : "native",
    idx: null, meta: null, row: null, rows: [], placed: [], unplaced: 0,
    cons: new Set(), ev: new Set(), excludeTraining: false,
    text: P.get("q") || "", pos: P.get("pos") ? Number(P.get("pos")) : null,
    sort: { key: "ev", dir: "asc" }, page: 0, open: new Set(),
    metric: "patho", wave: null, mol: null,
  };

  const TRACK = {
    native: {
      name: "Our model only", file: "native", license: "CC BY 4.0",
      note: n => `Scores from our model alone. Free for any use, including commercial, if you credit the source (CC BY 4.0).`,
    },
    routed: {
      name: "With AlphaMissense", file: "routed", license: "CC BY-NC-SA 4.0",
      note: n => `Our model combined with AlphaMissense for the ${fmt(n)} variants here that AlphaMissense covers, and our model alone for the rest. Non-commercial use only (CC BY-NC-SA 4.0).`,
    },
  };

  function setURL() {
    const u = new URLSearchParams();
    u.set("g", GENE);
    if (S.track === "routed") u.set("t", "routed");
    if (S.text) u.set("q", S.text);
    if (S.pos) u.set("pos", S.pos);
    history.replaceState(null, "", `gene.html?${u}`);
  }

  function decode(shard, idx) {
    const pre = shard.hgvs_prefixes, cs = idx.clinsig, rv = idx.review, routed = shard.track === "routed";
    return shard.rows.map(r => {
      const o = {
        vid: r[0], rs: r[1], hp: r[2] >= 0 ? pre[r[2]] : "", hs: r[3], cons: r[4], pv: r[5],
        clinsig: r[6] == null ? "" : cs[r[6]], review: r[7] == null ? "" : rv[r[7]],
        train: r[8] === 1, ev: r[9], score: r[10], lr: r[11], act: r[12], mech: r[13],
        lis: r[14] === 1, pos: r[15], cls: r[16] === 1,
      };
      if (routed) { o.plus = r[17] === 1; o.hasAM = r[18] === 1; o.v2n = r[19]; o.v2p = r[20]; o.am = r[21]; }
      o.hgvs = o.hp + o.hs;
      o.mclip = o.mech == null ? null : Math.max(0, Math.min(1, o.mech));
      return o;
    });
  }

  const mechWord = m => m == null ? "" : m >= 0.6 ? "leans activity-specific" : m <= 0.4 ? "leans expression-driven" : "unclear";

  /* ---------- skeleton ---------- */
  function renderShell() {
    const m = S.meta, r = S.row;
    const nTM = TVMap.topology(m).nTM;
    document.title = `${GENE}: ${m.name || "variants"} | Transporter Variants`;
    app.innerHTML = `
      <header class="ghead wrap">
        <div class="fam">${esc(familyName(m.family))}, UniProt <a href="https://www.uniprot.org/uniprotkb/${esc(m.uniprot)}/entry" rel="noopener">${esc(m.uniprot)}</a></div>
        <h1>${esc(GENE)}</h1>
        <div class="pname">${esc(m.name)}</div>
        <div class="facts">
          <span><b>${fmt(r[4])}</b> ClinVar variants</span>
          <span><b>${fmt(r[5])}</b> missense</span>
          <span><b>${fmt(m.length)}</b> amino acids${nTM ? `, <b>${nTM}</b> membrane helices` : ""}</span>
          <span><b>${fmt(r[8])}</b> seen in training</span>
          <span><a href="https://www.ncbi.nlm.nih.gov/clinvar/?term=${encodeURIComponent(GENE)}%5Bgene%5D" rel="noopener">ClinVar</a>, <a href="https://alphafold.ebi.ac.uk/entry/${esc(m.uniprot)}" rel="noopener">AlphaFold DB</a></span>
        </div>
        <div class="trackbar" role="group" aria-label="Data version">
          <div class="seg" id="trackseg">
            <button type="button" data-t="native">Our model only<small>CC BY 4.0</small></button>
            <button type="button" data-t="routed">With AlphaMissense<small>non-commercial</small></button>
          </div>
          <p class="tnote" id="tnote"></p>
        </div>
        <p class="evbar-h">How this gene's ${fmt(r[4])} variants split across evidence levels</p>
        <div class="evbar" id="evbar" aria-hidden="true"></div>
        <div class="evlegend" id="evlegend"></div>
      </header>

      <section class="panel wrap" aria-labelledby="map-h">
        <div class="panel-head">
          <h2 id="map-h">Protein map</h2>
          <p>The protein chain drawn across the cell membrane. Red stretches are where ClinVar missense variants tend to be damaging, blue where they tend to be tolerated, and grey where there are too few variants to say. <a href="help.html#map">How to read the map</a></p>
        </div>
        <div class="mapbox">
          <div class="map-controls">
            <label for="metricseg">Colour by</label>
            <div class="seg" id="metricseg">
              <button type="button" data-m="patho">Pathogenicity</button>
              <button type="button" data-m="act">Activity</button>
              <button type="button" data-m="mech">Mechanism<span class="tag-exp">exploratory</span></button>
              <button type="button" data-m="topo">Membrane topology</button>
              <button type="button" data-m="plddt">AlphaFold confidence</button>
            </div>
            <label class="check"><input type="checkbox" id="showv" checked> Show variants</label>
            <button type="button" class="btn ghost" id="unzoom" hidden>Show whole protein</button>
          </div>
          <div id="metricnote"></div>
          <div id="wave"></div>
          <div class="map-foot">
            <span id="legend"></span>
            <span id="symlegend"></span>
          </div>
          <div class="map-split">
            <div class="mol" id="mol">
              <div class="viewer"></div>
              <div class="molmsg"></div>
              <button type="button" class="btn ghost molreset" hidden>Reset view</button>
              <div class="molhint">Drag to rotate. Click a residue to select it.<span class="more"> Click the structure, then scroll to zoom. Colours are averaged over residues within about 10 Å of each other.</span></div>
            </div>
            <aside class="resinfo" id="resinfo" aria-live="polite"></aside>
          </div>
        </div>
        <p class="small muted" id="mapnote" style="margin-top:8px">Drag across the map to zoom in. Double-click it, or press “Show whole protein”, to zoom out.</p>
      </section>

      <section class="panel wrap" aria-labelledby="tab-h" style="margin-bottom:70px">
        <div class="panel-head">
          <h2 id="tab-h">Variants</h2>
          <p>Use the evidence level when classifying a variant. The raw score is there for reference; a cut-off of 0.5 on it calls too many variants pathogenic (<a href="help.html#score">why</a>). Click a row for details.</p>
        </div>
        <div class="filters" id="filters"></div>
        <div class="tablewrap"><table class="vt" id="vt"></table></div>
        <div class="pager" id="pager"></div>
      </section>`;

    app.querySelector("#trackseg").addEventListener("click", e => {
      const b = e.target.closest("button[data-t]");
      if (b && b.dataset.t !== S.track) switchTrack(b.dataset.t);
    });
    app.querySelector("#metricseg").addEventListener("click", e => {
      const b = e.target.closest("button[data-m]");
      if (b) setMetric(b.dataset.m);
    });
    app.querySelector("#showv").addEventListener("change", e => {
      S.wave.showVariants = e.target.checked; S.wave.render();
      if (S.mol) S.mol.paint();
    });
    app.querySelector("#unzoom").addEventListener("click", () => S.wave.zoom(null));
  }

  function familyName(f) {
    return { "ion channel": "Ion channel", SLC: "Solute carrier (SLC)", ABC: "ABC transporter", "P-type ATPase": "P-type ATPase",
      CFTR: "CFTR", "connexin/pannexin": "Connexin or pannexin", "V-type ATPase": "V-type ATPase", aquaporin: "Aquaporin",
      unassigned: "Family not assigned" }[f] || f;
  }

  /* ---------- track ---------- */
  async function loadTrack() {
    const shard = await getJSON(`data/${S.track}/${GENE}.json`);
    S.rows = decode(shard, S.idx);
    const seq = S.meta.sequence || "";
    S.placed = []; S.unplaced = 0;
    for (const o of S.rows) {
      if (o.cons !== 0 || !o.pos) continue;
      if (seq[o.pos - 1] === o.pv[0]) S.placed.push(o); else S.unplaced++;
    }
  }

  function renderTrackBar() {
    app.querySelectorAll("#trackseg button").forEach(b => b.setAttribute("aria-pressed", String(b.dataset.t === S.track)));
    const nPlus = S.rows.filter(o => o.plus).length;
    app.querySelector("#tnote").textContent = TRACK[S.track].note(nPlus);
    const counts = new Array(8).fill(0);
    for (const o of S.rows) counts[o.ev]++;
    const n = S.rows.length;
    app.querySelector("#evbar").innerHTML = counts.map((c, i) => c ? `<i style="width:${(c / n * 100).toFixed(3)}%;background:var(--${TV.EV_KEYS[i]})"></i>` : "").join("");
    app.querySelector("#evlegend").innerHTML = counts.map((c, i) => c ? `<span class="ev" style="--c:var(--${TV.EV_KEYS[i]})">${EV_LABEL[i]} ${fmt(c)}</span>` : "").join("");
  }

  async function switchTrack(t) {
    S.track = t; S.page = 0; S.open.clear();
    setURL();
    await loadTrack();
    renderTrackBar();
    S.wave.setVariants(S.placed);
    S.wave.render();
    if (S.mol && S.mol.viewer) {
      S.mol.spatial = {};
      if (S.mol.mapped) {
        for (const [k, m] of Object.entries(TVMap.METRICS)) {
          const pts = S.placed.map(v => [v.pos, m.value(v)]).filter(x => x[1] != null && !Number.isNaN(x[1]));
          S.mol.spatial[k] = TVMap.spaceProfile(S.mol.ca, pts);
        }
      }
      S.mol.paint();
    }
    renderResInfo();
    renderFilters();
    renderTable();
  }

  /* ---------- map ---------- */
  function setMetric(m) {
    S.metric = m;
    app.querySelectorAll("#metricseg button").forEach(b => b.setAttribute("aria-pressed", String(b.dataset.m === m)));
    S.wave.setMetric(m);
    if (S.mol) S.mol.setMetric(m);
    app.querySelector("#legend").innerHTML = TVMap.legendRamp(m);
    const note = app.querySelector("#metricnote");
    if (m === "mech") {
      note.innerHTML = `<div class="notice" style="margin:12px 14px 0"><b>Exploratory.</b> Mechanism estimates whether a damaging variant stops the protein working once it is in place (activity-specific, near 1) or stops it getting there at all (expression-driven, near 0). It is only given for confident pathogenic calls. On genes without lab data, single predictions are unreliable: on SLC19A2 and SLC19A3 two independent models gave almost unrelated answers. Use it to choose variants to test in the lab. <a href="help.html#activity">More</a></div>`;
    } else if (m === "act") {
      note.innerHTML = `<div class="notice" style="margin:12px 14px 0;background:var(--sunk)">Activity is a rough estimate of how much transport or enzyme activity the protein keeps with each missense change, from 0 (none) to 1 (normal). <a href="help.html#activity">More</a></div>`;
    } else note.innerHTML = "";
    renderResInfo();
  }

  function renderSymLegend() {
    const f = S.meta.features, has = t => f.some(x => x[0] === t);
    const parts = [];
    parts.push(`<span class="scale"><svg width="14" height="14" aria-hidden="true"><circle cx="7" cy="7" r="4.5" fill="var(--ps)"/></svg> ClinVar missense variant, coloured by evidence</span>`);
    if (has("Binding site")) parts.push(`<span class="scale"><svg width="14" height="14" aria-hidden="true"><circle cx="7" cy="7" r="5" class="fs bind"/></svg> Binding site</span>`);
    if (has("Active site")) parts.push(`<span class="scale"><svg width="14" height="14" aria-hidden="true"><rect x="3" y="3" width="8" height="8" transform="rotate(45 7 7)" class="fs act"/></svg> Active site</span>`);
    if (has("Site")) parts.push(`<span class="scale"><svg width="14" height="14" aria-hidden="true"><rect x="2.5" y="2.5" width="9" height="9" class="fs site"/></svg> Functional site</span>`);
    app.querySelector("#symlegend").innerHTML = parts.join(" ");
  }

  function mapNotes() {
    const m = S.meta, notes = [];
    const topo = S.wave.topo;
    if (!topo.hasTM) notes.push(topo.annotated ? "UniProt lists no membrane-crossing helix for this protein, so the chain is drawn on one side." : "UniProt has no membrane topology for this protein, so the chain is drawn as a plain line.");
    else if (topo.inferred) notes.push("UniProt does not say which side some loops are on. They are drawn alternating across each helix.");
    if (m.uniprot_sequence_matches === false) notes.push("UniProt has updated this sequence since the model was run, so its domains and sites are not drawn.");
    const nCls = S.rows.filter(o => o.cls).length;
    if (nCls) notes.push(`${fmt(nCls)} missense variant${nCls > 1 ? "s" : ""} could not be placed on this sequence, so ${nCls > 1 ? "they are" : "it is"} not on the map. ${nCls > 1 ? "They have" : "It has"} a class-level score only; look for the “Class-level score” flag in the table.`);
    notes.push("Drag across the map to zoom in. Double-click it, or press “Show whole protein”, to zoom out.");
    app.querySelector("#mapnote").textContent = notes.join(" ");
  }

  function regionList(list) {
    if (!list.length) return "";
    return `<ul class="hot">${list.map(h => `<li data-pos="${h.at}" data-a="${h.a}" data-b="${h.b}"><span>${esc(S.meta.sequence[h.a - 1])}${h.a}–${esc(S.meta.sequence[h.b - 1])}${h.b}<br><span class="small muted">${esc(h.region)}</span></span><span class="small" style="text-align:right">${esc(TVMap.binName(h.v))}</span></li>`).join("")}</ul>`;
  }

  function renderResInfo() {
    const box = app.querySelector("#resinfo");
    if (!box || !S.wave) return;
    const p = S.pos;
    if (!p) {
      const prof = S.wave.profiles.patho;
      const hot = TVMap.hotspots(S.meta.length, prof, S.wave.topo, 5, 18, 1.0, true).filter(h => h.v >= Math.log10(4.33));
      const cool = TVMap.hotspots(S.meta.length, prof, S.wave.topo, 3, 18, 1.0, false).filter(h => h.v <= -Math.log10(4.33));
      box.innerHTML = S.placed.length ? `
        <h3>Regions with the most damaging predictions</h3>
        <p class="empty small">Based on ${S.track === "routed" ? "the With AlphaMissense version" : "our model"}. Click a region to see its variants.</p>
        ${hot.length ? regionList(hot) : `<p class="small muted">No region averages moderate pathogenic evidence or stronger.</p>`}
        ${cool.length ? `<h3 style="margin-top:18px">Regions where changes are tolerated</h3>${regionList(cool)}` : ""}
        <p class="small muted" style="margin-top:14px">These regions reflect where ClinVar variants have been reported. A part of the protein with few reported variants can still be important.</p>`
        : `<h3>No missense variants to map</h3><p class="empty">None of this gene's ClinVar missense variants could be placed on its sequence. Every variant is still listed in the table below.</p>`;
      box.querySelectorAll("li[data-pos]").forEach(li => li.addEventListener("click", () =>
        selectPos(Number(li.dataset.pos), [Number(li.dataset.a), Number(li.dataset.b)])));
      return;
    }
    const seq = S.meta.sequence, aa = seq[p - 1] || "?";
    const R = S.range;
    const vs = R ? S.placed.filter(v => v.pos >= R[0] && v.pos <= R[1]).sort((a, b) => a.pos - b.pos || a.ev - b.ev)
      : (S.wave.byPos.get(p) || []);
    const feats = S.meta.features.filter(f => f[1] <= p && p <= f[2] && !["Topological domain", "Transmembrane"].includes(f[0]));
    const rows = [];
    rows.push(`<dt>Region</dt><dd>${esc(S.wave.topo.regionOf[p])}</dd>`);
    for (const f of feats) rows.push(`<dt>${esc(f[0])}</dt><dd>${esc(f[3])}</dd>`);
    for (const k of ["patho", "act", "mech"]) {
      const v = S.wave.profiles[k].val[p];
      if (!Number.isNaN(v)) rows.push(`<dt>${TVMap.METRICS[k].short} nearby${k === "mech" ? " <span class=\"tag-exp\">exploratory</span>" : ""}</dt><dd>${esc(TVMap.METRICS[k].describe(v))}</dd>`);
    }
    if (S.mol && S.mol.mapped && S.metric in TVMap.METRICS) {
      const v3 = S.mol.regional(p);
      if (v3 !== undefined) rows.push(`<dt>${TVMap.METRICS[S.metric].short} nearby in 3D</dt><dd>${esc(TVMap.METRICS[S.metric].describe(v3))}</dd>`);
    }
    if (S.mol && S.mol.plddt && S.mol.plddt[p] != null) rows.push(`<dt>AlphaFold confidence</dt><dd>${S.mol.plddt[p].toFixed(0)} of 100${S.mol.plddt[p] < 50 ? ", very low, likely disordered" : ""}</dd>`);
    box.innerHTML = `
      <h3>${esc(aa)}${p}</h3>
      ${S.part ? `<p class="parthint"><span class="sw-part"></span>${S.mol && S.mol.mapped ? "Highlighted on the map and the structure" : "Highlighted on the map"}: ${esc(S.part.name)} (${esc(seq[S.part.range[0] - 1])}${S.part.range[0]}–${esc(seq[S.part.range[1] - 1])}${S.part.range[1]})</p>` : ""}
      <dl>${rows.join("")}</dl>
      <h3>${R ? `Variants in ${esc(seq[R[0] - 1])}${R[0]}–${esc(seq[R[1] - 1])}${R[1]}` : "Variants at this residue"}</h3>
      ${vs.length ? `<ul>${vs.map(v => `<li><span class="mono">${esc(v.pv)}</span>${evChip(v.ev)}</li>`).join("")}</ul>` : `<p class="empty small">No ClinVar missense variant ${R ? "in this region" : `at ${esc(aa)}${p}`}.</p>`}
      <p style="margin-top:14px"><button type="button" class="btn ghost" id="clearpos">Clear selection</button></p>`;
    box.querySelector("#clearpos").addEventListener("click", () => selectPos(null));
  }

  // A click picks a residue; the part of the protein it sits in (its helix or loop) is highlighted
  // on the map and in 3D. A click on a listed region picks that region as the part.
  function partFor(p, range) {
    if (!p) return null;
    if (range) return { range, name: `${S.wave.topo.regionOf[p]}, region ${range[0]}–${range[1]}` };
    return S.wave.partAround(p);
  }
  function selectPos(p, range = null, origin = null) {
    S.pos = p; S.range = p ? range : null; S.page = 0;
    S.part = partFor(p, range);
    const r = S.part ? S.part.range : null;
    S.wave.select(p, r);
    // a click on the structure itself keeps the camera where the visitor put it
    if (S.mol) S.mol.select(p, r, origin !== "mol");
    // when the 3D view sits below the map (narrow screens), bring it into view after a map click
    if (origin === "map" && p && S.mol && S.mol.viewer) {
      const box = app.querySelector("#mol").getBoundingClientRect();
      if (box.top > innerHeight * 0.7 || box.bottom < 0) app.querySelector("#mol").scrollIntoView({ behavior: "smooth", block: "center" });
    }
    renderResInfo();
    renderFilters();
    renderTable();
    setURL();
  }

  /* ---------- table ---------- */
  function textMatch(o, q) {
    if (!q) return true;
    const Q = q.trim();
    if (/^vid:\d+$/i.test(Q)) return String(o.vid) === Q.slice(4);
    if (/^rs\d+$/i.test(Q)) return o.rs != null && `rs${o.rs}` === Q.toLowerCase();
    const pk = proteinKey(Q);
    if (pk && o.pv === pk) return true;
    const ql = Q.toLowerCase();
    return o.hgvs.toLowerCase().includes(ql) || (o.pv && o.pv.toLowerCase() === ql) || String(o.vid) === Q;
  }
  function filtered() {
    return S.rows.filter(o =>
      (!S.cons.size || S.cons.has(o.cons)) && (!S.ev.size || S.ev.has(o.ev)) &&
      (!S.excludeTraining || !o.train) && textMatch(o, S.text) &&
      (!S.pos || o.cons === 0 && (S.range ? o.pos >= S.range[0] && o.pos <= S.range[1] : o.pos === S.pos)));
  }

  function renderFilters() {
    const box = app.querySelector("#filters");
    const consC = new Map(), evC = new Array(8).fill(0);
    for (const o of S.rows) { consC.set(o.cons, (consC.get(o.cons) || 0) + 1); evC[o.ev]++; }
    const consKeys = [...consC.keys()].sort((a, b) => consC.get(b) - consC.get(a));
    const nTrain = S.rows.filter(o => o.train).length;
    box.innerHTML = `
      <div class="row"><span class="lbl">Variant type</span><div class="chips">${consKeys.map(k =>
        `<button type="button" class="chip" data-cons="${k}" aria-pressed="${S.cons.has(k)}">${esc(CONS_LABEL[S.idx.cons[k]])}<span class="ct">${fmt(consC.get(k))}</span></button>`).join("")}</div></div>
      <div class="row"><span class="lbl">Evidence level</span><div class="chips">${evC.map((c, i) => c ?
        `<button type="button" class="chip" data-ev="${i}" aria-pressed="${S.ev.has(i)}"><span class="dot" style="--c:var(--${TV.EV_KEYS[i]})"></span>${EV_LABEL[i]}<span class="ct">${fmt(c)}</span></button>` : "").join("")}</div></div>
      <div class="row">
        <span class="lbl">Find</span>
        <input type="search" id="tq" placeholder="Find a variant, e.g. c.515G>A or p.Gly172Asp" value="${esc(S.text)}" aria-label="Filter variants by HGVS, protein change or rsID">
        ${nTrain ? `<label class="check"><input type="checkbox" id="extrain" ${S.excludeTraining ? "checked" : ""}> Hide the ${fmt(nTrain)} variants seen in training</label>`
          : `<span class="small muted">No variant in this gene was used in training.</span>`}
        ${S.pos ? `<span class="active-filter">${S.range ? `Residues ${esc(S.meta.sequence[S.range[0] - 1])}${S.range[0]}–${esc(S.meta.sequence[S.range[1] - 1])}${S.range[1]}` : `Residue ${esc((S.meta.sequence || "")[S.pos - 1] || "")}${S.pos}`}<button type="button" id="xpos" aria-label="Clear residue filter">×</button></span>` : ""}
        <button type="button" class="btn ghost" id="dlgene" style="margin-left:auto">Download this gene (CSV, ${TRACK[S.track].license})</button>
      </div>`;
    box.querySelectorAll("[data-cons]").forEach(b => b.addEventListener("click", () => {
      const k = Number(b.dataset.cons); S.cons.has(k) ? S.cons.delete(k) : S.cons.add(k); S.page = 0; renderFilters(); renderTable();
    }));
    box.querySelectorAll("[data-ev]").forEach(b => b.addEventListener("click", () => {
      const k = Number(b.dataset.ev); S.ev.has(k) ? S.ev.delete(k) : S.ev.add(k); S.page = 0; renderFilters(); renderTable();
    }));
    const tq = box.querySelector("#tq");
    let tmr;
    tq.addEventListener("input", () => { clearTimeout(tmr); tmr = setTimeout(() => { S.text = tq.value; S.page = 0; renderTable(); setURL(); }, 160); });
    const ex = box.querySelector("#extrain");
    if (ex) ex.addEventListener("change", e => { S.excludeTraining = e.target.checked; S.page = 0; renderTable(); });
    const xp = box.querySelector("#xpos");
    if (xp) xp.addEventListener("click", () => selectPos(null));
    box.querySelector("#dlgene").addEventListener("click", downloadGene);
  }

  const SORTS = {
    hgvs: o => [o.pos ?? 1e9, o.hs], cons: o => [o.cons], ev: o => [o.ev, -(o.lr ?? 0)], lr: o => [o.lr ?? -1],
    score: o => [o.score ?? -1], act: o => [o.act ?? 9], mech: o => [o.mclip ?? 9], clin: o => [o.clinsig],
  };
  function cmp(a, b) { for (let i = 0; i < a.length; i++) { if (a[i] < b[i]) return -1; if (a[i] > b[i]) return 1; } return 0; }

  function head() {
    const routed = S.track === "routed";
    const th = (k, label, sub = "", cls = "") => {
      const dir = S.sort.key === k ? S.sort.dir : null;
      return `<th class="${cls}" ${dir ? `aria-sort="${dir === "asc" ? "ascending" : "descending"}"` : ""}><button type="button" data-sort="${k}" ${dir ? `data-dir="${dir}"` : ""}>${label}</button>${sub ? `<span class="sub">${sub}</span>` : ""}</th>`;
    };
    return `<thead><tr>
      ${th("hgvs", "Variant", "HGVS, rsID and flags")}
      ${th("cons", "Type")}
      ${th("ev", "Evidence level", "use this")}
      ${th("lr", "LR", "likelihood ratio", "num")}
      ${th("score", "Raw score", "for reference", "num")}
      ${th("act", "Activity", "0 none, 1 normal", "num")}
      ${th("mech", "Mechanism", "exploratory", "num")}
      ${th("clin", "ClinVar", "classification")}
    </tr></thead>`;
  }

  function flags(o) {
    const b = [];
    if (o.train) b.push(`<span class="badge train" title="The model saw this variant during fitting. Exclude it from any accuracy estimate.">Seen in training</span>`);
    if (o.cls) b.push(`<span class="badge noam" title="ClinVar's protein change could not be placed on the UniProt sequence, so this missense variant was scored by the consequence-class model, which knows it is missense but not where or what. Its evidence is weak.">Class-level score</span>`);
    if (S.track === "routed") {
      b.push(o.hasAM ? `<span class="badge am" title="AlphaMissense covers this variant; the score combines it with our model.">AlphaMissense</span>`
        : `<span class="badge noam" title="AlphaMissense does not cover this variant; the native score is served.">No AlphaMissense</span>`);
    }
    if (o.lis && (o.act != null || o.mech != null)) b.push(`<span class="badge lis" title="This protein has a deep mutational scan that trained the activity and mechanism layers, so those values are in-sample.">Layers in-sample</span>`);
    return b.join("");
  }

  function detail(o) {
    const routed = S.track === "routed";
    const pred = [
      `<dt>Version</dt><dd>${routed ? (o.plus ? "With AlphaMissense: our model combined with AlphaMissense" : "With AlphaMissense: our model alone (AlphaMissense does not cover this variant)") : "Our model only"}</dd>`,
      `<dt>Evidence</dt><dd>${evChip(o.ev, true)}</dd>`,
      `<dt>Likelihood ratio</dt><dd>${fmtLR(o.lr)}</dd>`,
      `<dt>Raw score</dt><dd>${o.score ?? ""}</dd>`,
    ];
    if (routed) {
      pred.push(`<dt>Our model alone</dt><dd>${o.v2n ?? ""}</dd>`);
      pred.push(`<dt>Our model with AlphaMissense</dt><dd>${o.v2p ?? ""}${o.hasAM ? "" : ` <span class="badge noam">not meaningful here</span>`}</dd>`);
      pred.push(`<dt>AlphaMissense</dt><dd>${o.am ?? "not covered"}</dd>`);
    }
    const lay = [];
    if (o.act != null) lay.push(`<dt>Activity</dt><dd>${o.act.toFixed(3)}</dd>`);
    if (o.mech != null) lay.push(`<dt>Mechanism <span class="tag-exp">exploratory</span></dt><dd>${o.mclip.toFixed(2)}, ${mechWord(o.mclip)}${o.mech !== o.mclip ? ` (model output ${o.mech.toFixed(3)}, shown clipped to the 0 to 1 range)` : ""}</dd>`);
    if (o.lis && lay.length) lay.push(`<dt>In-sample</dt><dd>This protein's deep mutational scan trained these layers.</dd>`);
    const notes = [];
    if (o.cls) notes.push("<b>Class-level score.</b> This missense change could not be placed on the UniProt sequence the model uses. It changes more than one DNA base, or ClinVar numbers it on a different transcript. It was scored by a simpler model that only knows the variant is missense, so its evidence is weak and it has no activity or mechanism value.");
    if (o.ev === 0) notes.push('<b>Count this as strong.</b> On genes held out from calibration, variants at this level were 99.2% pathogenic, which is the strength of the strong level. <a href="help.html#levels">Details</a>');
    if (o.ev === 4) notes.push('<b>Indeterminate is not neutral here.</b> On held-out genes, 76% of variants at this level were pathogenic. <a href="help.html#levels">Details</a>');
    if (o.train) notes.push("<b>Seen in training.</b> The model was trained on this variant's ClinVar classification, so this prediction is not an independent test. Leave it out of any accuracy estimate.");
    if (routed && !o.hasAM) notes.push("AlphaMissense does not cover this variant, so the score shown is our model's alone. Ignore the “combined with AlphaMissense” value for this row.");
    if (o.mech != null) notes.push("Mechanism is exploratory. Single predictions on genes without lab data are unreliable.");
    return `<div class="det-grid">
      <div><h4>Identifiers</h4><dl>
        <dt>HGVS</dt><dd class="mono">${esc(o.hgvs)}</dd>
        ${o.pv ? `<dt>Protein</dt><dd class="mono">${esc(o.pv)}</dd>` : ""}
        <dt>ClinVar</dt><dd><a href="https://www.ncbi.nlm.nih.gov/clinvar/variation/${o.vid}/" rel="noopener">VariationID ${o.vid}</a></dd>
        <dt>dbSNP</dt><dd>${o.rs ? `<a href="https://www.ncbi.nlm.nih.gov/snp/rs${o.rs}" rel="noopener">rs${o.rs}</a>` : "none"}</dd>
        <dt>Assertion</dt><dd>${esc(o.clinsig) || "none"}</dd>
        <dt>Review</dt><dd>${esc(o.review) || "none"}</dd>
      </dl></div>
      <div><h4>Prediction</h4><dl>${pred.join("")}</dl></div>
      ${lay.length ? `<div><h4>Activity and mechanism</h4><dl>${lay.join("")}</dl></div>` : ""}
    </div>${notes.length ? `<p class="det-note">${notes.join(" ")}</p>` : ""}`;
  }

  function renderTable() {
    const rows = filtered();
    const f = SORTS[S.sort.key], dir = S.sort.dir === "asc" ? 1 : -1;
    rows.sort((a, b) => dir * cmp(f(a), f(b)) || cmp(SORTS.hgvs(a), SORTS.hgvs(b)));
    const pages = Math.max(1, Math.ceil(rows.length / PAGE));
    S.page = Math.min(S.page, pages - 1);
    const slice = rows.slice(S.page * PAGE, S.page * PAGE + PAGE);
    const tx = mostCommonPrefix();
    const hl = S.text ? new Set(rows.filter(o => textMatch(o, S.text)).map(o => o.vid)) : new Set();
    const body = slice.map(o => {
      const open = S.open.has(o.vid);
      const full = o.hp === tx ? o.hs : o.hgvs;
      const m = full.match(/^(.*?)\s*\((p\.[^)]*)\)\s*$/);
      const hg = m ? m[1] : full, pp = m ? m[2] : "";
      return `<tr class="v${open ? " open" : ""}${S.text && hl.has(o.vid) && rows.length <= 3 ? " hl" : ""}" data-vid="${o.vid}" tabindex="0" aria-expanded="${open}">
        <td class="vcell"><div class="hg">${esc(hg)}</div>${pp ? `<div class="pp">${esc(pp)}</div>` : ""}<div class="rs">${o.rs ? `rs${o.rs}` : ""}</div><div>${flags(o)}</div></td>
        <td>${esc(CONS_LABEL[S.idx.cons[o.cons]])}</td>
        <td>${evChip(o.ev, true)}</td>
        <td class="num">${fmtLR(o.lr)}</td>
        <td class="num sc">${o.score == null ? "" : o.score.toFixed(3)}</td>
        <td class="num">${o.act == null ? "" : o.act.toFixed(2)}</td>
        <td class="num exp">${o.mclip == null ? "" : o.mclip.toFixed(2)}</td>
        <td class="small">${esc(o.clinsig)}</td>
      </tr>${open ? `<tr class="det"><td colspan="8">${detail(o)}</td></tr>` : ""}`;
    }).join("");
    const t = app.querySelector("#vt");
    t.innerHTML = head() + `<tbody>${body || `<tr><td colspan="8" style="padding:28px 12px" class="muted">No variants match these filters. Clear a filter above to see more.</td></tr>`}</tbody>` +
      (tx ? `<caption class="visually-hidden">Variants in ${esc(GENE)}; HGVS shown relative to ${esc(tx)}</caption>` : "");
    t.querySelectorAll("button[data-sort]").forEach(b => b.addEventListener("click", e => {
      e.stopPropagation();
      const k = b.dataset.sort;
      S.sort = S.sort.key === k ? { key: k, dir: S.sort.dir === "asc" ? "desc" : "asc" } : { key: k, dir: ["lr", "score"].includes(k) ? "desc" : "asc" };
      renderTable();
    }));
    t.querySelectorAll("tr.v").forEach(tr => {
      const toggle = () => { const v = Number(tr.dataset.vid); S.open.has(v) ? S.open.delete(v) : S.open.add(v); renderTable(); };
      tr.addEventListener("click", e => { if (!e.target.closest("a")) toggle(); });
      tr.addEventListener("keydown", e => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); toggle(); } });
    });
    const pg = app.querySelector("#pager");
    pg.innerHTML = `<span>${fmt(rows.length)} of ${fmt(S.rows.length)} variants${tx ? `. HGVS is on ${esc(tx.replace(/:$/, ""))} unless shown in full.` : ""}</span>
      <span class="pg"><button type="button" data-p="-1" ${S.page === 0 ? "disabled" : ""}>Previous</button>
      <span>Page ${S.page + 1} of ${pages}</span>
      <button type="button" data-p="1" ${S.page >= pages - 1 ? "disabled" : ""}>Next</button></span>`;
    pg.querySelectorAll("button[data-p]").forEach(b => b.addEventListener("click", () => { S.page += Number(b.dataset.p); renderTable(); app.querySelector("#tab-h").scrollIntoView({ block: "start" }); }));
    if (S.text && rows.length === 1 && !S.open.size) { S.open.add(rows[0].vid); renderTable(); }
  }
  function mostCommonPrefix() {
    const c = new Map();
    for (const o of S.rows) c.set(o.hp, (c.get(o.hp) || 0) + 1);
    let best = "", n = 0;
    for (const [k, v] of c) if (v > n) { best = k; n = v; }
    return best;
  }

  function downloadGene() {
    const routed = S.track === "routed";
    const cols = ["gene", "variation_id", "rsid", "hgvs", "consequence", "protein_variant", "clinsig", "review_status",
      "in_training", "class_level_score", "evidence", "lr", "score", ...(routed ? ["served_by", "has_alphamissense", "v2_native", "v2_plus", "alphamissense"] : []),
      "activity_pred", "mechanism_pred", "layers_in_sample"];
    const q = v => v == null ? "" : /[",\n]/.test(String(v)) ? `"${String(v).replace(/"/g, '""')}"` : String(v);
    const lines = [cols.join(",")];
    for (const o of S.rows) {
      lines.push([GENE, o.vid, o.rs ? `rs${o.rs}` : "", o.hgvs, S.idx.cons[o.cons], o.pv, o.clinsig, o.review,
        o.train ? "TRUE" : "FALSE", o.cls ? "TRUE" : "FALSE", S.idx.ev[o.ev], o.lr, o.score,
        ...(routed ? [o.plus ? "plus" : "native", o.hasAM ? "TRUE" : "FALSE", o.v2n, o.v2p, o.am] : []),
        o.act, o.mech, o.act == null && o.mech == null ? "" : (o.lis ? "TRUE" : "FALSE")].map(q).join(","));
    }
    const blob = new Blob([lines.join("\n") + "\n"], { type: "text/csv" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = `${GENE}_${S.track}_${routed ? "CC-BY-NC-SA-4.0" : "CC-BY-4.0"}.csv`;
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(a.href), 2000);
  }

  /* ---------- boot ---------- */
  async function main() {
    if (!GENE) { location.replace("index.html"); return; }
    try {
      S.idx = await TV.index();
    } catch (e) {
      app.innerHTML = `<div class="wrap loading">The gene index could not be loaded (${esc(e.message)}). Reload to try again.</div>`;
      return;
    }
    S.row = S.idx.genes.find(r => r[0] === GENE) || S.idx.genes.find(r => r[0].toUpperCase() === GENE.toUpperCase());
    if (!S.row) {
      document.title = "Gene not found | Transporter Variants";
      app.innerHTML = `<div class="wrap loading"><h1 style="font-size:2rem">${esc(GENE)} is not in this release</h1><p class="muted" style="margin-top:12px">The release covers 783 transporter, ion-channel, ABC-transporter and ATPase genes. Search for another gene above, or <a href="index.html#browse">browse the full list</a>.</p></div>`;
      return;
    }
    if (S.row[0] !== GENE) { location.replace(TV.geneUrl(S.row[0])); return; }
    try {
      S.meta = await getJSON(`data/protein/${GENE}.json`);
      await loadTrack();
    } catch (e) {
      app.innerHTML = `<div class="wrap loading">The data for ${esc(GENE)} could not be loaded (${esc(e.message)}). Reload to try again.</div>`;
      return;
    }
    renderShell();
    renderTrackBar();

    S.wave = new TVMap.Wave(app.querySelector("#wave"), {
      meta: S.meta, variants: S.placed,
      onSelect: p => selectPos(S.pos === p ? null : p, null, "map"),
      onZoom: z => { app.querySelector("#unzoom").hidden = !z; },
    });
    renderSymLegend();
    mapNotes();
    setMetric("patho");
    if (S.pos) { S.part = partFor(S.pos, null); S.wave.select(S.pos, S.part && S.part.range); }
    renderResInfo();
    renderFilters();
    renderTable();

    if (S.meta.uniprot && S.meta.sequence) {
      S.mol = new TVMap.Mol(app.querySelector("#mol"), { meta: S.meta, wave: S.wave, onSelect: p => selectPos(p, null, "mol") });
      const go = () => S.mol.load().then(() => { S.mol.setMetric(S.metric); if (S.pos) S.mol.select(S.pos, S.part && S.part.range); renderResInfo(); });
      if ("IntersectionObserver" in window) {
        const io = new IntersectionObserver(es => { if (es.some(e => e.isIntersecting)) { io.disconnect(); go(); } }, { rootMargin: "200px" });
        io.observe(app.querySelector("#mol"));
      } else go();
    }
    if (S.text) setTimeout(() => app.querySelector("#tab-h").scrollIntoView({ block: "start" }), 50);
  }
  main();
})();
