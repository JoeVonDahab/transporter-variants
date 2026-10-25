/* Landing page: hero map, coverage, calibration, gene browser, downloads. */
"use strict";

(async () => {
  const { esc, fmt, getJSON, EV_LABEL, EV_KEYS, CONS_LABEL } = TV;
  let idx;
  try { idx = await TV.index(); } catch (e) {
    document.getElementById("genes").innerHTML = `<p class="muted" style="padding:16px">The gene index could not be loaded (${esc(e.message)}). Reload to try again.</p>`;
    return;
  }
  const S = idx.summary;

  /* ---------- coverage ---------- */
  const cons = S.by_consequence.slice().sort((a, b) => b[1] - a[1]);
  // neutral greys and lipid tones: the crimson-to-teal scale is reserved for evidence
  const palette = ["var(--ink-2)", "#8c9aa1", "#a9b4b9", "#c4cdd0", "var(--lipid-head)", "#ddb95c", "var(--lipid)", "#efd894",
    "#d6c7a1", "#bfae86", "#a7976f", "#8f8159", "#776b45"];
  const tot = S.variants;
  document.getElementById("cov-nonmis").textContent = fmt(tot - S.missense);
  document.getElementById("cov-am").textContent = fmt(S.has_alphamissense);
  document.getElementById("cov").innerHTML = `
    <div class="cov-bar" role="img" aria-label="Variants by consequence class; AlphaMissense covers ${fmt(S.has_alphamissense)} of ${fmt(tot)}, all missense">
      ${cons.map(([c, n, am], i) => `<div style="flex:${n};background:${palette[i]}" title="${esc(CONS_LABEL[c])}: ${fmt(n)}">${am ? `<span class="am" style="width:${(am / n * 100).toFixed(2)}%"></span>` : ""}</div>`).join("")}
    </div>
    <div class="cov-legend">${cons.map(([c, n], i) => `<div><span class="sw" style="background:${palette[i]}"></span>${esc(CONS_LABEL[c])} <b>${fmt(n)}</b></div>`).join("")}</div>`;

  /* ---------- calibration ---------- */
  const cal = S.calibration;
  const order = ["pathogenic very strong", "pathogenic strong", "pathogenic moderate", "pathogenic supporting",
    "indeterminate", "benign supporting", "benign moderate", "benign strong"];
  const rows = order.map(b => cal.rows.find(r => r.bin === b)).filter(Boolean);
  const unscored = cal.rows.find(r => r.bin === "unscored");
  document.getElementById("caltable").innerHTML = `
    <thead><tr><th>Evidence bin</th><th class="n">Variants</th><th>Observed pathogenic or likely pathogenic</th><th></th></tr></thead>
    <tbody>${rows.map(r => {
      const i = order.indexOf(r.bin);
      return `<tr><td>${TV.evChip(i)}</td><td class="n">${fmt(r.n)}</td>
        <td class="barcell"><div class="track"><div class="fill" style="width:${r.observed}%;background:var(--${EV_KEYS[i]})"></div><div class="base" style="left:${cal.base_rate}%"></div></div></td>
        <td class="pct">${r.observed.toFixed(1)}%</td></tr>`;
    }).join("")}</tbody>
    <caption>${fmt(cal.n)} labelled out-of-sample variants across ${fmt(cal.genes)} genes. The dashed line is the base rate, ${cal.base_rate}% pathogenic. The top two bins invert by 0.3 points, which is saturation above 99% rather than miscalibration.${unscored ? ` ${fmt(unscored.n)} variants were unscored in this check.` : ""}</caption>`;

  if (S.uniprot_release) document.getElementById("up-rel").textContent = `, release ${S.uniprot_release}`;
  const MONTHS = ["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"];
  const [yy, mo, dd] = S.released.split("-").map(Number);
  document.querySelectorAll(".relver").forEach(e => { e.textContent = S.version; });
  document.querySelectorAll(".reldate").forEach(e => { e.textContent = `${dd} ${MONTHS[mo - 1]} ${yy}`; });
  if (S.class_level_missense != null) document.getElementById("cls-n").textContent = fmt(S.class_level_missense);

  /* ---------- downloads ---------- */
  const mb = b => (b / 1e6).toFixed(1) + " MB";
  const d = S.downloads;
  document.getElementById("dl-native").textContent = `Download native track, ${mb(d.native.bytes)}, CC BY 4.0`;
  document.getElementById("dl-routed").textContent = `Download routed track, ${mb(d.routed.bytes)}, CC BY-NC-SA 4.0`;
  document.getElementById("sha-native").textContent = `SHA-256 ${d.native.sha256}`;
  document.getElementById("sha-routed").textContent = `SHA-256 ${d.routed.sha256}`;

  /* ---------- browse ---------- */
  const FAM = { "ion channel": "Ion channels", SLC: "Solute carriers", ABC: "ABC transporters", "P-type ATPase": "P-type ATPases",
    CFTR: "CFTR", "connexin/pannexin": "Connexins and pannexins", "V-type ATPase": "V-type ATPases", aquaporin: "Aquaporins", unassigned: "Unassigned" };
  const fams = S.families.map(f => f[0]);
  let fam = null, filt = "", showAll = false;
  const COLLAPSED = 40;
  const chips = document.getElementById("famchips");
  function drawChips() {
    chips.innerHTML = `<button type="button" class="chip" data-f="" aria-pressed="${fam === null}">All<span class="ct">${fmt(idx.genes.length)}</span></button>` +
      S.families.map(([f, g]) => `<button type="button" class="chip" data-f="${esc(f)}" aria-pressed="${fam === f}">${esc(FAM[f] || f)}<span class="ct">${fmt(g)}</span></button>`).join("");
  }
  const famRank = new Map(fams.map((f, i) => [f, i]));
  function drawGenes() {
    const F = filt.trim().toUpperCase();
    const list = idx.genes.filter(g => (!fam || g[1] === fam) && (!F || g[0].toUpperCase().includes(F) || g[3].toUpperCase().includes(F)));
    list.sort((a, b) => a[0].localeCompare(b[0], "en", { numeric: true }));
    const cut = !showAll && !F && !fam && list.length > COLLAPSED;
    const shown = cut ? list.slice(0, COLLAPSED) : list;
    const more = document.getElementById("gmore");
    more.hidden = !cut;
    more.textContent = `Show all ${fmt(list.length)} genes`;
    document.getElementById("genes").innerHTML = shown.length ? shown.map(g => {
      const bins = g[6], n = g[4];
      return `<a href="${TV.geneUrl(g[0])}"><span class="n">${fmt(n)}</span><span class="g">${esc(g[0])}</span><span class="p" title="${esc(g[3])}">${esc(g[3] || " ")}</span>
        <span class="minibar">${bins.map((c, i) => c ? `<i style="width:${(c / n * 100).toFixed(2)}%;background:var(--${EV_KEYS[i]})"></i>` : "").join("")}</span></a>`;
    }).join("") : `<p class="muted" style="background:var(--surface);padding:16px;margin:0">No gene matches “${esc(filt)}”.</p>`;
  }
  chips.addEventListener("click", e => {
    const b = e.target.closest("button[data-f]");
    if (!b) return;
    fam = b.dataset.f || null; drawChips(); drawGenes();
  });
  document.getElementById("gfilter").addEventListener("input", e => { filt = e.target.value; drawGenes(); });
  document.getElementById("gmore").addEventListener("click", () => { showAll = true; drawGenes(); });
  drawChips(); drawGenes();

  /* ---------- hero map ---------- */
  const G = "SLC19A2";
  try {
    const [meta, shard] = await Promise.all([getJSON(`data/protein/${G}.json`), getJSON(`data/native/${G}.json`)]);
    const vs = [];
    for (const r of shard.rows) {
      const pos = r[15], pv = r[5];
      if (r[4] !== 0 || !pos || meta.sequence[pos - 1] !== pv[0]) continue;
      vs.push({ pos, pv, ev: r[9], lr: r[11], act: r[12], mech: r[13] });
    }
    const host = document.getElementById("hero-wave");
    const w = new TVMap.Wave(host, { meta, variants: vs, compact: true, onSelect: p => { location.href = `gene.html?g=${G}&pos=${p}`; } });
    w.render();
    const nTM = w.topo.nTM;
    document.getElementById("hero-cap").innerHTML = `${esc(G)}, thiamine transporter 1, threaded through the membrane by its ${nTM} helices. Each dot is one of its ${fmt(vs.length)} ClinVar missense variants, coloured by evidence; the chain is coloured by the average prediction nearby, crimson where variants tend to be damaging and blue where they are tolerated. ${meta.features.some(f => f[0] === "Site") ? "Squares mark residues UniProt lists as functional sites. " : ""}Hover to read a residue, or click to open it on the <a href="gene.html?g=${G}">${esc(G)} page</a>.`;
    host.insertAdjacentHTML("beforeend", `<div class="map-foot">${TVMap.legendRamp("patho")}</div>`);
  } catch (e) {
    document.getElementById("hero-wave").hidden = true;
  }

  /* ---------- example variants for the search hint ---------- */
  const ex = document.getElementById("ex-variants");
  ex.innerHTML = `<a href="#" data-q="rs113993960">rs113993960</a> or <a href="#" data-q="SLC19A2 p.Gly172Asp">SLC19A2 p.Gly172Asp</a>`;
  ex.addEventListener("click", e => {
    const a = e.target.closest("a[data-q]");
    if (!a) return;
    e.preventDefault();
    const input = document.getElementById("q-hero");
    input.value = a.dataset.q;
    input.focus();
    input.dispatchEvent(new KeyboardEvent("keydown", { key: "Enter", bubbles: true }));
  });
})();
