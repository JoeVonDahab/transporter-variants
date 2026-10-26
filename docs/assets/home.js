/* Home page: key numbers, hero map, gene list. */
"use strict";

(async () => {
  const { esc, fmt, getJSON, EV_KEYS } = TV;
  let idx;
  try { idx = await TV.index(); } catch (e) {
    document.getElementById("genelist").innerHTML = `<p class="muted" style="padding:16px">The gene list could not be loaded (${esc(e.message)}). Reload the page to try again.</p>`;
    return;
  }
  const S = idx.summary;
  TV.fillRelease(S);
  document.getElementById("f-variants").textContent = fmt(S.variants);
  document.getElementById("f-genes").textContent = fmt(S.genes);
  document.getElementById("f-activity").textContent = fmt(S.activity);

  /* ---------- gene list ---------- */
  const FAM = { "ion channel": "Ion channels", SLC: "Solute carriers", ABC: "ABC transporters", "P-type ATPase": "P-type ATPases",
    CFTR: "CFTR", "connexin/pannexin": "Connexins and pannexins", "V-type ATPase": "V-type ATPases", aquaporin: "Aquaporins", unassigned: "Other" };
  let fam = null, filt = "", showAll = false;
  const COLLAPSED = 40;
  const chips = document.getElementById("famchips");
  function drawChips() {
    chips.innerHTML = `<button type="button" class="chip" data-f="" aria-pressed="${fam === null}">All<span class="ct">${fmt(idx.genes.length)}</span></button>` +
      S.families.map(([f, g]) => `<button type="button" class="chip" data-f="${esc(f)}" aria-pressed="${fam === f}">${esc(FAM[f] || f)}<span class="ct">${fmt(g)}</span></button>`).join("");
  }
  function drawGenes() {
    const F = filt.trim().toUpperCase();
    const list = idx.genes.filter(g => (!fam || g[1] === fam) && (!F || g[0].toUpperCase().includes(F) || g[3].toUpperCase().includes(F)));
    list.sort((a, b) => a[0].localeCompare(b[0], "en", { numeric: true }));
    const cut = !showAll && !F && !fam && list.length > COLLAPSED;
    const shown = cut ? list.slice(0, COLLAPSED) : list;
    const more = document.getElementById("gmore");
    more.hidden = !cut;
    more.textContent = `Show all ${fmt(list.length)} genes`;
    document.getElementById("genelist").innerHTML = shown.length ? shown.map(g => {
      const bins = g[6], n = g[4];
      return `<a href="${TV.geneUrl(g[0])}"><span class="n">${fmt(n)} variants</span><span class="g">${esc(g[0])}</span><span class="p" title="${esc(g[3])}">${esc(g[3] || " ")}</span>
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
    document.getElementById("hero-cap").innerHTML = `${esc(G)}, the thiamine transporter, drawn across the cell membrane with its ${w.topo.nTM} helices. Red stretches are where variants tend to be damaging, blue where they tend to be tolerated. Each dot is one of its ${fmt(vs.length)} ClinVar missense variants. Point at the chain to read a position, or click to open it on the <a href="gene.html?g=${G}">${esc(G)} page</a>.`;
    host.insertAdjacentHTML("beforeend", `<div class="map-foot">${TVMap.legendRamp("patho")}</div>`);
  } catch (e) {
    document.getElementById("hero-wave").hidden = true;
  }
  // the map and gene list arrive after the browser has jumped to #genes; jump again now they are in place
  if (location.hash) {
    const t = document.getElementById(location.hash.slice(1));
    if (t) t.scrollIntoView();
  }

  /* ---------- example searches ---------- */
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
