/* Shared helpers: data loading, the evidence scale, formatting and search. */
"use strict";

const TV = (() => {
  const EV_KEYS = ["pvs", "ps", "pm", "pp", "ind", "bp", "bm", "bs"];
  const EV_LABEL = ["Pathogenic very strong", "Pathogenic strong", "Pathogenic moderate",
    "Pathogenic supporting", "Indeterminate", "Benign supporting", "Benign moderate", "Benign strong"];
  const CONS_LABEL = {
    missense: "Missense", synonymous: "Synonymous", intronic: "Intronic", splice_region: "Splice region",
    splice_canonical: "Canonical splice", frameshift: "Frameshift", nonsense: "Nonsense",
    inframe_indel: "In-frame indel", utr3: "3′ UTR", utr5: "5′ UTR", start_loss: "Start loss",
    stop_loss: "Stop loss", other: "Other",
  };
  const AA3 = { Ala: "A", Arg: "R", Asn: "N", Asp: "D", Cys: "C", Gln: "Q", Glu: "E", Gly: "G", His: "H",
    Ile: "I", Leu: "L", Lys: "K", Met: "M", Phe: "F", Pro: "P", Ser: "S", Thr: "T", Trp: "W", Tyr: "Y", Val: "V" };

  const cache = new Map();
  function getJSON(url) {
    if (!cache.has(url)) {
      cache.set(url, fetch(url).then(r => {
        if (!r.ok) throw new Error(`${url}: HTTP ${r.status}`);
        return r.json();
      }).catch(e => { cache.delete(url); throw e; }));
    }
    return cache.get(url);
  }
  const index = () => getJSON("data/index.json");

  function cssVar(name) { return getComputedStyle(document.documentElement).getPropertyValue(name).trim(); }
  function evColor(i) { return cssVar(`--${EV_KEYS[i]}`); }
  const evVar = i => `var(--${EV_KEYS[i]})`;

  const fmt = n => n == null ? "" : n.toLocaleString("en-US");
  function fmtLR(x) {
    if (x == null) return "";
    if (x >= 100) return Math.round(x).toLocaleString("en-US");
    if (x >= 10) return x.toFixed(1);
    if (x >= 1) return x.toFixed(2);
    return x.toFixed(3);
  }
  function esc(s) {
    return String(s ?? "").replace(/[&<>"']/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
  }
  function el(tag, attrs = {}, html = "") {
    const e = document.createElement(tag);
    for (const [k, v] of Object.entries(attrs)) {
      if (k === "class") e.className = v; else if (v !== false && v != null) e.setAttribute(k, v === true ? "" : v);
    }
    if (html) e.innerHTML = html;
    return e;
  }
  function evChip(i, strong) {
    return `<span class="ev${strong ? " strong-label" : ""}" style="--c:${evVar(i)}">${EV_LABEL[i]}</span>`;
  }

  // "p.Gly111Arg" / "Gly111Arg" / "G111R" -> "G111R"
  function proteinKey(q) {
    let s = q.trim().replace(/^p\./i, "").replace(/[()]/g, "");
    const m3 = s.match(/^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$/);
    if (m3 && AA3[m3[1]] && AA3[m3[3]]) return AA3[m3[1]] + m3[2] + AA3[m3[3]];
    const m1 = s.match(/^([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWY])$/i);
    if (m1) return (m1[1] + m1[2] + m1[3]).toUpperCase();
    return null;
  }

  // release number and date in footers and citations
  function fillRelease(S) {
    const M = ["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"];
    const [y, m, d] = String(S.released).split("-").map(Number);
    document.querySelectorAll(".relver").forEach(e => { e.textContent = S.version; });
    if (y && m && d) document.querySelectorAll(".reldate").forEach(e => { e.textContent = `${d} ${M[m - 1]} ${y}`; });
  }

  const geneUrl = (g, q, extra = "") => `gene.html?g=${encodeURIComponent(g)}${q ? `&q=${encodeURIComponent(q)}` : ""}${extra}`;

  /* ---------- search ----------
     Gene symbol or protein name: matched against the index.
     rs123 / ClinVar VariationID: looked up in 100 small shards.
     HGVS: the gene comes from "(GENE)" or from the transcript accession. */
  async function resolve(raw) {
    const q = raw.trim();
    if (!q) return { kind: "none" };
    const idx = await index();
    const genes = idx.genes;
    const bySym = new Map(genes.map(g => [g[0].toUpperCase(), g[0]]));

    let m = q.match(/^rs(\d+)$/i);
    if (m) {
      const shard = await getJSON(`data/lookup/rs/${String(Number(m[1]) % 100).padStart(2, "0")}.json`);
      const hit = shard[String(Number(m[1]))];
      return hit ? { kind: "variant", genes: hit, q: `rs${Number(m[1])}` } : { kind: "miss", what: `rs${m[1]}` };
    }
    m = q.match(/^(?:clinvar:|vcv0*)?(\d{2,9})$/i);
    if (m) {
      const id = String(Number(m[1]));
      const shard = await getJSON(`data/lookup/vid/${String(Number(id) % 100).padStart(2, "0")}.json`);
      return shard[id] ? { kind: "variant", genes: [shard[id]], q: `vid:${id}` } : { kind: "miss", what: `ClinVar variation ${id}` };
    }
    m = q.match(/^([A-Z]{2}_\d+)(?:\.\d+)?(?:\(([^)]+)\))?:(.+)$/i);
    if (m) {
      let g = m[2] && bySym.get(m[2].toUpperCase());
      if (!g) {
        const tx = m[1].toUpperCase();
        const hit = genes.find(r => r[9].includes(tx));
        if (hit) g = hit[0];
      }
      return g ? { kind: "variant", genes: [g], q: m[3].trim() } : { kind: "miss", what: m[1] };
    }
    // "SLC19A2 G111R", "SLC19A2:c.68G>C", "SLC19A2 p.Gly111Arg"
    m = q.match(/^([A-Za-z0-9-]+)[\s:]+(.+)$/);
    if (m && bySym.has(m[1].toUpperCase())) return { kind: "variant", genes: [bySym.get(m[1].toUpperCase())], q: m[2].trim() };
    if (bySym.has(q.toUpperCase())) return { kind: "gene", gene: bySym.get(q.toUpperCase()) };
    return { kind: "suggest" };
  }

  function suggestions(idx, q) {
    const Q = q.trim().toUpperCase();
    if (!Q) return [];
    const out = [];
    for (const g of idx.genes) {
      const s = g[0].toUpperCase();
      let score = -1;
      if (s === Q) score = 0; else if (s.startsWith(Q)) score = 1;
      else if (Q.length >= 3 && g[3].toUpperCase().includes(Q)) score = 3;
      else if (Q.length >= 3 && s.includes(Q)) score = 2;
      if (score >= 0) out.push([score, s.length, g]);
    }
    out.sort((a, b) => a[0] - b[0] || a[1] - b[1] || a[2][0].localeCompare(b[2][0]));
    return out.slice(0, 8).map(x => x[2]);
  }

  function attachSearch(root) {
    const input = root.querySelector("input");
    const list = el("ul", { class: "suggest", role: "listbox", id: `${input.id}-list` });
    list.hidden = true;
    root.appendChild(list);
    input.setAttribute("role", "combobox");
    input.setAttribute("aria-controls", list.id);
    input.setAttribute("aria-expanded", "false");
    input.setAttribute("aria-autocomplete", "list");
    let items = [], sel = -1;

    function show(html) {
      list.innerHTML = html;
      list.hidden = !html;
      input.setAttribute("aria-expanded", String(!!html));
      items = [...list.querySelectorAll("li[data-href]")];
      sel = -1;
    }
    function mark(i) {
      items.forEach((li, k) => li.setAttribute("aria-selected", String(k === i)));
      sel = i;
    }
    async function update() {
      const q = input.value.trim();
      if (!q) return show("");
      const idx = await index();
      const looksVariant = /^rs\d+$/i.test(q) || /^\d{2,9}$/.test(q) || /:/.test(q) || /\s/.test(q);
      const rows = [];
      if (looksVariant) rows.push(`<li data-go="1" data-href=""><span class="sym">Look up</span><span class="nm">${esc(q)}</span><span class="hint">Enter</span></li>`);
      for (const g of suggestions(idx, q)) {
        rows.push(`<li role="option" data-href="${geneUrl(g[0])}"><span class="sym">${esc(g[0])}</span><span class="nm">${esc(g[3])}</span><span class="hint">${fmt(g[4])} variants</span></li>`);
      }
      if (!rows.length) rows.push(`<li class="msg">No gene matches “${esc(q)}”. Try a gene symbol, an rsID such as rs113993960, a ClinVar variation ID, or an HGVS name.</li>`);
      show(rows.join(""));
    }
    async function go(useSel) {
      if (useSel && sel >= 0 && items[sel].dataset.href) { location.href = items[sel].dataset.href; return; }
      const r = await resolve(input.value).catch(() => ({ kind: "error" }));
      if (r.kind === "gene") location.href = geneUrl(r.gene);
      else if (r.kind === "variant" && r.genes.length === 1) location.href = geneUrl(r.genes[0], r.q);
      else if (r.kind === "variant") {
        show(r.genes.map(g => `<li data-href="${geneUrl(g, r.q)}"><span class="sym">${esc(g)}</span><span class="nm">${esc(r.q)} is listed in this gene</span></li>`).join(""));
      } else if (r.kind === "miss") show(`<li class="msg">${esc(r.what)} is not among the 349,056 variants in this release.</li>`);
      else if (r.kind === "error") show(`<li class="msg">Search could not load its index. Check the connection and try again.</li>`);
      else {
        const idx = await index();
        const s = suggestions(idx, input.value);
        if (s.length) location.href = geneUrl(s[0][0]);
        else show(`<li class="msg">No gene or variant matches “${esc(input.value)}”.</li>`);
      }
    }
    input.addEventListener("input", update);
    input.addEventListener("focus", () => { if (input.value) update(); });
    input.addEventListener("keydown", e => {
      if (e.key === "ArrowDown" && items.length) { e.preventDefault(); mark(Math.min(items.length - 1, sel + 1)); }
      else if (e.key === "ArrowUp" && items.length) { e.preventDefault(); mark(Math.max(0, sel - 1)); }
      else if (e.key === "Enter") { e.preventDefault(); go(true); }
      else if (e.key === "Escape") show("");
    });
    list.addEventListener("mousedown", e => {
      const li = e.target.closest("li");
      if (!li || li.classList.contains("msg")) return;
      e.preventDefault();
      if (li.dataset.go) go(false); else if (li.dataset.href) location.href = li.dataset.href;
    });
    document.addEventListener("click", e => { if (!root.contains(e.target)) show(""); });
  }

  document.addEventListener("DOMContentLoaded", () => {
    document.querySelectorAll("[data-search]").forEach(attachSearch);
  });

  return { EV_KEYS, EV_LABEL, CONS_LABEL, getJSON, index, cssVar, evColor, evVar, fmt, fmtLR, esc, el, evChip, proteinKey, geneUrl, fillRelease };
})();
