/* "How to use" and "Download and cite" pages: fill in numbers from the index. */
"use strict";

// highlight the section being read in the side contents
document.addEventListener("DOMContentLoaded", () => {
  const links = [...document.querySelectorAll(".toc a[href^='#']")];
  if (!links.length || !("IntersectionObserver" in window)) return;
  const byId = new Map(links.map(a => [a.getAttribute("href").slice(1), a]));
  const io = new IntersectionObserver(entries => {
    for (const e of entries) {
      if (!e.isIntersecting) continue;
      links.forEach(a => a.classList.remove("on"));
      const a = byId.get(e.target.id);
      if (a) a.classList.add("on");
    }
  }, { rootMargin: "-20% 0px -70% 0px" });
  document.querySelectorAll(".doc section[id]").forEach(s => io.observe(s));
  links.forEach(a => a.addEventListener("click", () => {
    links.forEach(x => x.classList.remove("on"));
    a.classList.add("on");
  }));
});

(async () => {
  const { esc, fmt, EV_LABEL, EV_KEYS, CONS_LABEL } = TV;
  let idx;
  try { idx = await TV.index(); } catch (e) { return; }
  const S = idx.summary;
  TV.fillRelease(S);

  const set = (sel, v) => document.querySelectorAll(sel).forEach(e => { e.textContent = v; });
  set(".n-am", fmt(S.has_alphamissense));
  set(".n-train", fmt(S.in_training));
  if (S.class_level_missense != null) set(".n-cls", fmt(S.class_level_missense));

  /* ---------- evidence levels on held-out genes ---------- */
  const cal = document.getElementById("caltable");
  if (cal) {
    const C = S.calibration;
    const order = EV_LABEL.map(s => s.toLowerCase());
    const rows = order.map(b => C.rows.find(r => r.bin === b)).filter(Boolean);
    const priorOdds = C.base_rate / (100 - C.base_rate);
    // which level the observed rate actually corresponds to, at the same thresholds
    const levelOf = p => {
      const lr = (Math.min(p, 99.95) / (100 - Math.min(p, 99.95))) / priorOdds;
      if (lr >= 350) return 0; if (lr >= 18.7) return 1; if (lr >= 4.33) return 2; if (lr >= 2.08) return 3;
      if (lr > 1 / 2.08) return 4; if (lr > 1 / 4.33) return 5; if (lr > 1 / 18.7) return 6; return 7;
    };
    cal.innerHTML = `
      <thead><tr><th>Evidence level</th><th class="n">Variants tested</th><th class="barcell"></th><th class="pct">Pathogenic</th><th>How it performed</th></tr></thead>
      <tbody>${rows.map(r => {
        const i = order.indexOf(r.bin), got = levelOf(r.observed);
        const perf = got === i ? `<span class="muted">As labelled</span>`
          : `<b class="warn-t">Like ${EV_LABEL[got].toLowerCase()}</b>`;
        return `<tr><td>${TV.evChip(i)}</td><td class="n">${fmt(r.n)}</td>
          <td class="barcell"><div class="track"><div class="fill" style="width:${r.observed}%;background:var(--${EV_KEYS[i]})"></div><div class="base" style="left:${C.base_rate}%"></div></div></td>
          <td class="pct">${r.observed.toFixed(1)}%</td><td>${perf}</td></tr>`;
      }).join("")}</tbody>
      <caption>${fmt(C.n)} labelled variants in ${fmt(C.genes)} test genes. The dashed line marks the overall share of pathogenic variants in the test set, ${C.base_rate}%. "How it performed" converts each share back into a likelihood ratio and reads it against the same thresholds.</caption>`;
  }

  /* ---------- AlphaMissense coverage bar ---------- */
  const cov = document.getElementById("cov");
  if (cov) {
    const cons = S.by_consequence.slice().sort((a, b) => b[1] - a[1]);
    const palette = ["var(--ink-2)", "#8c9aa1", "#a9b4b9", "#c4cdd0", "var(--lipid-head)", "#ddb95c", "var(--lipid)", "#efd894",
      "#d6c7a1", "#bfae86", "#a7976f", "#8f8159", "#776b45"];
    cov.innerHTML = `
      <div class="cov-bar" role="img" aria-label="Variants by type; AlphaMissense covers ${fmt(S.has_alphamissense)} of ${fmt(S.variants)}, all missense">
        ${cons.map(([c, n, am], i) => `<div style="flex:${n};background:${palette[i]}" title="${esc(CONS_LABEL[c])}: ${fmt(n)}">${am ? `<span class="am" style="width:${(am / n * 100).toFixed(2)}%"></span>` : ""}</div>`).join("")}
      </div>
      <div class="cov-legend">${cons.map(([c, n], i) => `<div><span class="sw" style="background:${palette[i]}"></span>${esc(CONS_LABEL[c])} <b>${fmt(n)}</b></div>`).join("")}</div>`;
  }

  /* ---------- downloads ---------- */
  const d = S.downloads;
  const mb = b => (b / 1e6).toFixed(1) + " MB";
  const dn = document.getElementById("dl-native"), dr = document.getElementById("dl-routed");
  if (dn && d) {
    dn.textContent = `Download, ${mb(d.native.bytes)}`;
    dr.textContent = `Download, ${mb(d.routed.bytes)}`;
    document.getElementById("sha-native").textContent = `SHA-256 ${d.native.sha256}`;
    document.getElementById("sha-routed").textContent = `SHA-256 ${d.routed.sha256}`;
  }
  const up = document.getElementById("up-rel");
  if (up && S.uniprot_release) up.textContent = ` release ${S.uniprot_release}`;
})();
