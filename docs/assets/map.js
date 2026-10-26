/* The protein map: a membrane-topology trace of the chain, coloured by a regional
   summary of the predictions, plus the AlphaFold structure coloured the same way.

   The regional values are kernel averages of the predictions for the ClinVar
   missense variants that exist at nearby residues. They are not a saturation
   scan: residues with no scored variant nearby are left grey. */
"use strict";

const TVMap = (() => {
  const { cssVar, esc, EV_LABEL, fmtLR } = TV;

  /* ---------- colour ---------- */
  function hex(c) {
    c = c.trim();
    if (c.startsWith("#")) {
      if (c.length === 4) c = "#" + [...c.slice(1)].map(x => x + x).join("");
      return [parseInt(c.slice(1, 3), 16), parseInt(c.slice(3, 5), 16), parseInt(c.slice(5, 7), 16)];
    }
    const m = c.match(/\d+(\.\d+)?/g);
    return m ? m.slice(0, 3).map(Number) : [128, 128, 128];
  }
  const toHex = rgb => "#" + rgb.map(v => Math.round(Math.max(0, Math.min(255, v))).toString(16).padStart(2, "0")).join("");
  function ramp(stops) {
    const s = stops.map(([x, c]) => [x, hex(c)]);
    return v => {
      if (v <= s[0][0]) return toHex(s[0][1]);
      for (let k = 1; k < s.length; k++) {
        if (v <= s[k][0]) {
          const t = (v - s[k - 1][0]) / (s[k][0] - s[k - 1][0]);
          return toHex(s[k - 1][1].map((a, j) => a + (s[k][1][j] - a) * t));
        }
      }
      return toHex(s[s.length - 1][1]);
    };
  }

  // log10 LR thresholds behind the ACMG combining rules (Tavtigian 2018)
  const T = { vs: Math.log10(350), s: Math.log10(18.7), m: Math.log10(4.33), p: Math.log10(2.08) };
  const METRICS = {
    patho: {
      label: "Predicted pathogenicity", short: "Pathogenicity",
      value: v => v.lr == null ? null : Math.max(-2.6, Math.min(2.6, Math.log10(v.lr))),
      stops: () => [[-T.vs, cssVar("--bs")], [-T.s, cssVar("--bm")], [-T.p, cssVar("--bp")], [0, cssVar("--ind")],
        [T.p, cssVar("--pp")], [T.m, cssVar("--pm")], [T.s, cssVar("--ps")], [T.vs, cssVar("--pvs")]],
      ends: ["benign", "pathogenic"],
      describe: x => `${cap(binName(x))}, LR ${fmtLR(10 ** x)}`,
    },
    act: {
      label: "Predicted activity", short: "Activity",
      value: v => v.act,
      stops: () => [[0.2, cssVar("--act-lo")], [0.5, cssVar("--act-mid")], [0.8, cssVar("--act-hi")]],
      ends: ["activity lost", "normal activity"],
      describe: x => `${x.toFixed(2)} (0 none, 1 normal)`,
    },
    mech: {
      label: "Mechanism", short: "Mechanism",
      exploratory: true,
      value: v => v.mech == null ? null : Math.max(0, Math.min(1, v.mech)),
      stops: () => [[0.1, cssVar("--mech-expr")], [0.5, cssVar("--ind")], [0.9, cssVar("--mech-act")]],
      ends: ["expression-driven", "activity-specific"],
      describe: x => `${x.toFixed(2)}, ${x >= 0.6 ? "leans activity-specific" : x <= 0.4 ? "leans expression-driven" : "unclear"}`,
    },
  };
  const cap = s => s.charAt(0).toUpperCase() + s.slice(1);
  function binName(x) {
    const lr = 10 ** x;
    if (lr >= 350) return EV_LABEL[0].toLowerCase();
    if (lr >= 18.7) return EV_LABEL[1].toLowerCase();
    if (lr >= 4.33) return EV_LABEL[2].toLowerCase();
    if (lr >= 2.08) return EV_LABEL[3].toLowerCase();
    if (lr > 1 / 2.08) return "indeterminate";
    if (lr > 1 / 4.33) return EV_LABEL[5].toLowerCase();
    if (lr > 1 / 18.7) return EV_LABEL[6].toLowerCase();
    return EV_LABEL[7].toLowerCase();
  }
  const PLDDT = v => v >= 90 ? "#0053d6" : v >= 70 ? "#65cbf3" : v >= 50 ? "#ffdb13" : "#ff7d45";

  /* ---------- topology ---------- */
  const OUT_RE = /extracellular|lumenal|luminal|periplasmic|vesicular|exoplasmic|intermembrane/i;
  const IN_RE = /cytoplasmic|matrix|nuclear|stromal/i;

  function topology(meta) {
    const L = meta.length;
    const st = new Array(L + 1).fill("u");
    let hasTM = false, annotated = false;
    for (const [type, a, b, d] of meta.features) {
      if (type === "Transmembrane" || type === "Intramembrane") {
        for (let i = a; i <= b && i <= L; i++) st[i] = type === "Transmembrane" ? "t" : "m";
        if (type === "Transmembrane") hasTM = true;
      }
    }
    for (const [type, a, b, d] of meta.features) {
      if (type !== "Topological domain") continue;
      const side = OUT_RE.test(d) ? "o" : IN_RE.test(d) ? "i" : null;
      if (!side) continue;
      annotated = true;
      for (let i = a; i <= b && i <= L; i++) if (st[i] === "u") st[i] = side;
    }
    // runs
    const runs = [];
    for (let i = 1; i <= L; i++) {
      const s = st[i];
      if (runs.length && runs[runs.length - 1].s === s) runs[runs.length - 1].b = i;
      else runs.push({ s, a: i, b: i });
    }
    let inferred = false;
    if (!hasTM) {
      for (const r of runs) if (r.s !== "m") r.s = r.s === "o" ? "o" : "i";
    } else {
      // unknown non-membrane runs: flip side across each helix from the nearest known run
      const isMem = r => r.s === "t";
      for (let k = 0; k < runs.length; k++) {
        if (runs[k].s !== "u") continue;
        inferred = true;
        let side = null;
        for (let j = k - 1, flips = 0; j >= 0; j--) {
          if (isMem(runs[j])) flips++;
          else if (runs[j].s === "o" || runs[j].s === "i") { side = flips % 2 ? (runs[j].s === "o" ? "i" : "o") : runs[j].s; break; }
        }
        if (!side) {
          for (let j = k + 1, flips = 0; j < runs.length; j++) {
            if (isMem(runs[j])) flips++;
            else if (runs[j].s === "o" || runs[j].s === "i") { side = flips % 2 ? (runs[j].s === "o" ? "i" : "o") : runs[j].s; break; }
          }
        }
        if (!side) {
          let flips = 0;
          for (let j = 0; j < k; j++) if (isMem(runs[j])) flips++;
          side = flips % 2 ? "o" : "i";
        }
        runs[k].s = side;
      }
    }
    // merge adjacent equal runs created by inference
    const merged = [];
    for (const r of runs) {
      if (merged.length && merged[merged.length - 1].s === r.s) merged[merged.length - 1].b = r.b;
      else merged.push({ ...r });
    }
    for (const r of merged) for (let i = r.a; i <= r.b; i++) st[i] = r.s;
    // region names
    let tm = 0;
    const regionOf = new Array(L + 1);
    merged.forEach((r, k) => {
      let name;
      if (r.s === "t") name = `Transmembrane helix ${++tm}`;
      else if (r.s === "m") name = "Intramembrane segment";
      else {
        const where = r.s === "o" ? "Outside (extracellular or lumenal)" : "Cytoplasmic";
        name = !hasTM ? (annotated ? where : "No membrane topology annotated")
          : k === 0 ? `N-terminus, ${where.toLowerCase()}`
          : k === merged.length - 1 ? `C-terminus, ${where.toLowerCase()}`
          : `${r.s === "o" ? "Outside" : "Cytoplasmic"} loop after helix ${tm}`;
      }
      r.name = name; r.tm = r.s === "t" ? tm : null;
      for (let i = r.a; i <= r.b; i++) regionOf[i] = name;
    });
    return { st, runs: merged, hasTM, annotated, inferred, regionOf, nTM: tm };
  }

  /* ---------- regional averages ---------- */
  function seqProfile(L, pts, sigma = 4, minW = 0.35) {
    const val = new Float64Array(L + 1).fill(NaN), w = new Float64Array(L + 1);
    const num = new Float64Array(L + 1);
    const R = Math.ceil(sigma * 3);
    for (const [p, v] of pts) {
      for (let i = Math.max(1, p - R); i <= Math.min(L, p + R); i++) {
        const k = Math.exp(-((i - p) ** 2) / (2 * sigma * sigma));
        w[i] += k; num[i] += k * v;
      }
    }
    for (let i = 1; i <= L; i++) if (w[i] >= minW) val[i] = num[i] / w[i];
    return { val, w };
  }
  function spaceProfile(ca, pts, sigma = 5, minW = 0.35) {
    // ca: Map resi -> [x,y,z]; pts: [pos, value]
    const out = new Map();
    const P = pts.map(([p, v]) => [ca.get(p), v]).filter(x => x[0]);
    const cut2 = (sigma * 2.6) ** 2;
    for (const [r, c] of ca) {
      let w = 0, n = 0;
      for (const [q, v] of P) {
        const d2 = (c[0] - q[0]) ** 2 + (c[1] - q[1]) ** 2 + (c[2] - q[2]) ** 2;
        if (d2 > cut2) continue;
        const k = Math.exp(-d2 / (2 * sigma * sigma));
        w += k; n += k * v;
      }
      if (w >= minW) out.set(r, n / w);
    }
    return out;
  }
  function hotspots(L, prof, topo, n = 5, sep = 18, minW = 1.0, highest = true) {
    const cand = [];
    for (let i = 1; i <= L; i++) {
      const v = prof.val[i];
      if (Number.isNaN(v) || prof.w[i] < minW) continue;
      const left = i > 1 ? prof.val[i - 1] : -Infinity, right = i < L ? prof.val[i + 1] : -Infinity;
      const isPeak = highest ? v >= (Number.isNaN(left) ? -Infinity : left) && v >= (Number.isNaN(right) ? -Infinity : right)
        : v <= (Number.isNaN(left) ? Infinity : left) && v <= (Number.isNaN(right) ? Infinity : right);
      if (isPeak) cand.push([i, v]);
    }
    cand.sort((a, b) => highest ? b[1] - a[1] : a[1] - b[1]);
    const picked = [];
    for (const [i, v] of cand) {
      if (picked.every(p => Math.abs(p.at - i) >= sep)) {
        let a = i, b = i;
        const lim = highest ? v - 0.3 : v + 0.3;
        const ok = x => !Number.isNaN(prof.val[x]) && (highest ? prof.val[x] >= lim : prof.val[x] <= lim);
        while (a > 1 && ok(a - 1) && i - a < 25) a--;
        while (b < L && ok(b + 1) && b - i < 25) b++;
        picked.push({ at: i, a, b, v, region: topo.regionOf[i] });
      }
      if (picked.length >= n) break;
    }
    return picked;
  }

  /* ---------- the 2D trace ---------- */
  class Wave {
    constructor(host, { meta, variants, compact = false, onSelect = null, onZoom = null }) {
      this.host = host; this.meta = meta; this.compact = compact;
      this.onSelect = onSelect; this.onZoom = onZoom;
      this.L = meta.length;
      this.topo = topology(meta);
      this.metric = "patho";
      this.showVariants = true;
      this.sel = null;
      this.z = [1, this.L];
      this.plddt = null;
      this.W = 1200;
      this.yT = 118; this.yB = 176;
      this.hasLane = !compact && meta.features.some(f => f[0] === "Domain" || f[0] === "Motif");
      this.H = this.hasLane ? 366 : 328;
      this.X0 = 74; this.X1 = this.W - 14;
      this.layoutY();
      this.setVariants(variants);
      host.classList.add("wave");
      this.tip = TV.el("div", { class: "tip", role: "status" });
      this.svgHost = TV.el("div");
      host.append(this.svgHost, this.tip);
      // on small screens the map is wider than the page and scrolls sideways; say so
      host.insertAdjacentElement("afterend", TV.el("p", { class: "swipe-hint" }, "Swipe sideways to see the whole protein."));
      this.bind();
    }
    setVariants(vs) {
      this.vars = vs;
      this.byPos = new Map();
      for (const v of vs) {
        if (!this.byPos.has(v.pos)) this.byPos.set(v.pos, []);
        this.byPos.get(v.pos).push(v);
      }
      for (const arr of this.byPos.values()) arr.sort((a, b) => a.ev - b.ev);
      this.profiles = {};
      for (const [k, m] of Object.entries(METRICS)) {
        const pts = vs.map(v => [v.pos, m.value(v)]).filter(x => x[1] != null && !Number.isNaN(x[1]));
        this.profiles[k] = seqProfile(this.L, pts);
      }
    }
    layoutY() {
      const { st, runs } = this.topo, yT = this.yT, yB = this.yB, L = this.L;
      const y = new Float64Array(L + 2);
      const sideY = s => s === "o" ? yT : yB;
      const h = n => Math.min(92, 12 + n * 0.85);
      runs.forEach((r, k) => {
        const n = r.b - r.a + 1;
        const prev = runs[k - 1], next = runs[k + 1];
        if (r.s === "t") {
          let from = prev ? (prev.s === "o" ? yT : yB) : (next && next.s === "o" ? yB : yT);
          if (prev && prev.s === "m") from = y[r.a - 1] < (yT + yB) / 2 ? yT : yB;
          const to = from === yT ? yB : yT;
          for (let i = r.a; i <= r.b; i++) y[i] = from + (to - from) * ((i - r.a + 0.5) / n);
        } else if (r.s === "m") {
          const from = prev && prev.s === "o" ? yT : yB;
          const depth = (yB - yT) * 0.62;
          for (let i = r.a; i <= r.b; i++) {
            const t = (i - r.a + 0.5) / n;
            y[i] = from + (from === yT ? 1 : -1) * depth * Math.sin(Math.PI * t);
          }
        } else if (!this.topo.hasTM) {
          for (let i = r.a; i <= r.b; i++) y[i] = r.s === "o" ? yT - 40 : yB + 40;
        } else {
          const base = sideY(r.s), dir = r.s === "o" ? -1 : 1, H = h(n);
          for (let i = r.a; i <= r.b; i++) {
            const t = (i - r.a + 0.5) / n;
            // termini leave the membrane and level off; internal loops arch between helices
            let f;
            const R = Math.min(n, 40);
            if (k === 0) f = Math.sin(Math.PI / 2 * Math.min(1, (r.b - i + 0.5) / R));
            else if (k === runs.length - 1) f = Math.sin(Math.PI / 2 * Math.min(1, (i - r.a + 0.5) / R));
            else f = Math.sin(Math.PI * t);
            y[i] = base + dir * H * f;
          }
        }
      });
      this.y = y;
    }
    x(i) { return this.X0 + (i - this.z[0]) / Math.max(1, this.z[1] - this.z[0]) * (this.X1 - this.X0); }
    resAt(px) { return Math.round(this.z[0] + (px - this.X0) / (this.X1 - this.X0) * (this.z[1] - this.z[0])); }
    setMetric(m) { this.metric = m; this.render(); }
    setPlddt(arr) { this.plddt = arr; if (this.metric === "plddt") this.render(); }
    select(pos, part = null) { this.sel = pos; this.part = pos ? part : null; this.render(); }
    // the stretch of chain a position belongs to: its helix or loop, or a window if the loop is long
    partAround(p) {
      const run = this.topo.runs.find(r => r.a <= p && p <= r.b);
      if (!run) return { range: [Math.max(1, p - 10), Math.min(this.L, p + 10)], name: `Residues near ${p}` };
      if (run.b - run.a + 1 <= 45) return { range: [run.a, run.b], name: run.name };
      const a = Math.max(run.a, p - 12), b = Math.min(run.b, p + 12);
      return { range: [a, b], name: `Part of the ${run.name.charAt(0).toLowerCase()}${run.name.slice(1)}` };
    }
    zoom(a, b) {
      if (a == null) this.z = [1, this.L];
      else {
        a = Math.max(1, Math.floor(a)); b = Math.min(this.L, Math.ceil(b));
        if (b - a < 12) { const c = (a + b) / 2; a = Math.max(1, Math.round(c - 6)); b = Math.min(this.L, a + 12); }
        this.z = [a, b];
      }
      this.render();
      if (this.onZoom) this.onZoom(this.z[0] !== 1 || this.z[1] !== this.L);
    }
    colorFn() {
      const m = this.metric;
      if (m === "topo") {
        const c = { o: cssVar("--bm"), i: cssVar("--ink-3"), t: cssVar("--lipid-head"), m: cssVar("--lipid-head") };
        return i => c[this.topo.st[i]] || cssVar("--rule");
      }
      if (m === "plddt") {
        return i => this.plddt && this.plddt[i] != null ? PLDDT(this.plddt[i]) : cssVar("--rule");
      }
      const r = ramp(METRICS[m].stops()), p = this.profiles[m].val, none = cssVar("--rule");
      return i => Number.isNaN(p[i]) ? none : r(p[i]);
    }
    render() {
      const { W, H, yT, yB, L } = this;
      const [z0, z1] = this.z;
      const lo = Math.max(1, Math.floor(z0)), hi = Math.min(L, Math.ceil(z1));
      const pxPer = (this.X1 - this.X0) / Math.max(1, z1 - z0);
      const out = [];
      out.push(`<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="Membrane topology map of ${esc(this.meta.gene)}, residues ${lo} to ${hi}">`);
      out.push(`<defs><clipPath id="clip-${this.meta.gene}"><rect x="${this.X0 - 6}" y="0" width="${this.X1 - this.X0 + 12}" height="${H}"/></clipPath></defs>`);
      // membrane
      if (this.topo.hasTM) {
        const lip = cssVar("--lipid"), head = cssVar("--lipid-head");
        out.push(`<rect x="${this.X0 - 8}" y="${yT}" width="${this.X1 - this.X0 + 16}" height="${yB - yT}" fill="${lip}" opacity="0.42"/>`);
        let heads = "";
        for (let x = this.X0 - 5; x <= this.X1 + 5; x += 7.2) heads += `<circle cx="${x.toFixed(1)}" cy="${yT + 3}" r="3"/><circle cx="${x.toFixed(1)}" cy="${yB - 3}" r="3"/>`;
        out.push(`<g fill="${head}" opacity="0.8">${heads}</g>`);
        out.push(`<text x="6" y="${yT - 10}" class="zl">Outside</text><text x="6" y="${yB + 20}" class="zl">Cytoplasm</text><text x="6" y="${(yT + yB) / 2 + 4}" class="zl">Membrane</text>`);
      }
      out.push(`<g clip-path="url(#clip-${this.meta.gene})">`);
      // the selected part, shaded so it matches what is highlighted in 3D
      if (this.part && this.part[1] >= lo && this.part[0] <= hi) {
        const xa = this.x(Math.max(lo, this.part[0]) - 0.5), xb = this.x(Math.min(hi, this.part[1]) + 0.5);
        const hb = (this.hasLane ? H - 52 : H - 30) - 2;
        out.push(`<rect x="${xa.toFixed(1)}" y="2" width="${Math.max(2, xb - xa).toFixed(1)}" height="${hb}" class="partband"/>`);
        out.push(`<line x1="${xa.toFixed(1)}" x2="${xa.toFixed(1)}" y1="2" y2="${hb + 2}" class="partedge"/><line x1="${xb.toFixed(1)}" x2="${xb.toFixed(1)}" y1="2" y2="${hb + 2}" class="partedge"/>`);
      }
      // helix outlines and numbers
      if (this.topo.hasTM) {
        for (const r of this.topo.runs) {
          if (r.s !== "t" || r.b < lo || r.a > hi) continue;
          const xa = this.x(r.a - 0.5), xb = this.x(r.b + 0.5), w = xb - xa;
          out.push(`<rect x="${xa.toFixed(1)}" y="${yT + 7}" width="${w.toFixed(1)}" height="${yB - yT - 14}" rx="${Math.min(8, w / 2).toFixed(1)}" class="helix"/>`);
          if (w >= 22) out.push(`<text x="${(xa + w / 2).toFixed(1)}" y="${yT + 21}" class="tmn">${r.tm}</text>`);
        }
      }
      // the chain, coloured residue by residue
      const col = this.colorFn();
      const sw = Math.max(3.2, Math.min(9, pxPer * 0.9));
      let under = "", segs = "";
      let prev = null;
      for (let i = lo; i <= hi; i++) {
        const xi = this.x(i), yi = this.y[i];
        if (prev) {
          under += `L${xi.toFixed(1)} ${yi.toFixed(1)}`;
        } else under += `M${xi.toFixed(1)} ${yi.toFixed(1)}`;
        if (i < hi) {
          const xn = this.x(i + 1), yn = this.y[i + 1];
          segs += `<line x1="${xi.toFixed(1)}" y1="${yi.toFixed(1)}" x2="${xn.toFixed(1)}" y2="${yn.toFixed(1)}" stroke="${col(i)}"/>`;
        }
        prev = i;
      }
      out.push(`<path d="${under}" fill="none" stroke="var(--surface)" stroke-width="${(sw + 3).toFixed(1)}" stroke-linejoin="round" stroke-linecap="round"/>`);
      out.push(`<g stroke-width="${sw.toFixed(1)}" stroke-linecap="round">${segs}</g>`);
      // UniProt functional sites
      let sites = "";
      for (const [type, a, b, d] of this.meta.features) {
        if (!["Binding site", "Active site", "Site"].includes(type)) continue;
        for (let i = a; i <= b; i++) {
          if (i < lo || i > hi) continue;
          const xi = this.x(i).toFixed(1), yi = this.y[i].toFixed(1);
          const r = Math.max(4.5, sw * 0.9);
          if (type === "Active site") sites += `<rect x="${(+xi - r * 0.75).toFixed(1)}" y="${(+yi - r * 0.75).toFixed(1)}" width="${(r * 1.5).toFixed(1)}" height="${(r * 1.5).toFixed(1)}" transform="rotate(45 ${xi} ${yi})" class="fs act"/>`;
          else if (type === "Binding site") sites += `<circle cx="${xi}" cy="${yi}" r="${r.toFixed(1)}" class="fs bind"/>`;
          else sites += `<rect x="${(+xi - r * 0.8).toFixed(1)}" y="${(+yi - r * 0.8).toFixed(1)}" width="${(r * 1.6).toFixed(1)}" height="${(r * 1.6).toFixed(1)}" class="fs site"/>`;
        }
      }
      out.push(`<g>${sites}</g>`);
      // variants, stacked outward from the membrane
      if (this.showVariants) {
        // dense proteins get smaller, shorter stacks so the regional colour stays readable;
        // zooming in restores full stacks
        const dr = Math.max(1.7, Math.min(4.6, pxPer * 0.45));
        const maxStack = pxPer < 1 ? 2 : pxPer < 2 ? 3 : 6;
        const dots = [];
        for (const [p, arr] of this.byPos) {
          if (p < lo || p > hi) continue;
          const side = this.topo.st[p];
          const yi = this.y[p];
          const dir = side === "o" ? -1 : side === "i" ? 1 : (yi < (yT + yB) / 2 ? -1 : 1);
          const n = Math.min(arr.length, maxStack);
          for (let k = 0; k < n; k++) {
            const v = arr[k];
            dots.push([v.ev, `<circle cx="${this.x(p).toFixed(1)}" cy="${(yi + dir * (sw / 2 + dr + 1.5 + k * (dr * 2 + 0.6))).toFixed(1)}" r="${dr.toFixed(1)}" fill="var(--${TV.EV_KEYS[v.ev]})"/>`]);
          }
        }
        dots.sort((a, b) => b[0] - a[0]);
        out.push(`<g class="vdots">${dots.map(d => d[1]).join("")}</g>`);
      }
      // selection
      if (this.sel && this.sel >= lo && this.sel <= hi) {
        const xs = this.x(this.sel).toFixed(1);
        out.push(`<line x1="${xs}" x2="${xs}" y1="4" y2="${this.hasLane ? H - 52 : H - 30}" class="selline"/>`);
        out.push(`<circle cx="${xs}" cy="${this.y[this.sel].toFixed(1)}" r="${(sw + 4).toFixed(1)}" class="selring"/>`);
      }
      out.push(`</g>`);
      // domains and motifs
      let axisY = H - 18;
      if (this.hasLane) {
        const dy = H - 54;
        let lane = "";
        for (const [type, a, b, d] of this.meta.features) {
          if (type !== "Domain" && type !== "Motif") continue;
          if (b < lo || a > hi) continue;
          const xa = Math.max(this.X0, this.x(a - 0.5)), xb = Math.min(this.X1, this.x(b + 0.5)), w = xb - xa;
          if (w <= 0) continue;
          const isM = type === "Motif";
          lane += `<g><title>${esc(type)}: ${esc(d)} (${a}–${b})</title><rect x="${xa.toFixed(1)}" y="${isM ? dy + 13 : dy}" width="${Math.max(2, w).toFixed(1)}" height="${isM ? 6 : 18}" rx="2" class="${isM ? "motif" : "dom"}"/>`;
          if (!isM && w > 60) lane += `<text x="${(xa + 6).toFixed(1)}" y="${dy + 13}" class="domt">${esc(d.length * 6.2 > w - 10 ? d.slice(0, Math.max(3, Math.floor((w - 16) / 6.2))) + "…" : d)}</text>`;
          lane += `</g>`;
        }
        if (lane) out.push(`<text x="6" y="${dy + 13}" class="zl">Domains</text>${lane}`);
      }
      // axis
      const span = z1 - z0;
      const step = [1, 2, 5, 10, 20, 25, 50, 100, 200, 250, 500, 1000].find(s => span / s <= 9) || 1000;
      let ax = `<line x1="${this.X0}" x2="${this.X1}" y1="${axisY - 9}" y2="${axisY - 9}" class="axis"/>`;
      for (let t = Math.ceil(lo / step) * step; t <= hi; t += step) {
        if (t < 1) continue;
        const xt = this.x(t).toFixed(1);
        ax += `<line x1="${xt}" x2="${xt}" y1="${axisY - 9}" y2="${axisY - 5}" class="axis"/><text x="${xt}" y="${axisY + 6}" class="tick">${t}</text>`;
      }
      out.push(ax + `<text x="6" y="${axisY + 6}" class="zl">Residue</text>`);
      out.push(`<rect class="brush" x="0" y="0" width="0" height="${H}" visibility="hidden"/>`);
      out.push(`</svg>`);
      this.svgHost.innerHTML = out.join("");
      this.svg = this.svgHost.querySelector("svg");
      this.brush = this.svg.querySelector(".brush");
    }
    toLocal(e) {
      const pt = this.svg.createSVGPoint();
      pt.x = e.clientX; pt.y = e.clientY;
      return pt.matrixTransform(this.svg.getScreenCTM().inverse());
    }
    tipFor(p) {
      const seq = this.meta.sequence;
      const lines = [`<b>${seq[p - 1]}${p}</b>, ${esc(this.topo.regionOf[p])}`];
      for (const [k, m] of Object.entries(METRICS)) {
        const v = this.profiles[k].val[p];
        if (!Number.isNaN(v) && (k === this.metric || k === "patho")) lines.push(`${m.short}${m.exploratory ? " (exploratory)" : ""} nearby: ${m.describe(v)}`);
      }
      if (this.metric === "plddt" && this.plddt && this.plddt[p] != null) lines.push(`AlphaFold confidence (pLDDT) ${this.plddt[p].toFixed(0)}`);
      const feats = this.meta.features.filter(f => ["Binding site", "Active site", "Site", "Motif"].includes(f[0]) && f[1] <= p && p <= f[2]);
      for (const f of feats) lines.push(`${esc(f[0])}: ${esc(f[3])}`);
      const vs = this.byPos.get(p);
      if (vs) lines.push(vs.slice(0, 4).map(v => `${esc(v.pv)} ${EV_LABEL[v.ev].toLowerCase()}`).join("<br>") + (vs.length > 4 ? `<br>and ${vs.length - 4} more` : ""));
      else if (Number.isNaN(this.profiles.patho.val[p])) lines.push("Too few variants nearby to say");
      return lines.join("<br>");
    }
    bind() {
      let down = null;
      const host = this.host;
      host.addEventListener("pointermove", e => {
        if (!this.svg) return;
        const pt = this.toLocal(e);
        if (down && down.pointerType === "mouse") {
          const a = Math.min(down.x, pt.x), b = Math.max(down.x, pt.x);
          if (b - a > 6) {
            this.brush.setAttribute("x", Math.max(this.X0, a)); this.brush.setAttribute("width", Math.min(this.X1, b) - Math.max(this.X0, a));
            this.brush.setAttribute("visibility", "visible");
          }
          this.tip.style.display = "none";
          return;
        }
        if (pt.x < this.X0 - 4 || pt.x > this.X1 + 4) { this.tip.style.display = "none"; return; }
        const p = Math.max(1, Math.min(this.L, this.resAt(pt.x)));
        this.tip.innerHTML = this.tipFor(p);
        this.tip.style.display = "block";
        const r = host.getBoundingClientRect();
        let left = e.clientX - r.left + 14;
        if (left + 300 > r.width) left = Math.max(4, e.clientX - r.left - 314);
        this.tip.style.left = left + host.scrollLeft + "px";
        this.tip.style.top = Math.max(4, e.clientY - r.top - 10) + "px";
      });
      host.addEventListener("pointerleave", () => { this.tip.style.display = "none"; });
      host.addEventListener("pointerdown", e => {
        if (!this.svg || e.button !== 0) return;
        const pt = this.toLocal(e);
        down = { x: pt.x, pointerType: e.pointerType };
      });
      window.addEventListener("pointerup", e => {
        if (!down) return;
        const pt = this.toLocal(e);
        const d = down; down = null;
        this.brush.setAttribute("visibility", "hidden");
        if (d.pointerType === "mouse" && Math.abs(pt.x - d.x) > 8 && !this.compact) {
          const a = this.resAt(Math.max(this.X0, Math.min(d.x, pt.x))), b = this.resAt(Math.min(this.X1, Math.max(d.x, pt.x)));
          this.zoom(a, b);
        } else if (pt.x >= this.X0 - 4 && pt.x <= this.X1 + 4 && host.contains(e.target)) {
          const p = Math.max(1, Math.min(this.L, this.resAt(pt.x)));
          if (this.onSelect) this.onSelect(p);
        }
      });
      host.addEventListener("dblclick", () => { if (!this.compact) this.zoom(null); });
    }
  }

  /* ---------- the 3D structure ---------- */
  let molLib = null;
  function load3Dmol() {
    if (window.$3Dmol) return Promise.resolve();
    if (!molLib) molLib = new Promise((res, rej) => {
      const s = document.createElement("script");
      s.src = "assets/vendor/3Dmol-min.js";
      s.onload = res; s.onerror = () => rej(new Error("3Dmol.js failed to load"));
      document.head.appendChild(s);
    });
    return molLib;
  }

  class Mol {
    constructor(host, { meta, wave, onSelect }) {
      this.host = host; this.meta = meta; this.wave = wave; this.onSelect = onSelect;
      this.msg = host.querySelector(".molmsg");
      this.metric = "patho"; this.sel = null; this.mapped = false;
    }
    say(t) {
      this.msg.innerHTML = t; this.msg.hidden = !t;
      this.host.querySelector(".molhint").hidden = !!t && !this.viewer;
    }
    async load() {
      const acc = this.meta.uniprot;
      this.say("Loading the AlphaFold model…");
      let entries;
      try {
        const r = await fetch(`https://alphafold.ebi.ac.uk/api/prediction/${encodeURIComponent(acc)}`);
        if (r.status === 404 || r.status === 400) {
          this.say(this.meta.length > 2700
            ? `AlphaFold DB has no single-chain model for ${esc(this.meta.gene)}; it publishes whole models only for proteins up to 2,700 residues. The map above still covers the full sequence.`
            : `AlphaFold DB has no model for UniProt ${esc(acc)}. The map above still covers the full sequence.`);
          return;
        }
        if (!r.ok) throw new Error(`HTTP ${r.status}`);
        entries = await r.json();
      } catch (e) {
        this.say(`The AlphaFold model could not be fetched from alphafold.ebi.ac.uk (${esc(e.message)}). Reload to try again.`);
        return;
      }
      const ent = entries.find(x => x.uniprotAccession === acc) || entries[0];
      if (!ent || !ent.pdbUrl) { this.say(`AlphaFold DB returned no coordinates for ${esc(acc)}.`); return; }
      this.entry = ent;
      this.mapped = ent.sequence === this.meta.sequence ||
        (ent.uniprotStart === 1 && ent.sequence && this.meta.sequence.startsWith(ent.sequence) && ent.sequence.length === this.meta.length);
      let pdb;
      try {
        const [txt] = await Promise.all([fetch(ent.pdbUrl).then(r => { if (!r.ok) throw new Error(`HTTP ${r.status}`); return r.text(); }), load3Dmol()]);
        pdb = txt;
      } catch (e) {
        this.say(`The AlphaFold model could not be loaded (${esc(e.message)}). Reload to try again.`);
        return;
      }
      this.say("");
      const el = this.host.querySelector(".viewer");
      this.viewer = $3Dmol.createViewer(el, { backgroundColor: cssVar("--paper"), antialias: true, id: `mol-${this.meta.gene}` });
      // Let the mouse wheel scroll the page until the visitor clicks the structure (or holds Ctrl),
      // so nobody gets stuck zooming the model while scrolling past it.
      this.active = false;
      el.addEventListener("wheel", e => { if (!this.active && !e.ctrlKey) e.stopImmediatePropagation(); }, { capture: true });
      el.addEventListener("pointerdown", () => { this.active = true; });
      el.addEventListener("mouseleave", () => { this.active = false; });
      this.model = this.viewer.addModel(pdb, "pdb");
      this.ca = new Map();
      this.plddt = [];
      for (const a of this.model.selectedAtoms({ atom: "CA" })) {
        this.ca.set(a.resi, [a.x, a.y, a.z]);
        this.plddt[a.resi] = a.b;
      }
      this.wave.setPlddt(this.plddt);
      this.spatial = {};
      if (this.mapped) {
        for (const [k, m] of Object.entries(METRICS)) {
          const pts = this.wave.vars.map(v => [v.pos, m.value(v)]).filter(x => x[1] != null && !Number.isNaN(x[1]));
          this.spatial[k] = spaceProfile(this.ca, pts);
        }
      } else {
        this.say("");
        this.host.querySelector(".molhint").textContent = "This AlphaFold model uses a different sequence version, so only its own confidence is shown; variants are not placed on it.";
      }
      this.viewer.setHoverable({ atom: "CA" }, true, (atom, viewer) => {
        if (this.hoverLabel) viewer.removeLabel(this.hoverLabel);
        const aa = this.meta.sequence[atom.resi - 1] || "";
        this.hoverLabel = viewer.addLabel(`${aa}${atom.resi}`, { position: atom, backgroundColor: cssVar("--ink"), fontColor: cssVar("--paper"), fontSize: 12, backgroundOpacity: 0.9, inFront: true });
      }, (atom, viewer) => { if (this.hoverLabel) { viewer.removeLabel(this.hoverLabel); this.hoverLabel = null; } });
      this.viewer.setClickable({}, true, atom => { if (this.mapped && this.onSelect) this.onSelect(atom.resi); });
      this.paint();
      // Keep the whole protein framed while the page settles, until the visitor moves it themselves
      // or selects a part. Fitting only once left the model off-centre when the layout shifted.
      this.userMoved = false;
      el.addEventListener("pointerdown", () => { this.userMoved = true; });
      const fit = () => {
        this.viewer.resize();
        if (el.clientWidth > 0 && !this.userMoved && !this.part) this.frameAll(0);
        this.viewer.render();
      };
      requestAnimationFrame(fit);
      new ResizeObserver(fit).observe(el);
      const reset = this.host.querySelector(".molreset");
      if (reset) {
        reset.hidden = false;
        reset.addEventListener("click", e => { e.stopPropagation(); this.resetView(); });
      }
    }
    // frame the confidently predicted core, so long disordered tails do not shrink the protein
    frameAll(ms = 500) {
      if (!this.core) {
        const good = [];
        for (const [r] of this.ca) if ((this.plddt[r] ?? 0) >= 70) good.push(r);
        this.core = good.length >= this.ca.size * 0.4 ? { resi: good } : {};
      }
      clearTimeout(this._zt);
      this.viewer.zoomTo(this.core, ms);
      if (ms) this._zt = setTimeout(() => this.viewer.zoom(1.1, 250), ms + 20);
      else this.viewer.zoom(1.1);
    }
    resetView() {
      this.userMoved = false;
      if (this.part) this.framePart(this.part); else this.frameAll(500);
    }
    framePart(part) {
      const sel = { resi: `${part[0]}-${part[1]}` };
      this.viewer.zoomTo(sel, 600);
      // pull back a little so the rest of the protein stays in view around the highlighted part
      clearTimeout(this._zt);
      this._zt = setTimeout(() => { this.viewer.zoom(0.6, 350); }, 620);
    }
    setMetric(m) { this.metric = m; this.paint(); }
    colorFor() {
      const none = cssVar("--rule");
      const m = this.mapped ? this.metric : "plddt";
      if (m === "plddt") return resi => PLDDT(this.plddt[resi] ?? 0);
      if (m === "topo") {
        const c = { o: cssVar("--bm"), i: cssVar("--ink-3"), t: cssVar("--lipid-head"), m: cssVar("--lipid-head") };
        return resi => c[this.wave.topo.st[resi]] || none;
      }
      const r = ramp(METRICS[m].stops()), sp = this.spatial[m];
      return resi => sp.has(resi) ? r(sp.get(resi)) : none;
    }
    paint() {
      if (!this.viewer) return;
      const f = this.colorFor();
      // with a part selected, it keeps its colours and the rest of the protein turns pale grey
      // with a part selected, the part is drawn in one bold highlight colour (the same as its band on
      // the map) and the rest of the protein turns pale grey, so it is obvious where the part is
      const part = this.part, dim = cssVar("--dim3d") || "#d3d9dc", hl = cssVar("--hl") || "#e8a200";
      const inPart = r => !part || (r >= part[0] && r <= part[1]);
      this.viewer.setStyle({}, { cartoon: { colorfunc: atom => part ? dim : f(atom.resi), thickness: 0.4 } });
      if (part) this.viewer.setStyle({ resi: `${part[0]}-${part[1]}` }, { cartoon: { color: hl, thickness: 1.0, width: 1.4 } });
      this.viewer.removeAllShapes();
      if (this.mapped && this.wave.showVariants) {
        const dense = this.wave.byPos.size > 400;
        for (const [p, arr] of this.wave.byPos) {
          if (!inPart(p)) continue;
          const c = this.ca.get(p);
          if (!c) continue;
          this.viewer.addSphere({ center: { x: c[0], y: c[1], z: c[2] }, radius: (dense && !part ? 0.45 : 0.7) + Math.min(arr.length, 4) * 0.1, color: cssVar(`--${TV.EV_KEYS[arr[0].ev]}`) });
        }
      }
      if (this.sel && this.ca.get(this.sel)) {
        this.viewer.setStyle({ resi: this.sel }, part
          ? { cartoon: { color: hl, thickness: 1.0, width: 1.4 }, stick: { colorscheme: "default", radius: 0.3 } }
          : { cartoon: { colorfunc: atom => f(atom.resi), thickness: 0.4 }, stick: { colorscheme: "default", radius: 0.25 } });
        const c = this.ca.get(this.sel);
        this.viewer.addSphere({ center: { x: c[0], y: c[1], z: c[2] }, radius: 2.8, color: cssVar("--ink"), opacity: 0.3 });
      }
      this.viewer.render();
    }
    select(p, part = null, focus = true) {
      this.sel = p;
      this.part = p ? part : null;
      this.paint();
      if (!focus || !this.viewer) return;
      if (this.part) this.framePart(this.part);
      else if (p && this.ca.get(p)) this.viewer.center({ resi: p }, 500);
      else { this.userMoved = false; this.frameAll(500); }
    }
    regional(p) { return this.spatial && this.spatial[this.metric] ? this.spatial[this.metric].get(p) : undefined; }
  }

  function legendRamp(metric) {
    if (metric === "plddt") {
      return `<span class="scale">AlphaFold confidence <span class="ramp" style="background:linear-gradient(90deg,#ff7d45 0 25%,#ffdb13 25% 50%,#65cbf3 50% 75%,#0053d6 75%)"></span> very low to very high</span>`;
    }
    if (metric === "topo") {
      return `<span class="scale"><span class="ev" style="--c:var(--bm)">Outside</span> <span class="ev" style="--c:var(--lipid-head)">Membrane</span> <span class="ev" style="--c:var(--ink-3)">Cytoplasm</span></span>`;
    }
    const m = METRICS[metric], st = m.stops();
    const lo = st[0][0], hi = st[st.length - 1][0];
    const g = st.map(([x, c]) => `${c} ${((x - lo) / (hi - lo) * 100).toFixed(0)}%`).join(",");
    return `<span class="scale">${esc(m.ends[0])} <span class="ramp" style="background:linear-gradient(90deg,${g})"></span> ${esc(m.ends[1])}</span><span class="scale"><span class="ramp" style="width:26px;background:var(--rule)"></span> too few variants nearby</span>`;
  }

  return { Wave, Mol, METRICS, topology, hotspots, legendRamp, binName, spaceProfile };
})();
