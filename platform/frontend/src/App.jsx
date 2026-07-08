import React, { useEffect, useState } from "react";

const api = (path, opts) =>
  fetch(path, { headers: { "Content-Type": "application/json" }, ...opts }).then((r) => {
    if (!r.ok) throw new Error(`${r.status}`);
    return r.json();
  });

function EnergyChart({ edRec, test }) {
  // bar chart of energy_total: ED (edRec, exact) vs QMC (test, with error bar)
  const items = [
    { label: edRec.code, o: edRec.observables.energy_total, color: "var(--ref)" },
    { label: test.code, o: test.observables.energy_total, color: "var(--accent)" },
  ];
  const vals = items.map((i) => i.o.value);
  const lo = Math.min(...vals) * 1.001, hi = Math.max(...vals) * 0.999;
  const W = 460, H = 180, pad = 44;
  const x = (i) => pad + 40 + i * 150;
  const y = (v) => H - pad - ((v - lo) / (hi - lo || 1)) * (H - 2 * pad);
  return (
    <svg viewBox={`0 0 ${W} ${H}`} width="100%">
      <line x1={pad} y1={H - pad} x2={W - 10} y2={H - pad} stroke="var(--line)" />
      {items.map((it, i) => {
        const v = it.o.value, e = it.o.error || 0;
        return (
          <g key={i}>
            <rect x={x(i) - 26} y={y(v)} width="52" height={H - pad - y(v)} fill={it.color} opacity="0.85" rx="3" />
            {e > 0 && (
              <line x1={x(i)} y1={y(v - e)} x2={x(i)} y2={y(v + e)} stroke="var(--text)" strokeWidth="2" />
            )}
            <text x={x(i)} y={H - pad + 16} textAnchor="middle">{it.label}</text>
            <text x={x(i)} y={y(v) - 6} textAnchor="middle" fill="var(--text)">{v.toFixed(4)}</text>
          </g>
        );
      })}
      <text x={pad - 6} y={y(lo)} textAnchor="end">{lo.toFixed(3)}</text>
      <text x={pad - 6} y={y(hi)} textAnchor="end">{hi.toFixed(3)}</text>
    </svg>
  );
}

function SweepChart({ res }) {
  // line: ED curve; points with error bars: QMC
  const pts = res.points;
  const xs = pts.map((p) => p.x);
  const ys = pts.flatMap((p) => [p.ed, p.qmc].filter((y) => y != null));
  const xlo = Math.min(...xs), xhi = Math.max(...xs);
  const ylo = Math.min(...ys), yhi = Math.max(...ys);
  const W = 520, H = 260, pad = 48;
  const X = (x) => pad + ((x - xlo) / (xhi - xlo || 1)) * (W - pad - 16);
  const Y = (y) => H - pad - ((y - ylo) / (yhi - ylo || 1)) * (H - 2 * pad);
  const edPath = pts.map((p, i) => `${i ? "L" : "M"}${X(p.x)},${Y(p.ed)}`).join(" ");
  return (
    <svg viewBox={`0 0 ${W} ${H}`} width="100%">
      <line x1={pad} y1={H - pad} x2={W - 10} y2={H - pad} stroke="var(--line)" />
      <line x1={pad} y1={pad} x2={pad} y2={H - pad} stroke="var(--line)" />
      <path d={edPath} fill="none" stroke="var(--ref)" strokeWidth="2" />
      {pts.map((p, i) => (
        <g key={i}>
          <circle cx={X(p.x)} cy={Y(p.ed)} r="3" fill="var(--ref)" />
          {p.qmc != null && (
            <>
              {p.qmc_err > 0 && (
                <line x1={X(p.x)} y1={Y(p.qmc - p.qmc_err)} x2={X(p.x)} y2={Y(p.qmc + p.qmc_err)} stroke="var(--accent)" strokeWidth="2" />
              )}
              <circle cx={X(p.x)} cy={Y(p.qmc)} r="3.5" fill="var(--accent)" />
            </>
          )}
          <text x={X(p.x)} y={H - pad + 15} textAnchor="middle">{p.x}</text>
        </g>
      ))}
      <text x={pad - 6} y={Y(ylo)} textAnchor="end">{ylo.toFixed(2)}</text>
      <text x={pad - 6} y={Y(yhi)} textAnchor="end">{yhi.toFixed(2)}</text>
      <text x={(W) / 2} y={H - 6} textAnchor="middle">{res.var}</text>
    </svg>
  );
}

function Heatmap({ m, title }) {
  // diverging blue(neg)–white(0)–coral(pos) heatmap of an n×n matrix
  const n = m.length;
  const amax = Math.max(...m.flat().map((x) => Math.abs(x)), 1e-9);
  const cell = Math.max(10, Math.min(22, Math.floor(160 / n)));
  const col = (x) => {
    const t = x / amax;
    return t >= 0 ? `rgba(216,90,48,${t})` : `rgba(55,138,221,${-t})`;
  };
  return (
    <div style={{ textAlign: "center" }}>
      <div style={{ fontSize: 11, color: "var(--muted)", marginBottom: 4 }}>{title}</div>
      <svg width={n * cell} height={n * cell} style={{ border: "1px solid var(--line)" }}>
        {m.map((row, i) => row.map((x, j) => (
          <rect key={i + "_" + j} x={j * cell} y={i * cell} width={cell} height={cell}
            fill={col(x)} stroke="var(--line)" strokeWidth="0.3" />
        )))}
      </svg>
    </div>
  );
}

export default function App() {
  const [runs, setRuns] = useState([]);
  const [ed, setEd] = useState({
    model: "hubbard", lx: 4, ly: 4, nup: 1, ndn: 1, t: 1, U: 3,
    t1: -1, t2: -1, t3: -1, t4: -1, uxy: 0, v: 0,
    nwalkers: 200, nequil: 120, nmeas: 250, bp: 0,
  });
  const [sweep, setSweep] = useState({ var: "U", values: "0,0.5,1,1.5,2", run_qmc: true });
  const [sweepRes, setSweepRes] = useState(null);
  const [tab, setTab] = useState("compare");   // compare | sweep | runs | speed | corr
  const [benchRes, setBenchRes] = useState(null);
  const [corrRes, setCorrRes] = useState(null);
  // d-wave pairing-vertex benchmark (its own params; default = NON-DEGENERATE 4x2 3+3)
  const [pair, setPair] = useState({
    lx: 4, ly: 2, nup: 3, ndn: 3, U: 4, t1: 0.3, tam: 0.2, etas: "0.0,0.5",
    dt: 0.02, bp: 28, nblocks: 30, npop: 12, nw: 60,
  });
  const [pairRes, setPairRes] = useState(null);
  const [ref, setRef] = useState(null);     // ED record
  const [sel, setSel] = useState(null);     // selected run summary
  const [test, setTest] = useState(null);   // selected run / live QMC record
  const [cmp, setCmp] = useState(null);     // compare result
  const [busy, setBusy] = useState(false);
  const [err, setErr] = useState(null);

  const loadRuns = () => api("/api/runs").then(setRuns).catch((e) => setErr(String(e)));
  useEffect(() => { loadRuns(); }, []);

  const computeED = async () => {
    setBusy(true); setErr(null);
    try { setRef(await api("/api/ed", { method: "POST", body: JSON.stringify(ed) })); setCmp(null); }
    catch (e) { setErr("ED failed: " + e); } finally { setBusy(false); }
  };

  const runQMC = async () => {
    setBusy(true); setErr(null);
    try {
      const rec = await api("/api/qmc", { method: "POST", body: JSON.stringify(ed) });
      setTest(rec); setSel({ code: rec.code, run_id: "live" }); setCmp(null); setTab("compare");
    } catch (e) { setErr("QMC run failed: " + e); } finally { setBusy(false); }
  };

  const runSweep = async () => {
    setBusy(true); setErr(null); setSweepRes(null);
    try {
      const values = sweep.values.split(",").map((s) => Number(s.trim())).filter((x) => !isNaN(x));
      const body = { base: ed, var: sweep.var, values, run_qmc: sweep.run_qmc };
      setSweepRes(await api("/api/sweep", { method: "POST", body: JSON.stringify(body) }));
    } catch (e) { setErr("sweep failed: " + e); } finally { setBusy(false); }
  };

  const pick = async (r) => {
    setSel(r); setCmp(null);
    setTest(await api(`/api/runs/${r.code}/${r.run_id}`));
    setTab("compare");
  };

  const runBench = async () => {
    setBusy(true); setErr(null); setBenchRes(null);
    try {
      setBenchRes(await api("/api/benchmark", { method: "POST", body: JSON.stringify({
        lx: ed.lx, ly: ed.ly, nup: ed.nup, ndn: ed.ndn, uxx: ed.U }) }));
    } catch (e) { setErr("benchmark failed: " + e); } finally { setBusy(false); }
  };

  const runCorr = async () => {
    setBusy(true); setErr(null); setCorrRes(null);
    try {
      setCorrRes(await api("/api/correlations", { method: "POST", body: JSON.stringify({
        lx: ed.lx, ly: ed.ly, nup: ed.nup, ndn: ed.ndn, U: ed.U, uxy: ed.uxy, v: ed.v,
        t1: ed.t1, t2: ed.t2, t3: ed.t3, t4: ed.t4, bp: ed.bp > 0 ? ed.bp : 12 }) }));
    } catch (e) { setErr("correlations failed: " + e); } finally { setBusy(false); }
  };

  const runPairing = async () => {
    setBusy(true); setErr(null); setPairRes(null);
    try { setPairRes(await api("/api/pairing", { method: "POST", body: JSON.stringify(pair) })); }
    catch (e) { setErr("pairing benchmark failed: " + e); } finally { setBusy(false); }
  };
  const pnum = (k, step = 1) => (
    <div><label>{k}</label>
      <input type="number" value={pair[k]} step={step}
        onChange={(e) => setPair({ ...pair, [k]: Number(e.target.value) })} /></div>
  );

  const doCompare = async () => {
    if (!ref || !test) return;
    setBusy(true); setErr(null);
    try { setCmp(await api("/api/compare", { method: "POST", body: JSON.stringify({ ref, test, ztol: 3 }) })); }
    catch (e) { setErr("compare failed: " + e); } finally { setBusy(false); }
  };

  const num = (k) => (
    <div><label>{k}</label>
      <input type="number" value={ed[k]} step={["t","U","t1","t2","t3","t4","uxy","v"].includes(k) ? 0.1 : 1}
        onChange={(e) => setEd({ ...ed, [k]: Number(e.target.value) })} /></div>
  );

  const tabs = [
    ["compare", "Compare"],
    ["sweep", "Sweep"],
    ["corr", "Correlations"],
    ["pairing", "Pairing vertex"],
    ["speed", "Speed"],
    ["runs", "Runs"],
  ];

  return (
    <div className="app">
      <h1>QMC Validation Platform</h1>
      <p className="sub">Exact diagonalization vs CPQMC — energies, error bars, σ-distance.</p>
      {err && <p className="err">{err}</p>}
      <div className="grid">
        {/* ---- left sidebar: shared parameters ---- */}
        <div className="panel">
          <h2>Model &amp; parameters</h2>
          <label>model</label>
          <select value={ed.model} onChange={(e) => setEd({ ...ed, model: e.target.value })}>
            <option value="hubbard">single-band Hubbard</option>
            <option value="altermagnet">two-orbital altermagnet</option>
          </select>
          <div className="row2">{num("lx")}{num("ly")}</div>
          <div className="row2">{num("nup")}{num("ndn")}</div>
          {ed.model === "altermagnet" ? (
            <>
              <div className="row2">{num("t1")}{num("t2")}</div>
              <div className="row2">{num("t3")}{num("t4")}</div>
              <div className="row2">{num("U")}{num("uxy")}</div>
              <div className="row2">{num("v")}{num("nwalkers")}</div>
            </>
          ) : (
            <>
              <div className="row2">{num("t")}{num("U")}</div>
              <div className="row2">{num("nwalkers")}{num("nmeas")}</div>
            </>
          )}
          <div className="row2">{num("bp")}<div><label>estimator</label>
            <input value={ed.bp > 0 ? "back-prop" : "mixed"} disabled style={{ opacity: .7 }} /></div></div>
          <div style={{ display: "flex", gap: 8 }}>
            <button onClick={computeED} disabled={busy}>{busy ? "…" : "Compute ED"}</button>
            <button className="secondary" onClick={runQMC} disabled={busy}>{busy ? "…" : "Run QMC ▶"}</button>
          </div>
          <div style={{ marginTop: 12, display: "flex", flexDirection: "column", gap: 4 }}>
            <span className="status">
              <span className={"sdot " + (ref ? "ok" : "")} />ED reference
              {ref ? `: ${ref.observables.energy_total.value.toFixed(4)}${ref.ed_basis_dim ? ` (${ref.ed_basis_dim} states)` : ""}` : " — not computed"}
            </span>
            <span className="status">
              <span className={"sdot " + (test ? "ok" : "")} />QMC
              {test ? `: ${test.observables.energy_total.value.toFixed(4)} (${test.estimator || test.code})` : " — not run"}
            </span>
          </div>
        </div>

        {/* ---- right: tabbed pane ---- */}
        <div>
          <div className="tabbar">
            {tabs.map(([id, label]) => (
              <button key={id} className={"tab" + (tab === id ? " active" : "")} onClick={() => setTab(id)}>
                {label}{id === "runs" && runs.length > 0 ? ` (${runs.length})` : ""}
              </button>
            ))}
          </div>

          {tab === "compare" && (
            <div className="panel tabpanel">
              <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center" }}>
                <h2 style={{ margin: 0 }}>ED vs {sel ? sel.code : "QMC"}</h2>
                <button onClick={doCompare} disabled={!ref || !test || busy} style={{ width: "auto", margin: 0, padding: "7px 16px" }}>
                  Compare
                </button>
              </div>
              {!ref || !test ? (
                <p className="muted" style={{ marginTop: 14 }}>
                  Set parameters on the left, then <b>Compute ED</b> and <b>Run QMC</b> (or pick a saved run in the Runs tab). Then Compare.
                </p>
              ) : (
                <>
                  {test.note && <p style={{ color: "var(--fail)", fontSize: 12, marginTop: 12 }}>⚠ {test.note}</p>}
                  <div className="legend" style={{ marginTop: 14 }}>
                    <span><span className="dot" style={{ background: "var(--ref)" }} />{ref.code} (exact)</span>
                    <span><span className="dot" style={{ background: "var(--accent)" }} />{test.code}</span>
                  </div>
                  <EnergyChart edRec={ref} test={test} />
                </>
              )}
              {cmp && (
                <>
                  <table>
                    <thead><tr><th>observable</th><th>ED</th><th>QMC</th><th>Δ</th><th>σ</th><th>z</th><th>verdict</th></tr></thead>
                    <tbody>
                      {cmp.rows.map((r) => (
                        <tr key={r.obs}>
                          <td>{r.obs}</td><td>{r.ref.toFixed(5)}</td><td>{r.test.toFixed(5)}</td>
                          <td>{r.delta.toFixed(5)}</td><td>{r.sigma.toFixed(5)}</td>
                          <td>{isFinite(r.z) ? r.z.toFixed(2) : "∞"}</td>
                          <td><span className={"badge " + r.verdict}>{r.verdict}</span></td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                  <p className="overall">Overall: <span className={"badge " + cmp.overall}>{cmp.overall}</span>
                    <span className="muted"> (worst z = {cmp.worst_z.toFixed(2)}, ztol = {cmp.ztol})</span></p>
                </>
              )}
            </div>
          )}

          {tab === "sweep" && (
            <div className="panel tabpanel">
              <h2>E vs interaction (ED curve + QMC points)</h2>
              <div style={{ display: "flex", gap: 12, alignItems: "flex-end", flexWrap: "wrap" }}>
                <div style={{ width: 110 }}><label>sweep variable</label>
                  <select value={sweep.var} onChange={(e) => setSweep({ ...sweep, var: e.target.value })}>
                    <option value="U">U / uxx</option>
                    <option value="uxy">uxy</option>
                    <option value="v">v</option>
                  </select></div>
                <div style={{ flex: 1, minWidth: 160 }}><label>values (comma-separated)</label>
                  <input value={sweep.values} onChange={(e) => setSweep({ ...sweep, values: e.target.value })} /></div>
                <label style={{ display: "flex", alignItems: "center", gap: 6, margin: "0 0 7px" }}>
                  <input type="checkbox" style={{ width: "auto" }} checked={sweep.run_qmc}
                    onChange={(e) => setSweep({ ...sweep, run_qmc: e.target.checked })} /> run QMC
                </label>
                <button onClick={runSweep} disabled={busy} style={{ width: "auto", padding: "8px 18px" }}>
                  {busy ? "running…" : "Run sweep"}
                </button>
              </div>
              <p className="muted" style={{ fontSize: 12, marginTop: 6 }}>
                Uses the parameters on the left as the base (set bp&gt;0 for unbiased QMC). QMC points are slow; ED-only is instant.
              </p>
              {sweepRes ? (
                <>
                  <div className="legend">
                    <span><span className="dot" style={{ background: "var(--ref)" }} />ED (exact)</span>
                    <span><span className="dot" style={{ background: "var(--accent)" }} />QMC</span>
                  </div>
                  <SweepChart res={sweepRes} />
                </>
              ) : <p className="muted" style={{ marginTop: 14 }}>Pick a variable and values, then Run sweep.</p>}
            </div>
          )}

          {tab === "corr" && (
            <div className="panel tabpanel">
              <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center" }}>
                <h2 style={{ margin: 0 }}>Equal-time GF + spin/charge correlations</h2>
                <button onClick={runCorr} disabled={busy} style={{ width: "auto", margin: 0, padding: "7px 16px" }}>
                  {busy ? "running…" : "Benchmark ED vs CPQMC"}
                </button>
              </div>
              <p className="muted" style={{ fontSize: 12, margin: "6px 0 0" }}>
                Compares the Green's function G↑ᵢⱼ=⟨c†ᵢcⱼ⟩, charge ⟨nᵢnⱼ⟩ and spin ⟨SᶻᵢSᶻⱼ⟩
                element-by-element (two-orbital, uses the left parameters). ~30s.
              </p>
              {corrRes && (
                <>
                  <table>
                    <thead><tr><th>observable</th><th>max |ED−QMC|</th><th>mean</th><th>ED range</th><th>verdict</th></tr></thead>
                    <tbody>
                      {corrRes.summary.map((s) => (
                        <tr key={s.obs}>
                          <td>{s.obs}</td><td>{s.max_dev.toFixed(4)}</td><td>{s.mean_dev.toFixed(4)}</td>
                          <td>[{s.ed_range[0].toFixed(2)}, {s.ed_range[1].toFixed(2)}]</td>
                          <td><span className={"badge " + s.verdict}>{s.verdict}</span></td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                  <p className="muted" style={{ fontSize: 12, margin: "12px 0 4px" }}>Green's function G↑ᵢⱼ ({corrRes.nsites}×{corrRes.nsites}):</p>
                  <div style={{ display: "flex", gap: 18, flexWrap: "wrap" }}>
                    <Heatmap m={corrRes.ed.green_up} title="ED" />
                    <Heatmap m={corrRes.qmc.green_up} title="CPQMC" />
                  </div>
                </>
              )}
              {!corrRes && <p className="muted" style={{ marginTop: 14 }}>Click Benchmark to measure and compare the correlations.</p>}
            </div>
          )}

          {tab === "pairing" && (
            <div className="panel tabpanel">
              <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center" }}>
                <h2 style={{ margin: 0 }}>d-wave pairing vertex vs full-Fock ED</h2>
                <button onClick={runPairing} disabled={busy} style={{ width: "auto", margin: 0, padding: "7px 16px" }}>
                  {busy ? "running…" : "Benchmark vs ED ▶"}
                </button>
              </div>
              <p className="muted" style={{ fontSize: 12, margin: "6px 0 10px" }}>
                Back-propagated AGP/BCS d-wave pairing <b>vertex</b> (and energy) vs the exact
                full-Fock ED, on a <b>non-degenerate</b> cluster — where the constrained-path
                CP-AFQMC reproduces both, unlike the frustrated half-filled (degenerate) point.
                η=0 is the free-electron bra; η&gt;0 the pairing (AGP) bra. ~1–2 min.
              </p>
              <div className="paramgrid" style={{ display: "grid", gridTemplateColumns: "repeat(6, 1fr)", gap: 8 }}>
                {pnum("lx")}{pnum("ly")}{pnum("nup")}{pnum("ndn")}{pnum("U", 0.5)}{pnum("t1", 0.1)}
                {pnum("tam", 0.1)}{pnum("bp")}{pnum("nblocks")}{pnum("nw")}{pnum("npop")}
                <div><label>etas</label>
                  <input value={pair.etas} onChange={(e) => setPair({ ...pair, etas: e.target.value })} /></div>
              </div>
              {pairRes && (
                <>
                  <p className="overall" style={{ marginTop: 14 }}>
                    {pairRes.params.lx}×{pairRes.params.ly}, {pairRes.params.nup}+{pairRes.params.ndn},
                    U={pairRes.params.U}, t1={pairRes.params.t1}, tam={pairRes.params.tam}
                    {" — "}
                    <span className={"badge " + (pairRes.degenerate ? "FAIL" : "PASS")}>
                      {pairRes.degenerate ? "DEGENERATE" : "non-degenerate"}
                    </span>{" "}
                    ⟨sign⟩={pairRes.sign?.toFixed(3)}
                  </p>
                  <p className="muted" style={{ fontSize: 12, margin: "2px 0 8px" }}>
                    full-Fock ED: energy <b>{pairRes.ed.energy?.toFixed(4)}</b>,
                    d-wave vertex <b>{pairRes.ed.S_vertex?.toFixed(3)}</b>
                  </p>
                  <table>
                    <thead><tr>
                      <th>bra</th><th>energy</th><th>z(E)</th><th>energy</th>
                      <th>d-vertex</th><th>frac of ED</th>
                    </tr></thead>
                    <tbody>
                      {pairRes.rows.map((r) => (
                        <tr key={r.eta}>
                          <td>{r.bra}</td>
                          <td>{r.energy.toFixed(4)} ± {r.energy_err.toFixed(4)}</td>
                          <td>{r.z != null ? r.z.toFixed(2) : "—"}</td>
                          <td><span className={"badge " + (r.energy_pass ? "PASS" : "FAIL")}>
                            {r.energy_pass ? "PASS" : "FAIL"}</span></td>
                          <td>{r.vertex.toFixed(3)} ± {r.vertex_err.toFixed(3)}</td>
                          <td>{r.vertex_frac != null ? (100 * r.vertex_frac).toFixed(0) + "%" : "—"}</td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                  <p className="muted" style={{ fontSize: 12, marginTop: 10 }}>
                    Energy PASS = within z&lt;3 of ED. The d-wave-vertex “frac of ED” is the
                    constrained-path estimate as a fraction of the exact vertex; on a
                    non-degenerate cluster the free bra (η=0) already lands near 100%.
                  </p>
                </>
              )}
              {!pairRes && <p className="muted" style={{ marginTop: 14 }}>
                Set parameters (default is the non-degenerate 4×2 3+3 benchmark) and click Benchmark.
              </p>}
            </div>
          )}

          {tab === "speed" && (
            <div className="panel tabpanel">
              <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center" }}>
                <h2 style={{ margin: 0 }}>Speed — Fortran vs Python CPQMC</h2>
                <button onClick={runBench} disabled={busy} style={{ width: "auto", margin: 0, padding: "7px 16px" }}>
                  {busy ? "running…" : "Run benchmark"}
                </button>
              </div>
              <p className="muted" style={{ fontSize: 12, margin: "6px 0 0" }}>
                Times both implementations and reports throughput (walker-updates/s). Fortran is
                compiled + MPI; Python is interpreted numpy. Small runs include MPI startup overhead.
              </p>
              {benchRes && (
                <>
                  <table>
                    <thead><tr><th>implementation</th><th>wall (s)</th><th>walkers</th><th>walker-steps</th><th>throughput /s</th></tr></thead>
                    <tbody>
                      <tr><td>Fortran (code/src)</td><td>{benchRes.fortran.wall_sec}</td><td>{benchRes.fortran.walkers}</td>
                        <td>{benchRes.fortran.walker_steps.toLocaleString()}</td><td>{benchRes.fortran.throughput.toLocaleString()}</td></tr>
                      <tr><td>Python (pyqmc)</td><td>{benchRes.python.wall_sec}</td><td>{benchRes.python.walkers}</td>
                        <td>{benchRes.python.walker_steps.toLocaleString()}</td><td>{benchRes.python.throughput.toLocaleString()}</td></tr>
                    </tbody>
                  </table>
                  <p className="overall">Fortran is <span className="badge PASS">{benchRes.speedup_fortran_over_python}×</span> the Python throughput.</p>
                </>
              )}
              {!benchRes && <p className="muted" style={{ marginTop: 14 }}>Click Run benchmark (uses the left lattice/filling; src-small must be built).</p>}
            </div>
          )}

          {tab === "runs" && (
            <div className="panel tabpanel">
              <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center" }}>
                <h2 style={{ margin: 0 }}>Saved QMC runs</h2>
                <button className="secondary" style={{ width: "auto", margin: 0, padding: "5px 12px" }} onClick={loadRuns}>↻ refresh</button>
              </div>
              <p className="muted" style={{ fontSize: 12, margin: "6px 0 4px" }}>
                Runs under results/&lt;version&gt;/. Click one to load it as the QMC side, then go to Compare.
              </p>
              <div className="runlist" style={{ maxHeight: 420 }}>
                {runs.map((r) => (
                  <div key={r.code + r.run_id} className={"run" + (sel && sel.run_id === r.run_id && sel.code === r.code ? " sel" : "")} onClick={() => pick(r)}>
                    <div><div className="code">{r.code}</div><div className="meta">{r.run_id}</div></div>
                    <div style={{ textAlign: "right" }}>
                      <div>{r.energy_total != null ? r.energy_total.toFixed(4) : "—"}</div>
                      <div className="meta">{r.error != null ? "±" + r.error.toFixed(4) : ""}</div>
                    </div>
                  </div>
                ))}
                {runs.length === 0 && <p className="muted">No runs yet — use Run QMC, or tools/run.sh.</p>}
              </div>
            </div>
          )}
        </div>
      </div>
    </div>
  );
}
