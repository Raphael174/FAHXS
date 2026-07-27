"""
TransientDashboard — self-contained dynamic (time-resolved) HTML dashboard.

Companion to the static `HXDashboard` (data_plotting.py). Consumes the
`time_series` dict produced by `transient_solver._build_time_series` and emits a
single self-contained .html file (embedded JSON + inline Canvas rendering, no
external assets, no server) that the user opens in a browser. A time scrubber +
play control animates the axial fields; the wall-temperature x–t heatmap gives
the headline thermal-shock view.

Design follows the dataviz skill: categorical hues assigned by role in fixed
order, a single blue sequential ramp for the heatmap, theme-aware light/dark, a
legend always present with direct labels. No dual axes.

@ author : Raphaël Aubry  (transient extension)
"""
import json
from pathlib import Path

import numpy as np


def _jsonable(obj):
    """Recursively convert numpy arrays/scalars to plain Python for json.dumps."""
    if isinstance(obj, dict):
        return {k: _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return [_jsonable(v) for v in obj.tolist()]
    if isinstance(obj, (np.floating, np.integer)):
        return float(obj)
    return obj


class TransientDashboard:
    def __init__(self, time_series, meta=None):
        self.ts = time_series
        self.meta = meta or {}

    def _payload(self):
        ts = self.ts
        return dict(
            t=_jsonable(ts["t"]),
            x=_jsonable(ts["x"]),
            fields={k: _jsonable(v) for k, v in ts["fields"].items()},
            scalars={k: _jsonable(v) for k, v in ts["scalars"].items()},
            meta=self.meta,
        )

    def to_html(self, path=None):
        if path is None:
            path = Path(__file__).resolve().parent.parent / "transient_dashboard.html"
        path = Path(path)
        data_json = json.dumps(self._payload())
        html = _HTML_TEMPLATE.replace("/*__DATA__*/", data_json)
        path.write_text(html, encoding="utf-8")
        return str(path)


_HTML_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Combustor-HX — Transient Dashboard</title>
<style>
  :root{
    --plane:#f9f9f7; --surface:#fcfcfb; --ink:#0b0b0b; --ink2:#52514e; --muted:#898781;
    --grid:#e1e0d9; --axis:#c3c2b7; --border:rgba(11,11,11,0.10);
    --s-gas:#e34948; --s-wg:#eb6834; --s-wc:#1baf7a; --s-c:#2a78d6; --s-ghost:#898781;
    --warn:#fab219; --good:#0ca30c;
    --seq0:#cde2fb; --seq1:#86b6ef; --seq2:#3987e5; --seq3:#1c5cab; --seq4:#0d366b;
  }
  @media (prefers-color-scheme: dark){
    :root{
      --plane:#0d0d0d; --surface:#1a1a19; --ink:#fff; --ink2:#c3c2b7; --muted:#898781;
      --grid:#2c2c2a; --axis:#383835; --border:rgba(255,255,255,0.10);
      --s-gas:#e66767; --s-wg:#d95926; --s-wc:#199e70; --s-c:#3987e5; --s-ghost:#898781;
    }
  }
  :root[data-theme=light]{ --plane:#f9f9f7; --surface:#fcfcfb; --ink:#0b0b0b; --ink2:#52514e;
    --grid:#e1e0d9; --axis:#c3c2b7; --border:rgba(11,11,11,0.10);
    --s-gas:#e34948; --s-wg:#eb6834; --s-wc:#1baf7a; --s-c:#2a78d6; }
  :root[data-theme=dark]{ --plane:#0d0d0d; --surface:#1a1a19; --ink:#fff; --ink2:#c3c2b7;
    --grid:#2c2c2a; --axis:#383835; --border:rgba(255,255,255,0.10);
    --s-gas:#e66767; --s-wg:#d95926; --s-wc:#199e70; --s-c:#3987e5; }
  *{box-sizing:border-box}
  body{margin:0;background:var(--plane);color:var(--ink);
    font-family:system-ui,-apple-system,"Segoe UI",sans-serif;font-size:14px}
  header{padding:16px 22px;border-bottom:1px solid var(--border);display:flex;
    align-items:baseline;gap:16px;flex-wrap:wrap}
  h1{font-size:17px;margin:0;font-weight:650}
  .sub{color:var(--ink2);font-size:13px}
  .theme-btn{margin-left:auto;background:var(--surface);border:1px solid var(--border);
    color:var(--ink2);border-radius:8px;padding:6px 12px;cursor:pointer;font-size:13px}
  .controls{display:flex;align-items:center;gap:14px;padding:12px 22px;flex-wrap:wrap;
    border-bottom:1px solid var(--border)}
  .play{background:var(--s-c);color:#fff;border:none;border-radius:8px;padding:8px 16px;
    cursor:pointer;font-size:14px;font-weight:600;min-width:80px}
  input[type=range]{flex:1;min-width:200px;accent-color:var(--s-c)}
  .tnow{font-variant-numeric:tabular-nums;color:var(--ink);min-width:110px;font-weight:600}
  .tiles{display:flex;gap:10px;flex-wrap:wrap;padding:14px 22px}
  .tile{background:var(--surface);border:1px solid var(--border);border-radius:10px;
    padding:10px 14px;min-width:120px}
  .tile .lbl{color:var(--muted);font-size:11px;text-transform:uppercase;letter-spacing:.04em}
  .tile .val{font-size:20px;font-weight:650;font-variant-numeric:tabular-nums;margin-top:2px}
  .grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(420px,1fr));gap:16px;
    padding:0 22px 28px}
  .panel{background:var(--surface);border:1px solid var(--border);border-radius:12px;padding:14px}
  .panel h2{font-size:13px;margin:0 0 4px;font-weight:600}
  .panel .cap{color:var(--muted);font-size:11px;margin:0 0 10px}
  .legend{display:flex;gap:14px;flex-wrap:wrap;margin:8px 0 0;font-size:12px;color:var(--ink2)}
  .legend span{display:inline-flex;align-items:center;gap:6px}
  .swatch{width:11px;height:11px;border-radius:3px;display:inline-block}
  canvas{width:100%;height:auto;display:block;border-radius:6px;overflow-x:auto}
  .note{color:var(--warn);font-size:11px;margin-top:6px}
</style>
</head>
<body>
<header>
  <h1>Combustor-HX — Transient Dashboard</h1>
  <span class="sub" id="metaLine"></span>
  <button class="theme-btn" id="themeBtn">Toggle theme</button>
</header>

<div class="controls">
  <button class="play" id="playBtn">▶ Play</button>
  <input type="range" id="tslider" min="0" value="0" step="1">
  <span class="tnow" id="tnow">t = 0.00 s</span>
</div>

<div class="tiles" id="tiles"></div>

<div class="grid">
  <div class="panel">
    <h2>Axial temperature profiles</h2>
    <p class="cap">Along the coil, at the selected time. Ghost = steady solution.</p>
    <canvas id="cProfiles"></canvas>
    <div class="legend">
      <span><i class="swatch" style="background:var(--s-gas)"></i>T_gas</span>
      <span><i class="swatch" style="background:var(--s-wg)"></i>T_wall,hot</span>
      <span><i class="swatch" style="background:var(--s-wc)"></i>T_wall,cold</span>
      <span><i class="swatch" style="background:var(--s-c)"></i>T_He</span>
      <span><i class="swatch" style="background:var(--s-ghost)"></i>steady ref</span>
    </div>
  </div>

  <div class="panel">
    <h2>Wall temperature — space × time</h2>
    <p class="cap">Mean wall T̄(x,t). Vertical cursor marks the selected time.</p>
    <canvas id="cHeat"></canvas>
    <div class="legend" id="heatScale"></div>
  </div>

  <div class="panel">
    <h2>Outlet temperature histories</h2>
    <p class="cap">He outlet flagged unreliable (shaded) while residence &gt; τ_wall (early ramp).</p>
    <canvas id="cOutlet"></canvas>
    <div class="legend">
      <span><i class="swatch" style="background:var(--s-c)"></i>He outlet</span>
      <span><i class="swatch" style="background:var(--s-gas)"></i>gas outlet</span>
    </div>
  </div>

  <div class="panel">
    <h2>Heat duty &amp; peak wall ΔT</h2>
    <p class="cap">Total transferred power and the max hot-to-cold wall gradient.</p>
    <canvas id="cPower"></canvas>
    <div class="legend">
      <span><i class="swatch" style="background:var(--s-wc)"></i>Q transferred [kW]</span>
      <span><i class="swatch" style="background:var(--s-wg)"></i>peak wall ΔT [K]</span>
    </div>
  </div>

  <div class="panel">
    <h2>Boundary schedules</h2>
    <p class="cap">Prescribed mass-flow ramps driving the transient.</p>
    <canvas id="cBC"></canvas>
    <div class="legend">
      <span><i class="swatch" style="background:var(--s-c)"></i>He ṁ [g/s]</span>
      <span><i class="swatch" style="background:var(--s-gas)"></i>gas ṁ [g/s]</span>
    </div>
  </div>
</div>

<script>
const DATA = /*__DATA__*/;
const T = DATA.t, X = DATA.x, F = DATA.fields, S = DATA.scalars;
const NT = T.length, NX = X.length;
let idx = 0, playing = false, raf = null;

const css = k => getComputedStyle(document.documentElement).getPropertyValue(k).trim();
function dpr(){ return Math.min(window.devicePixelRatio||1, 2); }
function setupCanvas(cv, hCss){
  const w = cv.clientWidth || 400, h = hCss;
  const r = dpr(); cv.width = w*r; cv.height = h*r; cv.style.height = h+'px';
  const ctx = cv.getContext('2d'); ctx.setTransform(r,0,0,r,0,0);
  return {ctx, w, h};
}
const PAD = {l:48, r:12, t:10, b:26};
function scaler(dmin,dmax,rmin,rmax){ const s=(rmax-rmin)/((dmax-dmin)||1); return v=>rmin+(v-dmin)*s; }
function niceExtent(arrs){ let lo=Infinity,hi=-Infinity; for(const a of arrs) for(const v of a){ if(v<lo)lo=v; if(v>hi)hi=v;} if(lo===hi){lo-=1;hi+=1;} const pad=(hi-lo)*0.06; return [lo-pad,hi+pad]; }

function axes(ctx,w,h,xdom,ydom,xlab,ylab){
  const sx=scaler(xdom[0],xdom[1],PAD.l,w-PAD.r), sy=scaler(ydom[0],ydom[1],h-PAD.b,PAD.t);
  ctx.strokeStyle=css('--grid'); ctx.fillStyle=css('--muted'); ctx.lineWidth=1;
  ctx.font='11px system-ui'; ctx.textAlign='right'; ctx.textBaseline='middle';
  for(let i=0;i<=4;i++){ const yv=ydom[0]+(ydom[1]-ydom[0])*i/4, yy=sy(yv);
    ctx.beginPath();ctx.moveTo(PAD.l,yy);ctx.lineTo(w-PAD.r,yy);ctx.stroke();
    ctx.fillText(yv.toFixed(yv>=100?0:1),PAD.l-5,yy); }
  ctx.textAlign='center'; ctx.textBaseline='top';
  for(let i=0;i<=4;i++){ const xv=xdom[0]+(xdom[1]-xdom[0])*i/4, xx=sx(xv);
    ctx.fillText(xv.toFixed(xv>=100?0:(xv>=10?1:2)),xx,h-PAD.b+5); }
  ctx.fillStyle=css('--muted'); ctx.textAlign='left';
  ctx.fillText(ylab,PAD.l-42,2); ctx.textAlign='right'; ctx.fillText(xlab,w-PAD.r,h-12);
  return {sx,sy};
}
function line(ctx,xs,ys,sx,sy,color,wid){ ctx.strokeStyle=color;ctx.lineWidth=wid||2;
  ctx.lineJoin='round';ctx.beginPath();
  for(let i=0;i<xs.length;i++){ const px=sx(xs[i]),py=sy(ys[i]); i?ctx.lineTo(px,py):ctx.moveTo(px,py);} ctx.stroke(); }

// ---- profiles panel (animated) ----
const cvP=document.getElementById('cProfiles');
let ydomP=null;
function drawProfiles(){
  const {ctx,w,h}=setupCanvas(cvP,240);
  ctx.clearRect(0,0,w,h);
  if(!ydomP) ydomP=niceExtent([F.T_g.flat(),F.T_c.flat(),F.T_wg.flat(),F.T_wc.flat()]);
  const {sx,sy}=axes(ctx,w,h,[X[0],X[NX-1]],ydomP,'coil arc-length x [m]','T [K]');
  if(DATA.meta.Tbar_steady) line(ctx,X,DATA.meta.Tbar_steady,sx,sy,css('--s-ghost'),1.5);
  line(ctx,X,F.T_g[idx],sx,sy,css('--s-gas'),2);
  line(ctx,X,F.T_wg[idx],sx,sy,css('--s-wg'),2);
  line(ctx,X,F.T_wc[idx],sx,sy,css('--s-wc'),2);
  line(ctx,X,F.T_c[idx],sx,sy,css('--s-c'),2);
}

// ---- wall heatmap (static + time cursor) ----
const cvH=document.getElementById('cHeat');
let heatDom=null;
function seqColor(f){ // f in [0,1] -> blue ramp
  const stops=[css('--seq0'),css('--seq1'),css('--seq2'),css('--seq3'),css('--seq4')];
  const p=Math.max(0,Math.min(1,f))*(stops.length-1); const i=Math.floor(p),t=p-i;
  const a=hex(stops[i]),b=hex(stops[Math.min(i+1,stops.length-1)]);
  return `rgb(${Math.round(a[0]+(b[0]-a[0])*t)},${Math.round(a[1]+(b[1]-a[1])*t)},${Math.round(a[2]+(b[2]-a[2])*t)})`;
}
function hex(c){ c=c.trim(); if(c[0]==='#'){const n=parseInt(c.slice(1),16);return [(n>>16)&255,(n>>8)&255,n&255];}
  const m=c.match(/\d+/g); return m?[+m[0],+m[1],+m[2]]:[0,0,0]; }
function drawHeat(){
  const {ctx,w,h}=setupCanvas(cvH,240);
  ctx.clearRect(0,0,w,h);
  if(!heatDom){ let lo=Infinity,hi=-Infinity; for(const row of F.Tbar) for(const v of row){if(v<lo)lo=v;if(v>hi)hi=v;} heatDom=[lo,hi]; }
  const x0=PAD.l,y0=PAD.t,pw=w-PAD.l-PAD.r,ph=h-PAD.t-PAD.b;
  const cw=pw/NX, chh=ph/NT;
  for(let j=0;j<NT;j++) for(let i=0;i<NX;i++){
    const f=(F.Tbar[j][i]-heatDom[0])/((heatDom[1]-heatDom[0])||1);
    ctx.fillStyle=seqColor(f);
    ctx.fillRect(x0+i*cw, y0+ph-(j+1)*chh, Math.ceil(cw)+0.5, Math.ceil(chh)+0.5);
  }
  // time cursor
  const yy=y0+ph-(idx+0.5)*chh;
  ctx.strokeStyle=css('--ink'); ctx.lineWidth=1.5; ctx.beginPath();ctx.moveTo(x0,yy);ctx.lineTo(x0+pw,yy);ctx.stroke();
  ctx.fillStyle=css('--muted');ctx.font='11px system-ui';ctx.textAlign='right';ctx.textBaseline='middle';
  ctx.fillText('t',x0-5,y0+ph/2); ctx.textAlign='center';ctx.textBaseline='top';
  ctx.fillText('x [m]',x0+pw/2,h-PAD.b+6);
  // scale legend
  const sc=document.getElementById('heatScale');
  sc.innerHTML=`<span>${heatDom[0].toFixed(0)} K</span>
    <span style="flex:1;height:10px;border-radius:3px;min-width:80px;background:linear-gradient(90deg,${css('--seq0')},${css('--seq2')},${css('--seq4')})"></span>
    <span>${heatDom[1].toFixed(0)} K</span>`;
}

// ---- generic history panel with time cursor + optional shading ----
function drawHistory(cv, seriesList, shadeMask){
  const {ctx,w,h}=setupCanvas(cv,220);
  ctx.clearRect(0,0,w,h);
  const ydom=niceExtent(seriesList.map(s=>s.y));
  const {sx,sy}=axes(ctx,w,h,[T[0],T[NT-1]],ydom,'t [s]',seriesList[0].unit||'');
  if(shadeMask){ ctx.fillStyle='rgba(250,178,25,0.14)';
    for(let j=0;j<NT;j++){ if(shadeMask[j]){ const x1=sx(T[j]),x2=sx(T[Math.min(j+1,NT-1)]);
      ctx.fillRect(x1,PAD.t,(x2-x1)||1,h-PAD.t-PAD.b);} } }
  for(const s of seriesList) line(ctx,T,s.y,sx,sy,s.color,2);
  const cx=sx(T[idx]); ctx.strokeStyle=css('--axis');ctx.lineWidth=1;
  ctx.setLineDash([3,3]);ctx.beginPath();ctx.moveTo(cx,PAD.t);ctx.lineTo(cx,h-PAD.b);ctx.stroke();ctx.setLineDash([]);
}
function drawOutlet(){ drawHistory(document.getElementById('cOutlet'),
  [{y:S.T_c_out,color:css('--s-c'),unit:'T [K]'},{y:S.T_g_out,color:css('--s-gas')}],
  S.He_outlet_reliable.map(v=>v===0)); }
function drawPower(){ drawHistory(document.getElementById('cPower'),
  [{y:S.Q_cold_kW,color:css('--s-wc'),unit:'kW / K'},{y:S.dT_wall_max,color:css('--s-wg')}]); }
function drawBC(){ drawHistory(document.getElementById('cBC'),
  [{y:S.mdot_c.map(v=>v*1000),color:css('--s-c'),unit:'ṁ [g/s]'},{y:S.mdot_g.map(v=>v*1000),color:css('--s-gas')}]); }

// ---- tiles ----
function tiles(){
  const el=document.getElementById('tiles');
  const items=[
    ['He outlet', S.T_c_out[idx].toFixed(0)+' K'],
    ['Gas outlet', S.T_g_out[idx].toFixed(0)+' K'],
    ['Q transferred', S.Q_cold_kW[idx].toFixed(0)+' kW'],
    ['Peak wall ΔT', S.dT_wall_max[idx].toFixed(1)+' K'],
    ['He ṁ', (S.mdot_c[idx]*1000).toFixed(1)+' g/s'],
    ['Wall T̄ max', Math.max(...F.Tbar[idx]).toFixed(0)+' K'],
  ];
  el.innerHTML=items.map(([l,v])=>`<div class="tile"><div class="lbl">${l}</div><div class="val">${v}</div></div>`).join('');
}

function redrawAll(){ drawProfiles(); drawHeat(); drawOutlet(); drawPower(); drawBC(); tiles();
  document.getElementById('tnow').textContent='t = '+T[idx].toFixed(2)+' s'; }

// ---- controls ----
const slider=document.getElementById('tslider'); slider.max=NT-1;
slider.addEventListener('input',()=>{ idx=+slider.value; redrawAll(); });
const playBtn=document.getElementById('playBtn');
function tick(){ if(!playing)return; idx=(idx+1)%NT; slider.value=idx; redrawAll();
  raf=setTimeout(()=>requestAnimationFrame(tick),60); }
playBtn.addEventListener('click',()=>{ playing=!playing;
  playBtn.textContent=playing?'❚❚ Pause':'▶ Play'; if(playing) tick(); else clearTimeout(raf); });
document.getElementById('themeBtn').addEventListener('click',()=>{
  const r=document.documentElement; const cur=r.getAttribute('data-theme');
  const next=cur==='dark'?'light':(cur==='light'?'dark':(matchMedia('(prefers-color-scheme: dark)').matches?'light':'dark'));
  r.setAttribute('data-theme',next); heatDom=null; ydomP=null; redrawAll(); });
window.addEventListener('resize',()=>{ heatDom=null; redrawAll(); });

document.getElementById('metaLine').textContent =
  (DATA.meta.config||'')+(DATA.meta.material?('  ·  '+DATA.meta.material):'')+
  '   ·   '+NT+' snapshots · '+NX+' nodes · t_end '+T[NT-1].toFixed(0)+' s';
redrawAll();
</script>
</body>
</html>
"""
