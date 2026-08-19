import base64
STYLE = open('report_style.css').read()
def img(p): return base64.b64encode(open('plots/'+p,'rb').read()).decode()
def img2(p): return base64.b64encode(open(p,'rb').read()).decode()

def block(txt, start):
    out=[]; on=False
    for ln in txt.splitlines():
        if not on and ln.strip().startswith(start): on=True; continue
        if on:
            if not ln.strip(): break
            out.append(ln)
    return out

L1=open('monitoring_ep-1.log').read(); L4=open('monitoring_ep4.log').read()

def pair_table(start, hdr):
    a=block(L1,start); b=block(L4,start)
    m={r.split()[0]: r.split() for r in b}
    body=''
    for r in a:
        c=r.split(); d=m.get(c[0],['','','','',''])
        body+=('<tr><td class="val">'+c[0]+'</td>'
               +''.join(f'<td class="val">{v}</td>' for v in c[1:5])
               +''.join(f'<td class="val">{v}</td>' for v in d[1:5])+'</tr>\n')
    return ('<div class="wrap"><table><thead><tr>'+''.join(f'<th>{h}</th>' for h in hdr)+
            f'</tr></thead><tbody>\n{body}</tbody></table></div>')

HDR=['x','y mean (&minus;1)','y width (&minus;1)','z mean (&minus;1)','z width (&minus;1)',
         'y mean (4)','y width (4)','z mean (4)','z width (4)']
PT  = pair_table('x     mean_y', HDR)   # first block in the log is pT
# eta block: find by section
def sect(txt, name):
    lines=txt.splitlines(); out=[]; on=False
    for i,ln in enumerate(lines):
        if ln.startswith('### DCA vs '+name): on=True; continue
        if on:
            if ln.startswith('###'): break
            if not ln.strip(): continue
            if ln.strip().startswith('x '): continue
            out.append(ln)
    return out
def pair2(name, hdr):
    a=sect(L1,name); b=sect(L4,name)
    m={r.split()[0]: r.split() for r in b}
    body=''
    for r in a:
        c=r.split(); d=m.get(c[0],['','','','',''])
        body+=('<tr><td class="val">'+c[0]+'</td>'
               +''.join(f'<td class="val">{v}</td>' for v in c[1:5])
               +''.join(f'<td class="val">{v}</td>' for v in d[1:5])+'</tr>\n')
    return ('<div class="wrap"><table><thead><tr>'+''.join(f'<th>{h}</th>' for h in hdr)+
            f'</tr></thead><tbody>\n{body}</tbody></table></div>')
PT  = pair2('pT', HDR)
ETA = pair2('eta', HDR)

HTML = f"""<title>Macro-Driven Monitoring</title>
<style>{STYLE}</style>

<div class="page">

<header class="mast">
  <p class="eyebrow">ALICE ITS2 &middot; its2-o2-decoupling &middot; run 901 &middot; 20 000 events</p>
  <h1>Monitoring on the Supplied Macros</h1>
  <p class="stand">DCA_y and DCA_z against p<sub>T</sub>, &eta;, z and &phi;; transverse and
  longitudinal residuals against &phi; and z, layer by layer. All produced by the macros you
  supplied, run unchanged, against the cluster-only fit &mdash; baseline and epoch 4 side by
  side.</p>
  <dl class="meta">
    <div><dt>Driver</dt><dd>tools/monitoring/make_monitoring_plots.sh</dd></div>
    <div><dt>Tree</dt><dd>TrkVtxer, 165 096 trk</dd></div>
    <div><dt>Fit</dt><dd>supplied fitter, scale = 2</dd></div>
    <div><dt>DCA_y at 1&ndash;1.25 GeV/c</dt><dd>36.4 &rarr; 35.9 &micro;m</dd></div>
    <div><dt>L4 transverse vs &phi;</dt><dd>&plusmn;50 &micro;m modulation</dd></div>
    <div><dt>Dropped</dt><dd>ratio-to-baseline panels</dd></div>
  </dl>
</header>

<section>
  <h2>How the macros are wired</h2>
  <p>The macros address the tree by the global names <code>TrkVtxer</code> and
  <code>ResMonitor</code>. <code>make_monitoring_plots.sh</code> points <strong>both</strong> at the
  <code>TrkVtxer</code> tree &mdash; the one <code>TrackVertexQualityEstimator</code> writes from
  <code>TrackerFit</code>, the fit that uses the clusters alone. The <code>ResMonitor</code> tree,
  written by the refit that carries the primary vertex as an eighth measured point, is not read at
  all, for either the DCA or the residuals.</p>
  <div class="wrap"><table>
    <thead><tr><th>Macro</th><th>Output</th></tr></thead>
    <tbody>
      <tr><td class="mono">check_vertex_pT_plots_trkvtxer_color</td><td>DCA_y, DCA_z vs p<sub>T</sub></td></tr>
      <tr><td class="mono">check_vertex_eta_plots_trkvtxer_color</td><td>DCA_y, DCA_z vs &eta;</td></tr>
      <tr><td class="mono">check_vertex_z_plots_trkvtxer_color</td><td>DCA_y, DCA_z vs z</td></tr>
      <tr><td class="mono">check_vertex_phi_plots_trkvtxer_color</td><td>DCA_y, DCA_z vs &phi;</td></tr>
      <tr><td class="mono">check_track_residualsHB_FitGaus</td><td>&Delta;s2 (XY) and &Delta;s1 (Z) vs &phi; and z, per layer, HB0/HB1 split</td></tr>
    </tbody>
  </table></div>
  <p>The p<sub>T</sub> macro is your <code>check_vertex_p_default_...</code> with the draw variable
  changed from <code>p</code> to <code>pT</code> and the axis titles relabelled; nothing else
  differs. Two things are done from the driver rather than by editing your files: the HB0/HB1
  overlays get separate colours and a key, and bins whose parent projection holds fewer than 200
  tracks are blanked &mdash; the fitter runs on every x bin including empty ones past the detector
  edge, and those came back with error bars larger than the pad.</p>
  <div class="note"><span class="tag">Dropped as asked</span>
  <p><code>trk_res_s1_ratio</code> and <code>trk_res_s2_ratio</code> are gone, and
  <code>plot_epoch_trends.C</code> no longer builds them.</p></div>
</section>

<section>
  <h2>DCA against p<sub>T</sub></h2>
  <figure><img src="data:image/png;base64,{img('dca_vs_pT_ep-1.png')}"
    alt="Four panels: DCA_y and DCA_z width and mean against transverse momentum at the baseline.">
  <figcaption><b>epoch &minus;1</b> Top row width, bottom row mean; DCA_y left, DCA_z
  right.</figcaption></figure>
  <figure><img src="data:image/png;base64,{img('dca_vs_pT_ep4.png')}"
    alt="The same four panels after five epochs of alignment.">
  <figcaption><b>epoch 4</b> Same axes.</figcaption></figure>
  {PT}
  <p>The width follows 1/p across the whole spectrum &mdash; 81.8&nbsp;&micro;m in the lowest bin
  down to 14&nbsp;&micro;m above 4&nbsp;GeV/c. Five epochs move it by well under a micron
  everywhere; the visible change is in the <code>DCA_y</code> mean, which comes down from
  +5.3 to +4.5&nbsp;&micro;m in the low bins.</p>
</section>

<section>
  <h2>DCA against &eta;</h2>
  <figure><img src="data:image/png;base64,{img('dca_vs_eta_ep-1.png')}"
    alt="DCA_y and DCA_z width and mean against pseudorapidity at the baseline.">
  <figcaption><b>epoch &minus;1</b> 20 bins over &minus;1&hellip;1.</figcaption></figure>
  <figure><img src="data:image/png;base64,{img('dca_vs_eta_ep4.png')}"
    alt="The same, after five epochs.">
  <figcaption><b>epoch 4</b></figcaption></figure>
  {ETA}
  <p><code>DCA_z</code>'s mean sweeps linearly across &eta; and crosses zero near mid-rapidity
  &mdash; a small longitudinal offset between the &chi;<sup>2</sup> vertex and where the tracks
  converge. Its width rises at the edges as 1/sin&theta;, as it must.</p>
</section>

<section>
  <h2>DCA against z and &phi;</h2>
  <figure><img src="data:image/png;base64,{img('dca_vs_z_ep-1.png')}"
    alt="DCA against the vertex z coordinate at the baseline.">
  <figcaption><b>vs z, epoch &minus;1</b> 30 bins over &plusmn;15&nbsp;cm, with the macro's own
  |z|&nbsp;&lt;&nbsp;10 cut applied.</figcaption></figure>
  <figure><img src="data:image/png;base64,{img('dca_vs_z_ep4.png')}"
    alt="DCA against the vertex z coordinate after five epochs.">
  <figcaption><b>vs z, epoch 4</b></figcaption></figure>
  <figure><img src="data:image/png;base64,{img('dca_vs_phi_ep-1.png')}"
    alt="DCA against azimuth at the baseline, 120 bins.">
  <figcaption><b>vs &phi;, epoch &minus;1</b> 120 bins. <code>DCA_y</code>'s mean carries a
  one-period modulation of roughly &plusmn;10&nbsp;&micro;m on top of its offset.</figcaption></figure>
  <figure><img src="data:image/png;base64,{img('dca_vs_phi_ep4.png')}"
    alt="DCA against azimuth after five epochs.">
  <figcaption><b>vs &phi;, epoch 4</b></figcaption></figure>
</section>

<section>
  <h2>Transverse residual, layer by layer</h2>
  <figure><img src="data:image/png;base64,{img('residualXY_vs_phi_ep-1.png')}"
    alt="Transverse residual mean against azimuth for seven layers at the baseline, showing a strong V-shaped modulation on the outer layers.">
  <figcaption><b>&Delta;s2 (XY) vs &phi;, epoch &minus;1</b> Inner barrel on the top row, outer on
  the bottom. Layers 3 and 4 swing &plusmn;50&nbsp;&micro;m and layers 5 and 6 &plusmn;100&nbsp;&micro;m
  with a clear one-period shape; the inner barrel stays inside
  &plusmn;5&nbsp;&micro;m.</figcaption></figure>
  <figure><img src="data:image/png;base64,{img('residualXY_vs_phi_ep4.png')}"
    alt="The same, after five epochs.">
  <figcaption><b>&Delta;s2 (XY) vs &phi;, epoch 4</b> The modulation is still there, essentially
  unchanged in amplitude. Five epochs shifted the layer means but did not flatten this
  &phi; structure.</figcaption></figure>
  <figure><img src="data:image/png;base64,{img('residualXY_vs_z_ep-1.png')}"
    alt="Transverse residual mean against z for seven layers, split by half-barrel, at the baseline.">
  <figcaption><b>&Delta;s2 (XY) vs z, epoch &minus;1</b> Black both half-barrels, blue HB0, orange
  HB1. Flat on the inner barrel; layers 5 and 6 sit 20&ndash;40&nbsp;&micro;m low across the whole
  z range with a shallow bow.</figcaption></figure>
  <figure><img src="data:image/png;base64,{img('residualXY_vs_z_ep4.png')}"
    alt="The same, after five epochs.">
  <figcaption><b>&Delta;s2 (XY) vs z, epoch 4</b></figcaption></figure>
</section>

<section>
  <h2>Longitudinal residual, layer by layer</h2>
  <figure><img src="data:image/png;base64,{img('residualZ_vs_phi_ep-1.png')}"
    alt="Longitudinal residual mean against azimuth for seven layers at the baseline, largely flat.">
  <figcaption><b>&Delta;s1 (Z) vs &phi;, epoch &minus;1</b> Flat to within 2&nbsp;&micro;m on the
  inner barrel and 10&nbsp;&micro;m on L3&ndash;L4. L5 and L6 carry a &minus;30&nbsp;&micro;m offset
  that grows toward &phi; &gt; 1 &mdash; the only longitudinal structure worth
  chasing.</figcaption></figure>
  <figure><img src="data:image/png;base64,{img('residualZ_vs_phi_ep4.png')}"
    alt="The same, after five epochs.">
  <figcaption><b>&Delta;s1 (Z) vs &phi;, epoch 4</b> Unchanged, consistent with the longitudinal
  width and mean not moving over the run.</figcaption></figure>
  <figure><img src="data:image/png;base64,{img('residualZ_vs_z_ep-1.png')}"
    alt="Longitudinal residual mean against z for seven layers, split by half-barrel.">
  <figcaption><b>&Delta;s1 (Z) vs z, epoch &minus;1</b> HB0 and HB1 separate visibly on L5 and L6
  &mdash; the two half-barrels are not at the same longitudinal position.</figcaption></figure>
  <figure><img src="data:image/png;base64,{img('residualZ_vs_z_ep4.png')}"
    alt="The same, after five epochs.">
  <figcaption><b>&Delta;s1 (Z) vs z, epoch 4</b></figcaption></figure>
</section>

<section>
  <h2>Epoch summary</h2>
  <p>The per-epoch views from <code>plot_epoch_trends.C</code>, kept alongside. Transverse before longitudinal;
  the ratio panels are removed.</p>
  <figure><img src="data:image/png;base64,{img2('residual_trend_vs_epoch.png')}"
    alt="Six panels of residual width and mean against epoch, transverse on top and longitudinal below.">
  <figcaption><b>residual vs epoch</b> Transverse narrows up to 2.3&nbsp;% and its mean marches;
  longitudinal is flat inside the fit error.</figcaption></figure>
  <figure><img src="data:image/png;base64,{img2('dca_trend_vs_epoch.png')}"
    alt="DCA mean and width against epoch.">
  <figcaption><b>DCA vs epoch</b> <code>DCA_y</code> mean +4.84 &rarr; +4.07&nbsp;&micro;m.</figcaption></figure>
  <figure><img src="data:image/png;base64,{img2('dca_width_vs_pT_by_epoch.png')}"
    alt="DCA width against transverse momentum with one curve per epoch.">
  <figcaption><b>DCA width vs p<sub>T</sub>, per epoch</b> All six curves on the same
  envelope.</figcaption></figure>
</section>

<footer>
  <p>Reproduce with <code>./make_monitoring_plots.sh &lt;epoch&gt; &lt;outdir&gt;</code> from
  the repository root; macros in <code>tools/monitoring/user_macros/</code>. Run 901, 20 000
  events, nEPOCH 5, nTrackMax 8, B = &minus;0.500673&nbsp;T. Fit scale 2, i.e. the
  &plusmn;2&sigma; pass of the supplied fitter after its half-maximum seed.</p>
</footer>

</div>
"""
open('report_macros.html','w').write(HTML)
print('wrote', len(HTML))
