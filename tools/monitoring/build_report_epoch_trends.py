import base64
STYLE = open('report_style.css').read()
def img(p): return base64.b64encode(open(p,'rb').read()).decode()
log = open('plot_epoch_trends.log').read()

def block(start):
    out=[]; on=False
    for ln in log.splitlines():
        if not on and ln.strip().startswith(start): on=True; continue
        if on:
            if not ln.strip(): break
            out.append(ln)
    return out

def tbl(rows, hdr, mono0=False):
    body=''
    for ln in rows:
        c=[x for x in ln.replace('|',' ').split()]
        body+='<tr>'+''.join(
            f'<td class="{"mono b" if (i==0 and mono0) else "val"}">{v}</td>' for i,v in enumerate(c))+'</tr>\n'
    return ('<div class="wrap"><table><thead><tr>'+''.join(f'<th>{h}</th>' for h in hdr)+
            f'</tr></thead><tbody>\n{body}</tbody></table></div>')

EPH = ['&minus;1','0','1','2','3','4']
wh = ['layer']+[f'&Delta;s2 trv e{e}' for e in EPH]+[f'&Delta;s1 lng e{e}' for e in EPH]
rows_all = log.splitlines()
i1 = [i for i,l in enumerate(rows_all) if l.startswith('lay      s2 e-1')]
W   = tbl(rows_all[i1[0]+1:i1[0]+8], wh, True)
M   = tbl(rows_all[i1[1]+1:i1[1]+8], wh, True)
CH  = tbl(block('lay     w2 trv'), ['layer','&Delta;s2 transverse width [%]','&Delta;s1 longitudinal width [%]','&Delta;s2 transverse mean shift [&micro;m]','&Delta;s1 longitudinal mean shift [&micro;m]'], True)
DCA = tbl(block('epoch     mean_y'), ['epoch','DCA_y mean','DCA_y width','DCA_z mean','DCA_z width'])
PT  = tbl(block('pT [GeV/c]       e-1'), ['p_{T} [GeV/c]']+[f'epoch {e}' for e in EPH])

HTML = f"""<title>Epoch-by-Epoch Monitoring</title>
<style>{STYLE}</style>

<div class="page">

<header class="mast">
  <p class="eyebrow">ALICE ITS2 &middot; its2-o2-decoupling &middot; run 901 &middot; 20 000 events</p>
  <h1>What Five Epochs Actually Moved</h1>
  <p class="stand">DCA and residuals for every epoch of the run, all from the cluster-only fit
  (<code>TrkVtxer</code>). Transverse is plotted before longitudinal throughout. The changes are
  fractions of a percent, so the overlaid distributions sit on top of each other &mdash; the
  per-epoch fits are where the movement is legible.</p>
  <dl class="meta">
    <div><dt>Epochs</dt><dd>&minus;1, 0, 1, 2, 3, 4</dd></div>
    <div><dt>Macro</dt><dd>monitor/plot_epoch_trends.C</dd></div>
    <div><dt>Transverse width change</dt><dd>up to &minus;2.3 %</dd></div>
    <div><dt>Longitudinal width change</dt><dd>&lt; 0.06 %, within error</dd></div>
    <div><dt>Transverse mean, L6</dt><dd>&minus;28.5 &rarr; &minus;19.1 &micro;m</dd></div>
    <div><dt>DCA_y width</dt><dd>42.93 &rarr; 42.60 &micro;m</dd></div>
  </dl>
</header>

<section>
  <h2>Which coordinate is which, and what moved</h2>
  <p>The two sensor-plane residuals are not interchangeable, and the code says which is which at
  the vertex point: <code>Residual_s1 = proj_GZc &minus; meas_GZc</code> is longitudinal, while
  <code>Residual_s2 = &plusmn;&radic;(&Delta;GX&sup2;+&Delta;GY&sup2;)</code> is transverse. So
  <strong>s2 pairs with DCA_y and s1 with DCA_z</strong>, and transverse leads in every panel and
  table below.</p>
  <div class="note"><span class="tag">The correction is transverse, and it is coherent</span>
  <p>Across five epochs the transverse width narrows by up to 2.3&nbsp;% and its mean marches
  monotonically &mdash; L6 from &minus;28.5 to &minus;19.1&nbsp;&micro;m, L5 from &minus;22.3 to
  &minus;14.2, L4 from &minus;7.1 to &minus;0.5. The longitudinal width changes by at most
  0.06&nbsp;% and its mean is stable to 0.03&nbsp;&micro;m, both inside the fit error.</p>
  <p>That lines up with the DCA: <code>DCA_y</code>'s mean walks from +4.84 to
  +4.07&nbsp;&micro;m while <code>DCA_z</code> sits at zero throughout. Both measurements say the
  same thing &mdash; this run is correcting a transverse mode and leaving the longitudinal one
  alone. Read together they are consistent, not a sign that one gradient path is dead.</p></div>
  {CH}
  <p>Change from epoch &minus;1 to epoch 4, transverse first. Widths in percent, mean shifts in
  &micro;m.</p>
</section>

<section>
  <h2>Residual width and mean against epoch</h2>
  <figure><img src="data:image/png;base64,{img('residual_trend_vs_epoch.png')}"
    alt="Six panels: residual width absolute and relative, and mean, against epoch for each layer, s1 on top and s2 below.">
  <figcaption><b>vs epoch</b> Top row transverse (&Delta;s2), bottom row longitudinal
  (&Delta;s1). Left: absolute width, log scale, seven layers spanning 2&ndash;415&nbsp;&micro;m.
  Middle: the same relative to epoch &minus;1, with fit errors &mdash; the transverse points move
  away from zero, the longitudinal ones do not. Right: fitted mean.</figcaption></figure>
  {W}
  <p>Fitted widths in &micro;m. Below, the fitted means.</p>
  {M}
</section>

<section>
  <h2>Distributions, every epoch overlaid</h2>
  <figure><img src="data:image/png;base64,{img('residual_XY_by_epoch.png')}"
    alt="Transverse residual distributions per layer with all six epochs overlaid, plus the DCA_y distribution.">
  <figcaption><b>&Delta;s2, transverse</b> Blue is the baseline, red is epoch 4. Seven layers plus
  <code>DCA_y</code> in the last pad &mdash; the cluster-only fit has no vertex point, so there is
  no VTX residual to draw there, and the transverse impact parameter is what that pad would have
  shown. Ranges are set per layer at about 4&sigma;.</figcaption></figure>
  <figure><img src="data:image/png;base64,{img('residual_Z_by_epoch.png')}"
    alt="Longitudinal residual distributions per layer with all six epochs overlaid, plus the DCA_z distribution.">
  <figcaption><b>&Delta;s1, longitudinal</b> Same layout, paired with
  <code>DCA_z</code>.</figcaption></figure>
</section>

<section>
  <h2>DCA against epoch</h2>
  <figure><img src="data:image/png;base64,{img('dca_trend_vs_epoch.png')}"
    alt="DCA mean and width against epoch for the y and z components.">
  <figcaption><b>vs epoch</b> Transverse plotted first in both pads. <code>DCA_y</code>'s mean
  falls steadily from +4.84 to +4.07&nbsp;&micro;m &mdash; a real 0.8&nbsp;&micro;m walk toward
  zero. Its width improves 0.33&nbsp;&micro;m. <code>DCA_z</code> starts centred and stays
  there.</figcaption></figure>
  {DCA}
  <figure><img src="data:image/png;base64,{img('dca_width_vs_pT_by_epoch.png')}"
    alt="DCA width against transverse momentum with one curve per epoch, all lying on the multiple-scattering curve.">
  <figcaption><b>width vs p<sub>T</sub>, per epoch</b> Transverse left, longitudinal right.
  All six curves lie on the same
  1/p envelope: 73&nbsp;&micro;m in the lowest bin, ~31&nbsp;&micro;m at 1&ndash;1.5&nbsp;GeV/c,
  10&ndash;16&nbsp;&micro;m above 4&nbsp;GeV/c. Five epochs of alignment do not move a
  scattering-dominated resolution, and should not.</figcaption></figure>
  {PT}
  <p>The gain is concentrated where scattering is not dominant: 31.80&nbsp;&rarr;&nbsp;31.14&nbsp;&micro;m
  at 1&ndash;1.5&nbsp;GeV/c, 26.61&nbsp;&rarr;&nbsp;26.01 at 1.5&ndash;2, against
  73.43&nbsp;&rarr;&nbsp;73.23 in the lowest bin where the material term swamps everything.</p>
</section>

<section>
  <h2>Correction to the widths I quoted before</h2>
  <div class="note warn"><span class="tag">Range-clipped and under-binned</span>
  <p>The earlier tables used one histogram range for all layers &mdash; &plusmn;100&nbsp;&micro;m
  for the inner barrel and &plusmn;500 for the outer. At L6 that truncates at 1.5&sigma;, and at L0
  a 3.3&nbsp;&micro;m bin cannot sample a 2&nbsp;&micro;m peak. Both distort the fit. The ranges
  are now set per layer at roughly 4&sigma;, 60 bins each.</p></div>
  <div class="wrap"><table>
    <thead><tr><th>layer</th><th>&Delta;s1 width as quoted</th><th>corrected</th><th>why</th></tr></thead>
    <tbody>
      <tr><td class="mono b">L0</td><td class="val">1.85</td><td class="val">2.08</td>
          <td>&plusmn;100&nbsp;&micro;m range, 3.3&nbsp;&micro;m bins &mdash; one bin across the peak</td></tr>
      <tr><td class="mono b">L5</td><td class="val">274.66</td><td class="val">293.44</td>
          <td>&plusmn;500&nbsp;&micro;m range clipped the core at 1.8&sigma;</td></tr>
      <tr><td class="mono b">L6</td><td class="val">335.28</td><td class="val">412.41</td>
          <td>same, clipped at 1.5&sigma;</td></tr>
    </tbody>
  </table></div>
  <p>Verified by scanning the range: L6 gives 335 &micro;m at &plusmn;500 and then 410, 412, 406,
  408&nbsp;&micro;m at &plusmn;1000 through &plusmn;4000 &mdash; a plateau. L5 gives 275 then
  295, 293, 291, 299. L0 gives 2.08 at &plusmn;5 through &plusmn;20 and only breaks at
  &plusmn;100. The new settings sit inside the plateau on every layer. The inner-barrel and DCA
  conclusions are unaffected; the outer-barrel widths are about 20&nbsp;% larger than I
  said.</p>
  <p>L0's 2.08&nbsp;&micro;m is still a biased residual &mdash; the layer is one of seven points in
  a five-parameter fit. It tracks an alignment step honestly; it is not a resolution.</p>
</section>

<footer>
  <p>Transverse before longitudinal in every panel and table. Run 901, 20 000 events, nEPOCH 5,
  nTrackMax 8, B = &minus;0.500673&nbsp;T. Tracks flagged
  <code>used</code> by the vertexer, from the <code>TrkVtxer</code> tree written by
  <code>TrackVertexQualityEstimator</code>. Widths and means from an iterative &plusmn;2&sigma;
  Gaussian fit seeded on the half-maximum points. Macro: <code>monitor/plot_epoch_trends.C</code>.</p>
</footer>

</div>
"""
open('report_epochs.html','w').write(HTML)
print('wrote', len(HTML))
