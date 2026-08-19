import base64, re, io, sys

STYLE = open('report_style.css').read()

def img(p):
    return base64.b64encode(open(p,'rb').read()).decode()

log = open('compare_dca_reference.log').read()

def block(start):
    out=[]; on=False
    for ln in log.splitlines():
        if ln.strip().startswith(start): on=True; continue
        if on:
            if not ln.strip(): break
            out.append(ln)
    return out

def table_from(tag, hdr):
    rows = block(tag)
    body=''
    for ln in rows:
        c = ln.split()
        body += '<tr>'+''.join(f'<td class="val">{v}</td>' for v in c)+'</tr>\n'
    th=''.join(f'<th>{h}</th>' for h in hdr)
    return f'<div class="wrap"><table><thead><tr>{th}</tr></thead><tbody>\n{body}</tbody></table></div>'

pt_e = table_from('[v_{est}] pT',  ['pT [GeV/c]','DCA_y mean','DCA_y width','DCA_z mean','DCA_z width'])
pt_r = table_from('[v_{reco}] pT', ['pT [GeV/c]','DCA_y mean','DCA_y width','DCA_z mean','DCA_z width'])

res = block('lay   |')
resrows=''
for ln in res:
    c=[x for x in ln.replace('|',' ').split()]
    resrows += '<tr><td class="mono b">'+c[0]+'</td>'+''.join(f'<td class="val">{v}</td>' for v in c[1:])+'</tr>\n'

HTML = f"""<title>Vertex Reference and DCA</title>
<style>{STYLE}</style>

<div class="page">

<header class="mast">
  <p class="eyebrow">ALICE ITS2 &middot; its2-o2-decoupling &middot; run 901 &middot; 20 000 events</p>
  <h1>The DCA Was Measured Against Its Own Fit</h1>
  <p class="stand">The impact parameter the module reports is taken against <code>v_est</code>, the
  adaptive vertex fitted from the same eight tracks &mdash; which also enters the circle fit as an
  eighth measured point. The width therefore collapses onto the sigma assigned to that point and
  loses its 1/p dependence. Measured against <code>v_reco</code> instead, the multiple-scattering
  behaviour reappears.</p>
  <dl class="meta">
    <div><dt>Reference in code</dt><dd>BeamCenter = fvertex_TRKF</dd></div>
    <div><dt>Call site</dt><dd>YMultiLayerPerceptron.cxx:5665</dd></div>
    <div><dt>DCA_y width vs v_est</dt><dd>3.19 &micro;m, flat</dd></div>
    <div><dt>DCA_y width vs v_reco</dt><dd>28.1 &rarr; 14.2 &micro;m</dd></div>
    <div><dt>Sigma_MEAS[7]</dt><dd>4.74 &micro;m</dd></div>
    <div><dt>Sigma_MSC[7]</dt><dd>3.32 &micro;m</dd></div>
  </dl>
</header>

<section>
  <h2>Which vertex is which</h2>
  <div class="wrap"><table>
    <thead><tr><th>Quantity</th><th>Variable in code</th><th>Branch</th><th>Role</th></tr></thead>
    <tbody>
      <tr><td class="mono b">v_reco</td><td class="mono">BeamPos[0..2] = tv1,tv2,tv3</td>
          <td class="mono">vtxX/Y/Z</td><td>read from the input file, recorded only</td></tr>
      <tr><td class="mono b">v_est</td><td class="mono">BeamCenter = fvertex_TRKF</td>
          <td class="mono">vtxevtX/Y/Z</td><td><strong>used as the DCA reference</strong></td></tr>
      <tr><td class="mono b">v_fit</td><td class="mono">Fitpar[2],Fitpar[3],parz</td>
          <td class="mono">vtxfitX/Y/Z</td><td>the track's own extrapolation</td></tr>
    </tbody>
  </table></div>
  <div class="note warn"><span class="tag">The circularity is double</span>
  <p><code>UpdateVertexByAlignment</code> builds <code>fvertex_TRKF</code> from the pairwise closest
  approach of the very tracks in the group. That same point is then appended to the circle fit as
  layer index 7, carrying <code>Sigma_MEAS[7] = 4.74&nbsp;&micro;m</code> and
  <code>Sigma_MSC[7] = 3.32&nbsp;&micro;m</code>. So the quantity called DCA is the fit residual of
  a constrained point, and it lands on the sigma it was given &mdash; the measured 3.19 / 3.33&nbsp;&micro;m
  are the same numbers as the VTX row of the residual table below, to two decimals.</p></div>
</section>

<section>
  <h2>Same tracks, two references</h2>
  <figure><img src="data:image/png;base64,{img('dca_reference.png')}"
    alt="DCA_y width against transverse momentum for two vertex references: flat near 3 micron against v_est, falling from 28 to 14 micron against v_reco.">
  <figcaption><b>The whole point</b> Grey: against <code>v_est</code>, flat at ~3&nbsp;&micro;m
  across the whole spectrum. Blue: against <code>v_reco</code>, 28.1&nbsp;&micro;m in the lowest pT
  bin falling to 14&ndash;17&nbsp;&micro;m above 4&nbsp;GeV/c. Only the blue curve behaves like an
  impact-parameter resolution.</figcaption></figure>
  <div class="tiles">
    <div class="warn"><span class="n">3.19</span><span class="k">DCA_y width vs v_est [&micro;m]</span></div>
    <div class="good"><span class="n">24.39</span><span class="k">DCA_y width vs v_reco [&micro;m]</span></div>
    <div class="good"><span class="n">28.1</span><span class="k">at 0&ndash;0.5 GeV/c [&micro;m]</span></div>
    <div class="good"><span class="n">21.5</span><span class="k">at 1&ndash;1.5 GeV/c [&micro;m]</span></div>
    <div><span class="n">14.2</span><span class="k">at 5.5&ndash;6 GeV/c [&micro;m]</span></div>
  </div>
  <p>The expectation was roughly 30&nbsp;&micro;m for a 1&nbsp;GeV/c particle from multiple
  scattering alone. Against <code>v_reco</code> the measurement gives 21.5&nbsp;&micro;m at
  1&ndash;1.5&nbsp;GeV/c and 24.7&nbsp;&micro;m at 0.5&ndash;1&nbsp;GeV/c, with the falling 1/p
  trend intact and a floor near 14&ndash;17&nbsp;&micro;m where scattering stops dominating. That is
  the right order and the right shape.</p>
</section>

<section>
  <h2>DCA against pT &mdash; against v_est</h2>
  <figure><img src="data:image/png;base64,{img('dca_vs_pT_vest.png')}"
    alt="DCA mean and width against transverse momentum, measured against the estimated vertex.">
  <figcaption><b>vs pT, v_est</b> Width flat at 2.7&ndash;3.8&nbsp;&micro;m. No momentum
  dependence at all &mdash; the signature of a constrained fit residual.</figcaption></figure>
  {pt_e}
</section>

<section>
  <h2>DCA against pT &mdash; against v_reco</h2>
  <figure><img src="data:image/png;base64,{img('dca_vs_pT_vreco.png')}"
    alt="DCA mean and width against transverse momentum, measured against the reconstructed vertex.">
  <figcaption><b>vs pT, v_reco</b> The width falls monotonically with momentum. The
  <code>DCA_z</code> mean sits at &minus;1808&nbsp;&micro;m in every bin &mdash; that is the known
  vertex-z discrepancy, rigid and unrelated to momentum.</figcaption></figure>
  {pt_r}
</section>

<section>
  <h2>eta, z and phi</h2>
  <figure><img src="data:image/png;base64,{img('dca_vs_eta_vreco.png')}"
    alt="DCA mean and width against pseudorapidity, against the reconstructed vertex.">
  <figcaption><b>vs &eta;</b> 20 bins over &minus;1&hellip;1. The width is flat at
  23&ndash;25&nbsp;&micro;m through the acceptance and rises at the edges where the track crosses
  more material. <code>DCA_z</code> drifts from &minus;1814 to &minus;1805&nbsp;&micro;m across the
  range, a 9&nbsp;&micro;m tilt on top of the 1.8&nbsp;mm offset.</figcaption></figure>
  <figure><img src="data:image/png;base64,{img('dca_vs_zvtx_vreco.png')}"
    alt="DCA mean and width against the reconstructed vertex z, against the reconstructed vertex.">
  <figcaption><b>vs z</b> 60 bins over &plusmn;15&nbsp;cm. Flat in the populated region; the
  outermost bins carry few tracks and their fits wander.</figcaption></figure>
  <figure><img src="data:image/png;base64,{img('dca_vs_phi_vreco.png')}"
    alt="DCA mean and width against azimuth, against the reconstructed vertex.">
  <figcaption><b>vs &phi;</b> 12 bins. Width flat at 22&ndash;26&nbsp;&micro;m. The
  1&nbsp;&micro;m sinusoidal modulation in <code>DCA_y</code> that was visible against
  <code>v_est</code> is buried under the 24&nbsp;&micro;m width here &mdash; the sharp reference is
  the more sensitive probe of that particular weak mode, even though it is not an impact
  parameter.</figcaption></figure>
  <figure><img src="data:image/png;base64,{img('dca_vs_phi_vest.png')}"
    alt="DCA mean and width against azimuth, against the estimated vertex.">
  <figcaption><b>vs &phi;, v_est</b> Kept for that reason: the coherent &phi; modulation of the
  mean is what the vertex constraint exists to suppress, and it reads more cleanly here.</figcaption></figure>
</section>

<section>
  <h2>Residuals and cost, unchanged</h2>
  <p>Nothing in this analysis touched the module, so the training is the same run reported before.
  Widths from the same iterative fit, in &micro;m.</p>
  <div class="wrap"><table>
    <thead><tr><th>layer</th><th>&Delta;s1 mean (&minus;1)</th><th>&Delta;s1 mean (4)</th>
    <th>&Delta;s1 width (&minus;1)</th><th>&Delta;s1 width (4)</th><th>&Delta;s2 mean (&minus;1)</th>
    <th>&Delta;s2 mean (4)</th><th>&Delta;s2 width (&minus;1)</th><th>&Delta;s2 width (4)</th></tr></thead>
    <tbody>{resrows}</tbody>
  </table></div>
  <p>The <code>VTX</code> row is the giveaway: 3.33 and 3.18&nbsp;&micro;m, the same values the
  DCA-against-<code>v_est</code> fit returns. It is one and the same residual.</p>
  <figure><img src="data:image/png;base64,{img('residual_s1.png')}"
    alt="Delta-s1 residual distributions for seven layers and the vertex point, 60 bins each.">
  <figcaption><b>Residual s1</b> 60 bins, &plusmn;100&nbsp;&micro;m for IB and VTX,
  &plusmn;500&nbsp;&micro;m for OB. Epoch &minus;1 solid, epoch 4 dashed, fitted Gaussian in
  black.</figcaption></figure>
  <figure><img src="data:image/png;base64,{img('residual_s2.png')}"
    alt="Delta-s2 residual distributions for seven layers and the vertex point.">
  <figcaption><b>Residual s2</b> Same layout for the second sensor-plane coordinate.</figcaption></figure>
  <figure><img src="data:image/png;base64,{img('cost_epoch.png')}"
    alt="Cost per degree of freedom against epoch for the training and test samples, both falling monotonically.">
  <figcaption><b>Cost vs epoch</b> Training &minus;0.848&nbsp;%, test &minus;0.677&nbsp;%, both
  monotonic, no epoch-guard backtracks.</figcaption></figure>
</section>

<section>
  <h2>What this does and does not settle</h2>
  <ul>
    <li><strong>Settled.</strong> The reported DCA width was small because the reference was the
    module's own vertex fit. Against an independent reference the width is 24&nbsp;&micro;m overall
    with the expected 1/p shape.</li>
    <li><strong>Not changed.</strong> The module still measures against <code>v_est</code>. That
    value feeds the <code>RANGE_IMPACTPARAMS_R</code> = 0.2&nbsp;cm selection at
    <code>:5666</code>, so switching it inside the module is a physics-selection change, not a
    monitoring change. Everything here was recomputed offline from
    <code>curvX/curvY/curvR</code>.</li>
    <li><strong>Still open.</strong> The 1.81&nbsp;mm gap between <code>v_est,z</code> and
    <code>tv3</code>. An independent check from the stored cluster positions alone put it at
    1.783&nbsp;mm, so the clusters agree with <code>v_est</code>, not with <code>tv3</code>. That
    is an input question, not a module defect.</li>
    <li><strong>Still unverified.</strong> <code>Angle2Alpha</code> and <code>kB2C</code> in
    <code>YO2Compat.h</code> &mdash; and <code>kB2C</code> now sets the seed radius as well as the
    impact parameter.</li>
  </ul>
  <div class="note"><span class="tag">If the module should report the physical DCA</span>
  <p>Two ways. Measure against <code>v_reco</code> &mdash; one-line change at <code>:5665</code>,
  but it inherits the 1.8&nbsp;mm z offset until that is resolved. Or leave one track out of the
  vertex estimate when measuring it, which keeps <code>v_est</code>'s accuracy and removes the
  circularity, at the cost of eight vertex fits per group instead of one.</p></div>
</section>

<footer>
  <p>Run 901, 20 000 events, nEPOCH 5, nTrackMax 8, B = &minus;0.500673&nbsp;T. Geometry cache
  verified against O<sup>2</sup> <code>getMatrixL2G</code> for all 24 120 chips. Widths from an
  iterative &plusmn;2&sigma; Gaussian fit seeded on the half-maximum points. Macro:
  <code>monitor/compare_dca_reference.C</code>.</p>
</footer>

</div>
"""
open('report_dca_ref.html','w').write(HTML)
print('wrote', len(HTML))
