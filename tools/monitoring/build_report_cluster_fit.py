import base64
STYLE = open('report_style.css').read()
def img(p): return base64.b64encode(open(p,'rb').read()).decode()
log = open('compare_track_fits.log').read()
ulog = open('supplied_macros_scan.log').read()

def block(txt, start, stop_blank=True):
    out=[]; on=False
    for ln in txt.splitlines():
        if not on and ln.strip().startswith(start): on=True; continue
        if on:
            if not ln.strip(): break
            out.append(ln)
    return out

def tbl(rows, hdr, first_mono=False):
    body=''
    for ln in rows:
        c=[x for x in ln.replace('|',' ').split()]
        cells=''
        for i,v in enumerate(c):
            cls='mono b' if (i==0 and first_mono) else 'val'
            cells+=f'<td class="{cls}">{v}</td>'
        body+=f'<tr>{cells}</tr>\n'
    th=''.join(f'<th>{h}</th>' for h in hdr)
    return f'<div class="wrap"><table><thead><tr>{th}</tr></thead><tbody>\n{body}</tbody></table></div>'

cmp_rows = block(log,'lay   |  s1 cluster')
res_rows = block(log,'lay   |    m1(-1)')
pt_rows  = block(log,'pT         mean_y')
eta_rows = block(log,'eta        mean_y')
phi_rows = block(log,'phi        mean_y')
up_rows  = block(ulog,'x     mean_y')

CMP = tbl(cmp_rows, ['layer','&Delta;s1 cluster-only','&Delta;s1 with vertex','&Delta;s2 cluster-only','&Delta;s2 with vertex'], True)
RES = tbl(res_rows, ['layer','&Delta;s1 mean (&minus;1)','&Delta;s1 mean (4)','&Delta;s1 width (&minus;1)','&Delta;s1 width (4)',
                     '&Delta;s2 mean (&minus;1)','&Delta;s2 mean (4)','&Delta;s2 width (&minus;1)','&Delta;s2 width (4)'], True)
PT  = tbl(pt_rows,  ['p_{T} [GeV/c]','DCA_y mean','DCA_y width','DCA_z mean','DCA_z width'])
ETA = tbl(eta_rows, ['&eta;','DCA_y mean','DCA_y width','DCA_z mean','DCA_z width'])
PHI = tbl(phi_rows, ['&phi; [rad]','DCA_y mean','DCA_y width','DCA_z mean','DCA_z width'])
UP  = tbl(up_rows,  ['p [GeV/c]','DCA_y mean','DCA_y width','DCA_z mean','DCA_z width'])

HTML = f"""<title>Cluster-Only Fit Monitoring</title>
<style>{STYLE}</style>

<div class="page">

<header class="mast">
  <p class="eyebrow">ALICE ITS2 &middot; its2-o2-decoupling &middot; run 901 &middot; 20 000 events</p>
  <h1>Monitoring From the Cluster-Only Fit</h1>
  <p class="stand">The module fits every track twice. The impact parameter and the residuals
  belong to the first fit, from clusters alone; the second exists only to supply the alignment
  gradient and carries the primary vertex as a fitted point. Reading the second one collapsed the
  DCA width onto the sigma that fit assigns to the vertex. Everything here is redone on the
  <code>TrkVtxer</code> tree.</p>
  <dl class="meta">
    <div><dt>Tree</dt><dd>TrkVtxer (165 096 trk)</dd></div>
    <div><dt>Fit</dt><dd>TrackerFit &mdash; clusters only</dd></div>
    <div><dt>DCA_y width</dt><dd>42.9 &micro;m integrated</dd></div>
    <div><dt>at 0&ndash;0.5 GeV/c</dt><dd>73.4 &micro;m</dd></div>
    <div><dt>at 1&ndash;1.5 GeV/c</dt><dd>31.8 &micro;m</dd></div>
    <div><dt>at 5.5&ndash;6 GeV/c</dt><dd>10.3 &micro;m</dd></div>
  </dl>
</header>

<section>
  <h2>Two fits, and which one owns which number</h2>
  <div class="wrap"><table>
    <thead><tr><th>Fit</th><th>Where</th><th>Tree</th><th>What it is for</th></tr></thead>
    <tbody>
      <tr><td class="mono b">TrackerFit</td><td class="mono">UpdateVertexByAlignment :4649</td>
          <td class="mono">TrkVtxer</td><td>clusters only, 7 points &mdash; <strong>the physical trajectory</strong></td></tr>
      <tr><td class="mono b">GetCost_Beam_CircleFit</td><td class="mono">:5259&hellip;5665</td>
          <td class="mono">ResMonitor</td><td>clusters + primary vertex as an 8th point &mdash; the alignment gradient</td></tr>
    </tbody>
  </table></div>
  <p><code>TrackVertexQualityEstimator</code> receives <code>TRKF_fparXY</code> /
  <code>TRKF_fparZR</code> and computes both the residuals and <code>fip[0..1]</code> from them,
  against the &chi;<sup>2</sup> vertex <code>mVertex_Chi2[2]</code>. The vertex is not a point in
  that fit, so the impact parameter is a real distance of closest approach.</p>
  <div class="note warn"><span class="tag">What reading the wrong tree cost</span>
  <p>In <code>ResMonitor</code> the vertex enters the fit carrying
  <code>Sigma_MEAS[7] = 4.74&nbsp;&micro;m</code> and <code>Sigma_MSC[7] = 3.32&nbsp;&micro;m</code>.
  Its <code>fip</code> is therefore the residual of a constrained point and lands on that sigma
  &mdash; 3.19&nbsp;&micro;m, flat in p<sub>T</sub>. Every layer residual is inflated the same way,
  because the extra constraint pulls the trajectory off the clusters.</p></div>
  {CMP}
  <p>Widths in &micro;m at epoch &minus;1. The outer barrel is worst hit: L6 reads
  335&nbsp;&micro;m from the clusters and 503&nbsp;&micro;m once the vertex is in the fit. The
  <code>VTX</code> row is zero on the cluster side because there is no vertex point to have a
  residual against.</p>
</section>

<section>
  <h2>DCA against p<sub>T</sub></h2>
  <figure><img src="data:image/png;base64,{img('dca_clusterfit_vs_pT.png')}"
    alt="DCA mean and width against transverse momentum from the cluster-only fit, width falling from 73 to 10 micron.">
  <figcaption><b>vs p<sub>T</sub></b> 0.5&nbsp;GeV/c bins. The width falls monotonically from
  73.4&nbsp;&micro;m in the lowest bin to 10&ndash;16&nbsp;&micro;m above 4&nbsp;GeV/c.</figcaption></figure>
  {PT}
  <p>This is the multiple-scattering behaviour that was missing. At 1&ndash;1.5&nbsp;GeV/c the
  transverse width is 31.8&nbsp;&micro;m, matching the ~30&nbsp;&micro;m expected for a
  1&nbsp;GeV/c particle. <code>DCA_z</code> tracks it at 34.4&nbsp;&micro;m in the same bin and its
  mean sits at zero throughout &mdash; the 1.8&nbsp;mm offset seen elsewhere is between
  <code>tv3</code> and the &chi;<sup>2</sup> vertex, and does not enter here because this DCA is
  measured against the latter.</p>
  <figure><img src="data:image/png;base64,{img('supplied_dca_vs_p.png')}"
    alt="Four-panel DCA width and mean against total momentum produced by the supplied macro.">
  <figcaption><b>Supplied macro</b> <code>check_vertex_p_default_plots_trkvtxer_color(4,21,2.0)</code>
  run unmodified on the same tree, binned in total p rather than p<sub>T</sub>. Same
  behaviour: 84.3&nbsp;&micro;m at 0.3&ndash;0.5&nbsp;GeV/c down to 15.2&nbsp;&micro;m above
  6&nbsp;GeV/c.</figcaption></figure>
  {UP}
</section>

<section>
  <h2>&eta;, z and &phi;</h2>
  <figure><img src="data:image/png;base64,{img('dca_clusterfit_vs_eta.png')}"
    alt="DCA mean and width against pseudorapidity from the cluster-only fit.">
  <figcaption><b>vs &eta;</b> 20 bins over &minus;1&hellip;1. <code>DCA_y</code> width is flat at
  42&ndash;45&nbsp;&micro;m. <code>DCA_z</code> width rises from 36&nbsp;&micro;m at &eta;&nbsp;=&nbsp;0
  to 58&ndash;60&nbsp;&micro;m at the edges, the 1/sin&theta; projection of the same
  scattering.</figcaption></figure>
  {ETA}
  <p>The <code>DCA_z</code> mean sweeps linearly from &minus;10.2&nbsp;&micro;m at
  &eta;&nbsp;=&nbsp;&minus;0.75 to +9.6&nbsp;&micro;m at +0.75, crossing zero at mid-rapidity.
  That is the signature of a small longitudinal offset between the &chi;<sup>2</sup> vertex and
  where the tracks actually converge &mdash; a 20&nbsp;&micro;m swing, not the 1.8&nbsp;mm one.</p>
  <figure><img src="data:image/png;base64,{img('dca_clusterfit_vs_phi.png')}"
    alt="DCA mean and width against azimuth, showing a clear one-period sinusoid in the y mean.">
  <figcaption><b>vs &phi;</b> 12 bins. The width is flat at 39&ndash;46&nbsp;&micro;m, but
  <code>DCA_y</code>'s mean traces a clean one-period sinusoid from &minus;8.7&nbsp;&micro;m near
  &phi;&nbsp;=&nbsp;0.3 to +17.4&nbsp;&micro;m near &phi;&nbsp;=&nbsp;5.5.</figcaption></figure>
  {PHI}
  <div class="note"><span class="tag">The &phi; modulation is real and it is bigger than it looked</span>
  <p>Against the constrained fit this same modulation measured about 1&nbsp;&micro;m. From the
  clusters it is a 13&nbsp;&micro;m amplitude, with <code>DCA_z</code> showing a second sinusoid of
  about 6&nbsp;&micro;m in quadrature. A one-period &phi; dependence of the transverse impact
  parameter is what a transverse displacement between the assumed and the true beam position
  produces. The constrained fit was absorbing it into the vertex point.</p></div>
  <figure><img src="data:image/png;base64,{img('dca_clusterfit_vs_zvtx.png')}"
    alt="DCA mean and width against the vertex z position.">
  <figcaption><b>vs z</b> 60 bins over &plusmn;15&nbsp;cm, against <code>tv3</code>. Flat where the
  sample is populated.</figcaption></figure>
  <figure><img src="data:image/png;base64,{img('supplied_dca_vs_eta.png')}"
    alt="Four-panel DCA plots against pseudorapidity from the supplied macro.">
  <figcaption><b>Supplied macro, &eta;</b> Same tree, same conclusion, the supplied fitting
  convention (half-maximum window, then &plusmn;2&sigma;).</figcaption></figure>
  <figure><img src="data:image/png;base64,{img('supplied_dca_vs_phi.png')}"
    alt="Four-panel DCA plots against azimuth from the supplied macro, 120 bins.">
  <figcaption><b>Supplied macro, &phi;</b> 120 bins resolve the sinusoid directly.</figcaption></figure>
</section>

<section>
  <h2>Residuals</h2>
  {RES}
  <figure><img src="data:image/png;base64,{img('residual_Z_distributions.png')}"
    alt="Delta-s1 residual distributions per layer from the cluster-only fit.">
  <figcaption><b>Residual s1</b> 60 bins, &plusmn;100&nbsp;&micro;m for the inner barrel,
  &plusmn;500&nbsp;&micro;m for the outer. Epoch &minus;1 solid, epoch 4 dashed, fitted Gaussian in
  black.</figcaption></figure>
  <figure><img src="data:image/png;base64,{img('residual_XY_distributions.png')}"
    alt="Delta-s2 residual distributions per layer from the cluster-only fit.">
  <figcaption><b>Residual s2</b> Same layout for the second sensor-plane coordinate.</figcaption></figure>
  <div class="note warn"><span class="tag">Read L0 with care</span>
  <p>1.85&nbsp;&micro;m on layer 0 is below the ALPIDE intrinsic resolution, which it must be:
  these are biased residuals. The layer is one of seven points in a five-parameter fit, so the
  trajectory follows it. The number is fine for tracking an alignment step epoch to epoch, but it
  is not a resolution. Unbiased residuals need the layer excluded from the fit, which
  <code>TrackerFit</code> does not currently offer.</p></div>
  <p>Layers 5 and 6 carry means of &minus;12 and &minus;18&nbsp;&micro;m in &Delta;s1 and
  &minus;20 and &minus;25&nbsp;&micro;m in &Delta;s2 at epoch &minus;1. Those are coherent outer-barrel
  offsets, and five epochs move them by under 2&nbsp;&micro;m &mdash; consistent with the
  &minus;0.85&nbsp;% cost change over the same run.</p>
</section>

<section>
  <h2>Where this leaves things</h2>
  <ul>
    <li><strong>Corrected.</strong> The 3&nbsp;&micro;m DCA width reported earlier came from
    <code>ResMonitor</code>, the fit that contains the vertex. From the cluster-only fit it is
    42.9&nbsp;&micro;m integrated, 31.8&nbsp;&micro;m at 1&ndash;1.5&nbsp;GeV/c, with a clean 1/p
    falloff. The same correction applies to every residual width quoted from that run &mdash; the
    outer-barrel numbers were inflated by roughly 1.5&times;.</li>
    <li><strong>New.</strong> The &phi; modulation of <code>DCA_y</code> is 13&nbsp;&micro;m, not
    1&nbsp;&micro;m. Worth acting on rather than watching.</li>
    <li><strong>Unchanged.</strong> The 1.81&nbsp;mm gap between <code>tv3</code> and the
    &chi;<sup>2</sup> vertex z. It does not appear in these plots because the DCA is measured
    against the latter, but it is still there in <code>vtxevtZ &minus; vtxZ</code>.</li>
    <li><strong>Still unverified.</strong> <code>Angle2Alpha</code> and <code>kB2C</code> in
    <code>YO2Compat.h</code>. <code>kB2C</code> sets <code>track_tpc_R</code> in
    <code>TrackVertexQualityEstimator</code> too, so it feeds these numbers.</li>
  </ul>
</section>

<footer>
  <p>Run 901, 20 000 events, nEPOCH 5, nTrackMax 8, B = &minus;0.500673&nbsp;T. Tracks flagged
  <code>used</code> by the vertexer. Widths from an iterative &plusmn;2&sigma; Gaussian fit seeded
  on the half-maximum points; the supplied-macro panels use their own equivalent convention.
  Macros: <code>monitor/compare_track_fits.C</code>, plus the four
  <code>check_vertex_*_plots_trkvtxer_color.C</code> run unmodified.</p>
</footer>

</div>
"""
open('report_trkvtxer.html','w').write(HTML)
print('wrote', len(HTML))
