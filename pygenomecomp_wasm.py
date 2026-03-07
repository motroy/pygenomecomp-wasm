"""
PyGenomeComp WASM - Genome Comparison Tool
Self-contained Python module for Pyodide (browser WASM runtime).
Alignments are performed by MUMmer4/nucmer via biowasm (in JavaScript),
and the resulting delta files are parsed here for visualization.
No external Python dependencies required.
"""
import math
import json

# ============================================================
# FASTA Parser
# ============================================================

def parse_fasta(text):
    """Parse FASTA text. Returns list of {'name': str, 'seq': str}."""
    sequences = []
    current_name = None
    current_seq = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if current_name is not None:
                sequences.append({'name': current_name, 'seq': ''.join(current_seq).upper()})
            current_name = line[1:].split()[0]
            current_seq = []
        else:
            current_seq.append(line)
    if current_name is not None:
        sequences.append({'name': current_name, 'seq': ''.join(current_seq).upper()})
    return sequences


# ============================================================
# GFF3 Parser
# ============================================================

def parse_gff3(text, target_seq_id):
    """Parse GFF3 text for features on target_seq_id."""
    features = []
    for line in text.splitlines():
        if line.startswith('#') or not line.strip():
            continue
        parts = line.strip().split('\t')
        if len(parts) != 9:
            continue
        if parts[0] != target_seq_id:
            continue
        try:
            start = int(parts[3])
            end = int(parts[4])
        except ValueError:
            continue
        attr_str = parts[8]
        attrs = {}
        for attr in attr_str.split(';'):
            if '=' in attr:
                k, v = attr.split('=', 1)
                attrs[k.strip()] = v.strip()
        strand_map = {'.': None, '?': None, '+': 1, '-': -1}
        features.append({
            'seq_id': parts[0],
            'start': start,
            'end': end,
            'feature_type': parts[2],
            'strand': strand_map.get(parts[6]),
            'attributes': attrs,
            'id': attrs.get('ID') or attrs.get('Name'),
        })
    return features


# ============================================================
# lastz General Format Parser
# ============================================================

def parse_lastz_general(alignment_text, query_file_name, query_index,
                        min_identity=70.0, min_coverage=0.0, min_length=100):
    """
    Parse a lastz --format=general output into hit dicts compatible with generate_svg().

    Expected header (ignoring leading #):
    score name1 strand1 size1 zstart1 end1 name2 strand2 size2 zstart2 end2 identity idPct coverage covPct

    Returns list of hit dicts.
    """
    lines = alignment_text.splitlines()
    hits = []

    # We map columns by index based on header if present, else fallback to standard defaults
    header_map = {}

    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith('#'):
            # Parse header
            if 'score' in line.lower() and 'name1' in line.lower():
                parts = line.lstrip('#').strip().split()
                header_map = {col.lower(): idx for idx, col in enumerate(parts)}
            continue

        parts = line.split('\t')

        # Determine indices
        if header_map:
            idx_name1   = header_map.get('name1', 1)
            idx_strand1 = header_map.get('strand1', 2)
            idx_size1   = header_map.get('size1', 3)
            idx_zstart1 = header_map.get('zstart1', 4)
            idx_end1    = header_map.get('end1', 5)

            idx_name2   = header_map.get('name2', 6)
            idx_strand2 = header_map.get('strand2', 7)
            idx_size2   = header_map.get('size2', 8)
            idx_zstart2 = header_map.get('zstart2', 9)
            idx_end2    = header_map.get('end2', 10)

            idx_idpct   = header_map.get('idpct', 12)
            idx_covpct  = header_map.get('covpct', 14)
        else:
            # Fallback to defaults
            idx_name1, idx_strand1, idx_size1, idx_zstart1, idx_end1 = 1, 2, 3, 4, 5
            idx_name2, idx_strand2, idx_size2, idx_zstart2, idx_end2 = 6, 7, 8, 9, 10
            idx_idpct, idx_covpct = 12, 14

        if len(parts) <= max(idx_end1, idx_end2, idx_idpct, idx_covpct):
            continue

        try:
            ref_name   = parts[idx_name1]
            query_name = parts[idx_name2]

            sstart = int(parts[idx_zstart1]) + 1 # Convert 0-based zstart to 1-based
            send   = int(parts[idx_end1])

            qstart = int(parts[idx_zstart2]) + 1
            qend   = int(parts[idx_end2])

            if parts[idx_strand2] == '-':
                qstart, qend = qend, qstart

            id_pct_str = parts[idx_idpct].rstrip('%')
            pident = float(id_pct_str)

            cov_pct_str = parts[idx_covpct].rstrip('%')
            cov_pct = float(cov_pct_str)

            ref_len   = int(parts[idx_size1])
            query_len = int(parts[idx_size2])

        except (ValueError, IndexError):
            continue

        ref_aln_len = abs(send - sstart) + 1

        if ref_aln_len < min_length:
            continue

        if pident < min_identity:
            continue

        if cov_pct < min_coverage:
            continue

        hits.append({
            'qseqid':      query_name,
            'sseqid':      ref_name,
            'pident':      round(pident, 2),
            'coverage':    round(cov_pct, 2),
            'length':      ref_aln_len,
            'mismatch':    0, # Approximation, exact errors not typically needed for visualization
            'gapopen':     0,
            'qstart':      qstart,
            'qend':        qend,
            'sstart':      min(sstart, send),
            'send':        max(sstart, send),
            'evalue':      0.0,
            'bitscore':    round(ref_aln_len * pident / 100.0, 1),
            'qlen':        query_len,
            'slen':        ref_len,
            'query_index': query_index,
        })

    return hits


# ============================================================
# Insertion Site Detection
# ============================================================

def find_insertion_sites(hits, max_ref_gap=50):
    """Identifies insertion sites where query has sequence absent from reference."""
    if len(hits) < 2:
        return []
    sorted_hits = sorted(hits, key=lambda h: min(h['sstart'], h['send']))
    sites = []
    for i in range(len(sorted_hits) - 1):
        h1, h2 = sorted_hits[i], sorted_hits[i + 1]
        ref_end        = max(h1['sstart'], h1['send'])
        ref_start_next = min(h2['sstart'], h2['send'])
        ref_gap   = ref_start_next - ref_end - 1
        q_end        = max(h1['qstart'], h1['qend'])
        q_start_next = min(h2['qstart'], h2['qend'])
        query_gap = q_start_next - q_end - 1
        if ref_gap <= max_ref_gap and query_gap > 0 and h1['sseqid'] == h2['sseqid']:
            sites.append({
                'ref_position': ref_end + 0.5,
                'insertion_size': query_gap,
                'query_index': h1.get('query_index', 0),
            })
    return sites


# ============================================================
# SVG Plot Generator (pure Python, no svgwrite dependency)
# ============================================================

SVG_WIDTH, SVG_HEIGHT = 1000, 920
CENTER_X, CENTER_Y = 450, 470

ANNOTATION_RING_RADIUS_INNER = 130.0
ANNOTATION_RING_THICKNESS    = 18.0
ANNOTATION_RING_RADIUS_OUTER = ANNOTATION_RING_RADIUS_INNER + ANNOTATION_RING_THICKNESS
REFERENCE_RING_RADIUS        = ANNOTATION_RING_RADIUS_OUTER + 4.0
BLAST_HIT_BASE_RADIUS        = REFERENCE_RING_RADIUS + 4.0
BLAST_RING_THICKNESS         = 25.0
BLAST_RING_SPACING           = 2.0

RING_COLORS = [
    (227, 26,  28),   # Red
    (55,  126, 184),  # Blue
    (77,  175, 74),   # Green
    (152, 78,  163),  # Purple
    (255, 127, 0),    # Orange
    (166, 86,  40),   # Brown
    (247, 129, 191),  # Pink
    (153, 153, 153),  # Grey
]

ANNOTATION_COLORS = {
    'cds':  '#4682B4',
    'gene': '#3CB371',
    'trna': '#FF6347',
    'rrna': '#EE82EE',
}
DEFAULT_ANNOTATION_COLOR = '#808080'


def _esc(s):
    return str(s).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;').replace('"', '&quot;')

def _rgb_hex(r, g, b):
    return f'#{r:02x}{g:02x}{b:02x}'

def _ring_color(qi):
    return RING_COLORS[qi % len(RING_COLORS)]

def _hit_color(qi, identity):
    r, g, b = _ring_color(qi)
    identity = max(70.0, min(100.0, identity))
    frac = (identity - 70.0) / 30.0
    strength = 0.4 + 0.6 * frac
    return _rgb_hex(int(255-(255-r)*strength), int(255-(255-g)*strength), int(255-(255-b)*strength))

def _ring_bg_color(qi):
    r, g, b = _ring_color(qi)
    f = 0.10
    return _rgb_hex(int(255-(255-r)*f), int(255-(255-g)*f), int(255-(255-b)*f))

def _annotation_color(ftype):
    return ANNOTATION_COLORS.get(ftype.lower(), DEFAULT_ANNOTATION_COLOR)

def _arc_path(cx, cy, ri, ro, s, e):
    large = 1 if abs(e - s) > math.pi else 0
    xso = cx + ro * math.cos(s);  yso = cy + ro * math.sin(s)
    xeo = cx + ro * math.cos(e);  yeo = cy + ro * math.sin(e)
    xsi = cx + ri * math.cos(s);  ysi = cy + ri * math.sin(s)
    xei = cx + ri * math.cos(e);  yei = cy + ri * math.sin(e)
    return (f"M {xso:.2f} {yso:.2f} "
            f"A {ro:.2f} {ro:.2f} 0 {large} 1 {xeo:.2f} {yeo:.2f} "
            f"L {xei:.2f} {yei:.2f} "
            f"A {ri:.2f} {ri:.2f} 0 {large} 0 {xsi:.2f} {ysi:.2f} Z")

def _arc(cx, cy, ri, ro, s, e, fill, opacity=1.0, stroke='none', sw=0):
    # SVG cannot render an arc whose start and end points coincide (degenerate arc).
    # This happens for near-full-circle spans (e.g. a query hit covering the entire
    # reference), where er - sr approaches 2π.  Split such arcs into two halves,
    # mirroring the same technique used by _full_ring().
    if abs(e - s) > 2 * math.pi - 0.001:
        mid = (s + e) / 2
        return (_arc(cx, cy, ri, ro, s, mid, fill, opacity, stroke, sw) +
                _arc(cx, cy, ri, ro, mid, e, fill, opacity, stroke, sw))
    p = _arc_path(cx, cy, ri, ro, s, e)
    return f'<path d="{p}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}" opacity="{opacity:.2f}"/>\n'

def _full_ring(cx, cy, ri, ro, fill, opacity=1.0):
    return (_arc(cx, cy, ri, ro, -math.pi/2,     math.pi/2,     fill, opacity) +
            _arc(cx, cy, ri, ro,  math.pi/2, 3 * math.pi/2,     fill, opacity))


def generate_svg(blast_hits, annotations, reference_length, query_names,
                 reference_display_name, show_gene_names=False,
                 insertion_sites=None):
    """Generate a BRIG-style circular SVG and return it as a string."""
    if reference_length == 0:
        raise ValueError("Reference length cannot be zero.")

    n_queries       = len(query_names)
    insertion_sites = insertion_sites or []

    outermost_r = (BLAST_HIT_BASE_RADIUS
                   + n_queries * BLAST_RING_THICKNESS
                   + max(0, n_queries - 1) * BLAST_RING_SPACING) if n_queries > 0 else REFERENCE_RING_RADIUS + 10
    tick_inner = outermost_r + 6
    tick_outer = tick_inner + 10
    label_r    = tick_outer + 18

    out = []

    # --- SVG header + defs ---
    out.append(f'<svg xmlns="http://www.w3.org/2000/svg" '
               f'xmlns:xlink="http://www.w3.org/1999/xlink" '
               f'width="{SVG_WIDTH}" height="{SVG_HEIGHT}" '
               f'viewBox="0 0 {SVG_WIDTH} {SVG_HEIGHT}">\n')
    out.append('<style>text { font-family: Arial, Helvetica, sans-serif; }</style>\n')
    out.append('<defs>\n')
    for i in range(min(n_queries, len(RING_COLORS))):
        gid = f'ig{i}'
        stops = ''.join(
            f'  <stop offset="{frac}" stop-color="{_hit_color(i, ident)}"/>\n'
            for frac, ident in ((0, 70), (0.33, 80), (0.66, 90), (1, 100))
        )
        out.append(f'<linearGradient id="{gid}" x1="0" y1="0" x2="1" y2="0">\n{stops}</linearGradient>\n')
    out.append('</defs>\n')

    # --- Background ---
    out.append(f'<rect width="{SVG_WIDTH}" height="{SVG_HEIGHT}" fill="white"/>\n')
    out.append(f'<circle cx="{CENTER_X}" cy="{CENTER_Y}" r="{ANNOTATION_RING_RADIUS_INNER-2:.1f}" fill="white"/>\n')

    # --- Center labels ---
    disp = _esc(reference_display_name)
    out.append(f'<text x="{CENTER_X}" y="{CENTER_Y-12}" text-anchor="middle" '
               f'font-size="13" font-weight="bold" fill="#333">{disp}</text>\n')
    out.append(f'<text x="{CENTER_X}" y="{CENTER_Y+10}" text-anchor="middle" '
               f'font-size="11" fill="#666">{reference_length:,} bp</text>\n')

    # --- Annotation ring ---
    if annotations:
        out.append(_full_ring(CENTER_X, CENTER_Y,
                              ANNOTATION_RING_RADIUS_INNER, ANNOTATION_RING_RADIUS_OUTER,
                              '#f0f0f0'))
        for feat in annotations:
            ms = feat['start'] - 1
            me = feat['end']
            sr = (ms / reference_length) * 2 * math.pi - math.pi / 2
            er = (me / reference_length) * 2 * math.pi - math.pi / 2
            if abs(er - sr) < 1e-6:
                continue
            out.append(_arc(CENTER_X, CENTER_Y,
                            ANNOTATION_RING_RADIUS_INNER, ANNOTATION_RING_RADIUS_OUTER,
                            sr, er, _annotation_color(feat['feature_type']), 0.9))
            if show_gene_names and feat['feature_type'].lower() == 'gene':
                mid = (sr + er) / 2
                gname = (feat['attributes'].get('Name') or feat.get('id') or
                         feat['attributes'].get('gene') or '')
                if gname:
                    lr2 = outermost_r + 40
                    x1 = CENTER_X + (ANNOTATION_RING_RADIUS_OUTER + 2) * math.cos(mid)
                    y1 = CENTER_Y + (ANNOTATION_RING_RADIUS_OUTER + 2) * math.sin(mid)
                    x2 = CENTER_X + (lr2 - 4) * math.cos(mid)
                    y2 = CENTER_Y + (lr2 - 4) * math.sin(mid)
                    out.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
                               f'stroke="#aaa" stroke-width="0.7"/>\n')
                    gx = CENTER_X + (outermost_r + 45) * math.cos(mid)
                    gy = CENTER_Y + (outermost_r + 45) * math.sin(mid)
                    deg = math.degrees(mid)
                    rot = deg + 90 if -90 <= deg <= 90 else deg - 90
                    out.append(f'<text x="{gx:.1f}" y="{gy:.1f}" font-size="9" text-anchor="middle" '
                               f'fill="#222" transform="rotate({rot:.1f},{gx:.2f},{gy:.2f})">'
                               f'{_esc(gname)}</text>\n')

    # --- Reference backbone ---
    out.append(f'<circle cx="{CENTER_X}" cy="{CENTER_Y}" r="{REFERENCE_RING_RADIUS:.1f}" '
               f'fill="none" stroke="#333" stroke-width="2"/>\n')

    # --- Query ring backgrounds ---
    for i in range(n_queries):
        ir = BLAST_HIT_BASE_RADIUS + i * (BLAST_RING_THICKNESS + BLAST_RING_SPACING)
        or_ = ir + BLAST_RING_THICKNESS
        out.append(_full_ring(CENTER_X, CENTER_Y, ir, or_, _ring_bg_color(i)))
        out.append(f'<circle cx="{CENTER_X}" cy="{CENTER_Y}" r="{ir:.1f}" fill="none" stroke="#ddd" stroke-width="0.5"/>\n')
        out.append(f'<circle cx="{CENTER_X}" cy="{CENTER_Y}" r="{or_:.1f}" fill="none" stroke="#ddd" stroke-width="0.5"/>\n')

    # --- Alignment hit segments ---
    for hit in blast_hits:
        qi  = hit['query_index']
        ir  = BLAST_HIT_BASE_RADIUS + qi * (BLAST_RING_THICKNESS + BLAST_RING_SPACING)
        or_ = ir + BLAST_RING_THICKNESS
        ss  = min(hit['sstart'], hit['send'])
        se  = max(hit['sstart'], hit['send'])
        sr  = ((ss - 1) / reference_length) * 2 * math.pi - math.pi / 2
        er  = ((se - 1) / reference_length) * 2 * math.pi - math.pi / 2
        if abs(er - sr) < 1e-6:
            continue
        out.append(_arc(CENTER_X, CENTER_Y, ir, or_, sr, er,
                        _hit_color(qi, hit['pident'])))

    # --- Insertion site markers ---
    INS_COLOR = '#1a1a1a'
    INS_MIN_HW = 0.012
    for site in insertion_sites:
        qi  = site['query_index']
        ir  = BLAST_HIT_BASE_RADIUS + qi * (BLAST_RING_THICKNESS + BLAST_RING_SPACING)
        or_ = ir + BLAST_RING_THICKNESS
        mid = (site['ref_position'] / reference_length) * 2 * math.pi - math.pi / 2
        hw  = max(INS_MIN_HW, (1.0 / reference_length) * 2 * math.pi)
        out.append(_arc(CENTER_X, CENTER_Y, ir, or_, mid-hw, mid+hw, INS_COLOR, 0.85))
        x1 = CENTER_X + tick_inner * math.cos(mid)
        y1 = CENTER_Y + tick_inner * math.sin(mid)
        x2 = CENTER_X + tick_outer * math.cos(mid)
        y2 = CENTER_Y + tick_outer * math.sin(mid)
        out.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
                   f'stroke="{INS_COLOR}" stroke-width="1.5"/>\n')
        ins = site['insertion_size']
        slabel = (f'+{ins/1e6:.1f} Mb' if ins >= 1_000_000 else
                  f'+{ins/1e3:.1f} kb' if ins >= 1_000 else f'+{ins} bp')
        xlbl = CENTER_X + label_r * math.cos(mid)
        ylbl = CENTER_Y + label_r * math.sin(mid)
        deg  = math.degrees(mid)
        rot  = deg + 90 if -90 <= deg <= 90 else deg - 90
        out.append(f'<text x="{xlbl:.1f}" y="{ylbl:.1f}" font-size="9" text-anchor="middle" '
                   f'fill="{INS_COLOR}" transform="rotate({rot:.1f},{xlbl:.2f},{ylbl:.2f})">'
                   f'{_esc(slabel)}</text>\n')

    # --- Outermost ring border ---
    if n_queries > 0:
        fr = (BLAST_HIT_BASE_RADIUS
              + (n_queries - 1) * (BLAST_RING_THICKNESS + BLAST_RING_SPACING)
              + BLAST_RING_THICKNESS)
        out.append(f'<circle cx="{CENTER_X}" cy="{CENTER_Y}" r="{fr:.1f}" '
                   f'fill="none" stroke="#999" stroke-width="1"/>\n')

    # --- Scale ticks ---
    for i in range(12):
        angle  = (i / 12) * 2 * math.pi - math.pi / 2
        pos_bp = (i / 12) * reference_length
        if   reference_length > 2_000_000: lbl = f'{pos_bp/1e6:.1f} Mb'
        elif reference_length > 2000:      lbl = f'{pos_bp/1e3:.0f} kb'
        else:                              lbl = f'{int(pos_bp)} bp'
        x1 = CENTER_X + tick_inner * math.cos(angle)
        y1 = CENTER_Y + tick_inner * math.sin(angle)
        x2 = CENTER_X + tick_outer * math.cos(angle)
        y2 = CENTER_Y + tick_outer * math.sin(angle)
        xt = CENTER_X + label_r    * math.cos(angle)
        yt = CENTER_Y + label_r    * math.sin(angle)
        out.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
                   f'stroke="#555" stroke-width="1"/>\n')
        out.append(f'<text x="{xt:.1f}" y="{yt:.1f}" text-anchor="middle" '
                   f'font-size="10" fill="#555">{_esc(lbl)}</text>\n')

    # --- Legend ---
    lx  = SVG_WIDTH - 155
    ly  = 50
    bs  = 14
    ih  = 22
    out.append(f'<text x="{lx}" y="{ly}" font-weight="bold" font-size="14" fill="#333">Legend</text>\n')
    ly += ih + 5

    if n_queries > 0:
        out.append(f'<text x="{lx}" y="{ly}" font-weight="bold" font-size="11" fill="#555">Query Rings</text>\n')
        ly += ih - 2
        for i, name in enumerate(query_names):
            r, g, b = _ring_color(i)
            c = _rgb_hex(r, g, b)
            out.append(f'<rect x="{lx}" y="{ly}" width="{bs}" height="{bs}" fill="{c}" rx="2"/>\n')
            disp_name = name[:18] + '…' if len(name) > 19 else name
            out.append(f'<text x="{lx+bs+7}" y="{ly+bs-2}" font-size="10" fill="#333">'
                       f'{_esc(disp_name)}</text>\n')
            ly += ih
        ly += 6

        out.append(f'<text x="{lx}" y="{ly}" font-weight="bold" font-size="11" fill="#555">% Identity</text>\n')
        ly += ih - 2
        gw, gh = 100, 12
        for i in range(min(n_queries, len(RING_COLORS))):
            out.append(f'<rect x="{lx}" y="{ly}" width="{gw}" height="{gh}" '
                       f'fill="url(#ig{i})" rx="2"/>\n')
            ly += gh + 4
        out.append(f'<text x="{lx}" y="{ly}" font-size="9" fill="#666">70%</text>\n')
        out.append(f'<text x="{lx+gw-25}" y="{ly}" font-size="9" fill="#666">100%</text>\n')
        ly += ih

    if insertion_sites:
        out.append(f'<text x="{lx}" y="{ly}" font-weight="bold" font-size="11" fill="#555">Insertions</text>\n')
        ly += ih - 2
        out.append(f'<rect x="{lx}" y="{ly}" width="{bs}" height="{bs}" fill="#1a1a1a" rx="2"/>\n')
        out.append(f'<text x="{lx+bs+7}" y="{ly+bs-2}" font-size="10" fill="#333">Insertion site</text>\n')
        ly += ih + 4

    if annotations:
        out.append(f'<text x="{lx}" y="{ly}" font-weight="bold" font-size="11" fill="#555">Annotations</text>\n')
        ly += ih - 2
        for label, ftype in (('CDS','cds'),('Gene','gene'),('tRNA','trna'),('rRNA','rrna')):
            c = _annotation_color(ftype)
            out.append(f'<rect x="{lx}" y="{ly}" width="{bs}" height="{bs}" fill="{c}" rx="2"/>\n')
            out.append(f'<text x="{lx+bs+7}" y="{ly+bs-2}" font-size="10" fill="#333">{label}</text>\n')
            ly += ih

    out.append('</svg>')
    return ''.join(out)


# ============================================================
# Main Pipeline
# ============================================================

def run_comparison(reference_fasta_text, query_alignment_texts, query_file_names,
                   annotation_text=None,
                   min_identity=70.0, min_coverage=0.0, min_length=100,
                   show_gene_names=False,
                   progress_callback=None):
    """
    Run genome comparison pipeline using pre-computed lastz alignments.

    Alignments are performed in JavaScript via biowasm/Aioli before this
    function is called. The resulting lastz --format=general outputs are parsed here.

    Args:
        reference_fasta_text   : str       – reference genome in FASTA format
        query_alignment_texts  : list[str] – lastz outputs (one per query)
        query_file_names       : list[str] – display names for query files
        annotation_text        : str|None  – optional GFF3 annotation text
        min_identity           : float     – minimum % identity (default 70)
        min_coverage           : float     – minimum % coverage (default 0)
        min_length             : int       – minimum alignment length bp (default 100)
        show_gene_names        : bool      – annotate gene names on the plot
        progress_callback      : callable(step, total, message) | None

    Returns:
        str – SVG plot as a string
    """
    def _prog(step, total, msg):
        if progress_callback:
            progress_callback(step, total, msg)

    n_queries   = len(query_alignment_texts)
    total_steps = 3 + n_queries

    # 1. Parse reference (for name and length only — alignment done externally)
    _prog(0, total_steps, 'Parsing reference genome…')
    ref_records = parse_fasta(reference_fasta_text)
    if not ref_records:
        raise ValueError("No sequences found in reference FASTA.")
    ref_record = ref_records[0]
    ref_name   = ref_record['name']
    ref_len    = len(ref_record['seq'])
    _prog(1, total_steps, f'Reference: {ref_name} ({ref_len:,} bp)')

    # 2. Parse lastz alignment files into hit dicts
    all_hits       = []
    all_insertions = []

    for i, (alignment_text, qname) in enumerate(zip(query_alignment_texts, query_file_names)):
        _prog(2 + i, total_steps, f'Processing alignments for {qname} ({i+1}/{n_queries})…')
        query_hits = parse_lastz_general(
            alignment_text, qname, i, min_identity, min_coverage, min_length
        )
        all_hits.extend(query_hits)
        all_insertions.extend(find_insertion_sites(query_hits))

    # 3. Parse annotations
    _prog(2 + n_queries, total_steps, 'Parsing annotations…')
    annotations = None
    if annotation_text and annotation_text.strip():
        try:
            annotations = parse_gff3(annotation_text, ref_name)
        except Exception as e:
            print(f'Warning: annotation parse error: {e}')

    # 4. Generate plot
    _prog(3 + n_queries - 1, total_steps, 'Generating circular plot…')
    svg = generate_svg(
        blast_hits=all_hits,
        annotations=annotations,
        reference_length=ref_len,
        query_names=query_file_names,
        reference_display_name=ref_name,
        show_gene_names=show_gene_names,
        insertion_sites=all_insertions,
    )
    _prog(total_steps, total_steps, 'Done!')
    return svg
