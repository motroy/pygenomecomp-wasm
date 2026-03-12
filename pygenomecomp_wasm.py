"""
PyGenomeComp WASM - Genome Comparison Tool
Self-contained Python module for Pyodide (browser WASM runtime).
Alignments are performed by lastz via biowasm (in JavaScript),
and the resulting general format outputs are parsed here for visualization.
No external Python dependencies required.
"""
import math
import json
import re

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

SVG_WIDTH, SVG_HEIGHT = 1400, 1100
CENTER_X, CENTER_Y = 600, 550

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

# Linear alignment view colours
LINEAR_FWD_COLOR      = '#f4a460'   # sandy orange – forward alignment
LINEAR_INV_COLOR      = '#b8a4c9'   # muted purple – inverted alignment
LINEAR_BACKBONE_COLOR = '#c8d4e8'   # light steel-blue – sequence backbone


def _esc(s):
    return str(s).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;').replace('"', '&quot;')

def _rgb_hex(r, g, b):
    return f'#{r:02x}{g:02x}{b:02x}'

def _ring_color(qi):
    return RING_COLORS[qi % len(RING_COLORS)]

def _hit_color(qi, identity, min_identity=70.0):
    r, g, b = _ring_color(qi)
    identity = max(min_identity, min(100.0, identity))
    span = 100.0 - min_identity
    frac = (identity - min_identity) / span if span > 0 else 1.0
    strength = 0.4 + 0.6 * frac
    return _rgb_hex(int(255-(255-r)*strength), int(255-(255-g)*strength), int(255-(255-b)*strength))

def _linear_hit_color(qi, identity, inverted, min_identity=70.0):
    """Color linear/alignment ribbon using the query's ring colour, shaded by identity.

    Mirrors the per-genome colour scheme used in the circular view so colours are
    consistent across all three view modes.  Forward alignments use the full ring
    colour; inverted alignments use a lighter shade of the same colour so strand
    direction is still visible.
    """
    r, g, b = _ring_color(qi)
    identity = max(min_identity, min(100.0, identity))
    span = 100.0 - min_identity
    frac = (identity - min_identity) / span if span > 0 else 1.0
    if inverted:
        strength = 0.15 + 0.35 * frac   # 0.15 (very pale) → 0.50 (medium)
    else:
        strength = 0.30 + 0.70 * frac   # 0.30 (pale) → 1.00 (full colour)
    return _rgb_hex(
        int(255 - (255 - r) * strength),
        int(255 - (255 - g) * strength),
        int(255 - (255 - b) * strength),
    )

def _ring_bg_color(qi):
    r, g, b = _ring_color(qi)
    f = 0.10
    return _rgb_hex(int(255-(255-r)*f), int(255-(255-g)*f), int(255-(255-b)*f))

def _annotation_color(ftype):
    return ANNOTATION_COLORS.get(ftype.lower(), DEFAULT_ANNOTATION_COLOR)

def _is_resistance_gene(feat):
    """Check if an annotation feature is a resistance gene based on its product or attributes."""
    if feat.get('feature_type', '').lower() not in ('cds', 'gene'):
        return False
    attrs = feat.get('attributes', {})
    product = attrs.get('product', '').lower()
    name = attrs.get('Name', '').lower()
    gene_attr = attrs.get('gene', '').lower()
    note = attrs.get('note', '').lower()

    # Common resistance keywords or patterns
    # These often include 'bla', 'ndm', 'kpc', 'oxa', 'tet', 'sul', 'mcr', etc.
    res_keywords = [
        'resistance', 'beta-lactamase', 'carbapenemase', 'efflux',
        'antibiotic', 'antimicrobial', 'aminoglycoside', 'tetracycline',
        'macrolide', 'chloramphenicol', 'sulfonamide', 'trimethoprim',
        'quinolone', 'streptomycin', 'spectinomycin', 'erythromycin'
    ]

    # Also check if the gene name looks like a common AMR gene (e.g. blaNDM, tetA, sul1)
    short_gene_pattern = re.compile(r'^(bla|tet|sul|aad|aph|aac|erm|mef|msr|lnu|vat|vga|dfr|qnr|mcr|cat|flo|cml|ere)[a-zA-Z0-9\-\(\)\']*$', re.IGNORECASE)

    text_to_search = f"{product} {name} {gene_attr} {note}"

    if any(kw in text_to_search for kw in res_keywords):
        return True

    if name and short_gene_pattern.match(name):
        return True

    if gene_attr and short_gene_pattern.match(gene_attr):
        return True

    # Check ontology terms for AMR-related GO terms if present
    if 'GO:0046677' in attrs.get('Ontology_term', ''): # response to antibiotic
        return True

    return False

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
                 insertion_sites=None, show_insertions=True, min_identity=70.0,
                 show_res_genes=False):
    """Generate a BRIG-style circular SVG and return it as a string."""
    if reference_length == 0:
        raise ValueError("Reference length cannot be zero.")

    n_queries       = len(query_names)
    insertion_sites = insertion_sites or []

    # Dynamically scale ring thickness/spacing so all rings fit within the SVG canvas.
    # Available radius = distance from center to edge minus space for ticks and labels.
    if n_queries > 0:
        _tick_label_margin = 34  # pixels reserved outside rings for ticks + labels
        _available_r = min(CENTER_X, CENTER_Y) - BLAST_HIT_BASE_RADIUS - _tick_label_margin
        _total_default = (n_queries * BLAST_RING_THICKNESS
                          + max(0, n_queries - 1) * BLAST_RING_SPACING)
        if _total_default > _available_r:
            _scale = _available_r / _total_default
            blast_ring_thickness = max(2.0, BLAST_RING_THICKNESS * _scale)
            blast_ring_spacing   = max(0.3, BLAST_RING_SPACING * _scale)
        else:
            blast_ring_thickness = BLAST_RING_THICKNESS
            blast_ring_spacing   = BLAST_RING_SPACING
    else:
        blast_ring_thickness = BLAST_RING_THICKNESS
        blast_ring_spacing   = BLAST_RING_SPACING

    outermost_r = (BLAST_HIT_BASE_RADIUS
                   + n_queries * blast_ring_thickness
                   + max(0, n_queries - 1) * blast_ring_spacing) if n_queries > 0 else REFERENCE_RING_RADIUS + 10
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
    _grad_fracs = (0, 0.33, 0.66, 1)
    _span = 100.0 - min_identity
    for i in range(min(n_queries, len(RING_COLORS))):
        gid = f'ig{i}'
        stops = ''.join(
            f'  <stop offset="{frac}" stop-color="{_hit_color(i, min_identity + frac * _span, min_identity)}"/>\n'
            for frac in _grad_fracs
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
    _gene_label_queue = []   # (orig_angle, name) collected for de-overlapping
    if annotations:
        out.append(_full_ring(CENTER_X, CENTER_Y,
                              ANNOTATION_RING_RADIUS_INNER, ANNOTATION_RING_RADIUS_OUTER,
                              '#f0f0f0'))
        for feat in annotations:
            ftype = feat['feature_type'].lower()
            if ftype not in ('cds', 'gene', 'trna', 'rrna'):
                continue
            if show_res_genes and not _is_resistance_gene(feat):
                continue
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
                    _gene_label_queue.append((mid, gname))

    # --- Gene name labels (radial tiers, horizontal text) ---
    if _gene_label_queue:
        _GL_FONT   = 9
        _GL_CHARW  = _GL_FONT * 0.58   # estimated px per character
        _GL_FONT_H = _GL_FONT * 1.3    # label line height
        _GL_PAD    = 3                  # bbox padding
        _GL_BASE_R = label_r + 18      # start just outside scale-tick labels
        _GL_TIER_S = 14                # radial px per tier

        def _gl_bbox(angle, name, tier):
            r  = _GL_BASE_R + tier * _GL_TIER_S
            bx = CENTER_X + r * math.cos(angle)
            by = CENTER_Y + r * math.sin(angle)
            w  = len(name) * _GL_CHARW
            h  = _GL_FONT_H
            if math.cos(angle) >= 0:   # right half: text grows rightward
                return (bx - _GL_PAD, bx + w + _GL_PAD,
                        by - h / 2 - _GL_PAD, by + h / 2 + _GL_PAD)
            else:                       # left half: text grows leftward
                return (bx - w - _GL_PAD, bx + _GL_PAD,
                        by - h / 2 - _GL_PAD, by + h / 2 + _GL_PAD)

        def _gl_overlaps(a, b):
            return a[0] < b[1] and a[1] > b[0] and a[2] < b[3] and a[3] > b[2]

        _gl_sorted = sorted(_gene_label_queue, key=lambda t: t[0])
        _gl_placed = []   # list of (angle, name, tier, bbox)

        for (_gl_angle, _gl_name) in _gl_sorted:
            _gl_t = 0
            while _gl_t <= 40:
                _gl_bb = _gl_bbox(_gl_angle, _gl_name, _gl_t)
                if all(not _gl_overlaps(_gl_bb, _p[3]) for _p in _gl_placed):
                    break
                _gl_t += 1
            _gl_placed.append((_gl_angle, _gl_name, _gl_t,
                                _gl_bbox(_gl_angle, _gl_name, _gl_t)))

        for (_gl_angle, _gl_name, _gl_tier, _) in _gl_placed:
            _gl_r  = _GL_BASE_R + _gl_tier * _GL_TIER_S
            _gl_lx = CENTER_X + _gl_r * math.cos(_gl_angle)
            _gl_ly = CENTER_Y + _gl_r * math.sin(_gl_angle)
            # Leader line from outermost ring edge to label anchor
            _gl_sx = CENTER_X + tick_inner * math.cos(_gl_angle)
            _gl_sy = CENTER_Y + tick_inner * math.sin(_gl_angle)
            out.append(
                f'<line x1="{_gl_sx:.1f}" y1="{_gl_sy:.1f}" '
                f'x2="{_gl_lx:.1f}" y2="{_gl_ly:.1f}" '
                f'stroke="#aaa" stroke-width="0.7"/>\n'
            )
            _gl_anchor = 'start' if math.cos(_gl_angle) >= 0 else 'end'
            out.append(
                f'<text x="{_gl_lx:.1f}" y="{_gl_ly:.1f}" font-size="{_GL_FONT}" '
                f'text-anchor="{_gl_anchor}" dominant-baseline="middle" fill="#222">'
                f'{_esc(_gl_name)}</text>\n'
            )

    # --- Reference backbone ---
    out.append(f'<circle cx="{CENTER_X}" cy="{CENTER_Y}" r="{REFERENCE_RING_RADIUS:.1f}" '
               f'fill="none" stroke="#333" stroke-width="2"/>\n')

    # --- Query ring backgrounds ---
    for i in range(n_queries):
        ir = BLAST_HIT_BASE_RADIUS + i * (blast_ring_thickness + blast_ring_spacing)
        or_ = ir + blast_ring_thickness
        out.append(_full_ring(CENTER_X, CENTER_Y, ir, or_, _ring_bg_color(i)))
        out.append(f'<circle cx="{CENTER_X}" cy="{CENTER_Y}" r="{ir:.1f}" fill="none" stroke="#ddd" stroke-width="0.5"/>\n')
        out.append(f'<circle cx="{CENTER_X}" cy="{CENTER_Y}" r="{or_:.1f}" fill="none" stroke="#ddd" stroke-width="0.5"/>\n')

    # --- Alignment hit segments ---
    for hit in blast_hits:
        qi  = hit['query_index']
        ir  = BLAST_HIT_BASE_RADIUS + qi * (blast_ring_thickness + blast_ring_spacing)
        or_ = ir + blast_ring_thickness
        ss  = min(hit['sstart'], hit['send'])
        se  = max(hit['sstart'], hit['send'])
        sr  = ((ss - 1) / reference_length) * 2 * math.pi - math.pi / 2
        er  = ((se - 1) / reference_length) * 2 * math.pi - math.pi / 2
        if abs(er - sr) < 1e-6:
            continue
        out.append(_arc(CENTER_X, CENTER_Y, ir, or_, sr, er,
                        _hit_color(qi, hit['pident'], min_identity)))

    # --- Insertion site markers ---
    INS_COLOR  = '#1a1a1a'
    INS_MIN_HW = 0.012
    if show_insertions:
        for site in insertion_sites:
            qi  = site['query_index']
            ir  = BLAST_HIT_BASE_RADIUS + qi * (blast_ring_thickness + blast_ring_spacing)
            or_ = ir + blast_ring_thickness
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
              + (n_queries - 1) * (blast_ring_thickness + blast_ring_spacing)
              + blast_ring_thickness)
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
        out.append(f'<text x="{lx}" y="{ly}" font-size="9" fill="#666">{min_identity:.0f}%</text>\n')
        out.append(f'<text x="{lx+gw-25}" y="{ly}" font-size="9" fill="#666">100%</text>\n')
        ly += ih

    if show_insertions and insertion_sites:
        out.append(f'<text x="{lx}" y="{ly}" font-weight="bold" font-size="11" fill="#555">Insertions</text>\n')
        ly += ih - 2
        out.append(f'<rect x="{lx}" y="{ly}" width="{bs}" height="{bs}" fill="#1a1a1a" rx="2"/>\n')
        out.append(f'<text x="{lx+bs+7}" y="{ly+bs-2}" font-size="10" fill="#333">Insertion site</text>\n')
        ly += ih + 4

    if annotations:
        out.append(f'<text x="{lx}" y="{ly}" font-weight="bold" font-size="11" fill="#555">Annotations</text>\n')
        ly += ih - 2

        seen_types = set()
        for feat in annotations:
            ft = feat['feature_type'].lower()
            if ft in ('cds', 'gene', 'trna', 'rrna') and ft not in seen_types:
                if show_res_genes and not _is_resistance_gene(feat):
                    continue
                seen_types.add(ft)

        for label, ftype in (('CDS','cds'),('Gene','gene'),('tRNA','trna'),('rRNA','rrna')):
            if ftype in seen_types:
                c = _annotation_color(ftype)
                out.append(f'<rect x="{lx}" y="{ly}" width="{bs}" height="{bs}" fill="{c}" rx="2"/>\n')
                out.append(f'<text x="{lx+bs+7}" y="{ly+bs-2}" font-size="10" fill="#333">{label}</text>\n')
                ly += ih

    out.append('</svg>')
    return ''.join(out)


# ============================================================
# Linear Alignment SVG Generator
# ============================================================

def generate_linear_svg(blast_hits, annotations, reference_length, reference_name,
                        query_names, show_gene_names=False, min_identity=70.0,
                        insertion_sites=None, show_insertions=True, show_res_genes=False):
    """Generate a linear (Mauve-style) alignment view SVG and return it as a string.

    Sequences are stacked as horizontal tracks. Alignment ribbons connect aligned
    regions between the reference (top) and each query (below), coloured by each
    query's ring colour (matching the circular view): forward = full colour,
    inverted = lighter shade of the same colour.
    """
    insertion_sites = insertion_sites or []
    if reference_length == 0:
        raise ValueError("Reference length cannot be zero.")

    n_queries = len(query_names)
    n_tracks  = n_queries + 1  # reference + one per query

    # Compute a display length for each query from the max qlen seen in hits
    query_lengths = {}
    for hit in blast_hits:
        qi = hit['query_index']
        ql = hit['qlen']
        if qi not in query_lengths or ql > query_lengths[qi]:
            query_lengths[qi] = ql
    for i in range(n_queries):
        if i not in query_lengths:
            query_lengths[i] = reference_length  # fallback

    # Shared scale: all genomes use the same pixels-per-bp ratio so lengths are comparable
    max_len = max(reference_length, max(query_lengths.values(), default=reference_length))

    # ── Layout constants ──────────────────────────────────────────────
    W           = 1200
    LEFT        = 150   # reserved for sequence labels
    RIGHT       = 190   # wide right margin so the legend never overlaps content
    BOTTOM      = 70
    CONTENT_W   = W - LEFT - RIGHT
    TRACK_H     = 12    # backbone bar height (px)
    TRACK_SPACE = 90    # center-to-center vertical spacing (px)
    GENE_H      = 14    # gene-arrow height (px)
    RIBBON_OP   = 0.55  # ribbon opacity

    # ── Pre-compute gene-label tiers (horizontal stacking) ────────────
    _LBL_FONT  = 9
    _LBL_CHARW = _LBL_FONT * 0.58   # estimated px per character
    _LBL_HPAD  = 8                   # min horizontal gap between labels in same tier
    _TIER_H    = 14                  # vertical px per tier

    _lbl_predata = []   # (x_center, name)
    if show_gene_names and annotations:
        for _f in annotations:
            _ft = _f['feature_type'].lower()
            if _ft not in ('cds', 'gene'):
                continue
            if show_res_genes and not _is_resistance_gene(_f):
                continue
            _ax1 = LEFT + (_f['start'] - 1) / reference_length * CONTENT_W
            _ax2 = LEFT + _f['end']          / reference_length * CONTENT_W
            if _ax2 - _ax1 < 0.5:
                continue
            _gn = (_f['attributes'].get('Name') or _f.get('id') or
                   _f['attributes'].get('gene') or '')
            if _gn:
                _lbl_predata.append(((_ax1 + _ax2) / 2, _gn))
    _lbl_predata.sort(key=lambda t: t[0])

    _tier_right = {}   # tier_index → rightmost used x
    _lbl_tiers  = []   # (x_center, name, tier)
    for (_lx, _ln) in _lbl_predata:
        _hw = len(_ln) * _LBL_CHARW / 2
        _t  = 0
        while _t in _tier_right and _lx - _hw < _tier_right[_t] + _LBL_HPAD:
            _t += 1
        _tier_right[_t] = _lx + _hw
        _lbl_tiers.append((_lx, _ln, _t))

    _n_tiers = max((_e[2] for _e in _lbl_tiers), default=-1) + 1
    TOP      = 20 + _n_tiers * _TIER_H + 12   # base margin + label rows + gap

    H = TOP + n_tracks * TRACK_SPACE + BOTTOM

    # ── Coordinate helpers ────────────────────────────────────────────
    def ref_x(pos):
        """Reference position (1-based bp) → SVG x (shared scale)."""
        return LEFT + (pos - 1) / max_len * CONTENT_W

    def qx(pos, qi):
        """Query position (1-based bp) → SVG x (same shared scale)."""
        return LEFT + (pos - 1) / max_len * CONTENT_W

    def track_cy(i):
        """Centre-y of track i (0 = reference)."""
        return TOP + i * TRACK_SPACE

    out = []

    # ── SVG header ────────────────────────────────────────────────────
    out.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{W}" height="{H}" viewBox="0 0 {W} {H}">\n'
    )
    out.append('<style>text { font-family: Arial, Helvetica, sans-serif; }</style>\n')

    # ── Gradient defs for identity legend (one pair per query) ───────
    _grad_fracs = (0, 0.33, 0.66, 1)
    _span = 100.0 - min_identity
    out.append('<defs>\n')
    for _qi in range(n_queries):
        _fwd = ''.join(
            f'  <stop offset="{f}" stop-color="{_linear_hit_color(_qi, min_identity + f * _span, False, min_identity)}"/>\n'
            for f in _grad_fracs
        )
        _inv = ''.join(
            f'  <stop offset="{f}" stop-color="{_linear_hit_color(_qi, min_identity + f * _span, True, min_identity)}"/>\n'
            for f in _grad_fracs
        )
        out.append(f'<linearGradient id="lin_fwd_grad_{_qi}" x1="0" y1="0" x2="1" y2="0">\n{_fwd}</linearGradient>\n')
        out.append(f'<linearGradient id="lin_inv_grad_{_qi}" x1="0" y1="0" x2="1" y2="0">\n{_inv}</linearGradient>\n')
    out.append('</defs>\n')

    out.append(f'<rect width="{W}" height="{H}" fill="white"/>\n')

    # ── Alignment ribbons (drawn first, behind track backbones) ───────
    for hit in blast_hits:
        qi  = hit['query_index']
        rx1 = ref_x(hit['sstart'])
        rx2 = ref_x(hit['send'])
        ry  = track_cy(0) + TRACK_H / 2        # bottom edge of reference track

        qs  = min(hit['qstart'], hit['qend'])
        qe  = max(hit['qstart'], hit['qend'])
        qx1 = qx(qs, qi)
        qx2 = qx(qe, qi)
        qy  = track_cy(qi + 1) - TRACK_H / 2  # top edge of query track

        # Inverted when qstart > qend (coords were swapped by the parser for '-' strand)
        inverted = hit['qstart'] > hit['qend']
        fill     = _linear_hit_color(qi, hit['pident'], inverted, min_identity)

        # Skip ribbons too narrow to see
        if abs(rx2 - rx1) < 0.5 and abs(qx2 - qx1) < 0.5:
            continue

        mid_y = (ry + qy) / 2
        out.append(
            f'<path d="'
            f'M {rx1:.1f},{ry:.1f} L {rx2:.1f},{ry:.1f} '
            f'C {rx2:.1f},{mid_y:.1f} {qx2:.1f},{mid_y:.1f} {qx2:.1f},{qy:.1f} '
            f'L {qx1:.1f},{qy:.1f} '
            f'C {qx1:.1f},{mid_y:.1f} {rx1:.1f},{mid_y:.1f} {rx1:.1f},{ry:.1f} Z" '
            f'fill="{fill}" opacity="{RIBBON_OP:.2f}"/>\n'
        )

    # ── Track backbones (width proportional to actual genome length) ──
    ref_cy   = track_cy(0)
    ref_w    = reference_length / max_len * CONTENT_W
    out.append(
        f'<rect x="{LEFT:.1f}" y="{ref_cy - TRACK_H / 2:.1f}" '
        f'width="{ref_w:.1f}" height="{TRACK_H:.1f}" '
        f'fill="{LINEAR_BACKBONE_COLOR}" rx="3"/>\n'
    )
    for i in range(n_queries):
        cy  = track_cy(i + 1)
        qw  = query_lengths.get(i, reference_length) / max_len * CONTENT_W
        out.append(
            f'<rect x="{LEFT:.1f}" y="{cy - TRACK_H / 2:.1f}" '
            f'width="{qw:.1f}" height="{TRACK_H:.1f}" '
            f'fill="{LINEAR_BACKBONE_COLOR}" rx="3"/>\n'
        )

    # ── Insertion site markers on reference track ─────────────────────
    if show_insertions and insertion_sites:
        INS_TICK_H = 8
        for site in insertion_sites:
            _qi   = site['query_index']
            ix    = ref_x(site['ref_position'])
            _r, _g, _b = _ring_color(_qi)
            ins_c = _rgb_hex(_r, _g, _b)
            tick_top = ref_cy - TRACK_H / 2 - INS_TICK_H
            out.append(
                f'<line x1="{ix:.1f}" y1="{ref_cy - TRACK_H/2:.1f}" '
                f'x2="{ix:.1f}" y2="{tick_top:.1f}" '
                f'stroke="{ins_c}" stroke-width="1.5" opacity="0.85"/>\n'
            )
            ins = site['insertion_size']
            lbl = (f'+{ins/1e6:.1f}M' if ins >= 1_000_000 else
                   f'+{ins/1e3:.0f}k'  if ins >= 1_000     else f'+{ins}')
            out.append(
                f'<text x="{ix:.1f}" y="{tick_top - 2:.1f}" '
                f'font-size="6" text-anchor="middle" fill="{ins_c}" opacity="0.85">'
                f'{_esc(lbl)}</text>\n'
            )

    # ── Gene-annotation arrows on the reference track ─────────────────
    if annotations:
        for feat in annotations:
            ftype = feat['feature_type'].lower()
            if ftype not in ('cds', 'gene', 'trna', 'rrna'):
                continue
            if show_res_genes and not _is_resistance_gene(feat):
                continue
            ax1 = LEFT + (feat['start'] - 1) / reference_length * CONTENT_W
            ax2 = LEFT + feat['end']           / reference_length * CONTENT_W
            aw  = ax2 - ax1
            if aw < 0.5:
                continue
            color  = _annotation_color(ftype)
            strand = feat.get('strand') or 1
            half_h = GENE_H / 2
            head_w = min(half_h * 1.2, aw * 0.4)
            if strand == 1:
                body_end = max(ax1, ax2 - head_w)
                pts = (
                    f'{ax1:.1f},{ref_cy - half_h:.1f} '
                    f'{body_end:.1f},{ref_cy - half_h:.1f} '
                    f'{ax2:.1f},{ref_cy:.1f} '
                    f'{body_end:.1f},{ref_cy + half_h:.1f} '
                    f'{ax1:.1f},{ref_cy + half_h:.1f}'
                )
            else:
                body_start = min(ax2, ax1 + head_w)
                pts = (
                    f'{ax2:.1f},{ref_cy - half_h:.1f} '
                    f'{body_start:.1f},{ref_cy - half_h:.1f} '
                    f'{ax1:.1f},{ref_cy:.1f} '
                    f'{body_start:.1f},{ref_cy + half_h:.1f} '
                    f'{ax2:.1f},{ref_cy + half_h:.1f}'
                )
            out.append(f'<polygon points="{pts}" fill="{color}" opacity="0.9"/>\n')

    # ── Gene name labels – tiered horizontal rows ─────────────────────
    if _lbl_tiers:
        _gene_top_y = ref_cy - GENE_H / 2
        for (_lx2, _ln2, _lt) in _lbl_tiers:
            _ly = _gene_top_y - _TIER_H * (_lt + 1)
            out.append(
                f'<line x1="{_lx2:.1f}" y1="{_gene_top_y - 1:.1f}" '
                f'x2="{_lx2:.1f}" y2="{_ly + 2:.1f}" '
                f'stroke="#aaa" stroke-width="0.7"/>\n'
            )
            out.append(
                f'<text x="{_lx2:.1f}" y="{_ly:.1f}" '
                f'font-size="{_LBL_FONT}" text-anchor="middle" '
                f'dominant-baseline="auto" fill="#333">'
                f'{_esc(_ln2)}</text>\n'
            )

    # ── Track labels ──────────────────────────────────────────────────
    out.append(
        f'<text x="{LEFT - 8:.1f}" y="{ref_cy + 4:.1f}" '
        f'font-size="12" font-weight="bold" text-anchor="end" fill="#333">'
        f'{_esc(reference_name)}</text>\n'
    )
    for i, name in enumerate(query_names):
        cy        = track_cy(i + 1)
        disp_name = (name[:22] + '\u2026') if len(name) > 23 else name
        out.append(
            f'<text x="{LEFT - 8:.1f}" y="{cy + 4:.1f}" '
            f'font-size="12" font-weight="bold" text-anchor="end" fill="#333">'
            f'{_esc(disp_name)}</text>\n'
        )

    # ── Scale ticks (shared scale, spans max_len) ─────────────────────
    scale_y = ref_cy + TRACK_H / 2 + 15
    out.append(
        f'<line x1="{LEFT:.1f}" y1="{scale_y:.1f}" '
        f'x2="{LEFT + CONTENT_W:.1f}" y2="{scale_y:.1f}" '
        f'stroke="#bbb" stroke-width="0.5"/>\n'
    )
    for t in range(11):
        tx  = LEFT + t / 10 * CONTENT_W
        pos = t / 10 * max_len
        if   max_len > 2_000_000: lbl = f'{pos/1e6:.1f} Mb'
        elif max_len > 2000:      lbl = f'{pos/1e3:.0f} kb'
        else:                     lbl = f'{int(pos)} bp'
        out.append(
            f'<line x1="{tx:.1f}" y1="{scale_y:.1f}" '
            f'x2="{tx:.1f}" y2="{scale_y + 4:.1f}" '
            f'stroke="#888" stroke-width="0.8"/>\n'
        )
        anchor = 'start' if t == 0 else ('end' if t == 10 else 'middle')
        out.append(
            f'<text x="{tx:.1f}" y="{scale_y + 13:.1f}" '
            f'font-size="9" text-anchor="{anchor}" fill="#666">'
            f'{_esc(lbl)}</text>\n'
        )

    # ── Legend ─────────────────────────────────────────────────────────
    # Positioned in the right margin (x > LEFT + CONTENT_W) so it never
    # overlaps the track or ribbon content.
    lx  = LEFT + CONTENT_W + 12   # 12 px inside the right margin
    ly  = TOP
    bs  = 12
    ih  = 20
    gw  = 95    # gradient bar width
    gh  = 9     # gradient bar height

    solid_items = [('Sequence Backbone', LINEAR_BACKBONE_COLOR)]
    if annotations:
        ann_types = []
        seen = set()
        for feat in annotations:
            ft = feat['feature_type'].lower()
            if ft in ('cds', 'gene', 'trna', 'rrna') and ft not in seen:
                if show_res_genes and not _is_resistance_gene(feat):
                    continue
                ann_types.append(ft)
                seen.add(ft)
        for ft in ann_types:
            label = {'cds': 'CDS', 'gene': 'Gene', 'trna': 'tRNA', 'rrna': 'rRNA'}[ft]
            solid_items.append((label, _annotation_color(ft)))

    # Per-query section height: name row + fwd bar + inv bar + spacing
    _pqh = (ih - 4) + (gh + 3) + (gh + 3) + 6
    # Height: solid items + "% Identity" header + per-query sections + min%/100% row + padding
    lbg_h = (len(solid_items) * ih
              + ih                      # "% Identity" header
              + n_queries * _pqh        # per-query gradient rows
              + ih                      # min%/100% label row
              + 20)                     # top/bottom padding
    if show_insertions and insertion_sites:
        lbg_h += ih
    lbg_w = gw + 38

    out.append(
        f'<rect x="{lx - 8:.1f}" y="{ly - 8:.1f}" width="{lbg_w:.1f}" height="{lbg_h:.1f}" '
        f'fill="white" fill-opacity="0.95" stroke="#ddd" stroke-width="1" rx="5"/>\n'
    )
    for label, color in solid_items:
        out.append(f'<rect x="{lx:.1f}" y="{ly:.1f}" width="{bs}" height="{bs}" fill="{color}" rx="2"/>\n')
        out.append(
            f'<text x="{lx + bs + 5:.1f}" y="{ly + bs - 1:.1f}" '
            f'font-size="10" fill="#333">{_esc(label)}</text>\n'
        )
        ly += ih

    # Per-query identity gradient section
    out.append(
        f'<text x="{lx:.1f}" y="{ly + 10:.1f}" '
        f'font-size="10" font-weight="bold" fill="#555">% Identity</text>\n'
    )
    ly += ih
    for _qi, _qname in enumerate(query_names):
        _r, _g, _b = _ring_color(_qi)
        _ring_hex = _rgb_hex(_r, _g, _b)
        _disp = (_qname[:15] + '\u2026') if len(_qname) > 16 else _qname
        # Query colour chip + name
        out.append(f'<rect x="{lx:.1f}" y="{ly:.1f}" width="{bs}" height="{bs}" fill="{_ring_hex}" rx="2"/>\n')
        out.append(
            f'<text x="{lx + bs + 5:.1f}" y="{ly + bs - 1:.1f}" '
            f'font-size="9" fill="#333">{_esc(_disp)}</text>\n'
        )
        ly += ih - 4
        # Forward gradient bar
        out.append(f'<rect x="{lx:.1f}" y="{ly:.1f}" width="{gw}" height="{gh}" '
                   f'fill="url(#lin_fwd_grad_{_qi})" rx="2"/>\n')
        out.append(f'<text x="{lx + gw + 3:.1f}" y="{ly + gh:.1f}" '
                   f'font-size="8" fill="#555">Fwd</text>\n')
        ly += gh + 3
        # Inverted gradient bar (lighter shade)
        out.append(f'<rect x="{lx:.1f}" y="{ly:.1f}" width="{gw}" height="{gh}" '
                   f'fill="url(#lin_inv_grad_{_qi})" rx="2"/>\n')
        out.append(f'<text x="{lx + gw + 3:.1f}" y="{ly + gh:.1f}" '
                   f'font-size="8" fill="#555">Inv</text>\n')
        ly += gh + 9
    out.append(f'<text x="{lx:.1f}" y="{ly + 10:.1f}" font-size="9" fill="#666">{min_identity:.0f}%</text>\n')
    out.append(f'<text x="{lx + gw - 18:.1f}" y="{ly + 10:.1f}" font-size="9" fill="#666">100%</text>\n')
    ly += ih

    # Insertions legend entry
    if show_insertions and insertion_sites:
        out.append(f'<rect x="{lx:.1f}" y="{ly:.1f}" width="{bs}" height="{bs}" fill="#555" rx="2"/>\n')
        out.append(
            f'<text x="{lx + bs + 5:.1f}" y="{ly + bs - 1:.1f}" '
            f'font-size="10" fill="#333">Insertion site</text>\n'
        )

    out.append('</svg>')
    return ''.join(out)


# ============================================================
# Alignment View SVG Generator (pygenomeviz-style, per-section)
# ============================================================

def generate_alignment_svg(blast_hits, annotations, reference_length, reference_name,
                           query_names, show_gene_names=False, min_identity=70.0,
                           insertion_sites=None, show_insertions=True, show_res_genes=False):
    """Generate a pygenomeviz-style alignment SVG.

    Each query gets its own section with the reference track on top and the
    query track below, connected by ribbons. Ribbons never cross between
    different query sections, making the view much easier to interpret.
    Ribbon colours match the circular view's per-genome ring colours.
    """
    insertion_sites = insertion_sites or []
    if reference_length == 0:
        raise ValueError("Reference length cannot be zero.")

    n_queries = len(query_names)

    # Derive display length for each query from max qlen seen in hits
    query_lengths = {}
    for hit in blast_hits:
        qi = hit['query_index']
        ql = hit['qlen']
        if qi not in query_lengths or ql > query_lengths[qi]:
            query_lengths[qi] = ql
    for i in range(n_queries):
        if i not in query_lengths:
            query_lengths[i] = reference_length

    # Shared scale: all genomes use the same pixels-per-bp ratio
    max_len = max(reference_length, max(query_lengths.values(), default=reference_length))

    # ── Layout constants ──────────────────────────────────────────────
    W            = 1200
    LEFT         = 160      # reserved for sequence labels
    RIGHT        = 190      # wide right margin so the legend never overlaps content
    CONTENT_W    = W - LEFT - RIGHT
    TRACK_H      = 12       # backbone bar height (px)
    GENE_H       = 14       # gene arrow height (px)
    RIBBON_ZONE  = 75       # vertical space for ribbons between the two tracks
    PAD_BOT      = 28       # padding below query track in each section
    SECTION_SEP  = 4        # gap line between sections
    RIBBON_OP    = 0.55     # ribbon opacity
    SCALE_H      = 22       # height reserved for scale labels (first section only)

    # ── Pre-compute gene-label tiers (horizontal stacking) ────────────
    _LBL_FONT  = 9
    _LBL_CHARW = _LBL_FONT * 0.58
    _LBL_HPAD  = 8
    _TIER_H    = 14

    _lbl_predata = []
    if show_gene_names and annotations:
        for _f in annotations:
            _ft = _f['feature_type'].lower()
            if _ft not in ('cds', 'gene'):
                continue
            if show_res_genes and not _is_resistance_gene(_f):
                continue
            _ax1 = LEFT + (_f['start'] - 1) / reference_length * CONTENT_W
            _ax2 = LEFT + _f['end']          / reference_length * CONTENT_W
            if _ax2 - _ax1 < 0.5:
                continue
            _gn = (_f['attributes'].get('Name') or _f.get('id') or
                   _f['attributes'].get('gene') or '')
            if _gn:
                _lbl_predata.append(((_ax1 + _ax2) / 2, _gn))
    _lbl_predata.sort(key=lambda t: t[0])

    _tier_right = {}
    _lbl_tiers  = []
    for (_lx, _ln) in _lbl_predata:
        _hw = len(_ln) * _LBL_CHARW / 2
        _t  = 0
        while _t in _tier_right and _lx - _hw < _tier_right[_t] + _LBL_HPAD:
            _t += 1
        _tier_right[_t] = _lx + _hw
        _lbl_tiers.append((_lx, _ln, _t))

    _n_tiers = max((_e[2] for _e in _lbl_tiers), default=-1) + 1
    PAD_TOP   = 10 + _n_tiers * _TIER_H + 6   # enough room for all label rows

    section_h = PAD_TOP + TRACK_H + RIBBON_ZONE + TRACK_H + PAD_BOT
    TOP_HEADER = 50         # space at top for legend

    # First section also needs scale ticks below the reference track
    first_extra = SCALE_H
    H = (TOP_HEADER
         + (section_h + first_extra)     # first section is taller
         + max(0, n_queries - 1) * (section_h + SECTION_SEP)
         + 20)

    # ── Coordinate helpers (shared scale) ────────────────────────────
    def ref_x(pos):
        return LEFT + (pos - 1) / max_len * CONTENT_W

    def qx(pos, qi):
        return LEFT + (pos - 1) / max_len * CONTENT_W

    def section_top(qi):
        if qi == 0:
            return TOP_HEADER
        return TOP_HEADER + (section_h + first_extra) + (qi - 1) * (section_h + SECTION_SEP)

    out = []

    # ── SVG header ────────────────────────────────────────────────────
    out.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{W}" height="{H}" viewBox="0 0 {W} {H}">\n'
    )
    out.append('<style>text { font-family: Arial, Helvetica, sans-serif; }</style>\n')

    # ── Gradient defs for identity legend (one pair per query) ───────
    _grad_fracs = (0, 0.33, 0.66, 1)
    _span = 100.0 - min_identity
    out.append('<defs>\n')
    for _qi in range(n_queries):
        _fwd = ''.join(
            f'  <stop offset="{f}" stop-color="{_linear_hit_color(_qi, min_identity + f * _span, False, min_identity)}"/>\n'
            for f in _grad_fracs
        )
        _inv = ''.join(
            f'  <stop offset="{f}" stop-color="{_linear_hit_color(_qi, min_identity + f * _span, True, min_identity)}"/>\n'
            for f in _grad_fracs
        )
        out.append(f'<linearGradient id="aln_fwd_grad_{_qi}" x1="0" y1="0" x2="1" y2="0">\n{_fwd}</linearGradient>\n')
        out.append(f'<linearGradient id="aln_inv_grad_{_qi}" x1="0" y1="0" x2="1" y2="0">\n{_inv}</linearGradient>\n')
    out.append('</defs>\n')

    out.append(f'<rect width="{W}" height="{H}" fill="white"/>\n')

    # Group hits by query index
    hits_by_query = {i: [] for i in range(n_queries)}
    for hit in blast_hits:
        qi = hit['query_index']
        if qi in hits_by_query:
            hits_by_query[qi].append(hit)

    # ── Draw each query section ───────────────────────────────────────
    for qi, qname in enumerate(query_names):
        stop   = section_top(qi)
        ref_cy = stop + PAD_TOP + TRACK_H / 2
        qry_cy = ref_cy + TRACK_H / 2 + RIBBON_ZONE + TRACK_H / 2

        ref_y_bot = ref_cy + TRACK_H / 2
        qry_y_top = qry_cy - TRACK_H / 2

        # Alternating section background for readability
        bg = '#f8fafc' if qi % 2 == 0 else '#f1f5f9'
        sh = (section_h + first_extra) if qi == 0 else section_h
        out.append(
            f'<rect x="0" y="{stop:.1f}" width="{W}" height="{sh:.1f}" '
            f'fill="{bg}"/>\n'
        )

        # Separator line
        if qi > 0:
            sep_y = stop - SECTION_SEP / 2
            out.append(
                f'<line x1="0" y1="{sep_y:.1f}" x2="{W}" y2="{sep_y:.1f}" '
                f'stroke="#cbd5e1" stroke-width="1"/>\n'
            )

        # ── Alignment ribbons (drawn first, behind backbones) ─────────
        for hit in hits_by_query[qi]:
            rx1 = ref_x(hit['sstart'])
            rx2 = ref_x(hit['send'])
            qs  = min(hit['qstart'], hit['qend'])
            qe  = max(hit['qstart'], hit['qend'])
            qx1 = qx(qs, qi)
            qx2 = qx(qe, qi)

            inverted = hit['qstart'] > hit['qend']
            fill     = _linear_hit_color(qi, hit['pident'], inverted, min_identity)

            if abs(rx2 - rx1) < 0.5 and abs(qx2 - qx1) < 0.5:
                continue

            mid_y = (ref_y_bot + qry_y_top) / 2
            out.append(
                f'<path d="'
                f'M {rx1:.1f},{ref_y_bot:.1f} L {rx2:.1f},{ref_y_bot:.1f} '
                f'C {rx2:.1f},{mid_y:.1f} {qx2:.1f},{mid_y:.1f} {qx2:.1f},{qry_y_top:.1f} '
                f'L {qx1:.1f},{qry_y_top:.1f} '
                f'C {qx1:.1f},{mid_y:.1f} {rx1:.1f},{mid_y:.1f} {rx1:.1f},{ref_y_bot:.1f} Z" '
                f'fill="{fill}" opacity="{RIBBON_OP:.2f}"/>\n'
            )

        # ── Reference backbone (width proportional to actual length) ──
        ref_w = reference_length / max_len * CONTENT_W
        out.append(
            f'<rect x="{LEFT:.1f}" y="{ref_cy - TRACK_H/2:.1f}" '
            f'width="{ref_w:.1f}" height="{TRACK_H:.1f}" '
            f'fill="{LINEAR_BACKBONE_COLOR}" rx="3"/>\n'
        )

        # ── Insertion site markers on this section's reference track ─────
        if show_insertions and insertion_sites:
            INS_TICK_H = 8
            for site in insertion_sites:
                if site['query_index'] != qi:
                    continue
                ix = ref_x(site['ref_position'])
                _r, _g, _b = _ring_color(qi)
                ins_c = _rgb_hex(_r, _g, _b)
                tick_top = ref_cy - TRACK_H / 2 - INS_TICK_H
                out.append(
                    f'<line x1="{ix:.1f}" y1="{ref_cy - TRACK_H/2:.1f}" '
                    f'x2="{ix:.1f}" y2="{tick_top:.1f}" '
                    f'stroke="{ins_c}" stroke-width="1.5" opacity="0.85"/>\n'
                )
                ins = site['insertion_size']
                lbl = (f'+{ins/1e6:.1f}M' if ins >= 1_000_000 else
                       f'+{ins/1e3:.0f}k'  if ins >= 1_000     else f'+{ins}')
                out.append(
                    f'<text x="{ix:.1f}" y="{tick_top - 2:.1f}" '
                    f'font-size="6" text-anchor="middle" fill="{ins_c}" opacity="0.85">'
                    f'{_esc(lbl)}</text>\n'
                )

        # ── Gene annotation arrows on reference ───────────────────────
        if annotations:
            for feat in annotations:
                ftype = feat['feature_type'].lower()
                if ftype not in ('cds', 'gene', 'trna', 'rrna'):
                    continue
                if show_res_genes and not _is_resistance_gene(feat):
                    continue
                ax1 = LEFT + (feat['start'] - 1) / reference_length * CONTENT_W
                ax2 = LEFT + feat['end']           / reference_length * CONTENT_W
                aw  = ax2 - ax1
                if aw < 0.5:
                    continue
                color  = _annotation_color(ftype)
                strand = feat.get('strand') or 1
                half_h = GENE_H / 2
                head_w = min(half_h * 1.2, aw * 0.4)
                if strand == 1:
                    body_end = max(ax1, ax2 - head_w)
                    pts = (
                        f'{ax1:.1f},{ref_cy - half_h:.1f} '
                        f'{body_end:.1f},{ref_cy - half_h:.1f} '
                        f'{ax2:.1f},{ref_cy:.1f} '
                        f'{body_end:.1f},{ref_cy + half_h:.1f} '
                        f'{ax1:.1f},{ref_cy + half_h:.1f}'
                    )
                else:
                    body_start = min(ax2, ax1 + head_w)
                    pts = (
                        f'{ax2:.1f},{ref_cy - half_h:.1f} '
                        f'{body_start:.1f},{ref_cy - half_h:.1f} '
                        f'{ax1:.1f},{ref_cy:.1f} '
                        f'{body_start:.1f},{ref_cy + half_h:.1f} '
                        f'{ax2:.1f},{ref_cy + half_h:.1f}'
                    )
                out.append(f'<polygon points="{pts}" fill="{color}" opacity="0.9"/>\n')

        # ── Gene name labels – tiered horizontal rows ──────────────────
        if _lbl_tiers:
            _gene_top_y = ref_cy - GENE_H / 2
            for (_lx2, _ln2, _lt) in _lbl_tiers:
                _ly = _gene_top_y - _TIER_H * (_lt + 1)
                out.append(
                    f'<line x1="{_lx2:.1f}" y1="{_gene_top_y - 1:.1f}" '
                    f'x2="{_lx2:.1f}" y2="{_ly + 2:.1f}" '
                    f'stroke="#aaa" stroke-width="0.7"/>\n'
                )
                out.append(
                    f'<text x="{_lx2:.1f}" y="{_ly:.1f}" '
                    f'font-size="{_LBL_FONT}" text-anchor="middle" '
                    f'dominant-baseline="auto" fill="#333">'
                    f'{_esc(_ln2)}</text>\n'
                )

        # ── Scale ticks (first section only, shared scale spans max_len)
        if qi == 0:
            scale_y = ref_cy + TRACK_H / 2 + 6
            out.append(
                f'<line x1="{LEFT:.1f}" y1="{scale_y:.1f}" '
                f'x2="{LEFT + CONTENT_W:.1f}" y2="{scale_y:.1f}" '
                f'stroke="#bbb" stroke-width="0.5"/>\n'
            )
            for t in range(11):
                tx  = LEFT + t / 10 * CONTENT_W
                pos = t / 10 * max_len
                if   max_len > 2_000_000: lbl = f'{pos/1e6:.1f} Mb'
                elif max_len > 2000:      lbl = f'{pos/1e3:.0f} kb'
                else:                     lbl = f'{int(pos)} bp'
                out.append(
                    f'<line x1="{tx:.1f}" y1="{scale_y:.1f}" '
                    f'x2="{tx:.1f}" y2="{scale_y + 4:.1f}" '
                    f'stroke="#888" stroke-width="0.8"/>\n'
                )
                anchor = 'start' if t == 0 else ('end' if t == 10 else 'middle')
                out.append(
                    f'<text x="{tx:.1f}" y="{scale_y + 14:.1f}" '
                    f'font-size="9" text-anchor="{anchor}" fill="#666">'
                    f'{_esc(lbl)}</text>\n'
                )

        # ── Query backbone (width proportional to actual length) ────────
        qw = query_lengths.get(qi, reference_length) / max_len * CONTENT_W
        out.append(
            f'<rect x="{LEFT:.1f}" y="{qry_cy - TRACK_H/2:.1f}" '
            f'width="{qw:.1f}" height="{TRACK_H:.1f}" '
            f'fill="{LINEAR_BACKBONE_COLOR}" rx="3"/>\n'
        )

        # ── Track labels ──────────────────────────────────────────────
        out.append(
            f'<text x="{LEFT - 8:.1f}" y="{ref_cy + 4:.1f}" '
            f'font-size="11" font-weight="bold" text-anchor="end" fill="#374151">'
            f'{_esc(reference_name)}</text>\n'
        )
        disp_name = (qname[:22] + '\u2026') if len(qname) > 23 else qname
        out.append(
            f'<text x="{LEFT - 8:.1f}" y="{qry_cy + 4:.1f}" '
            f'font-size="11" font-weight="bold" text-anchor="end" fill="#374151">'
            f'{_esc(disp_name)}</text>\n'
        )

    # ── Legend ─────────────────────────────────────────────────────────
    # Positioned in the right margin (x > LEFT + CONTENT_W) so it never
    # overlaps the track or ribbon content.
    lx  = LEFT + CONTENT_W + 12   # 12 px inside the right margin
    ly  = TOP_HEADER + 10
    bs  = 12
    ih  = 20
    gw  = 95    # gradient bar width
    gh  = 9     # gradient bar height

    solid_items = [('Sequence Backbone', LINEAR_BACKBONE_COLOR)]
    if annotations:
        ann_types = []
        seen = set()
        for feat in annotations:
            ft = feat['feature_type'].lower()
            if ft in ('cds', 'gene', 'trna', 'rrna') and ft not in seen:
                if show_res_genes and not _is_resistance_gene(feat):
                    continue
                ann_types.append(ft)
                seen.add(ft)
        for ft in ann_types:
            label = {'cds': 'CDS', 'gene': 'Gene', 'trna': 'tRNA', 'rrna': 'rRNA'}[ft]
            solid_items.append((label, _annotation_color(ft)))

    # Per-query section height: name row + fwd bar + inv bar + spacing
    _pqh = (ih - 4) + (gh + 3) + (gh + 3) + 6
    # Height: solid items + "% Identity" header + per-query sections + min%/100% row + padding
    lbg_h = (len(solid_items) * ih
              + ih                      # "% Identity" header
              + n_queries * _pqh        # per-query gradient rows
              + ih                      # min%/100% label row
              + 20)                     # top/bottom padding
    if show_insertions and insertion_sites:
        lbg_h += ih
    lbg_w = gw + 38

    out.append(
        f'<rect x="{lx - 8:.1f}" y="{ly - 8:.1f}" width="{lbg_w:.1f}" height="{lbg_h:.1f}" '
        f'fill="white" fill-opacity="0.95" stroke="#ddd" stroke-width="1" rx="5"/>\n'
    )
    for label, color in solid_items:
        out.append(f'<rect x="{lx:.1f}" y="{ly:.1f}" width="{bs}" height="{bs}" fill="{color}" rx="2"/>\n')
        out.append(
            f'<text x="{lx + bs + 5:.1f}" y="{ly + bs - 1:.1f}" '
            f'font-size="10" fill="#333">{_esc(label)}</text>\n'
        )
        ly += ih

    # Per-query identity gradient section
    out.append(
        f'<text x="{lx:.1f}" y="{ly + 10:.1f}" '
        f'font-size="10" font-weight="bold" fill="#555">% Identity</text>\n'
    )
    ly += ih
    for _qi, _qname in enumerate(query_names):
        _r, _g, _b = _ring_color(_qi)
        _ring_hex = _rgb_hex(_r, _g, _b)
        _disp = (_qname[:15] + '\u2026') if len(_qname) > 16 else _qname
        # Query colour chip + name
        out.append(f'<rect x="{lx:.1f}" y="{ly:.1f}" width="{bs}" height="{bs}" fill="{_ring_hex}" rx="2"/>\n')
        out.append(
            f'<text x="{lx + bs + 5:.1f}" y="{ly + bs - 1:.1f}" '
            f'font-size="9" fill="#333">{_esc(_disp)}</text>\n'
        )
        ly += ih - 4
        # Forward gradient bar
        out.append(f'<rect x="{lx:.1f}" y="{ly:.1f}" width="{gw}" height="{gh}" '
                   f'fill="url(#aln_fwd_grad_{_qi})" rx="2"/>\n')
        out.append(f'<text x="{lx + gw + 3:.1f}" y="{ly + gh:.1f}" '
                   f'font-size="8" fill="#555">Fwd</text>\n')
        ly += gh + 3
        # Inverted gradient bar (lighter shade)
        out.append(f'<rect x="{lx:.1f}" y="{ly:.1f}" width="{gw}" height="{gh}" '
                   f'fill="url(#aln_inv_grad_{_qi})" rx="2"/>\n')
        out.append(f'<text x="{lx + gw + 3:.1f}" y="{ly + gh:.1f}" '
                   f'font-size="8" fill="#555">Inv</text>\n')
        ly += gh + 9
    out.append(f'<text x="{lx:.1f}" y="{ly + 10:.1f}" font-size="9" fill="#666">{min_identity:.0f}%</text>\n')
    out.append(f'<text x="{lx + gw - 18:.1f}" y="{ly + 10:.1f}" font-size="9" fill="#666">100%</text>\n')
    ly += ih

    # Insertions legend entry
    if show_insertions and insertion_sites:
        out.append(f'<rect x="{lx:.1f}" y="{ly:.1f}" width="{bs}" height="{bs}" fill="#555" rx="2"/>\n')
        out.append(
            f'<text x="{lx + bs + 5:.1f}" y="{ly + bs - 1:.1f}" '
            f'font-size="10" fill="#333">Insertion site</text>\n'
        )

    out.append('</svg>')
    return ''.join(out)


# ============================================================
# Similarity-based Query Clustering
# ============================================================

def _coverage_vector(hits, ref_len, n_bins=500):
    """Binary coverage vector at n_bins resolution for a set of hits."""
    bins = [0] * n_bins
    for hit in hits:
        s = int((min(hit['sstart'], hit['send']) - 1) / ref_len * n_bins)
        e = int((max(hit['sstart'], hit['send']) - 1) / ref_len * n_bins) + 1
        for b in range(max(0, s), min(n_bins, e)):
            bins[b] = 1
    return bins


def _jaccard(a, b):
    inter = sum(x & y for x, y in zip(a, b))
    union = sum(x | y for x, y in zip(a, b))
    return inter / union if union > 0 else 0.0


def cluster_queries_by_similarity(blast_hits, query_names, ref_len):
    """Reorder queries using greedy nearest-neighbour on reference-coverage Jaccard similarity.

    Returns (new_names_list, old_index_to_new_index_dict).
    The query most similar to the reference (highest coverage) is placed first,
    and subsequent queries are chosen to be as similar as possible to the previous one.
    """
    n = len(query_names)
    if n <= 1:
        return list(query_names), {i: i for i in range(n)}

    hits_by_q = {i: [] for i in range(n)}
    for hit in blast_hits:
        hits_by_q[hit['query_index']].append(hit)

    vecs     = [_coverage_vector(hits_by_q[i], ref_len) for i in range(n)]
    coverage = [sum(v) for v in vecs]

    # Start with the query that covers the most reference
    start    = coverage.index(max(coverage))
    visited  = [False] * n
    order    = [start]
    visited[start] = True

    for _ in range(n - 1):
        last = order[-1]
        best_sim, best_j = -1, -1
        for j in range(n):
            if visited[j]:
                continue
            sim = _jaccard(vecs[last], vecs[j])
            if sim > best_sim:
                best_sim, best_j = sim, j
        order.append(best_j)
        visited[best_j] = True

    new_names   = [query_names[i] for i in order]
    old_to_new  = {old: new for new, old in enumerate(order)}
    return new_names, old_to_new


# ============================================================
# Main Pipeline
# ============================================================

def run_comparison(reference_fasta_text, query_alignment_texts, query_file_names,
                   annotation_text=None,
                   min_identity=70.0, min_coverage=0.0, min_length=100,
                   show_gene_names=False, show_insertions=True,
                   progress_callback=None, show_res_genes=False):
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

    # 3. Cluster queries by similarity so most-similar genomes are adjacent
    _prog(2 + n_queries, total_steps, 'Clustering queries by similarity…')
    ordered_names, idx_remap = cluster_queries_by_similarity(all_hits, query_file_names, ref_len)
    for hit in all_hits:
        hit['query_index'] = idx_remap[hit['query_index']]
    for site in all_insertions:
        site['query_index'] = idx_remap[site['query_index']]

    # 4. Parse annotations
    _prog(2 + n_queries, total_steps, 'Parsing annotations…')
    annotations = None
    if annotation_text and annotation_text.strip():
        try:
            annotations = parse_gff3(annotation_text, ref_name)
        except Exception as e:
            print(f'Warning: annotation parse error: {e}')

    # 5. Generate plots
    _prog(3 + n_queries - 1, total_steps, 'Generating plots…')
    circular_svg = generate_svg(
        blast_hits=all_hits,
        annotations=annotations,
        reference_length=ref_len,
        query_names=ordered_names,
        reference_display_name=ref_name,
        show_gene_names=show_gene_names,
        insertion_sites=all_insertions,
        show_insertions=show_insertions,
        min_identity=min_identity,
        show_res_genes=show_res_genes,
    )
    linear_svg = generate_linear_svg(
        blast_hits=all_hits,
        annotations=annotations,
        reference_length=ref_len,
        reference_name=ref_name,
        query_names=ordered_names,
        show_gene_names=show_gene_names,
        min_identity=min_identity,
        insertion_sites=all_insertions,
        show_insertions=show_insertions,
        show_res_genes=show_res_genes,
    )
    alignment_svg = generate_alignment_svg(
        blast_hits=all_hits,
        annotations=annotations,
        reference_length=ref_len,
        reference_name=ref_name,
        query_names=ordered_names,
        show_gene_names=show_gene_names,
        min_identity=min_identity,
        insertion_sites=all_insertions,
        show_insertions=show_insertions,
        show_res_genes=show_res_genes,
    )
    _prog(total_steps, total_steps, 'Done!')
    return circular_svg, linear_svg, alignment_svg
