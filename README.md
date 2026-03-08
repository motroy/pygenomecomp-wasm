# PyGenomeComp Web

A self-contained **WebAssembly** port of [pygenomecomp](https://github.com/motroy/pygenomecomp) — a circular genome comparison visualizer. Runs entirely in your browser; no server, no installation required.

**Live app →** [motroy.github.io/pygenomecomp-wasm](https://motroy.github.io/pygenomecomp-wasm) *(enabled after GitHub Pages is configured)*

---

## What it does

Upload one or more genome assemblies in FASTA format and the tool will:

1. Align each **query** genome against the **reference** using **lastz** (via [biowasm](https://biowasm.com))
2. Optionally overlay **GFF3 annotations** (CDS, gene, rRNA, tRNA features)
3. Render a **BRIG-style circular plot**, as well as **Linear** and **Alignment** views, showing sequence identity as coloured rings or ribbons
4. Let you **download** the SVGs for publication

All computation runs in your browser via [biowasm](https://biowasm.com) and [Pyodide](https://pyodide.org) — so your genomic data never leaves your machine.

---

## Key differences from the CLI version

| Feature | CLI (`pygenomecomp`) | Web WASM |
|---|---|---|
| Alignment engine | BLAST+ (external binary) | lastz (via biowasm) |
| Annotation formats | GFF3 + GenBank (BioPython) | GFF3 only (pure Python) |
| Output | Circular SVG file on disk | Circular, Linear, and Alignment SVGs rendered in browser + download |
| Dependencies | Python, BLAST+, biopython, pandas | None (self-contained) |

---

## Repository layout

```
index.html              # Single-page web app (HTML + JS)
pygenomecomp_wasm.py    # Python logic loaded by Pyodide at runtime
.github/
  workflows/
    deploy.yml          # GitHub Pages deployment workflow
```

---

## Local development

No build step needed — just serve the files over HTTP:

```bash
# Python built-in server
python -m http.server 8080
# then open http://localhost:8080
```

The Pyodide runtime is loaded from CDN on first page visit (~10 MB, cached by browser).

---

## GitHub Pages deployment

Push to `main` or `master` and the included GitHub Actions workflow deploys automatically.

To enable GitHub Pages in your repository:
1. Go to **Settings → Pages**
2. Under **Source**, select **GitHub Actions**
3. The next push will deploy the app

---

## License

MIT — same as the original pygenomecomp project.
