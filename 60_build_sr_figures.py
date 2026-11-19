# -*- coding: utf-8 -*-
"""Build submission-ready 300 dpi TIFF figures for Scientific Reports."""
import os
from PIL import Image, ImageDraw, ImageFont

BASE = r'C:\Users\19354\Doubao\chats\2026-11-18\new-chat\pulp_bioinfo'
RES = os.path.join(BASE, 'results')
DEEP = os.path.join(RES, 'deep')
SC = os.path.join(BASE, 'sc', 'results')
OUTDIR = os.path.join(RES, 'sr_figures')
os.makedirs(OUTDIR, exist_ok=True)

FONT = r'C:\Windows\Fonts\arialbd.ttf'

MAIN = {
    1: [('a', os.path.join(RES, 'volcano_GSE77459.png')),
        ('b', os.path.join(DEEP, 'GSEA_Hallmark_bar.png'))],
    2: [(None, os.path.join(RES, 'immune_heatmap.png'))],
    3: [('a', os.path.join(SC, 'GSE185222_UMAP_clusters.png')),
        ('b', os.path.join(SC, 'GSE185222_composition_bar.png'))],
    4: [('a', os.path.join(SC, 'EMTAB_composition_stage.png')),
        ('b', os.path.join(SC, 'EMTAB_endothelial_subtype_pct.png'))],
    5: [('a', os.path.join(SC, 'EMTAB_celltype_GSEA_heatmap.png')),
        ('b', os.path.join(SC, 'EMTAB_endothelial_marker_heatmap.png'))],
    6: [('a', os.path.join(DEEP, 'TF_activity_top20.png')),
        ('b', os.path.join(DEEP, 'PROGENy_bulk.png'))],
    7: [(None, os.path.join(DEEP, 'DSigDB_drug_direction_top40.png'))],
    8: [('a', os.path.join(SC, 'EMTAB_LR_focus_axes.png')),
        ('b', os.path.join(SC, 'EMTAB_pseudotime_umap.png')),
        ('c', os.path.join(SC, 'EMTAB_pseudotime_gene_trend.png'))],
}

SUPP = {
    1: [os.path.join(RES, 'top50DEG_heatmap.png')],
    2: [os.path.join(RES, 'GO_BP_up.png')],
    3: [os.path.join(RES, 'KEGG_up.png')],
    4: [os.path.join(DEEP, 'GSEA_GOBP_bar.png')],
    5: [os.path.join(RES, 'immune_boxplot.png')],
    6: [os.path.join(RES, 'ROC_curves.png')],
    7: [os.path.join(DEEP, 'GSEA_top3_curves.png')],
    8: [os.path.join(SC, 'GSE185222_keygene_dotplot_white.png')],
    9: [os.path.join(SC, 'GSE185222_LR_delta_top20.png')],
    10: [os.path.join(SC, 'GSE185222_trajectory.png')],
    11: [os.path.join(SC, 'GSE77459_deconv_barplot.png')],
    12: [os.path.join(DEEP, 'WGCNA_softthreshold.png')],
    13: [os.path.join(DEEP, 'WGCNA_dendrogram.png')],
    14: [os.path.join(DEEP, 'WGCNA_module_trait_heatmap.png')],
    15: [os.path.join(DEEP, 'WGCNA_eigengene_box.png')],
    16: [os.path.join(DEEP, 'ML_feature_votes.png')],
    17: [os.path.join(DEEP, 'ML_stability.png')],
    18: [os.path.join(DEEP, 'ML_LOOCV_ROC.png')],
    19: [os.path.join(RES, 'PBM_LPS_reversal_scatter.png')],
    20: [os.path.join(RES, 'PBM_key_genes_bar.png')],
    21: [os.path.join(SC, 'EMTAB_signature7_heatmap.png')],
    22: [os.path.join(SC, 'EMTAB_M1M2_stage.png')],
    23: [os.path.join(SC, 'EMTAB_GSEA_Endothelial_cells_3.png'),
         os.path.join(SC, 'EMTAB_GSEA_Fibroblast___Stromal_cells_1.png'),
         os.path.join(SC, 'EMTAB_GSEA_Blood___Immune_cells_1.png')],
    24: [os.path.join(SC, 'EMTAB_immune_M1M2_stage.png')],
    25: [os.path.join(SC, 'GSE185222_mac_M1M2_boxplot.png')],
    26: [os.path.join(SC, 'EMTAB_PROGENy_stage.png')],
    27: [os.path.join(SC, 'EMTAB_TF_activity_stages.png')],
    28: [os.path.join(SC, 'EMTAB_LR_pair_gain.png')],
}

MANIFEST = []


def load_rgb(path, target_w):
    im = Image.open(path)
    if im.mode == 'P':
        im = im.convert('RGBA')
    if im.mode in ('RGBA', 'LA'):
        # PNG 常为 RGBA：背景 RGB 为黑、alpha=0。直接 convert('RGB') 会丢掉 alpha
        # 并把整片背景变成纯黑，导致黑底黑字不可读。必须先合成到白底。
        rgba = im.convert('RGBA')
        bg = Image.new('RGB', rgba.size, (255, 255, 255))
        bg.paste(rgba, (0, 0), rgba)
        im = bg
    elif im.mode != 'RGB':
        im = im.convert('RGB')
    if im.width != target_w:
        h = int(round(im.height * target_w / im.width))
        im = im.resize((target_w, h), Image.LANCZOS)
    return im


def stack(paths, labels, out_path, min_w=1200, max_w=3000):
    """Single panel, or vertical stack for two, or 'a on top + b/c side by side' for three."""
    widths = []
    for p in paths:
        with Image.open(p) as im0:
            widths.append(im0.width)
    W = min(max(max(widths), min_w), max_w)
    gap = 26
    label_h = 70 if labels and any(labels) else 0

    if len(paths) <= 2:
        panels = [load_rgb(p, W) for p in paths]
        total_h = sum(p.height for p in panels) + label_h * len(panels) + gap * (len(panels) - 1)
        canvas = Image.new('RGB', (W, total_h), 'white')
        y = 0
        for i, im in enumerate(panels):
            if labels and labels[i]:
                ImageDraw.Draw(canvas).text((12, y + 4), f'({labels[i]})',
                                            font=ImageFont.truetype(FONT, 58), fill='black')
                y += label_h
            canvas.paste(im, (0, y))
            y += im.height + gap
    else:
        labs = list(labels) if labels else [None, None, None]
        lh = 70 if any(labs) else 0
        a = load_rgb(paths[0], W)
        half = (W - gap) // 2
        b = load_rgb(paths[1], half)
        c = load_rgb(paths[2], half)
        row2 = max(b.height, c.height)
        total_h = lh + a.height + gap + lh + row2
        canvas = Image.new('RGB', (W, total_h), 'white')
        dr = ImageDraw.Draw(canvas)
        fnt = ImageFont.truetype(FONT, 58)
        if labs[0]:
            dr.text((12, 4), f'({labs[0]})', font=fnt, fill='black')
        canvas.paste(a, (0, lh))
        y2 = lh + a.height + gap
        if labs[1]:
            dr.text((12, y2 + 4), f'({labs[1]})', font=fnt, fill='black')
        if labs[2]:
            dr.text((half + gap + 12, y2 + 4), f'({labs[2]})', font=fnt, fill='black')
        canvas.paste(b, (0, y2 + lh))
        canvas.paste(c, (half + gap, y2 + lh))

    max_h = 2760
    if canvas.height > max_h:
        s = max_h / canvas.height
        canvas = canvas.resize((int(canvas.width * s), max_h), Image.LANCZOS)
    canvas.save(out_path, format='TIFF', compression='tiff_lzw', dpi=(300, 300))
    MANIFEST.append((os.path.basename(out_path), canvas.width, canvas.height,
                     ';'.join(os.path.basename(p) for p in paths)))


for num, panels in MAIN.items():
    labels = [p[0] for p in panels]
    paths = [p[1] for p in panels]
    for p in paths:
        if not os.path.exists(p):
            print('MISSING:', p)
    out = os.path.join(OUTDIR, f'Figure{num}.tif')
    stack(paths, labels, out)
    print('SAVED:', out)

for num, paths in SUPP.items():
    for p in paths:
        if not os.path.exists(p):
            print('MISSING:', p)
    out = os.path.join(OUTDIR, f'SupplementaryFig{num}.tif')
    stack(paths, None, out, min_w=1200)
    print('SAVED:', out)

with open(os.path.join(OUTDIR, 'figure_manifest.txt'), 'w', encoding='utf-8') as f:
    f.write('file\twidth_px\theight_px\t300dpi_inches\tsource_pngs\n')
    for name, w, h, src in MANIFEST:
        f.write(f'{name}\t{w}\t{h}\t{w/300:.1f}x{h/300:.1f}\t{src}\n')
print('MANIFEST:', os.path.join(OUTDIR, 'figure_manifest.txt'))
