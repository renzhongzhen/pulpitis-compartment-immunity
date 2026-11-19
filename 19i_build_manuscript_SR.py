# -*- coding: utf-8 -*-
"""Scientific Reports formatted manuscript: structured as Intro-Results-Discussion-Methods, 8 composite figures, superscript Nature-style citations."""
import re
from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_ALIGN_VERTICAL
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

OUT = 'C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/results/Manuscript_pulpitis_ScientificReports_v1.docx'

doc = Document()
sec = doc.sections[0]
sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)
sec.top_margin = sec.bottom_margin = sec.left_margin = sec.right_margin = Cm(2.5)

st = doc.styles['Normal']
st.font.name = 'Times New Roman'
st.font.size = Pt(12)
st.font.color.rgb = RGBColor(0, 0, 0)
st.paragraph_format.line_spacing = 2.0
st.paragraph_format.space_before = Pt(0)
st.paragraph_format.space_after = Pt(0)
st.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.LEFT

hs1 = doc.styles['Heading 1']
hs1.font.name = 'Times New Roman'
hs1.font.size = Pt(12)
hs1.font.bold = True
hs1.font.color.rgb = RGBColor(0, 0, 0)
hs1.paragraph_format.space_before = Pt(14)
hs1.paragraph_format.space_after = Pt(6)
hs1.paragraph_format.keep_with_next = True
hs1.paragraph_format.line_spacing = 2.0
rpr = hs1.element.get_or_add_rPr()
rf = rpr.find(qn('w:rFonts'))
if rf is None:
    rf = OxmlElement('w:rFonts'); rpr.append(rf)
rf.set(qn('w:eastAsia'), '宋体')

CITE = re.compile(r'\[([\d,\s]+)\]')


def title(text):
    text = uk(text)
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_after = Pt(12); p.paragraph_format.space_before = Pt(6)
    p.paragraph_format.keep_with_next = True
    r = p.add_run(text); r.font.name = 'Times New Roman'; r.font.size = Pt(14); r.font.bold = True


def center(text, size=12, italic=False, space_after=6):
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_after = Pt(space_after); p.paragraph_format.line_spacing = 1.5
    r = p.add_run(text); r.font.name = 'Times New Roman'; r.font.size = Pt(size); r.font.italic = italic


def h1(text):
    return doc.add_heading(text, level=1)


def h2(text):
    text = uk(text)
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(10); p.paragraph_format.space_after = Pt(4)
    p.paragraph_format.keep_with_next = True; p.paragraph_format.line_spacing = 2.0
    r = p.add_run(text); r.font.name = 'Times New Roman'; r.font.size = Pt(12); r.font.bold = True


UK_MAP = [
    ('Normalization', 'Normalisation'), ('normalization', 'normalisation'),
    ('Normalized', 'Normalised'), ('normalized', 'normalised'),
    ('normalize', 'normalise'), ('normalizes', 'normalises'),
    ('Localization', 'Localisation'), ('localization', 'localisation'),
    ('Localized', 'Localised'), ('localized', 'localised'), ('localize', 'localise'),
    ('compartmentalized', 'compartmentalised'),
    ('emphasized', 'emphasised'),
    ('generalization', 'generalisation'), ('generalized', 'generalised'), ('generalizes', 'generalises'),
    ('polarization', 'polarisation'), ('polarized', 'polarised'),
    ('organized', 'organised'), ('organizes', 'organises'),
    ('summarized', 'summarised'),
    ('colored', 'coloured'),
    ('signaling', 'signalling'),
]


def uk(text):
    for a, b in UK_MAP:
        text = text.replace(a, b)
    return text


def body(text):
    text = uk(text)
    p = doc.add_paragraph()
    pos = 0
    for m in CITE.finditer(text):
        if m.start() > pos:
            r = p.add_run(text[pos:m.start()]); r.font.name = 'Times New Roman'; r.font.size = Pt(12)
        r = p.add_run('[' + m.group(1) + ']'); r.font.name = 'Times New Roman'; r.font.size = Pt(12)
        r.font.superscript = True
        pos = m.end()
    if pos < len(text):
        r = p.add_run(text[pos:]); r.font.name = 'Times New Roman'; r.font.size = Pt(12)
    return p


def table_caption(text):
    text = uk(text)
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.LEFT
    p.paragraph_format.space_before = Pt(10); p.paragraph_format.space_after = Pt(4)
    p.paragraph_format.keep_with_next = True; p.paragraph_format.line_spacing = 1.5
    r = p.add_run(text); r.font.name = 'Times New Roman'; r.font.size = Pt(10.5); r.font.bold = True


def _border(cell, tag, sz):
    tcPr = cell._tc.get_or_add_tcPr()
    b = tcPr.find(qn('w:tcBorders'))
    if b is None:
        b = OxmlElement('w:tcBorders'); tcPr.append(b)
    el = OxmlElement(tag); el.set(qn('w:val'), 'single'); el.set(qn('w:sz'), str(sz)); el.set(qn('w:color'), '000000')
    b.append(el)


def add_table(headers, rows, aligns):
    t = doc.add_table(rows=1 + len(rows), cols=len(headers))
    t.alignment = 1
    tblPr = t._tbl.tblPr
    borders = OxmlElement('w:tblBorders')
    for tag, sz in [('w:top', '8'), ('w:bottom', '8')]:
        el = OxmlElement(tag); el.set(qn('w:val'), 'single'); el.set(qn('w:sz'), sz); el.set(qn('w:color'), '000000')
        borders.append(el)
    for tag in ['w:left', 'w:right', 'w:insideH', 'w:insideV']:
        el = OxmlElement(tag); el.set(qn('w:val'), 'none'); borders.append(el)
    tblPr.append(borders)
    for i, hd in enumerate(headers):
        c = t.rows[0].cells[i]; c.vertical_alignment = WD_ALIGN_VERTICAL.CENTER
        _border(c, 'w:bottom', '4')
        p = c.paragraphs[0]; p.alignment = WD_ALIGN_PARAGRAPH.CENTER; p.paragraph_format.line_spacing = 1
        r = p.add_run(str(hd)); r.font.name = 'Times New Roman'; r.font.size = Pt(10); r.font.bold = True
    for ri, row in enumerate(rows):
        trPr = t.rows[ri + 1]._tr.get_or_add_trPr()
        trPr.append(OxmlElement('w:cantSplit'))
        for ci, val in enumerate(row):
            c = t.rows[ri + 1].cells[ci]; c.vertical_alignment = WD_ALIGN_VERTICAL.CENTER
            p = c.paragraphs[0]
            p.alignment = {'l': WD_ALIGN_PARAGRAPH.LEFT, 'c': WD_ALIGN_PARAGRAPH.CENTER, 'r': WD_ALIGN_PARAGRAPH.RIGHT}[aligns[ci]]
            p.paragraph_format.line_spacing = 1
            r = p.add_run(str(val)); r.font.name = 'Times New Roman'; r.font.size = Pt(10)
    t.rows[0]._tr.get_or_add_trPr().append(OxmlElement('w:tblHeader'))


def ref(text):
    p = doc.add_paragraph()
    p.paragraph_format.line_spacing = 1.5; p.paragraph_format.space_after = Pt(3)
    hind = OxmlElement('w:ind'); hind.set(qn('w:left'), '360'); hind.set(qn('w:hanging'), '360')
    p._p.get_or_add_pPr().append(hind)
    r = p.add_run(text); r.font.name = 'Times New Roman'; r.font.size = Pt(10.5)


# ================= Title block =================
title('Compartment-specific immunity in irreversible pulpitis centres on an inflammatory venule')

center('Zhongzhen Ren¹*', 12)
center('¹School of Stomatology, Shandong Medical and Pharmaceutical University, Yantai 264003, China', 10.5, italic=True)
center('*Corresponding author: Zhongzhen Ren, School of Stomatology, Shandong Medical and Pharmaceutical University, Yantai 264003, China', 10.5, italic=True)
center('E-mail: 1573741309@qq.com; ORCID: 0009-0008-6735-5224', 10.5, italic=True)

# ================= Abstract (no heading, <=150 words per Scientific Reports) =================
body('Pulpitis lacks a cellular definition. We integrated bulk microarray, single-cell RNA sequencing and an independent three-stage cohort, and re-analysed an 808 nm photobiomodulation (PBM) panel. Bulk analysis (six healthy, six inflamed pulps) identified 713 up- and 273 downregulated genes activating TNF-alpha/NF-kB, interferon and IL6-JAK-STAT3 programmes. Single-cell profiling (6,582 cells) showed B-cell appearance, plasma-cell expansion and fibroblast contraction, with IL6, CXCL10 and SELE localized to endothelium. In 12 patients across three stages (53,407 cells), endothelium expanded directionally from 20.8% to 48.0% while stroma contracted from 56.9% to 18.9%, and cell-type-resolved enrichment placed these programmes in endothelium and immune cells; the expanding cluster was ACKR1+SELE+ venular endothelium. NFKB1/RELA/SPI1 drove the programme, STAT1 and MYC rose across stages, and 3,375 drug signatures recovered NF-kB and JAK-STAT inhibitors. PBM reversed 32 of 72 panel targets while preserving ALPL and TWIST1. Pulpitis is a compartment-specific state organized around an inflammatory venule; all findings are hypotheses.')

# ================= Introduction (no heading per Scientific Reports) =================
body('Pulpitis, inflammation of the dental pulp, is commonly triggered by caries, trauma or restorative procedures. Its classification as reversible or irreversible still relies largely on patient symptoms and on cold, heat and electric pulp testing, which are subjective and of limited prognostic value, so a molecular and cellular definition would support objective diagnosis and reveal therapeutic targets.')
body('Bulk transcriptome studies have begun to map pulpitis. Galicia et al. showed upregulation of immune-response and cytokine-cytokine receptor genes and downregulation of dentinogenesis genes [1], long non-coding RNAs have been profiled in inflamed pulp [2], and a merged microarray analysis proposed diagnostic candidates including ICAM1 and hub genes such as PTPRC, TLR8, CD86, IL6, CCL2, MMP9, CXCL8, SPI1 and LCP2 [3]. These studies converge on an innate myeloid-dominated programme, but tissue averages cannot say which cells transcribe it or which features survive systematic selection.')
body('Single-cell RNA sequencing has begun to resolve pulpal cellular architecture. In a pilot comparison of sound and carious teeth, deep but not enamel caries altered immune composition, with plasma-cell enrichment and appearance of B cells and CD103+ dendritic cells [4], suggesting that progression involves a composition transition rather than uniform activation. These observations have not been connected to bulk biomarkers or to therapeutic modulation.')
body('Photobiomodulation (PBM) uses low-energy red or near-infrared light to modulate cell function, and 808 nm irradiation attenuates inflammation and modulates osteogenic and inflammatory gene expression in lipopolysaccharide (LPS)-challenged human dental pulp stem cells (DPSCs) [5]. PBM is attractive for vital pulp therapy, but its gene-level effects have not been placed in the cellular context of pulpitis.')
body('We therefore asked one question in five layers: is the pulpitis inflammatory programme uniform across the tissue, or is it owned by a specific compartment, driven by specific regulators and addressable by specific interventions? Bulk transcriptomics defined the programme, single-cell and three-stage cohorts assigned it to compartments, cell-type-resolved enrichment and endothelial subtyping identified the cellular interface, independent regulator, drug-signature and communication algorithms tested convergence, and robustness analyses delimited what survives. An 808 nm PBM panel was then re-analysed to connect the compartment map to a candidate intervention.')

# ================= Results =================
h1('Results')
h2('The bulk transcriptome defines the inflammatory programme')
body('To establish what the inflamed pulp transcribes before asking which cells produce it, six pulpitis and six normal samples were compared. There were 713 upregulated and 273 downregulated genes (Fig. 1a and Supplementary Figs. 1-4). Upregulated genes were enriched in leukocyte migration (adjusted P = 3.6e-46), chemotaxis (4.4e-45) and immune response-regulating cell surface receptor signalling (4.5e-42), and KEGG analysis highlighted Staphylococcus aureus infection (9.3e-23), cytokine-cytokine receptor interaction (2.1e-17), chemokine signalling (9.5e-16) and osteoclast differentiation (3.7e-09). Preranked GSEA then tested the same data without a threshold (Fig. 1b): allograft rejection (NES 3.26, adjusted P = 4.4e-50), TNF-alpha signalling via NF-kB (3.23, 3.1e-47), inflammatory response (3.14, 1.5e-41), interferon-gamma response (2.96, 8.8e-32), IL6-JAK-STAT3 signalling (2.80, 8.0e-18) and complement (2.62, 8.9e-19) were activated, whereas oxidative phosphorylation (-1.88), cytoplasmic translation, RNA processing and cilium assembly were suppressed. The co-activation of NF-kB, interferon and IL6-JAK-STAT3 arms, and the suppression of cilium assembly, which is required for odontoblast polarization, define the programme that the following layers must localize.')
h2('The immune infiltrate dominates the bulk signal')
body('To quantify which immune populations account for the bulk signal, 12 immune scores were computed by ssGSEA (Fig. 2 and Supplementary Fig. 5). All 12 increased in pulpitis (Wilcoxon, FDR < 0.01); because ssGSEA scores are relative within each sample, part of this uniform rise may reflect a global immune signal rather than cell-type-specific expansion. The macrophage balance shifted towards M1: the M1-minus-M2 difference moved from about -0.18 to -0.04, driven by IL1B (log2FC 3.36), CXCL10 (2.72), PTGS2 (2.44), IL6 (2.12), CXCL9 (1.83) and CD86 (1.77), whereas CD163, MRC1, MSR1, IL10 and ARG1 did not pass the fold-change threshold and only CCL18 rose among M2 markers (Supplementary Table 1). The top ten up- and downregulated genes reached an in-sample AUC of 1.00 (Table 2; Supplementary Figs. 6 and 7), values explicitly internal to this 12-sample cohort and not used as a diagnostic claim. Whether the infiltrate is polarized and which cells produce the endothelial chemokines are questions that bulk averaging cannot answer.')
h2('Single-cell profiling resolves composition and gene addresses')
body('To resolve which cells generate these signals, 6,582 cells from sound/enamel-caries and deep-caries teeth were integrated into 12 cell types plus a mixed/low-quality class (Fig. 3a and Supplementary Fig. 8). Deep caries differed mainly in immune and progenitor compartments (Fig. 3b; Table 3): B cells appeared (0 to 192 cells, 5.2%), plasma cells expanded from 0.07% to 2.7%, T/NK cells from 7.4% to 10.0% and macrophages/monocytes from 6.1% to 7.0%, while mature fibroblasts decreased from 27.8% to 15.1% and fibroblast/MSC progenitors expanded from 5.3% to 14.2%. Marker localization assigned the bulk signals to compartments (Supplementary Fig. 8 and Supplementary Table 2): IL1B, PTGS2, CD86 and MMP9 were highest in macrophages/monocytes, but IL6, CXCL10, CXCL9, SELE, SOCS3 and ICAM1 in endothelial cells, and ALPL, RUNX2, VEGFB, FGFR2 and TNFSF11 in fibroblast/MSC progenitors. Two observations qualified the bulk interpretation: the M1-minus-M2 module score within macrophages did not differ between conditions (0.14 versus 0.22, P = 0.27; Supplementary Fig. 25), and endothelial chemokines exceeded the myeloid signal within the endothelial compartment, pointing to the vasculature as an active participant.')
h2('Endothelial autocrine chemokine signalling and a progenitor trajectory')
body('To test whether the endothelium signals actively, secreted single-subunit ligand-receptor pairs were scored with permutation testing (Supplementary Fig. 9). Deep caries enhanced endothelial autocrine CXCL2-ACKR1, CCL2-ACKR1, CXCL8-ACKR1 and CXCL10-ACKR1 signalling, and acquired PTN-SDC2 crosstalk between odontoblasts and fibroblast/MSC progenitors, SPP1-CD44 between SPP1+ neural cells and pericytes and MDK-NCL within the MSC compartment. Trajectory inference placed MSC-like cells at the origin of a lineage progressing through fibroblasts to odontoblast-like cells (mean pseudotime 1.7, 6.0 and 11.8), with deep-caries cells at earlier pseudotime (7.6 versus 8.7; Supplementary Fig. 10); a second lineage ran from the same origin through fibroblasts to DMP1+/IBSP+ odontoblasts (mean pseudotime 1.7, 6.5 and 21.3). NNLS deconvolution of the bulk cohort (Supplementary Fig. 11 and Supplementary Table 3) reproduced the immune expansion (macrophages 6.1% to 13.1%, P = 0.045; T/NK 2.9% to 7.4%, P = 0.008; plasma cells 0.2% to 6.8%, P = 0.005) and odontoblast loss (46.0% to 23.1%, P = 0.005), but the robustness analysis showed that macrophage and T/NK estimates are signature-size sensitive, so they are not used as evidence.')
h2('Exploratory co-expression and machine-learning analyses')
body('These analyses are exploratory: 12 samples cannot support formal inference and their discrimination estimates are internally optimistic. With that stated, WGCNA of the 6,000 most variable genes (soft power 5, scale-free R2 = 0.91) yielded five modules (Supplementary Figs. 12-15 and Supplementary Table 4) dominated by an immune module (2,402 genes, r = 0.906 with pulpitis; hubs IGHM, IGKC, ZFP36, CXCR4, LCP1, FPR1) and a developmental module (2,028 genes, r = -0.901; SOX9, DLX6, CDH12, DPT). Feature selection on the 986 differentially expressed genes converged on seven genes selected by at least two of three algorithms (ADGRL4, CCL21, FOSB, ICAM2, KIT, PI3, SELE; Supplementary Figs. 16 and 17 and Supplementary Table 5), and a logistic model separated groups with leave-one-out and five-fold AUC of 1.00 (Supplementary Fig. 18), which must not be read as diagnostic accuracy. The value of this layer is the compartment split it shares with the single-cell data: ADGRL4, SELE and ICAM2 are endothelial, FOSB and KIT are fibroblast/MSC-progenitor genes, CCL21 was only weakly detected (highest in endothelium) and PI3 was not detected in single cells.')
h2('PBM selectively reverses LPS-induced changes in panel targets')
body('To connect the compartment map to a candidate intervention, the 808 nm panel was re-analysed after correction of the sample-to-group assignment against GEO sample titles (Table 5 and Supplementary Figs. 19 and 20). Of 72 targets with complete data, PBM reversed LPS-induced changes in 32 (44%), not different from chance (binomial P = 0.41), indicating selective rather than global reversal. The strongest effect was VDR: LPS induced a DeltaCt of -10.07 (about 1,074-fold) and PBM returned it to control (DeltaCt +10.44), matching the original report that PBM decreased LPS-induced VDR by 37.66-fold [5]. MMP8 (-2.32 to +2.63) and VEGFB (-1.12 to +2.00) followed, whereas ALPL (-2.74 to -2.31) and TWIST1 (-2.05 to -0.80) remained induced and RUNX2 and TNFSF11 recovered. PBM therefore suppresses resorption- and inflammation-associated genes while preserving the odontogenic programme, the profile desired for vital pulp therapy. The overlap with pulpitis genes is partial (SPP1 and ITGA2 link the datasets; VDR and MMP8 are not differentially expressed in tissue) and hypothesis-generating.')
h2('An independent three-stage cohort reproduces the compartment map')
body('To test whether the map generalizes, 12 patients sampled at three stages were analysed (Fig. 4, Table 4 and Supplementary Figs. 21 and 22). Relative to healthy pulp, the endothelial compartment expanded monotonically (20.8% to 40.8% to 48.0%), a directionally consistent but not formally significant shift (Kruskal-Wallis P = 0.058), mirrored by the endothelial-to-stromal ratio (0.38, 1.29 and 3.59; P = 0.058). Endothelial cluster 3 rose from 3.3% to 7.4% and 11.7% (P = 0.030; advanced versus healthy Wilcoxon P = 0.029) and endothelial cluster 2 from 5.1% to 14.1% and 13.8% (P = 0.105), while fibroblast/stromal cluster 1 contracted from 56.9% to 25.9% and 18.9% (P = 0.055) and blood/immune cells increased from 5.5% to 7.3% and 10.1%. The four signature genes rose across stages (SELE 8.71 to 9.64, ICAM2 3.15 to 4.80, ADGRL4 8.41 to 8.87, FOSB 9.73 to 10.07 log2 CPM) and KIT fell (8.19 to 7.77), with the endothelial assignment reproduced in the independent cohort (Table 4). Across two independent cohorts, the seven-gene signature therefore resolves to the same endothelial and progenitor compartments.')
h2('Inflammatory programmes are compartment-specific and the venule is their interface')
body('To determine which compartment owns the inflammatory programme, patient-level pseudobulk differential expression and GSEA were computed per cluster (Fig. 5a and Supplementary Fig. 23). Seven clusters contained significantly changed genes, led by fibroblast/stromal cluster 1 (68 genes), glial cluster 2 (21) and endothelial cluster 3 (11). Endothelial cluster 1 was the dominant pathway responder with 29 enriched sets, including TNF-alpha signalling via NF-kB (NES 2.41, adjusted P = 1.8e-13), interferon-gamma (2.02), inflammatory response (2.08) and interferon-alpha (2.36); endothelial cluster 2 carried interferon-gamma (2.22) and interferon-alpha (2.34); blood/immune cluster 1 interferon-gamma (2.11); and mural cells interferon-gamma (1.97) and IL6-JAK-STAT3 (1.96). In the opposite direction, fibroblast/stromal cluster 1 suppressed TNF-alpha signalling (-2.23, 2.6e-8), KRAS signalling (-2.01) and inflammatory response (-1.86). The inflammatory programme is therefore the property of the vasculature and immune cells, and the contracting stroma actively downregulates it.')
body('Subtype annotation identified which endothelium carries the programme (Fig. 5b). Endothelial cluster 3, the steepest-expanding cluster, expressed the highest ACKR1 (0.43), SELE (1.95), ICAM1 (1.13) and CCL2 (0.59), the combination characteristic of inflammatory post-capillary venular endothelium; cluster 1 expressed SELE (1.53) and CXCL8 (0.44) as cytokine-activated endothelium; cluster 2 carried angiogenic receptors FLT1 (2.80), KDR (1.73) and EFNB2 (1.13); and cluster 4 the fenestrated marker PLVAP (0.93). The expanding compartment is therefore a specific ACKR1+ venular identity, the same compartment that in GSE185222 displayed autocrine CXCL2/CCL2/CXCL8-ACKR1 signalling. Within the immune compartment, M1 and M2 scores rose in parallel (M1 3.25 to 5.51, M2 3.36 to 5.62 log2 CPM) with no monotonic M1-minus-M2 shift (Kruskal-Wallis P = 0.926; Supplementary Fig. 24), so bulk M1 skewing reflects infiltration volume rather than a polarization switch.')
h2('Robustness and a negative metabolic control')
body('To delimit what survives perturbation, four checks were run (Supplementary Tables 6, 7 and 10). In 20 downsampling iterations the advanced-minus-healthy proportion difference kept its sign in all iterations for 11 of 12 cell types, with fibroblast/stromal contraction (-38.3 percentage points, 95% CI -39.1 to -37.5) and endothelial expansion stable throughout. Leave-one-patient-out analysis preserved the direction of every shift: the advanced-minus-healthy endothelial difference stayed positive after all 12 omissions (+16.5 to +32.8 percentage points) and endothelial cluster 3 likewise (+5.0 to +9.9), although the total-endothelium P value ranged from 0.036 to 0.122, so the direction is robust but significance is power-limited (Supplementary Table 10). The immune WGCNA module was reproduced in all five parameter variants (Jaccard 0.43-0.46, r = 0.91-0.92), whereas the developmental module had a stable core and a parameter-sensitive periphery (Jaccard 0.06-0.17). Deconvolution with the top 100 instead of top 200 specificity genes retained B cells, plasma cells, odontoblasts and fibroblasts but reversed macrophages, lost T/NK significance and reversed Schwann cells, so only directionally stable estimates are emphasized. As a negative control, cell-type-resolved metabolic scores did not change across stages for any of nine Hallmark metabolic sets (smallest P = 0.055); because oxidative phosphorylation is highest in odontoblast-lineage cells, which are progressively lost, the bulk suppression of oxidative phosphorylation is compositional rather than metabolic reprogramming.')
h2('Independent algorithms converge on the NF-kB and interferon-STAT1 axis')
body('To test whether the enriched pathways have identifiable drivers, transcription-factor (TF) and pathway activities were inferred with independent algorithms (Fig. 6 and Supplementary Figs. 26 and 27). In GSE77459 the NF-kB family dominated (NFKB1 +4.20, adjusted P = 0.008; RELA +3.63; REL +1.79), followed by myeloid master regulators SPI1 (+3.50), CEBPB (+2.58), CEBPA (+2.15), ETS1 (+2.39) and ETS2 (+2.10), and the AP-1 and STAT arms JUN (+2.77), STAT3 (+2.52) and FOS (+2.03) (Table 6). PROGENy independently reported higher EGFR (+1.87), MAPK (+1.80), JAK-STAT (+1.75), TNF-alpha (+1.72), NF-kB (+1.65) and PI3K (+1.40) activities and lower androgen-receptor activity (-1.42). In the three-stage cohort, STAT1 (14.43 to 15.98, P = 0.0498) and MYC (14.40 to 15.67, P = 0.0498) rose monotonically and significantly, and NFKB1, RELA, STAT3, IRF1, SPI1, CEBPB, JUN and FOS rose monotonically without reaching significance, consistent with an intensifying interferon-STAT1 axis superimposed on the NF-kB programme.')
h2('A drug-signature screen recovers the same axis')
body('As a knowledge-independent test, 4,027 drug perturbation signatures were screened directionally against the pulpitis transcriptome, of which 3,375 satisfied the gene-set size filter (Fig. 7 and Supplementary Table 8). The strongest opposing signatures were rimexolone, a glucocorticoid (NES 2.60, adjusted P = 8.6e-11), parthenolide, an NF-kB inhibitor (2.43, 1.9e-12), CP-690334-01, the JAK inhibitor tofacitinib (2.30, 8.1e-07), progesterone (2.38), valproic acid (2.33) and geldanamycin (2.11). The aligned control list contained prostaglandins alprostadil (PGE1) and dinoprost (PGF2-alpha) and several HL60 myeloid-differentiation signatures (tretinoin, etoposide, mebendazole), whose overlap with the pulpitis programme reflects shared myeloid genes rather than a shared mechanism; the aligned list is used only as a specificity control and is not interpreted as evidence that these drugs aggravate pulpitis. That a screen containing no pathway information recovers NF-kB and JAK-STAT inhibitors and a glucocorticoid is a convergence result, not a therapeutic recommendation.')
h2('Communication shifts to the venule and the trajectory orders the events')
body('To place the compartment map in time and intercellular space, communication was quantified across stages and pseudotime fitted independently of grading (Fig. 8, Supplementary Fig. 28 and Supplementary Tables 9 and 11). Total communication strength declined with progression (4,844, 3,842 and 3,270), mirroring stromal contraction. Patient-level analysis with 10,000 label permutations confirmed that fibroblast-to-endothelial-3 signalling gained share beyond chance (0.0050, 0.0105 and 0.0113; permutation P = 0.0048, FDR = 0.038), whereas fibroblast-to-immune signalling did not change significantly (P = 0.19) and immune-to-endothelial-1 signalling lost share (P = 0.039) (Supplementary Table 11). Two axes linking endothelial cluster 3 and immune cells changed direction between pooled and patient-level estimates and are therefore not interpreted. A slingshot trajectory recovered a single lineage that strongly separated stages (Kruskal-Wallis P = 1.1e-189; median pseudotime 7.25, 14.68 and 15.28) while early and advanced cells overlapped along the continuum (Spearman rho = 0.128). Along the trajectory SELE and ACKR1 peaked in middle bins (z 1.80 and 2.51), CCL2, CXCL8, IL1B, CCL5, SPP1 and NFKB1 culminated in terminal bins, and COL1A1 and DCN declined overall, with a transient rebound of COL1A1 in the fifth pseudotime bin and of DCN in the eighth bin. The order of events is therefore stromal contraction, venule activation and inflammatory amplification.')

# ================= Discussion =================
h1('Discussion')
body('We asked one question in five layers and obtained one answer: irreversible pulpitis is not a uniformly inflamed tissue but a compartment-specific state organized around an ACKR1+ inflammatory venular endothelium, driven by NF-kB and interferon-STAT1 programmes, with active suppression of inflammatory programmes in the contracting stroma and loss of the odontogenic developmental programme. Against this map, 808 nm PBM selectively suppresses LPS-induced VDR, MMP8 and VEGFB while preserving the odontogenic programme in DPSCs.')
body('The layered design answers the standard objection that in silico studies accumulate methods without accumulating evidence. Each layer resolved what the previous layer could not: bulk averaging could not say which cells transcribe the chemokines; single-cell profiling of four samples could not say whether the pattern generalizes; the three-stage cohort could not say which regulators drive it; and no single algorithm could establish that regulators, pathways and drug signatures agree. Analyses that failed this test were reported as negative or downgraded rather than kept for completeness: cell-type-resolved metabolic scoring showed no change (smallest P = 0.055) and serves only to explain that bulk oxidative-phosphorylation suppression is compositional; PROGENy in the validation cohort added nothing and is not interpreted; deconvolution estimates that reversed under signature perturbation are excluded; and WGCNA and machine learning are labelled exploratory because 12 samples cannot support them. The seven-gene signature is therefore a compartment-localization result, not a diagnostic claim.')
body('The endothelial finding is the strongest convergence. Bulk enrichment, cell-type-resolved GSEA, subtype annotation, autocrine communication scoring and the stage-resolved ligand-receptor analysis independently place TNF-alpha/NF-kB and interferon programmes in the vasculature and identify the expanding compartment as ACKR1+ post-capillary venular endothelium expressing E-selectin, ICAM1 and CCL2. ACKR1 on venular endothelium presents chemokines to rolling leukocytes and modulates transmigration, so its expansion provides a concrete cellular route through which the bulk chemokine signal is generated and explains why the communication centre of gravity shifts towards fibroblast-to-venule signalling as disease advances. The suppression of inflammatory programmes in the contracting stroma is equally informative: it shows that the inflammatory transcriptome belongs to specific compartments rather than to the tissue as a whole, exactly the distinction bulk studies cannot make.')
body('These findings extend, rather than duplicate, the two studies that generated the single-cell data. Opasawatchai et al. annotated pulpal endothelium by SELE and PECAM1 and reported a CCL2-ACKR1 interaction between macrophages and endothelial cells in deep caries [4], and Ha et al. resolved venular, arteriolar and capillary endothelial subsets and described vascular remodelling across decay stages [6]. Neither study asked which compartment owns the tissue-level inflammatory programme. The present analysis addresses that question: the programmes detected in bulk pulpitis transcriptomes are carried predominantly by the vasculature and immune cells while the contracting stroma suppresses them, the assignment is reproduced across two independent cohorts and platforms, and two common bulk-level interpretations, namely M1 polarisation and metabolic suppression, are shown to be compositional rather than regulatory. The ACKR1+ venule is therefore not claimed as a newly discovered cell type, but as the compartment that integrates these independently reproduced observations.')
body('Two apparent discrepancies are resolved rather than left standing. First, the bulk M1 skewing of the discovery cohort is not reproduced within the immune compartment of the validation cohort, where M1 and M2 rise in parallel and the M1-minus-M2 difference does not shift monotonically; the bulk signal is infiltration volume, not polarization. Second, the loss of oxidative-phosphorylation and cilium-assembly programmes in bulk is compositional and developmental rather than a metabolic switch within surviving cells. Both resolutions follow from cell-type resolution converting bulk averages into statements about compartments.')
body('The regulatory and pharmacological layers reinforce the same axis. NF-kB family and myeloid master regulators drive the bulk programme, interferon-STAT1 activity intensifies along the gradient, PROGENy independently reports higher JAK-STAT, TNF-alpha and NF-kB activity, and a screen of 3,375 drug signatures containing no pathway information recovers NF-kB and JAK-STAT inhibitors and a glucocorticoid as opposing agents. The aligned drug list is deliberately not over-interpreted: it is dominated by HL60 myeloid-differentiation signatures whose overlap reflects shared myeloid genes, and serves only as a specificity control. For vital pulp therapy the axis is nonetheless actionable, because NF-kB and JAK-STAT are druggable and PBM has been reported to modulate NF-kB signalling in inflamed dental tissues [5], suggesting a convergence point between pharmacological and photonic intervention.')
body('The PBM re-analysis adds a candidate intervention. After correction of the sample assignment against GEO sample titles, the panel reproduces the original study and quantifies its selectivity: 44% of genes reversed, with VDR showing the most dramatic effect, a roughly 1,074-fold LPS induction returned to baseline. The suppression of MMP8 and VEGFB and the preservation of ALPL and TWIST1 sketch a profile in which PBM damps resorption- and inflammation-associated genes without ablating the odontogenic differentiation programme. Because responsive genes localize to fibroblast/MSC progenitors and odontoblast-like cells, two compartments remodelled in pulpitis, the panel provides a cellular rationale for testing PBM while progenitor compartments are still preserved.')
body('Limitations are stated at the same resolution as the findings. The bulk cohort has 12 samples and the single-cell discovery cohort 4; both are too small for the exploratory network and machine-learning analyses to support formal inference, and the perfect cross-validated discrimination reflects complete separation of two homogeneous groups and must not be read as diagnostic performance. The validation cohort has four patients per stage, one healthy sample contributes 52% of its cells, and the deposited object provides only superclass-level annotation, so its comparisons are directional and underpowered and it cannot serve as a classification test set. The single-cell data have limited depth, so MMP8, ARG1, PI3 and CCL21 were poorly or not detected, and no B cells were annotated in the sound/enamel-caries group, a difference that may partly reflect sampling and sequencing depth rather than true absence. Communication scores are co-expression approximations rather than measured protein interactions, the trajectory is inferred from anchored labels, and the pseudotime ordering is descriptive. The PBM panel has two replicates per group at a single dose and time point, and its supplemented group labels conflict with the GEO titles, so the analysis relies exclusively on deposited Ct values and GEO annotation. The robustness checks address only specific perturbations (composition downsampling, leave-one-patient-out, WGCNA parameter variants and deconvolution signature size); other analytical choices remain untested. Total endothelial expansion is a directional trend (P = 0.058) whose significance varies with patient omission (up to P = 0.122), so the venule conclusion rests on endothelial cluster 3 and on the permutation-supported communication gain rather than on the total compartment shift. Most importantly, no experimental validation is presented: every statement is a hypothesis generated from public data, and the drug-signature, regulator and communication results are computational predictions rather than measured activities.')

# ================= Methods =================
h1('Methods')
h2('Analytical strategy')
body('Analyses were arranged in five layers: (1) bulk differential expression, enrichment, ssGSEA and ROC define the inflammatory programme; (2) single-cell profiling and deconvolution identify the owning compartment; (3) cell-type-resolved pseudobulk analysis in the three-stage cohort tests generalization; (4) transcription-factor, pathway, drug-signature and communication inference test convergence; (5) parameter perturbation delimits what survives. Analyses that changed no conclusion are reported as negative controls, and network and machine-learning sections are labelled exploratory.')
h2('Data sources')
body('Four public datasets were used (Table 1) [6, 7]. GSE77459 (GPL17692) contains six healthy and six irreversible pulpitis samples [1]. GSE185222 provides scRNA-seq of pulp from sound/enamel-caries and deep-caries teeth (4 samples, 6,582 cells) [4]. GSE193448 provides an RT-PCR panel of DPSCs treated with LPS or LPS plus 808 nm laser irradiation [5]; raw Ct values were taken from the supplementary file. E-MTAB-14896 provides scRNA-seq of pulp from 12 patients across healthy, early and advanced stages (53,407 cells; Parse Evercode WT v2) [6]; the processed Seurat object was used as deposited, and patient stage was taken from its Health_state metadata (4 healthy, 4 early, 4 advanced). The eight ArrayExpress samples are sequencing sub-libraries of the 12 multiplexed patients (24 assays), not eight patients.')
h2('Bulk microarray processing')
body('The GSE77459 series matrix contains RMA-normalized log2 values; probes were mapped to gene symbols (27,279 genes). Differential expression used limma [8] with |log2FC| > 1 and Benjamini-Hochberg adjusted P < 0.05. GO biological process and KEGG over-representation enrichment of upregulated genes used clusterProfiler [9] with adjusted P < 0.05.')
h2('Immune infiltration and in-sample discrimination')
body('ssGSEA implemented in GSVA [10] scored 12 immune cell types from a curated marker set; group differences used the Wilcoxon rank-sum test with FDR correction. ROC analysis of candidate genes was performed within the same 12 samples and is reported strictly as in-sample discrimination.')
h2('Genome-wide gene set enrichment analysis')
body('All genes were ranked by limma log2 fold change and analysed with preranked GSEA [11] using the Hallmark collection (MSigDB 2024.1.Hs) [12] and GO biological process gene sets derived from org.Hs.eg.db and GO.db. Gene sets were restricted to 10-500 (Hallmark) or 15-500 (GO) genes, with exact P values and Benjamini-Hochberg adjustment.')
h2('Exploratory co-expression and machine-learning analyses')
body('Two analyses are exploratory because 12 samples cannot support formal inference. First, genes with mean log2 expression above 2.5 were retained and the 6,000 most variable were used to build a signed co-expression network with WGCNA [13]; the soft-thresholding power was the smallest value reaching scale-free topology R2 > 0.90, modules were detected with blockwiseModules (minimum size 30, merge height 0.30, deepSplit 3, signed TOM), and hubs were ranked by module membership. Second, the 986 differentially expressed genes formed a candidate pool: LASSO with leave-one-out cross-validation [14], random forest with 1,000 trees and SVM-RFE with a linear kernel [15]; genes selected by at least two methods formed a consensus set, a logistic model was fitted, and performance was assessed by leave-one-out and by 50 repetitions of stratified five-fold cross-validation with selection repeated inside every training fold. All discrimination estimates are internal and explicitly optimistic.')
h2('Single-cell RNA-seq processing')
body('Four GSE185222 matrices were merged on common genes into a Seurat object [16]. Highly variable genes (2,000), scaling and 30 principal components were computed, batch effects were removed with Harmony, and Harmony embeddings were used for graph construction (20 dimensions), Louvain clustering (resolution 0.5) and UMAP. Clusters were annotated into 12 cell types plus a mixed/low-quality class using canonical markers: fibroblasts (THY1, PDGFRA, COL1A1), fibroblast/MSC progenitors (THY1, PDGFRA, MCAM), odontoblasts (DMP1, IBSP), odontoblast-like cells (OMD), endothelial cells (PECAM1, SELE), pericytes/SMC (ACTA2, PDGFRB), Schwann cells (MPZ, S100B), SPP1+ neural cells (SPP1), T/NK cells (CD3D, NKG7), B cells (CD79A, MS4A1), plasma cells (MZB1, JCHAIN) and macrophages/monocytes (CD68, LYZ, CSF1R).')
h2('Composition, localization, communication and deconvolution')
body('Proportions were computed per condition and mean marker expression per cell type from log-normalized data. M1 and M2 module scores of macrophages used AddModuleScore. Secreted single-subunit ligand-receptor pairs were extracted from CellChatDB.human [17]; sender-receiver scores were the geometric mean of mean ligand expression in the sender and mean receptor expression in the receiver, and significance of the top 300 pairs in deep caries was assessed by 200 permutations of cell-type labels. A signature matrix was built from the top 200 specificity genes per cell type and bulk expression deconvolved sample by sample with non-negative least squares [18]. Slingshot [19] trajectory inference used Harmony dimensions with the MSC-like cluster as the starting population.')
h2('PBM panel re-analysis')
body('Ct values were taken from the supplementary file for 84 target genes, five housekeeping genes and one DNA control (90 wells). GEO sample titles (GSM5813588-GSM5813591) annotate samples 3-4 as LPS-induced and samples 5-6 as LPS plus laser-irradiated DPSCs, consistent with the original publication [5], so the supplemented Group 1/Group 2 labels were reassigned. Quality-control wells were removed and 72 genes passed replicate QC; the LPS effect was DeltaCt = Ct(LPS) - Ct(control) and PBM effect DeltaCt = Ct(LPS + laser) - Ct(LPS), positive values indicating downregulation. Genes were called reversed when the two effects opposed, and the reversal rate was compared with 50% by an exact binomial test. NFKB1 was excluded because its LPS replicates were discordant.')
h2('Three-stage validation cohort')
body('The deposited Seurat object (60,721 genes, 53,407 cells) was used with its original annotation, resolving five superclasses into 12 clusters. Patients were grouped by the deposited Health_state field into healthy (4 patients, 30,786 cells), early (4, 12,573) and advanced (4, 10,048). Patient-level pseudobulk profiles were built by summing raw counts and converting to log2 counts per million, and signature values, M1/M2 module scores and per-patient proportions were compared by Kruskal-Wallis test. Localization used mean log-normalized expression per cluster. One healthy sample contributed 27,810 cells (52%), so comparisons rest on four patients per group and are interpreted as directional.')
h2('Cell-type-resolved differential expression and GSEA')
body('For each cluster, patient-level pseudobulk profiles were built by summing raw counts, retaining genes with at least 5 counts in at least 3 patients and patients with at least 10 cells, and requiring at least two patients per stage. Differential expression between advanced and healthy stages used DESeq2 with poscounts size factors. Hallmark enrichment used preranked GSEA with fgsea (fgseaMultilevel, gene sets 10-500) on DESeq2 statistics collapsed to symbols.')
h2('Endothelial subtype annotation and immune M1/M2 scoring')
body('The four endothelial clusters were annotated with vessel-type markers: inflammatory venule (ACKR1, SELE, VCAM1, ICAM1, CCL2, CXCL8), capillary (CA4, RGCC, VWA1, PLVAP), arteriolar (GJA5, HEY1, SEMA3G, EFNB2, FBLN5), angiogenic tip (ESM1, KDR, DLL4, FLT1, PGF) and lymphatic (PROX1, LYVE1, PDPN, CCL21). Within the immune compartment, M1 and M2 scores were the mean z-scored expression of each marker set at cell level and mean log2 CPM at patient level, compared by Wilcoxon or Kruskal-Wallis tests.')
h2('Robustness analyses')
body('Four checks were performed. Composition stability: each patient was downsampled to 1,000 cells, the advanced-minus-healthy difference recomputed 20 times per cell type, recording the fraction of iterations sharing its sign. Leave-one-patient-out analysis repeated every proportion test after omitting each patient in turn. Co-expression stability: WGCNA was repeated under five parameter variants (power 5 or 7; minimum module size 30 or 50; deepSplit 2 or 3; 5,500 or 6,000 most variable genes), matched to the original modules by Jaccard overlap. Deconvolution stability: the NNLS signature was rebuilt with the top 100 instead of top 200 specificity genes.')
h2('Transcription-factor and pathway activity inference')
body('TF activities were inferred with decoupleR [20] using the univariate linear model and DoRothEA regulons restricted to confidence levels A and B [21] (118 TFs, 6,408 interactions, minimum regulon size 5), in GSE77459 (Wilcoxon with Benjamini-Hochberg correction) and in patient-level pseudobulk of E-MTAB-14896 (Kruskal-Wallis). PROGENy [22] inferred the activities of 14 signalling pathways (top 100 responsive genes, scaled) in GSE77459.')
h2('Direction-aware drug-signature reversal screen')
body('Drug perturbation signatures were obtained from DSigDB [23] (4,027 gene sets; 3,375 passed the 15-500 size filter). Preranked GSEA (gene sets 15-500) used the pulpitis fold-change ranking, and signatures were interpreted directionally: a drug whose down-regulated set was enriched at the pulpitis up-regulated end, or whose up-regulated set was enriched at the pulpitis down-regulated end, was classified as opposing; the converse as aligned. Only opposing drugs with adjusted P < 0.05 were interpreted, and aligned signatures served solely as a specificity control.')
h2('Ligand-receptor communication and pseudotime in the validation cohort')
body('Secreted single-subunit pairs were taken from CellChatDB v1 human [17] (1,199 pairs). Mean log2 CPM of each ligand and receptor was computed per cell type and stage; pairwise strength was the outer product of sender ligand and receiver receptor means, summed over pairs and normalized by stage total. The same calculation was repeated within each patient, and stage labels were permuted 10,000 times preserving group sizes to obtain empirical two-sided P values for advanced-minus-healthy share differences on eight prespecified axes, with Benjamini-Hochberg correction. Slingshot [19] used the embedded UMAP with Health_state classes as labels and Healthy as start; pseudotime was compared across stages by Kruskal-Wallis and gene dynamics summarized as z-scored means in ten bins.')
h2('Statistical analysis')
body('Analyses were performed in R 4.5.3 [24]; two-sided P < 0.05 was considered significant, with Benjamini-Hochberg correction where applicable.')

# ================= Data availability / statements =================
h1('Ethics approval and consent to participate')
body('This study analysed only publicly available, de-identified human transcriptomic datasets deposited in GEO and ArrayExpress. No new human tissue samples or animal experiments were generated, and no ethical approval or participant consent was therefore required.')
h1('Consent for publication')
body('Not applicable.')
h1('Supplementary information')
body('The online version contains supplementary material comprising Supplementary Figs. 1-28 and Supplementary Tables 1-11.')
h1('Data availability')
body('All datasets analysed are publicly available: GSE77459 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77459), GSE185222 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185222), GSE193448 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193448) and E-MTAB-14896 (https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14896).')
h1('Code availability')
body('The R scripts used for data processing, statistical analysis and figure generation are available from the corresponding author on reasonable request.')
h1('Acknowledgements')
body('Not applicable.')
h1('Author contributions')
body('Z.R. conceived and designed the study, performed all computational analyses, interpreted the results and wrote the manuscript.')
h1('Competing interests')
body('The author declares no competing interests.')
h1('Funding')
body('This research received no specific grant from any funding agency in the public, commercial or not-for-profit sectors.')

# ================= References (Nature style) =================
h1('References')
refs = [
    'Galicia, J. C., Henson, B. R., Parker, J. S. & Khan, A. A. Gene expression profile of pulpitis. Genes Immun. 17, 239–243 (2016).',
    'Huang, X. & Chen, K. Differential expression of long noncoding RNAs in normal and inflamed human dental pulp. J. Endod. 44, 62–72 (2018).',
    'Chen, M., Zeng, J., Yang, Y. & Wu, B. Diagnostic biomarker candidates for pulpitis revealed by bioinformatics analysis of merged microarray gene expression datasets. BMC Oral Health 20, 279 (2020).',
    'Opasawatchai, A. et al. Single-cell transcriptomic profiling of human dental pulp in sound and carious teeth: a pilot study. Front. Dent. Med. 2, 806294 (2021).',
    'Rocha, E. A. et al. Laser photobiomodulation 808 nm: effects on gene expression in inflammatory and osteogenic biomarkers in human dental pulp stem cells. Front. Pharmacol. 12, 782095 (2022).',
    'Ha, H. T. et al. Human atlas of tooth decay progression: identification of cellular mechanisms driving the switch from dental pulp repair toward irreversible pulpitis. Adv. Sci. 13, e10096 (2026).',
    'Barrett, T. et al. NCBI GEO: archive for functional genomics data sets—update. Nucleic Acids Res. 41, D991–D995 (2013).',
    'Ritchie, M. E. et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res. 43, e47 (2015).',
    'Wu, T. et al. clusterProfiler 4.0: a universal enrichment tool for interpreting omics data. Innovation 2, 100141 (2021).',
    'Hänzelmann, S., Castelo, R. & Guinney, J. GSVA: gene set variation analysis for microarray and RNA-seq data. BMC Bioinformatics 14, 7 (2013).',
    'Subramanian, A. et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proc. Natl. Acad. Sci. USA 102, 15545–15550 (2005).',
    'Liberzon, A. et al. The Molecular Signatures Database (MSigDB) hallmark gene set collection. Cell Syst. 1, 417–425 (2015).',
    'Langfelder, P. & Horvath, S. WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 9, 559 (2008).',
    'Friedman, J., Hastie, T. & Tibshirani, R. Regularization paths for generalized linear models via coordinate descent. J. Stat. Softw. 33, 1–22 (2010).',
    'Guyon, I., Weston, J., Barnhill, S. & Vapnik, V. Gene selection for cancer classification using support vector machines. Mach. Learn. 46, 389–422 (2002).',
    'Hao, Y. et al. Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nat. Biotechnol. 42, 293–304 (2024).',
    'Jin, S. et al. Inference and analysis of cell-cell communication using CellChat. Nat. Commun. 12, 1088 (2021).',
    'Mullen, K. M. & van Stokkum, I. H. M. nnls: the Lawson-Hanson algorithm for non-negative least squares (NNLS). R package version 1.6 (2012).',
    'Street, K. et al. Slingshot: cell lineage and pseudotime inference for single-cell transcriptomics. BMC Genomics 19, 477 (2018).',
    'Badia-i-Mompel, P. et al. decoupleR: ensemble of computational methods to infer biological activities from omics data. Bioinform. Adv. 2, vbac016 (2022).',
    'Garcia-Alonso, L. et al. Benchmark and integration of resources for the estimation of human transcription factor activities. Genome Res. 29, 1363–1375 (2019).',
    'Schubert, M. et al. Perturbation-response genes reveal signaling footprints in cancer gene expression. Nat. Commun. 9, 20 (2018).',
    'Yoo, M. et al. DSigDB: drug signatures database for gene set analysis. Bioinformatics 31, 3069–3071 (2015).',
    'R Core Team. R: a language and environment for statistical computing. R Foundation for Statistical Computing (2025).',
]
for i, r in enumerate(refs, 1):
    ref(f'{i}. {r}')

# ================= Figure legends =================
h1('Figure legends')
body('Figure 1. The bulk pulpitis transcriptome. (a) Volcano plot of differential expression between pulpitis and normal pulp (GSE77459; |log2FC| > 1, adjusted P < 0.05). (b) Preranked GSEA of Hallmark gene sets, showing activated inflammatory programmes and suppressed oxidative-phosphorylation, translation and cilium-assembly programmes.')
body('Figure 2. Immune infiltration of the inflamed pulp. Heatmap of ssGSEA scores for 12 immune cell types across individual samples (GSE77459).')
body('Figure 3. Single-cell atlas of human dental pulp. (a) UMAP of 6,582 cells coloured by cluster (left) and by sample (right) (GSE185222). (b) Cell-type composition of the two deep-caries and two sound/enamel-caries samples; the mixed/low-quality class is displayed but excluded from compartment comparisons.')
body('Figure 4. Compartment shifts in the independent three-stage cohort. (a) Cell-type composition in healthy, early-stage and advanced-stage pulp. (b) Patient-mean percentage of endothelial clusters 1-4 across stages, with individual patient values overlaid, showing the steepest expansion of cluster 3 (E-MTAB-14896).')
body('Figure 5. Compartment-specific inflammatory programmes and the inflammatory venule. (a) Hallmark GSEA normalized enrichment scores per annotated cell type (advanced versus healthy, patient-level pseudobulk). (b) Vessel-subtype marker expression across the four endothelial clusters, identifying cluster 3 as ACKR1+ venular endothelium (E-MTAB-14896).')
body('Figure 6. Convergent regulatory inference. (a) Transcription-factor activity differences between pulpitis and normal pulp (DoRothEA A+B regulons, decoupleR ULM; GSE77459). (b) PROGENy pathway activities in pulpitis versus normal pulp.')
body('Figure 7. Drug perturbation signatures opposing or aligned with the pulpitis transcriptome (DSigDB, direction-aware preranked GSEA). A signature is called opposing when its down-regulated set is enriched at the pulpitis up-regulated end (or its up-regulated set at the pulpitis down-regulated end), and aligned when same-direction gene sets coincide. Aligned signatures are shown as a specificity control only.')
body('Figure 8. Communication and event ordering in the three-stage cohort. (a) Normalized ligand-receptor communication shares on key cell-type axes across stages. (b) UMAP colored by slingshot pseudotime. (c) Z-scored dynamics of selected genes in ten pseudotime bins (E-MTAB-14896).')

# ================= Tables =================
table_caption('Table 1. Public datasets used in this study.')
add_table(['Dataset', 'Platform', 'Samples', 'Cells / features', 'Role'],
          [['GSE77459', 'Affymetrix Human Gene 2.1 ST (GPL17692)', '6 healthy pulp, 6 irreversible pulpitis', '27,279 genes', 'Bulk discovery'],
           ['GSE185222', 'Illumina HiSeq 4000', '2 sound/enamel-caries, 2 deep-caries teeth', '6,582 cells', 'Single-cell discovery'],
           ['GSE193448', 'RT-PCR panel, 808 nm PBM', 'Control, LPS, LPS + laser (2 wells each)', '90 wells (84 targets plus 6 controls)', 'PBM re-analysis'],
           ['E-MTAB-14896', 'Parse Evercode WT v2, NovaSeq X', '4 healthy, 4 early, 4 advanced patients', '53,407 cells, 60,721 genes', 'Three-stage validation']],
          ['l', 'l', 'l', 'c', 'l'])

table_caption('Table 2. Top candidate genes in bulk pulpitis (in-sample AUC = 1.00 for all).')
add_table(['Upregulated', 'log2FC', 'Downregulated', 'log2FC'],
          [['PI3', '4.74', 'NR1D1', '-2.47'],
           ['CCL21', '4.79', 'GAP43', '-2.06'],
           ['SELE', '3.34', 'CRISPLD1', '-1.95'],
           ['ACP5', '2.99', 'SPON1', '-1.93'],
           ['SOCS3', '2.77', 'MOXD1', '-1.74'],
           ['CH25H', '2.67', 'DBP', '-1.67'],
           ['ATP6V0D2', '2.55', 'NR1D2', '-1.59'],
           ['PTGS2', '2.44', 'RSPO2', '-1.57'],
           ['FOSB', '2.41', 'ECM2', '-1.26'],
           ['SLC37A2', '1.80', 'RN7SK', '-1.26']],
          ['l', 'c', 'l', 'c'])

table_caption('Table 3. Cell-type composition of human dental pulp (GSE185222). Percentages are rounded and may not sum to exactly 100.')
add_table(['Cell type', 'Sound/enamel caries (cells)', 'Percent', 'Deep caries (cells)', 'Percent'],
          [['Fibroblasts/MSC', '155', '5.3', '522', '14.2'],
           ['Fibroblasts', '808', '27.8', '556', '15.1'],
           ['Odontoblast-like (OMD+)', '900', '31.0', '939', '25.5'],
           ['Odontoblasts (DMP1+/IBSP+)', '22', '0.8', '142', '3.9'],
           ['Endothelial', '318', '10.9', '341', '9.3'],
           ['Pericytes/SMC', '219', '7.5', '138', '3.8'],
           ['Macrophages/Monocytes', '176', '6.1', '256', '7.0'],
           ['T/NK cells', '215', '7.4', '368', '10.0'],
           ['B cells', '0', '0', '192', '5.2'],
           ['Plasma cells', '2', '0.1', '99', '2.7'],
           ['Neural/Schwann (SPP1+)', '60', '2.1', '60', '1.6'],
           ['Schwann cells (MPZ+)', '29', '1.0', '39', '1.1'],
           ['Mixed/low-quality', '1', '0.0', '25', '0.7']],
          ['l', 'c', 'c', 'c', 'c'])

table_caption('Table 4. Compartments and signature genes in the three-stage cohort (E-MTAB-14896). P values are Kruskal-Wallis tests across stages on patient-level values.')
add_table(['Cell type', 'Healthy (%)', 'Early (%)', 'Advanced (%)', 'P', 'Signature genes maximal here'],
          [['Endothelial, all', '20.80', '40.77', '48.04', '0.058', 'SELE, ICAM2, ADGRL4'],
           ['Endothelial/stromal ratio', '0.38', '1.29', '3.59', '0.058', '-'],
           ['Endothelial 3', '3.27', '7.42', '11.70', '0.030', 'CCL2, ACKR1, CXCL10, IL6'],
           ['Endothelial 2', '5.09', '14.06', '13.82', '0.105', 'ADGRL4, ICAM2'],
           ['Fibroblast/stromal 1', '56.93', '25.94', '18.87', '0.055', 'ALPL, RUNX2, VDR'],
           ['Fibroblast/stromal 2', '6.46', '9.10', '9.82', '0.472', 'KIT, FGFR2, TWIST1'],
           ['Blood/immune, all', '5.54', '7.31', '10.07', '-', 'CXCL2, CXCL8, IL1B, CD163'],
           ['Blood/immune 1', '2.89', '4.19', '7.17', '0.298', 'FOSB, PTGS2, MRC1'],
           ['Mural', '2.96', '6.45', '6.33', '0.292', '-'],
           ['Glial 1', '1.01', '6.56', '2.96', '0.026', '-'],
           ['Glial 2', '0.80', '2.83', '2.45', '0.031', 'SPP1']],
          ['l', 'c', 'c', 'c', 'c', 'l'])

table_caption('Table 5. Representative genes of the 808 nm PBM panel (GSE193448) under the corrected assignment. Ct values are rounded to two decimals; DeltaCt values were computed from unrounded data.')
add_table(['Gene', 'Ct Control', 'Ct LPS', 'Ct LPS+Laser', 'DeltaCt LPS', 'DeltaCt PBM', 'Pattern'],
          [['VDR', '32.27', '22.20', '32.64', '-10.07', '+10.44', 'Complete suppression'],
           ['MMP8', '36.92', '34.60', '37.23', '-2.32', '+2.63', 'Suppressed'],
           ['VEGFB', '30.57', '29.45', '31.45', '-1.12', '+2.00', 'Suppressed'],
           ['BMP4', '31.22', '30.89', '31.59', '-0.33', '+0.70', 'Mildly suppressed'],
           ['MMP10', '29.19', '30.34', '28.29', '+1.16', '-2.06', 'Restored'],
           ['RUNX2', '29.77', '30.31', '29.17', '+0.55', '-1.15', 'Restored'],
           ['TNFSF11', '32.60', '33.66', '33.10', '+1.07', '-0.56', 'Partially restored'],
           ['CD36', '30.55', '34.35', '33.57', '+3.81', '-0.79', 'Partially restored'],
           ['ALPL', '36.51', '33.77', '31.45', '-2.74', '-2.31', 'Further enhanced'],
           ['TWIST1', '33.69', '31.64', '30.84', '-2.05', '-0.80', 'Maintained']],
          ['l', 'c', 'c', 'c', 'c', 'c', 'l'])

table_caption('Table 6. Transcription-factor activities (DoRothEA A+B, decoupleR ULM).')
add_table(['TF', 'Bulk difference (pulpitis minus normal)', 'Adjusted P', 'Stage trend (E-MTAB-14896)', 'Stage P'],
          [['NFKB1', '+4.20', '0.008', '7.78 to 8.76', '0.211'],
           ['RELA', '+3.63', '0.008', '7.20 to 8.00', '0.397'],
           ['SPI1', '+3.50', '0.008', '6.84 to 7.34', '0.397'],
           ['JUN', '+2.77', '0.008', '7.08 to 7.72', '0.276'],
           ['CEBPB', '+2.58', '0.008', '2.90 to 3.51', '0.211'],
           ['STAT3', '+2.52', '0.008', '6.76 to 7.86', '0.211'],
           ['ETS1', '+2.39', '0.008', '-', '-'],
           ['CEBPA', '+2.15', '0.008', '-', '-'],
           ['ETS2', '+2.10', '0.008', '-', '-'],
           ['FOS', '+2.03', '0.008', '5.65 to 6.31', '0.584'],
           ['STAT1', '-', '-', '14.43 to 15.98', '0.0498'],
           ['MYC', '-', '-', '14.40 to 15.67', '0.0498']],
          ['l', 'c', 'c', 'c', 'c'])

# ================= Footer page numbers (Arabic, centred; Scientific Reports) =================
footer = sec.footer
fp = footer.paragraphs[0]
fp.alignment = WD_ALIGN_PARAGRAPH.CENTER
fp.paragraph_format.line_spacing = 1.0
fr = fp.add_run()
fr.font.name = 'Times New Roman'
fr.font.size = Pt(10.5)
fb = OxmlElement('w:fldChar'); fb.set(qn('w:fldCharType'), 'begin')
fi = OxmlElement('w:instrText'); fi.set('{http://www.w3.org/XML/1998/namespace}space', 'preserve'); fi.text = ' PAGE '
fsep = OxmlElement('w:fldChar'); fsep.set(qn('w:fldCharType'), 'separate')
fe = OxmlElement('w:fldChar'); fe.set(qn('w:fldCharType'), 'end')
fr._r.append(fb); fr._r.append(fi); fr._r.append(fsep); fr._r.append(fe)

# continuous line numbering (Scientific Reports advises including it to help reviewers)
sectPr = doc.sections[0]._sectPr
for old in sectPr.findall(qn('w:lnNumType')):
    sectPr.remove(old)
ln = OxmlElement('w:lnNumType')
ln.set(qn('w:countBy'), '1')
ln.set(qn('w:restart'), 'continuous')
pgmar = sectPr.find(qn('w:pgMar'))
if pgmar is not None:
    pgmar.addnext(ln)
else:
    sectPr.append(ln)

doc.save(OUT)
print('SAVED:', OUT)
