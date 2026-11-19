# -*- coding: utf-8 -*-
"""Supplementary Information in Scientific Reports style: Supplementary Fig. 1-28 and Supplementary Table 1-9."""
from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_ALIGN_VERTICAL
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

OUT = 'C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/results/Supplementary_Information_ScientificReports.docx'

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
st.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY

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


def title(text):
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_after = Pt(18); p.paragraph_format.space_before = Pt(6)
    p.paragraph_format.keep_with_next = True
    r = p.add_run(text); r.font.name = 'Times New Roman'; r.font.size = Pt(18); r.font.bold = True


def body(text):
    return doc.add_paragraph(text)


def h1(text):
    return doc.add_heading(text, level=1)


def h2(text):
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(12); p.paragraph_format.space_after = Pt(3)
    p.paragraph_format.keep_with_next = True; p.paragraph_format.line_spacing = 2.0
    r = p.add_run(text); r.font.name = 'Times New Roman'; r.font.size = Pt(12); r.font.bold = True


def table_caption(text):
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
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
            r = p.add_run(str(val)); r.font.name = 'Times New Roman'; r.font.size = Pt(9.5)
    t.rows[0]._tr.get_or_add_trPr().append(OxmlElement('w:tblHeader'))


title('Supplementary Information')
body('Compartment-specific immunity in irreversible pulpitis centres on an inflammatory venule')
body('Zhongzhen Ren¹*')
body('¹School of Stomatology, Shandong Medical and Pharmaceutical University, Yantai 264003, China')
body('*Corresponding author: Zhongzhen Ren, School of Stomatology, Shandong Medical and Pharmaceutical University, Yantai 264003, China. E-mail: 1573741309@qq.com; ORCID: 0009-0008-6735-5224')
body('This Supplementary Information contains the analyses that support but do not drive the conclusions of the main text: unthresholded descriptive plots, the exploratory co-expression and machine-learning material, the photobiomodulation panel plots, the stage-resolved secondary analyses, the negative controls, and the numerical tables behind the main figures. Figures are cited in the main text as Supplementary Fig. 1 to 28 and tables as Supplementary Table 1 to 11.')

h1('Supplementary Figures')
h2('Bulk transcriptome (GSE77459)')
body('Supplementary Fig. 1. Heatmap of the 50 most significant differentially expressed genes between pulpitis and normal pulp (GSE77459). Rows are genes, columns are samples; colour is row-scaled log2 expression.')
body('Supplementary Fig. 2. Top enriched Gene Ontology biological process terms among upregulated genes in pulpitis (over-representation analysis, clusterProfiler). Bar length is gene count; colour is adjusted P.')
body('Supplementary Fig. 3. Top enriched KEGG pathways among upregulated genes in pulpitis (over-representation analysis).')
body('Supplementary Fig. 4. Top activated and suppressed Gene Ontology biological processes in pulpitis (preranked GSEA). Positive NES indicates activation in pulpitis.')
body('Supplementary Fig. 5. Boxplots of ssGSEA scores for 12 immune cell types in normal pulp and pulpitis, showing sample-level distributions and the heterogeneity of the infiltrate.')
body('Supplementary Fig. 6. ROC curves of the top 20 candidate genes within the discovery cohort (in-sample AUC = 1.00 for all genes; explicitly internal and optimistic).')
body('Supplementary Fig. 7. Enrichment plots for the three most strongly activated Hallmark gene sets in pulpitis (allograft rejection, TNF-alpha signalling via NF-kB and inflammatory response).')
h2('Single-cell discovery cohort (GSE185222)')
body('Supplementary Fig. 8. Dot plot of candidate, PBM-panel and M1/M2 genes across the 12 cell types of GSE185222, showing the compartment to which each bulk-level signal maps.')
body('Supplementary Fig. 9. Top 20 ligand-receptor communications increased in deep caries relative to sound/enamel caries (CellChatDB secreted single-subunit pairs, 200 permutations).')
body('Supplementary Fig. 10. Slingshot trajectory of the MSC-like to fibroblast to odontoblast-like lineage (left, cell types with fitted lineage; right, pseudotime) and marker dynamics along lineage 1.')
body('Supplementary Fig. 11. NNLS deconvolution of GSE77459 using the single-cell marker signature (top 200 specificity genes per cell type). Estimates for macrophages and T/NK cells are signature-size sensitive and are not used as evidence in the main text (see Supplementary Table 7).')
h2('Exploratory co-expression and machine-learning analyses')
body('Supplementary Fig. 12. Scale-free topology fit and mean connectivity across soft-thresholding powers for WGCNA (6,000 most variable genes).')
body('Supplementary Fig. 13. Gene dendrogram and module assignment in WGCNA (soft power 5, minimum module size 30, deepSplit 3).')
body('Supplementary Fig. 14. Module eigengene correlation with pulp status (normal versus pulpitis).')
body('Supplementary Fig. 15. Eigengene distributions of the top pulpitis-associated modules in normal and inflamed pulp (immune module M1 and developmental module M2).')
body('Supplementary Fig. 16. Genes selected by at least two of the three feature-selection algorithms (LASSO, random forest, SVM-RFE).')
body('Supplementary Fig. 17. Selection frequency of genes across 50 repetitions of five-fold cross-validation (LASSO stability selection).')
body('Supplementary Fig. 18. Leave-one-out cross-validated ROC curve of the seven-gene logistic model (AUC = 1.00; internal estimate only, not a diagnostic claim).')
h2('Photobiomodulation panel (GSE193448)')
body('Supplementary Fig. 19. Scatter plot of the LPS effect versus the PBM effect for 72 targets of the panel under the corrected sample-to-group assignment. Quadrants II and IV contain reversed genes (32 of 72, 44%; exact binomial P = 0.41 against 50%).')
body('Supplementary Fig. 20. LPS and PBM effects on representative panel genes; positive DeltaCt indicates downregulation. VDR is returned completely to control, MMP8 and VEGFB are suppressed, and ALPL and TWIST1 remain induced.')
h2('Three-stage cohort (E-MTAB-14896)')
body('Supplementary Fig. 21. Seven-gene signature expression across the three disease stages (patient-level pseudobulk log2 CPM). SELE, ICAM2, ADGRL4 and FOSB rise and KIT falls.')
body('Supplementary Fig. 22. M1 and M2 module scores across the three disease stages (patient-level pseudobulk). The two scores rise in parallel, without a monotonic M1-minus-M2 shift.')
body('Supplementary Fig. 23. Top Hallmark gene sets in the three key cell types: endothelial cluster 3, fibroblast/stromal cluster 1 and blood/immune cluster 1 (advanced versus healthy, patient-level pseudobulk).')
body('Supplementary Fig. 24. M1 and M2 module scores within the blood/immune compartment across disease stages (cell level). Cell-level tests are dominated by cell numbers and are treated as descriptive.')
body('Supplementary Fig. 25. M1 and M2 module scores within macrophages/monocytes of GSE185222 (cell level). The difference was not significant (P = 0.27), showing that the bulk M1 skewing is not reproduced at single-cell resolution in this cohort.')
body('Supplementary Fig. 26. PROGENy pathway activities across the three disease stages (patient-level pseudobulk, E-MTAB-14896). No pathway reached significance with four patients per stage.')
body('Supplementary Fig. 27. Transcription-factor activities across the three disease stages (patient-level pseudobulk; STAT1 and MYC rise monotonically, P = 0.0498).')
body('Supplementary Fig. 28. Pair-level ligand-receptor gains across disease stages for signalling from fibroblast/stromal 1 to endothelial 3 (E-MTAB-14896).')

h1('Supplementary Tables')
table_caption('Supplementary Table 1. Expression changes of M1/M2 macrophage marker genes in pulpitis (GSE77459).')
add_table(['Polarisation', 'Gene', 'log2FC', 'Adjusted P', 'Direction'],
          [['M1', 'IL1B', '3.36', '4.55e-03', 'Upregulated'],
           ['M1', 'CXCL10', '2.72', '8.50e-06', 'Upregulated'],
           ['M1', 'PTGS2', '2.44', '3.80e-08', 'Upregulated'],
           ['M1', 'IL6', '2.12', '6.46e-06', 'Upregulated'],
           ['M1', 'CXCL9', '1.83', '1.53e-03', 'Upregulated'],
           ['M1', 'CD86', '1.77', '1.27e-03', 'Upregulated'],
           ['M2', 'CCL18', '2.12', '4.38e-03', 'Upregulated'],
           ['M2', 'CD163', '0.75', '1.36e-02', 'Below FC threshold'],
           ['M2', 'MRC1', '0.40', '1.28e-01', 'Not significant'],
           ['M2', 'MSR1', '0.86', '1.64e-02', 'Below FC threshold'],
           ['M2', 'IL10', '0.70', '5.98e-02', 'Not significant'],
           ['M2', 'ARG1', '0.18', '6.07e-01', 'Not significant']],
          ['c', 'l', 'c', 'c', 'l'])

table_caption('Supplementary Table 2. Cell-type localisation of representative genes (GSE185222).')
add_table(['Gene', 'Top cell type', 'Mean expression'],
          [['IL1B', 'Macrophages/Monocytes', '1.54'],
           ['PTGS2', 'Macrophages/Monocytes', '0.59'],
           ['CD86', 'Macrophages/Monocytes', '0.38'],
           ['CD163', 'Macrophages/Monocytes', '0.89'],
           ['IL6', 'Endothelial', '0.74'],
           ['CXCL10', 'Endothelial', '0.35'],
           ['SELE', 'Endothelial', '1.63'],
           ['SOCS3', 'Endothelial', '2.24'],
           ['ICAM1', 'Endothelial', '1.39'],
           ['FOSB', 'Fibroblasts/MSC', '2.33'],
           ['ALPL', 'Fibroblasts/MSC', '0.63'],
           ['RUNX2', 'Fibroblasts/MSC', '0.42'],
           ['VEGFB', 'Fibroblasts/MSC', '0.10'],
           ['FGFR2', 'Fibroblasts/MSC', '0.06'],
           ['TNFSF11', 'Fibroblasts/MSC', '0.03'],
           ['BMP4', 'Odontoblasts', '0.31'],
           ['CD36', 'Odontoblasts', '0.69'],
           ['VDR', 'Odontoblast-like', '0.10'],
           ['TWIST1', 'Odontoblast-like', '0.22'],
           ['NR1D1', 'Odontoblast-like', '0.12'],
           ['SPP1', 'Neural/Schwann (SPP1+)', '1.25'],
           ['MMP9', 'Macrophages/Monocytes', '0.37']],
          ['l', 'l', 'c'])

table_caption('Supplementary Table 3. NNLS deconvolution of GSE77459 with the single-cell marker signature (top 200 genes per cell type).')
add_table(['Cell type', 'Normal (%)', 'Pulpitis (%)', 'P (Wilcoxon)'],
          [['Odontoblasts (DMP1+/IBSP+)', '46.0', '23.1', '0.005'],
           ['Odontoblast-like (OMD+)', '9.9', '0.2', '0.005'],
           ['Fibroblasts', '19.3', '30.0', '0.093'],
           ['Fibroblasts/MSC', '2.0', '0', '0.003'],
           ['Macrophages/Monocytes', '6.1', '13.1', '0.045'],
           ['T/NK cells', '2.9', '7.4', '0.008'],
           ['B cells', '7.5', '11.0', '0.005'],
           ['Plasma cells', '0.2', '6.8', '0.005'],
           ['Endothelial', '2.3', '2.8', '0.174'],
           ['Schwann cells (MPZ+)', '0.3', '0.7', '0.005']],
          ['l', 'c', 'c', 'c'])

table_caption('Supplementary Table 4. WGCNA modules and their association with pulpitis (exploratory).')
add_table(['Module', 'Genes', 'Correlation with pulpitis (r)', 'P', 'Dominant content'],
          [['M1', '2,402', '0.906', '4.8e-05', 'Immune / immunoglobulin / myeloid'],
           ['M2', '2,028', '-0.901', '6.3e-05', 'Developmental / odontogenic (SOX9, DLX6)'],
           ['M0', '1,343', '-0.310', '0.33', 'Unassigned'],
           ['M3', '134', '0.177', '0.58', 'Small module'],
           ['M4', '93', '-0.043', '0.90', 'Small module']],
          ['c', 'c', 'c', 'c', 'l'])

table_caption('Supplementary Table 5. Consensus feature-selection genes: selection, module membership and single-cell localisation (exploratory).')
add_table(['Gene', 'Methods', 'Module', 'kME', 'log2FC', 'Adjusted P', 'Top cell type (mean, % detected)'],
          [['ADGRL4', 'LASSO + SVM-RFE', 'M1', '0.93', '1.33', '3.8e-06', 'Endothelial (1.12, 62.1%)'],
           ['CCL21', 'LASSO + SVM-RFE', 'M1', '0.89', '4.79', '1.4e-08', 'Endothelial (0.01, 1.2%)'],
           ['FOSB', 'LASSO + SVM-RFE', 'M1', '0.90', '2.41', '1.6e-08', 'Fibroblasts/MSC (2.33, 91.0%)'],
           ['ICAM2', 'LASSO + SVM-RFE', 'M1', '0.93', '1.28', '2.0e-06', 'Endothelial (0.64, 42.2%)'],
           ['KIT', 'LASSO + SVM-RFE', 'M2', '0.93', '-1.04', '7.8e-06', 'Fibroblasts/MSC (0.51, 37.8%)'],
           ['PI3', 'LASSO + SVM-RFE', 'M1', '0.90', '4.74', '3.0e-11', 'Not detected'],
           ['SELE', 'LASSO + SVM-RFE', 'M1', '0.90', '3.34', '6.8e-08', 'Endothelial (1.63, 54.8%)']],
          ['l', 'l', 'c', 'c', 'c', 'c', 'l'])

table_caption('Supplementary Table 6. Robustness of the WGCNA immune and developmental modules across parameter variants. Module counts are pre-merge colour assignments; the main network described in the text retained five merged modules (M0-M4).')
add_table(['Variant', 'Modules', 'Immune-like Jaccard', 'Immune-like r', 'P', 'Developmental-like Jaccard', 'Developmental-like r', 'P'],
          [['Original (power 5, min 30, deep 3)', '21 (5 after merging)', '0.444', '0.921', '2.2e-05', '0.155', '-0.738', '0.006'],
           ['Power 7', '25', '0.439', '0.910', '4.0e-05', '0.121', '-0.798', '0.002'],
           ['Minimum size 50', '14', '0.451', '0.921', '2.2e-05', '0.165', '-0.638', '0.026'],
           ['deepSplit 2', '18', '0.463', '0.921', '2.2e-05', '0.062', '-0.280', '0.378'],
           ['Top 5,500 genes', '17', '0.425', '0.919', '2.4e-05', '0.169', '-0.789', '0.002']],
          ['l', 'c', 'c', 'c', 'c', 'c', 'c', 'c'])

table_caption('Supplementary Table 7. Deconvolution robustness: top 200 versus top 100 signature genes per cell type.')
add_table(['Cell type', 'Normal (%) top200', 'Pulpitis (%) top200', 'P', 'Normal (%) top100', 'Pulpitis (%) top100', 'P', 'Direction'],
          [['B cells', '7.5', '11.0', '0.005', '0.4', '4.3', '0.005', 'Stable up'],
           ['Plasma cells', '0.2', '6.8', '0.005', '0.1', '5.1', '0.004', 'Stable up'],
           ['Odontoblasts (DMP1+/IBSP+)', '46.0', '23.1', '0.005', '64.2', '26.7', '0.005', 'Stable down'],
           ['Fibroblasts', '19.3', '30.0', '0.093', '20.0', '51.5', '0.005', 'Stable up'],
           ['Macrophages/Monocytes', '6.1', '13.1', '0.045', '7.1', '4.7', '0.005', 'Reversed'],
           ['T/NK cells', '2.9', '7.4', '0.008', '3.5', '5.0', '0.298', 'Unstable'],
           ['Schwann cells (MPZ+)', '0.3', '0.7', '0.005', '1.5', '0.3', '0.005', 'Reversed'],
           ['Odontoblast-like (OMD+)', '9.9', '0.2', '0.005', '0', '0', '-', 'Zero in variant'],
           ['Fibroblasts/MSC', '2.0', '0', '0.003', '0', '0', '-', 'Zero in variant']],
          ['l', 'c', 'c', 'c', 'c', 'c', 'c', 'l'])

table_caption('Supplementary Table 8. Drug perturbation signatures versus the pulpitis transcriptome (DSigDB, direction-aware preranked GSEA).')
add_table(['Drug', 'Direction', 'NES', 'Adjusted P', 'Note'],
          [['Rimexolone', 'Opposing', '2.60', '8.6e-11', 'Glucocorticoid'],
           ['Parthenolide', 'Opposing', '2.43', '1.9e-12', 'NF-kB inhibitor'],
           ['Progesterone', 'Opposing', '2.38', '2.6e-07', 'Immunomodulator'],
           ['Valproic acid', 'Opposing', '2.33', '1.8e-10', 'HDAC inhibitor'],
           ['CP-690334-01', 'Opposing', '2.30', '8.1e-07', 'JAK inhibitor (tofacitinib code)'],
           ['Geldanamycin', 'Opposing', '2.11', '3.0e-04', 'HSP90 inhibitor'],
           ['Tretinoin', 'Aligned', '3.09', '1.8e-39', 'HL60 myeloid differentiation'],
           ['Clonidine', 'Aligned', '2.90', '1.9e-18', '-'],
           ['Alprostadil', 'Aligned', '2.87', '5.6e-38', 'Prostaglandin E1'],
           ['Dinoprost', 'Aligned', '2.77', '6.4e-15', 'Prostaglandin F2-alpha'],
           ['Mebendazole', 'Aligned', '3.18', '1.1e-30', 'HL60 myeloid differentiation'],
           ['Etoposide', 'Aligned', '3.13', '2.1e-30', 'HL60 myeloid differentiation']],
          ['l', 'l', 'c', 'c', 'l'])

table_caption('Supplementary Table 9. Normalised ligand-receptor communication share on key axes (E-MTAB-14896). Pooled estimates across patients; patient-level permutation tests are given in Supplementary Table 11.')
add_table(['Axis', 'Healthy', 'Early', 'Advanced', 'Advanced minus Healthy'],
          [['Fibroblast/stromal 1 to Endothelial 3', '0.00825', '0.01146', '0.01122', '+0.00297'],
           ['Fibroblast/stromal 1 to immune', '0.01285', '0.01396', '0.01288', '+0.00003'],
           ['Immune to Endothelial 3', '0.01544', '0.01385', '0.01465', '-0.00078'],
           ['Immune to Fibroblast/stromal 1', '0.01056', '0.00946', '0.00803', '-0.00253'],
           ['Endothelial 1 to immune', '0.01696', '0.01379', '0.01432', '-0.00263'],
           ['Endothelial 3 to immune', '0.01215', '0.00960', '0.00834', '-0.00381'],
           ['Immune to Endothelial 1', '0.01918', '0.01551', '0.01554', '-0.00365'],
           ['Immune to Glial 2', '0.01442', '0.01040', '0.00814', '-0.00628']],
          ['l', 'c', 'c', 'c', 'c'])

# ---- Supplementary Tables 10 and 11: generated directly from analysis CSVs ----
import csv, os
SCR = r'C:\Users\19354\Doubao\chats\2026-11-18\new-chat\pulp_bioinfo\sc\results'

def read_csv(name):
    with open(os.path.join(SCR, name), encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))

lopo = read_csv('EMTAB_LOPO_endo.csv')
table_caption('Supplementary Table 10. Leave-one-patient-out analysis of endothelial expansion (E-MTAB-14896). P values are Kruskal-Wallis tests across the three stages after omitting the indicated patient; differences are advanced minus healthy patient means in percentage points.')
rows10 = []
for r in lopo:
    rows10.append([r['removed'],
                   '%.3f' % float(r['p_endo']), '%.3f' % float(r['p_ratio']), '%.3f' % float(r['p_endo3']),
                   '%+.2f' % (float(r['endo_A']) - float(r['endo_H'])),
                   '%+.2f' % (float(r['endo3_A']) - float(r['endo3_H']))])
add_table(['Removed patient', 'P endothelium', 'P ratio', 'P Endo3', 'Endo diff (pp)', 'Endo3 diff (pp)'],
          rows10, ['l', 'c', 'c', 'c', 'c', 'c'])

perm = read_csv('EMTAB_LR_patient_permutation.csv')
table_caption('Supplementary Table 11. Patient-level ligand-receptor communication shares and permutation tests (E-MTAB-14896). Shares are means of four patients per stage. Delta is advanced minus healthy. Permutation P values are two-sided from 10,000 stage-label permutations preserving group sizes; FDR is Benjamini-Hochberg across the eight axes. KW, Kruskal-Wallis.')
rows11 = []
for r in perm:
    rows11.append([r['axis'],
                   '%.4f' % float(r['healthy']), '%.4f' % float(r['early']), '%.4f' % float(r['advanced']),
                   '%+.4f' % float(r['delta']),
                   '%.4f' % float(r['p_kw']), '%.4f' % float(r['p_perm_two']), '%.4f' % float(r['p_perm_fdr'])])
add_table(['Axis', 'Healthy', 'Early', 'Advanced', 'Delta', 'KW P', 'Permutation P', 'FDR'],
          rows11, ['l', 'c', 'c', 'c', 'c', 'c', 'c', 'c'])

doc.save(OUT)
print('SAVED:', OUT)
