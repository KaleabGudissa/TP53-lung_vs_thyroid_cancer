Inputs: cBioPortal TCGA PanCancer Atlas study bundles
- luad_tcga_pan_can_atlas_2018.tar.gz
- thca_tcga_pan_can_atlas_2018.tar.gz

Extracted tables analyzed:
- data_mutations.txt (or data_mutations_extended.txt)
- data_clinical_patient.txt

Definitions:
- patient_id = first 12 chars of Tumor_Sample_Barcode
- TP53 mutated = any non-silent Variant_Classification

Outputs:
- Table1: counts, %, Wilson 95% CI
- Table2: chi-square + Fisher, OR, RR (95% CI)
- Figure1: % TP53 (95% CI bars)
- Figure2: KM (log-rank)
- Table3: Cox PH (if possible)
