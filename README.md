# BIB-24-0876
Open access to output data and source codes of BIB paper 'A multi-omics study of brain tissue transcription and DNA methylation revealing the genetic pathogenesis of ADHD'

## Directory names
• TWAS_output
  • `eTWAS_[method]`--eTWAS results
    • `Expression_[method]_result_[tissue].csv`

  - `sTWAS_[method]`--sTWAS results
    - `Splicing_[method]_result_[tissue].csv`

- Mediation_output
  - `[Phenotype]_to_[Phenotype]`--causal inference results for five different regulatory relationship
    - `[Phenotype]_to_[Phenotype]_[method].csv`

## Options
method:\
`[IVW]`--Inverse-variance weighted method\
`[MR_Egger]`--MR-Egger method\
`[LDA]`--LDA MR-Egger method\
`[PMR]`--PMR-Egger method\
\
tissue:\
`[Amygdala][Anterior cingulate cortex][Caudate (basal ganglia)][Cerebellar hemisphere][Cerebellum][Cortex][Frontal cortex][Hippocampus][Hypothalamus][Nucleus accumbens (basal ganglia)][Pituitary][Putamen (basal ganglia)][Spinal cord (cervical c-1)][Substantia nigra]`\
\
Phenotype:\
`[Trait]`\
`[Expression]`\
`[Methylation]`\
`[Splicing]`
