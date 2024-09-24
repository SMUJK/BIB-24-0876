# BIB-24-0876
Open access to output data and source codes of paper 'A multi-omics study of brain tissue transcription and DNA methylation revealing the genetic pathogenesis of ADHD' published on journal _Briefings in Bioinformatics_.

## Directory names
• `/TWAS_output`

&nbsp;&nbsp;&nbsp;&nbsp;• `/eTWAS_<Method>`  &nbsp;--Directory of eTWAS results

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;• `/Expression_<Method>_result_<Tissue>.csv`

&nbsp;&nbsp;&nbsp;&nbsp;• `/sTWAS_<Method>`  &nbsp;--Directory of sTWAS results

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;• `/Splicing_<Method>_result_<Tissue>.csv`

• `/Mediation_output`

&nbsp;&nbsp;&nbsp;&nbsp;• `<Phenotype>_to_<Phenotype>`  &nbsp;--Causal inference results for five different regulatory relationship

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;• `<Phenotype>_to_<Phenotype>_<Method>.csv`
## Options
### \<Method\>
`[IVW]`  &nbsp;--Inverse-variance weighted method

`[MR_Egger]`  &nbsp;--MR-Egger method

`[LDA]`  &nbsp;--LDA MR-Egger method

`[PMR]`  &nbsp;--PMR-Egger method
### \<Tissue\>
`[Amygdala][Anterior cingulate cortex][Caudate (basal ganglia)][Cerebellar hemisphere][Cerebellum][Cortex][Frontal cortex][Hippocampus][Hypothalamus][Nucleus accumbens (basal ganglia)][Pituitary][Putamen (basal ganglia)][Spinal cord (cervical c-1)][Substantia nigra]`

### \<Phenotype\>
`[Trait]`  &nbsp;--ADHD

`[Expression]`  &nbsp;--Gene expression

`[Methylation]`  &nbsp;--DNA methylation

`[Splicing]`  &nbsp;--Alternative splicing event expression
