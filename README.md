# A multi-omics study of brain tissue transcription and DNA methylation revealing the genetic pathogenesis of ADHD
This paper is published on _Briefings in Bioinformatics_ journal. 
# Code availability
Please find all the source code on GitHub files.
# Output data
Please find the open access to output data of the paper at [https://console.cloud.google.com/storage/browser/output20240516](https://console.cloud.google.com/storage/browser/output20240516) under folder `/BIB-24-0876_output_data`.

## File names under folder /_BIB-24-0876_output_data_
• `/TWAS_output`

&nbsp;&nbsp;&nbsp;&nbsp;• `/eTWAS_<Method>`  &nbsp;--Directory of eTWAS results

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;• `/Expression_<Method>_result_<Tissue>.csv`  &nbsp;--eTWAS result of a certain tissue implemented by a certain method

&nbsp;&nbsp;&nbsp;&nbsp;• `/sTWAS_<Method>`  &nbsp;--Directory of sTWAS results

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;• `/Splicing_<Method>_result_<Tissue>.csv`  &nbsp;--sTWAS result of a certain tissue implemented by a certain method

• `/Mediation_output`

&nbsp;&nbsp;&nbsp;&nbsp;• `<Phenotype>_to_<Phenotype>`  &nbsp;--Causal inference results for five different regulatory relationship

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;• `<Phenotype>_to_<Phenotype>_<Method>.csv`  &nbsp;--Mediation analysis between phenotypes of two omics implemented by a certain method
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
