Paralogue Annotation utilizes information from evolutionarily related proteins, specifically paralogues, to help inform the clinical significance of missense variants associated with human diseases. The original methodology and implementation of [Paralogue Annotation on arrhythmia syndrome genes](https://www.cardiodb.org/paralogue_annotation/) was published [here](https://onlinelibrary.wiley.com/doi/full/10.1002/humu.22114) and [here](https://jmg.bmj.com/content/51/1/35). This web app extends Paralogue Annotation exome-wide, using paralogues defined by [Ensembl's gene trees](https://www.ensembl.org/Help/View?id=137) and pathogenic/likely pathogenic missense variants defined by [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/).
<br>
This web app is currently being built using [Shiny](http://shiny.rstudio.com), the source code is available at [https://github.com/ImperialCardioGenetics/Paralogue_Annotation_App](https://github.com/ImperialCardioGenetics/Paralogue_Annotation_App).
<br>
<br>
<br>
## Frequently Asked Questions (FAQ)
<br>
**Q. What genome build coordinates do my variants need to be in?**  
**A.** Currently only GRCh37 coordinates are supported. We recommend using [Ensembl's liftover service](https://www.ensembl.org/Homo_sapiens/Tools/AssemblyConverter) for coordinate conversions.
<br>
<br>
**Q. What are the Para_z scores?**  
**A.** The Para_z scores are a measure of paralogue conservation independently derived by [_Lal et al._ (2020)](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-020-00725-6). You may therefore find in your results that some Para_z scores do not agree with your expectations. This is because the paralogue alignments used to generate the scores are different to the alignments used here. The Para_z scores can thus be thought as a third-party confidence score of paralogue conservation across aligned positions.
<br>
<br>
**Q. How can we use the conservation of Ref/Alt alleles to filter out results?**  
**A.** That is not currently available in this version of the web app.
<br>
<br>
**Q. What paralogue alignments do you use here?**  
**A.** We utilize paralogue alignments at the protein level generated by Ensembl, which were obtained through [Compara](http://www.ensembl.org/info/docs/api/compara/compara_tutorial.html)
<br>
<br>
**Q. Why do the results for arrhythmia genes from the original [Paralogue Annotation](https://www.cardiodb.org/paralogue_annotation/) and here differ?**  
**A.** This is mainly because the original Paralogue Annotation utilized T-COFFEE for the alignments, whereas Ensembl's alignments are generated by CLUSTAL W instead. Furthermore, variants from [HGMD](http://www.hgmd.cf.ac.uk/ac/index.php) were used instead of [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/).
<br>
<br>
**Q. What formats do my input variants have to be in?**  
**A.** Currently variants have to be submitted using their chromosome, position, reference allele, and alternate allele using any delimiter in the format of "CHROM:POS:REF:ALT" with separate variants on newlines. Alternatively we also accept [VCF](https://www.internationalgenome.org/wiki/Analysis/vcf4.0/) as well.
<br>
<br>
<br>
For more details on specific methods, code of how Paralogue Annotation functions, or any other questions please email nyl112@ic.ac.uk
<br>
<br>
<br>
This web app is a **work in progress,** final version may differ.