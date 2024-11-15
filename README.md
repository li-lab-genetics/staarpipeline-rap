# staarpipeline (DNAnexus Platform App)

This is the source code for the staarpipeline app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://documentation.dnanexus.com.

### Applet Usage
The **staarpipeline** app can run single variant, gene-centric coding, gene-centric noncoding, ncRNA, sliding window, and dynamic window tests for biobank-scale whole-genome/whole-exome sequencing data. It will account for relatedness using a kinship/relatedness matrix and dynamically incorporates multiple functional annotations to empower rare variant (set) association analysis.

Specifically, this app will

1. Fit the null model. This is fitting your model with your outcome, adjustments and kinship/genetic relatedness matrix, but does not use the genotypes;

2. Take the null model object from the first step and run your association analysis, while dynamically incorporating multiple functional annotations to empower rare variant (set) association analysis using the STAAR method. The same null model can be used for single variant or aggregate tests.

The single variant or aggregate test results generated from this app can be summarized by the <a href="https://github.com/xihaoli/staarpipelinesummary_varset-rap">**staarpipelinesummary_varset**</a> and <a href="https://github.com/xihaoli/staarpipelinesummary_indvar-rap">**staarpipelinesummary_indvar**</a> apps.

Please see the <a href="https://tinyurl.com/staarpipelineapps">**user manual and tutorial**</a> for detailed usage of staarpipeline app.

### Cloning an Applet
To acquire the staarpipeline applet, you will need to compile this applet for your respective DNANexus project, by cloning the repository from github and `dx build` an APPLET into your own workspace.

1. Clone this github repo to some directory:

```commandline
git clone https://github.com/li-lab-genetics/staarpipeline-rap.git
```

This will create a folder named staarpipeline-rap, you can then:

2. Compile the source code:

```commandline
dx build -f staarpipeline-rap
```

the `-f` flag just tells DNANexus to overwrite older versions of the applet within the same project if it is already there.

You can then run the following to run this applet:

```commandline
dx run staarpipeline <options>
```

### Version
The current version is 0.9.7.1.

### Citation
Zilin Li*, Xihao Li*, Hufeng Zhou, Sheila M. Gaynor, Margaret Sunitha Selvaraj, Theodore Arapoglou, Corbin Quick, Yaowu Liu, Han Chen, Ryan Sun, Rounak Dey, Donna K. Arnett, Paul L. Auer, Lawrence F. Bielak, Joshua C. Bis, Thomas W. Blackwell, John Blangero, Eric Boerwinkle, Donald W. Bowden, Jennifer A. Brody, Brian E. Cade, Matthew P. Conomos, Adolfo Correa, L. Adrienne Cupples, Joanne E. Curran, Paul S. de Vries, Ravindranath Duggirala, Nora Franceschini, Barry I. Freedman, Harald H. H. Göring, Xiuqing Guo, Rita R. Kalyani, Charles Kooperberg, Brian G. Kral, Leslie A. Lange, Bridget M. Lin, Ani Manichaikul, Alisa K. Manning, Lisa W. Martin, Rasika A. Mathias, James B. Meigs, Braxton D. Mitchell, May E. Montasser, Alanna C. Morrison, Take Naseri, Jeffrey R. O’Connell, Nicholette D. Palmer, Patricia A. Peyser, Bruce M. Psaty, Laura M. Raffield, Susan Redline, Alexander P. Reiner, Muagututi’a Sefuiva Reupena, Kenneth M. Rice, Stephen S. Rich, Jennifer A. Smith, Kent D. Taylor, Margaret A. Taub, Ramachandran S. Vasan, Daniel E. Weeks, James G. Wilson, Lisa R. Yanek, Wei Zhao, NHLBI Trans-Omics for Precision Medicine (TOPMed) Consortium, TOPMed Lipids Working Group, Jerome I. Rotter, Cristen J. Willer, Pradeep Natarajan, Gina M. Peloso, & Xihong Lin. (2023). **A framework for detecting noncoding rare variant associations of large-scale whole-genome sequencing studies**. _Nature Methods_, _19_(12), 1599-1611. PMID: <a href="https://www.ncbi.nlm.nih.gov/pubmed/36303018">36303018</a>. PMCID: <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10008172/">PMC10008172</a>. DOI: <a href="https://doi.org/10.1038/s41592-022-01640-x">10.1038/s41592-022-01640-x</a>.

Xihao Li*, Zilin Li*, Hufeng Zhou, Sheila M. Gaynor, Yaowu Liu, Han Chen, Ryan Sun, Rounak Dey, Donna K. Arnett, Stella Aslibekyan, Christie M. Ballantyne, Lawrence F. Bielak, John Blangero, Eric Boerwinkle, Donald W. Bowden, Jai G. Broome, Matthew P. Conomos, Adolfo Correa, L. Adrienne Cupples, Joanne E. Curran, Barry I. Freedman, Xiuqing Guo, George Hindy, Marguerite R. Irvin, Sharon L. R. Kardia, Sekar Kathiresan, Alyna T. Khan, Charles L. Kooperberg, Cathy C. Laurie, X. Shirley Liu, Michael C. Mahaney, Ani W. Manichaikul, Lisa W. Martin, Rasika A. Mathias, Stephen T. McGarvey, Braxton D. Mitchell, May E. Montasser, Jill E. Moore, Alanna C. Morrison, Jeffrey R. O'Connell, Nicholette D. Palmer, Akhil Pampana, Juan M. Peralta, Patricia A. Peyser, Bruce M. Psaty, Susan Redline, Kenneth M. Rice, Stephen S. Rich, Jennifer A. Smith, Hemant K. Tiwari, Michael Y. Tsai, Ramachandran S. Vasan, Fei Fei Wang, Daniel E. Weeks, Zhiping Weng, James G. Wilson, Lisa R. Yanek, NHLBI Trans-Omics for Precision Medicine (TOPMed) Consortium, TOPMed Lipids Working Group, Benjamin M. Neale, Shamil R. Sunyaev, Gonçalo R. Abecasis, Jerome I. Rotter, Cristen J. Willer, Gina M. Peloso, Pradeep Natarajan, & Xihong Lin. (2020). **Dynamic incorporation of multiple in silico functional annotations empowers rare variant association analysis of large whole-genome sequencing studies at scale**. _Nature Genetics_, _52_(9), 969-983. PMID: <a href="https://www.ncbi.nlm.nih.gov/pubmed/32839606">32839606</a>. PMCID: <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7483769/">PMC7483769</a>. DOI: <a href="https://doi.org/10.1038/s41588-020-0676-4">10.1038/s41588-020-0676-4</a>.


