# staarpipeline (DNAnexus Platform App)

This is the source code for the staarpipeline app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://documentation.dnanexus.com/.

### Applet Usage
The staarpipeline app can run single variant, gene-centric coding, gene-centric noncoding, ncRNA, sliding window, and dynamic window tests for biobank-scale whole-genome/whole-exome sequencing data. It will account for relatedness using a kinship/relatedness matrix and dynamically incorporates multiple functional annotations to empower rare variant (set) association analysis.

Specifically, this app will

1. Fit the null model. This is fitting your model with your outcome, adjustments and kinship/genetic relatedness matrix, but does not use the genotypes;

2. Take the null model object from the first step and run your association analysis, while dynamically incorporating multiple functional annotations to empower rare variant (set) association analysis using the STAAR method. The same null model can be used for single variant or aggregate tests.


Please see the <a href="https://tinyurl.com/staarpipeline">**user manual and tutorial**</a> for detailed usage of staarpipeline app.

### Cloning an Applet
To acquire the staarpipeline applet, you will need to compile this applet for your respective DNANexus project, by cloning the repository from github and `dx build` an APPLET into your own workspace.

1. Clone this github repo to some directory:

```commandline
git clone https://github.com/xihaoli/staarpipeline-rap.git
```

This will create a folder named staarpipeline-rap, you can then:

2. Compile the source code:

```commandline
dx build -f staarpipeline-rap
```

the `-f` flag just tells DNANexus to overwrite older versions of the applet within the same project if it is already there.

You can then run the following to run this applet:

```commandline
dx run staarpipeline-rap <options>
```
