# staarpipeline (DNAnexus Platform App)

This is the source code for the staarpipeline app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://documentation.dnanexus.com/.

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