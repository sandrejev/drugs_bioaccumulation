source("exp012depletion_sankey.R")
source("exp2depletion.R")
source("exp3metabolomics.R")
source("exp3proteomics.R")
source("exp4spentmedia_metabolomics.R")
source("exp5concentration.R")
source("exp6transfers.R")
source("exp10bioaccumulation.R")
source("exp12nmr.R")
source("exp7spentmedia_Erectale_growth.R")
source("summary.R")
source("screen_vs_assay.R")
source("chemical_diversity.R")


# Create docker container
create.dockerfile = function()
{
  library(containerit)
  dfile = containerit::dockerfile(
    copy="script_dir", 
    cmd=containerit::Cmd(params="run_everything.R"), 
    entrypoint=containerit::Entrypoint("Rscript", params=list("run_everything.R")), 
    versioned_packages=T)
  write(dfile, file=file.path(getwd(), "Dockerfile"))
  
  system("docker build --file Dockerfile --tag sandrejev/drugs_bioaccumulation .", wait=T)
  system("docker save -o drugs_bioaccumulation.tar sandrejev/drugs_bioaccumulation:latest .", wait=T)
  system("docker push sandrejev/drugs_bioaccumulation", wait=T)
}

exp012depletion.sankey() # Done (Figure 1)
exp12nmr_analyze() # Figure 2b
exp10bioaccumulation.analyze() # Figure 2c
exp6transfers.analyze() # Figure 3a
exp5concentration.analyze() # Figure 3b
exp7spentmedia.Erectale_growth() # Figure 3c


exp2depletion.total2supernatant()
exp3metabolomics.analyze()
exp3proteomics.analyze()
exp4spentmedia.metabolomics()
exp5concentration.analyze()
exp6transfers.analyze()
summary.bugs_piechart()
summary.enzymatic_coverage()
screen_vs_assay()
chemical_diversity()
exp7spentmedia.Erectale_growth()

#
# Export to Metabolights
#