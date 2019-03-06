source("exp012depletion_sankey.R")
source("exp2depletion.R")
source("exp3metabolomics.R")
source("exp3proteomics.R")
source("exp4spentmedia_metabolomics.R")
source("exp5concentration.R")
source("exp6transfers.R")
source("summary.R")
source("screen_vs_assay.R")


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
}

exp012depletion.sankey()
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