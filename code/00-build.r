# Use renv to download R packages from the renv.lock in the project directory ####

# Disable automatic snapshots to avoid writing over the lock file ####
options(renv.config.auto.snapshot = F)

# Initialize renv, but don't let it search for dependencies ####
renv::init(bare = T)

# Begin to build libraries from the lock file ####
renv::restore(prompt = F)
