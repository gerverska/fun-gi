# Use renv to download R packages from renv.lock to the project directory ####

# Change the snapshot protocol ####
options(renv.settings.snapshot.type = 'all')

# Initialize renv, but don't let it search for dependencies ####
renv::init(bare = T)

# Begin to build libraries from the lock file ####
renv::restore(prompt = F)
