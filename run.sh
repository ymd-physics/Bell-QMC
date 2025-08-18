# ------------------------------------------------
#   Input params
# ------------------------------------------------
l=8                 # number of sites
beta=$((l*6))       # inverse temperature (large enough to approximate the ground state)
J=1                 # coupling strength of the Ising interaction
h=1                 # strength of the magnetic fields
num_thm=20000       # number of MC steps for thermalization
num_stat=50000      # number of MC samples for each binning data
num_bins=5          # number of bins
seed=2025           # seed for PRNG

# -------------------------------------
#   Initialize "/data" and compile 
# -------------------------------------
cargo build --release   
path="./data/l"$l"_beta"$beta"_J"$J"_h"$h
mkdir $path

# ---------------------------
#   Run the program
# ---------------------------
./target/release/bell_qmc_tfim_1d $l $beta $J $h $num_thm $num_stat $num_bins $path $seed