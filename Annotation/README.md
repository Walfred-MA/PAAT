The pipeline to annotate pangenome-alleles.

Example running command:

example command on HPC cluster: snakemake -k --cluster "sbatch --account=mchaisso_100 --time 50:00:00 --partition=qcb --time=200:00:00 {resources.slurm_extra} " --default-resources "mem_mb=3000" --jobs 50 --rerun-incomplete --nt --wait-for-files --latency-wait 100 --resources mem_gb=1000 --configfile example.json &
