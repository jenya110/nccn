This repo contains data and analysis pipeline for:
1) Comparison of multispecies conservation of non-consensus nucleotides with neutral control (L.py);
2) Patterns of selection pressure on CNs and NCNs (4nes.py);
3) Selection acting on the energy of binding between TF and TFBS (weight_problem.py);
4) Overlap with unknown sites (lenrat.py);

To run any of these scripts (except for weight_problem.py) for a bacterial family, 4nes.py needs to be run first as it generates files neccesary for the other scripts.
All the scripts and input files are already located for the scripts to run properly.
Within each bacterial family, each TF folder contains two input files: output.xlsx (with all the identified TFBSs for this TF) and clustered.xlsx (only orthologous TFBSs groupped by orthology).
