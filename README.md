# perlin_microfibrosis
Perlin Noise for the Creation of Microfibrotic Patterns
https://www.biorxiv.org/content/10.1101/668848v1


## Main Folder
Contains MATLAB code that includes the Perlin noise pattern generator,
as well as the SMC-ABC algorithm used to match to a provided target pattern.
Code for plotting results is also provided, but the histological targets,
used in the paper, are not made available as they are the property of the
original journal (see deJong reference in linked paper). Provided .mat files are
the final particle sets, and a set of seed data used by the Perlin noise
generator.

Chastefiles folder contains a Chaste project (qutemu) used for running the
electrophysiological simulations presented in the work.

Contact: b.lawson@qut.edu.au