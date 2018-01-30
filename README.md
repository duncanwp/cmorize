# cmorize
Python scripts for converting raw GCM output to CMOR (and CF) compliant NetCDF files

Example usage:

    ./cmorize_ham.py '/nerc/n02/n02/duncanwp/pdrmip-bcx10_fsst/Raw/pdrmip-bcx10_fsst_2*_.nc' /nerc/n02/n02/duncanwp/pdrmip-bcx10_fsst/Data/ECHAM6.3-HAM2.2_bcx10_fsst --experiment "BCx10 fixed SST run using ECHAM6.3-HAM2.2" --contact "Duncan Watson-Parris: duncan.watson-parris@physics.ox.ac.uk" --institute "University of Oxford" --time '200001-203101' --pdrmip --pdrmip_format