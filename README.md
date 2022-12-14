# Miniscope_cleaning_GUI
[Matlab] GUI to clean preprocessed Miniscope V3 data in the Willuhn lab at the Netherlands Institute for Neuroscience. 

This Matlab GUI requires CNMF-E output as input and enables the user to inspect all (potential) neurons.

Use it to inspect all neurons, delete false positives (background), and merge duplicates.

Options:
- Distance threshold: max distance between neurons
- Correlation threshold: min correlation between neurons
- Contour threshold: fraction of pixels to plot contour of neurons
- SNR: min signal-to-noise (based on OASIS)

GUI helps you to either delete or merge cells.

Delete (delete neuron spatially and temporal trace):
- use slider to look at contour of neuron (middle panel: Spatial location of cells) and temporal activity pattern (left panel: SNR- value, Cell number)
- use "Add to delete list" to add neuron to delete list
- use "Remove from the delete list" to remove neuron from delete list
- once happy, press "Delete" to delete cells on the delete list
- use "save files to" to save output file
- RESTART GUI 

Merge (average neuron pair spatially and temporal trace):
- use "Previous" and "Next" to look through all pairs of neurons that match your criteria (Distance threshold and Correlation threshold)
- either use "Add Lazy" to add current pair to merge list, or select neuron numbers (top right) and select "Add to merge list"
- use "Delete from merge list" to delete a pair of neurons from the merge list
- one happy, press "Check Merge" to check if a neuron is in multiple pairs
- notification "Uh oh. You have cells that are repeated in the list" means there are duplicates. Press "Multi merge!" to merge duplicates. Press "Check Merge" again, if you still get the error, manually delete the pairs from the merge list
- notification "Everything is awesome!!" means no duplicates
- press "Merge" to merge the pairs on the merge list
- use "save files to" to save output file
- RESTART GUI 

GUI written by Aishwarya Parthasarathy, PhD
