# Real-time monitoring of reaction mechanisms from spectroscopic data using hidden semi-Markov models for mode identification

Description of the datasets:
1. The spectroscopic data associated with a chemical reaction process subject to a decreasing temperature sequence can be found in the Decreasing_temperature_data.csv file. 

2. The spectroscopic data associated when a similar process is subject to a randomized sequence of temperature changes can be found in the Randomized_temperature_data.csv file.

In both the datasets the first column contains the temperature sequence, the second column contains the randomly sampled residence time durations at each of the temperatures, and the remaining columns are the spectral absorbances recorded at channels specified in the wavenumbers.csv file.

The numerical implementation of the HSMM follows the calculation of the conditional probabilities using logarithms, that are duly scaled by normalization to avoid instabilities during computation, as outlined in Mann, T. P.. “Numerically Stable Hidden Markov Model Implementation.” (2006).

Run the main.m code file to train the HSMM on spectroscopic data, followed by using it for state decoding via the Viterbi algorithm, and inferences via the calculation of posterior probabilities.
