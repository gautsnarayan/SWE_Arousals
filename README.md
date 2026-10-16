# SWE_Arousals
The study tests whether slow wave entrainment during sleep delivered by a smartwatch can increase measurements of sleep quality, alertness, and cognitive function


SWE Results Section:
  Write-up of the section where we validate the simulation to calculate predicted arousal rate.

fig 1,2,3:
  Graphs of the real vs predicted arousal rates

full_batch_h8_results (.csv, .mat):
  This is the full output variable needed to run all downstream analysis scripts.
  It contains all the participants/files and their real stim/arousal/arousal rate and predicted stim/arousal/arousal rate

h13_v4_101626:
  This is the script used to run all the analyses for pre-registered hypothesis H13. 
  The first pass of analysis will show how stim type (Audio or Haptic) interact with arousal rate to predict behavior outcomes
  The second pass of analysis will show how entrainment (vs sham) interact with arousal rate to predict behavior outcomes

main_effects_lme_h1_h7:
  This script looks at the main effects of
    i) Stim Type (Audio vs Haptic)
    ii) Auditory entrainment (Audio vs Sham)
    iii) Haptic entrainment (Haptic vs Sham)
    iv) Entrainment (Both groups vs Sham)
  on the behavior outcomes, using a linear mixed model.

year 2 study data.csv:
  This is the file with only behavioral outcomes that i have been using for correlations.
