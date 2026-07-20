# EyeCal
Code for the receptive fields' eye movement correction paper `Subspace reverse-correlation estimation of receptive fields during free viewing` published in the Journal of Vision. To use the code, 
1. Download the data from the figshare repository [https://doi.org/10.6084/m9.figshare.32257377](https://doi.org/10.6084/m9.figshare.32257377) and place it under ./data/ directory where the code is.
2. Open eye_cal file and specify the experiment name.
3. Run eye_cal

The experiment databases with the .db extensions include a Matlab table named "log" with the following structure:  

      frame       sign    kx     ky      sbxframe      xpos        ypos        area        speed           signal     
    __________    ____    ___    ___    __________    _______    _________    ______    ___________    _______________

             0     -1       0    -10     4.527e+05     2.3999     -0.13561    191.25              0    {10×505 double}
            12     -1      -4     -9    4.5818e+05     2.4495     -0.15445     190.9              0    {10×505 double}
            .      .        .      .      .              .            .          .                .        .  
            .      .        .      .      .              .            .          .                .        .
            .      .        .      .      .              .            .          .                .        .
            .      .        .      .      .              .            .          .                .        .
            
  In this data structure, each row is associated with a trial in the experiment.
  
  Referring to the paper, we have k,l, and s parameters that specify the stimuli. In the above data set, k="kx", l="ky", s="sign".
  
  "xpos" and "ypos" are the eye positions in the camera coordinates.
  
  "signal" specifies the neuronal responses. In each cell, there are 10 time bins (each trial is 200ms and is divided into intervals of 20ms) and N neurons. In the above example, we have 505 total neurons.
            
