## Log with changes in the program 

### Thursday, 26th March: First Upload

### Monday, 30th March 
1. Progress bar indicating the status of the program *MPC_progress.m*
2. Visualising 3D POD modes and saving each one in specified directory *podmodes3dim.m*
3. Modal analysis, allowing to visualise mode shape through time *modeanimation.m*

### Monday, 6th April
1. *preprocessdmdid* and *preprocessdmdval* functions created, in other to reuse same scaling factors in identification data in validation data scaling.  
2 *MPC_wakesteering*: state reconstruction is computed by projecting the low order states into the high order states using the left singular vectors of snapshot matrices.
3. Animation functions changed to allow handling of different data sets.  
4. *suptitle* function replaced by function *suplabel.m* added to *Functions* directory. *suptitle* function is part of the BioInformatics toolbox, which is now not necessary.
5. *evaluatetimevaryingerror.m* uses normalized root mean squared error (NRMSE) to evaluate states reconstruction deviations
6. Possibility of performing SVD on non linear observables (function of main states). *koopman* variable must be set to 1 in the beginning and function *koopmanextension.m* selects the non linear observables - functions of the state or other flow field components - to be included alongside main states;
7. Flow field is also reconstructed with validation data. Variable *xstatesvalid* is computed and compared with SOWFA's validation flow field (if available).
8. Program automatically saves 3 important results: the FIT for the two turbines for validation models, the reconstruction error for the states and the main results arising from the SVD. These allows to compare the variables for different methods used.

### Monday, 13th April

### Monday, 20th April

### Monday, 27th April

### Monday, 4th May

### Monday, 11th May

### Monday, 18th May

### Monday, 25th May

### Monday, 1st June

### Monday, 8th June

### Monday, 15th June

### Monday, 22th June

### Monday, 29th June

### Monday, 6th July

### Monday, 13th July

### Monday, 20th July

### Monday, 27th July