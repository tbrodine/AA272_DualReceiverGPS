# AA272_DualReceiverGPS
Final Project MATLAB code repository for AA 272 at Stanford

-driver.m calls the gnss data preprocessing scripts to obtain quantities of interest for the solver. This script saves these QOIs to a CSV file. \
-PreProcessGNSS.m calls the ProcessGnssMeasScript.m script from the open source GPS Tools repository from Google, then returns the QOIs. \
-Final_Project_Solver_1phone.m solves for the position solution using one phone's gnss data. \
-Final_Project_Solver_2phones.m solves for the position solution using 2 phones' gnss data. \
-The data directory contains CVS outputted by driver.m with the quantities of interest to be used in the solver. \
\
The GPS Tools library used to pre-process the raw GNSS measurements can be found at: https://github.com/google/gps-measurement-tools. It is not inlcuded in this repository. 
