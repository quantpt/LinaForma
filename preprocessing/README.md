# E3_loop_THERIAK
This code allows the user to run THERIAK in a loop and save the results to a CSV file. The pressure-temperature path can be changed in the PT_path.csv file.

### NOTE: YOU WILL NEED BASH INSTALLED ON YOUR COMPUTER TO RUN THIS CODE

The package is split into two parts:

## 1) E3a_create_grid.m
If the user wants to run THERIAK in a P-T grid, this code allows the user to create an appropriate PT_path.csv file. The dimensions of the grid can be chosen in the input section of the code.

## 2) E3b_automate_THERIAK.m
This code runs THERIAK at the P-T points given in the PT_path.csv file. 

### Data Inputs
1) The P-T path of interest must be specified in the PT_path.csv file.
2) Please create a folder (perhaps called THERIAK_automation) which contains all of the THERIAK-DOMINO programs, a THERIN file, and the database of interest.

### Once you have added the THERIAK-DOMINO programs, please make sure the theriak.ini settings look like this:
 <p align="center">
<img src="https://github.com/TMackay-Champion/LinaForma/blob/08dfa78d9dfcccdae68b7bf9e0d73bd3189681c3/images/E3_settings.png", width="80%">
</p>

### E3b script input
% ====== Run parameters ======\
**directory = '?'**\
This needs to be the full-path to the folder containing your THERIAK-DOMINO programs. Please use a forward slash ('/') to separate sub-directories.\
**database = '?'**\
This should match the name of the database file in the automation folder.\
**path = '?'**\
This should be a CSV file with the P-T path information. Grids can be created using script E3a.\
**mac = ?**\
Are you using a MAC or WINDOWS machine? Mac = 1, Windows = 0.


### E3b script output
The code will output a CSV file containing all of the THERIAK data.  This CSV file is appropriate for all further scripts in the LinaForma package. The CSV file will look like this:

 <p align="center">
<img src="https://github.com/TMackay-Champion/LinaForma/blob/b17f80121984bd8460ea58220c26b06a878bd627/images/E3_output.png", width="80%">
</p>
