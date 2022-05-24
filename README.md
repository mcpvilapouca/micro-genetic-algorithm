# micro-genetic-algorithm

Micro-Genetic algorithm for the optimization of constitutive models material parameters
-Runs in matlab, through Main.m
-Requires python3 and abaqus
-to run from the terminal in background use the following command:
nohup matlab -r -nosplash -nodesktop -nojvm < Main.m > output.txt &

General explanation:

It minimizes the least square differences between an experimental curve
and numerical results from ABAQUS.

Results from the example:


Notes on implementation:

-In the Main.m file it is necessary to change:
    - no. parameters to optimize
    - dimension of the population (tipically 5)
    - no. decimal cases of the parameters
    - domain of each parameter
    - random initial population (yes/no)
    - absolute convergence parameters
    - elitism option ('y'/'n') --> we recommend using elitism for better results
    - no. of times the micro-genetic algorithm should be run. Usually 1 time is enough.

-In the merit.m file change:
    - name of the file to be run in abaqus
    - name of the user-define subroutine to be applied

In each folder xr1 through xrn (n=dimension of the population) change:
    - param_orig.inp (to have the parameters of the required model. Note that
       the parameters that need to be changed must be equal to 1.0 ,i.e. C10=1.0, to be consistent
       with the file subs_param.py)
    - change the ABAQUS input file name in merit.m
    - change the umat name in merit.m (note: to run faster, compile the umat first using abaqus make library=your_umat.for. This creates a file: your_umat-std.o which can be used instead of umat.for in the merit.m function)
    - create your get_rpt.py to generate the appropriate report (to do this, run your input file once,
      open the viewer, create the necessarry x-y results then go to Report and generate a .rpt file with
      the name report.rpt. Then close the viewer and change the file abaqus.rpy to get_rpt.py, changing
      the extension name)
    - see the number of string lines the report.rpt has and change the get_num.py accordingly, to delete those lines
    - change the file subs_param.py to address your parameters
    - modify the experimental.txt to consider the required experimental curve (be aware that the number of
    experimental and numerical points has to be the same!)