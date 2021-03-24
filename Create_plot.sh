#This progran is created to run the scripts required to create the figure 3.2c
#The first 4 scripts process the data and the final creates the plot.
#There is no need to re-run the processing scripts if they have completed succesfully
#When running this script the path to the main directory is required as an arguement.
#. Create_plot.sh path_to_maindirectory

echo 'Program has started'
echo 'The main directory is'
echo $1
mkdir -p $1/SavedModData
cd $1
echo running Process_CMIP5.sh
. scripts/Process_CMIP5.sh $1
cd $1
echo running Process_CMIP6.sh
. scripts/Process_CMIP6.sh $1
cd $1
echo running Process_CESM.py
python3 scripts/Process_CESM.py $1
cd $1
echo running Process_volcs.py
python3 scripts/Process_volcs.py $1
cd $1
echo running plot_past1000.py
python3 scripts/plot_past1000.py $1
echo 'Program has finished'
