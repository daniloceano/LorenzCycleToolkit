
track_dir='../../SWSA-cyclones_energetic-analysis/tracks_LEC-format/intense/'
era5_dir='../../SWSA-cyclones_energetic-analysis/met_data/ERA5/DATA/'

for file in "$track_dir"*
do

    filename=$(basename "$file")
    id=$(echo "$filename" | sed 's/track_//')
    era5_file="${era5_dir}${id}_ERA5.nc"
    
    echo "\nFile ID: $id"
    echo "ERA5 data: $era5_file"
    
    cp $file ../inputs/track
    
    python lorenz-cycle.py $era5_file -t -r -g

done
