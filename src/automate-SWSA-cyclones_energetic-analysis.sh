
track_dir='../../SWSA-cyclones_energetic-analysis/tracks_LEC-format/intense/'
era5_dir='../../SWSA-cyclones_energetic-analysis/met_data/ERA5/DATA/'

for file in "$track_dir"*
do

    filename=$(basename "$file")
    id=$(echo "$filename" | sed 's/track_//')
   
    echo "$file"
    echo "$id"

    echo "File ID: $id"
    cp $file ../inputs/track
    
    era5_file="${era5_dir}${id}_ERA5.nc"
    python lorenz-cycle.py $era5_file -t -r -g    

done
