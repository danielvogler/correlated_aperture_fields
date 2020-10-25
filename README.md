# correlated_aperture_fields

generate spatially correlated aperture fields

### Usage
- Generate correlated aperture fields with settings in <code>generate_correlated_field.R</code>. This generates a txt-file with the aperture field as entries.
- Visualize results with: <code>python read_variogram.py name_of_generated_field_file.txt</code>

### Dependencies
- R libraries (tidyr, gstat)

### Credits 
Daniel Vogler



![Example image field](/images/aperture_field_example.png)
![Example image distribution](/images/aperture_distribution_example.png)
