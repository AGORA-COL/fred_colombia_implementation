This script aims to estimate the contact time series for all departments based on the contact time series and mobility trends of Bogotá. By leveraging a proportion between the two metrics, it extrapolates the contact time series for the rest of the departments using their respective mobility trends. Below is an explanation of the code, broken down into its main components.

### Load and Preprocess Data

1. **Community Survey Data for Bogotá (Department Code 11001)**:
   - The script starts by loading a CSV file containing community survey data for Bogotá. This dataset includes dates and the total number of contacts reported on those dates.

2. **Mobility Trends Data for Bogotá**:
   - Next, it loads mobility trends data for Bogotá, converting the 'date' column to datetime format.
   - The mobility trend values are adjusted by calculating an exponentially weighted moving average (EWMA) with a smoothing factor `alpha`, and adding 1 to the result. This adjustment likely aims to normalize the mobility trend data and make it more comparable across time.


### Estimation Process for Other Departments


**Combining and Estimating Contacts Data**:
   - The script calculates the ratio of the EWMA of the department's mobility trend to the EWMA of Bogotá's mobility trend. This ratio aims to capture the relative mobility between the department and Bogotá over time.
   - The contact estimation for the department is then calculated by multiplying this ratio by the total number of contacts from the Bogotá dataset. This step assumes that the contact dynamics in each department can be approximated by adjusting the Bogotá contact data based on the relative mobility trends.