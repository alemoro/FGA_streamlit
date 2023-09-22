import streamlit as st
import io
import numpy as np
import pandas as pd
from scipy.optimize import leastsq
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 'truetype'
import seaborn as sns
sns.set_theme(style='ticks', rc={'axes.spines.right':False, 'axes.spines.top':False})
from PIL import Image

icon = Image.open("./img/Neurospector.png")

st.set_page_config(
    page_title="FGA ELISA analysis",
    page_icon=icon,
    layout='wide',
)

# Start by defining the functions for the standards
def logistic4(x,a,b,c,d):
    """
    The logistic4 function calculates the value of a sigmoidal (S-shaped) curve
    defined by four parameters: a, b, c, and d.

    Parameters:
    - x (float): The input value at which to calculate the logistic function.
    - a (float): The upper asymptote (maximum value) of the logistic curve.
    - b (float): The Hill slope or steepness of the curve.
    - c (float): The x-value of the sigmoid's midpoint.
    - d (float): The lower asymptote (minimum value) of the logistic curve.

    Returns:
    float: The value of the 4-parameter logistic function at the input point x.
    """
    return ((a-d) / (1+((x/c)**b))) + d

def residual(p, y, x):
    """
    This function calculates the residuals by subtracting the predicted values of a
    4-parameter logistic curve from the observed data at corresponding input points.

    Parameters:
    - p (tuple of floats): A tuple containing the four logistic curve parameters (a, b, c, d).
    - y (array-like): Observed data values (dependent variable).
    - x (array-like): Input values corresponding to the observed data (independent variable).

    Returns:
    array-like: An array of residuals, representing the differences between observed data and
                the predicted values of the 4-parameter logistic curve.
    """
    a,b,c,d = p # Extract the logistic curve parameters
    # Calculate the residuals by subtracting predicted values from observed data
    err = y-logistic4(x,a,b,c,d)
    return err

def peval(x, p):
    """
    This function calculates the predicted values of a 4-parameter logistic curve
    using the given parameters `p` at specified input points `x`.

    Parameters:
    - x (array-like): Input values at which to evaluate the logistic curve.
    - p (tuple of floats): A tuple containing the four logistic curve parameters (a, b, c, d).

    Returns:
    array-like: An array of predicted values representing the 4-parameter logistic curve
                evaluated at the specified input points.
    """
    a, b, c, d = p  # Extract the logistic curve parameters
    # Calculate the predicted values of the 4-parameter logistic curve
    return logistic4(x, a, b, c, d)

def fit_concentrations(y,p):
    """
    This function applies a logistic function to fit concentrations (`x`) based on observed
    values (`y`) and a set of logistic function parameters (`p`).

    Parameters:
    - y (float): Observed value (dependent variable).
    - p (tuple of floats): A tuple containing the four logistic curve parameters (a, b, c, d).

    Returns:
    float: The fitted concentration value.
    """
    a,b,c,d = p
    # solving the logistic function per x
    x = c * ((a-d)/(y-d) - 1)**(1/b)
    return x

def calculate_standards():
    # first unpack the session states
    elisa_df = st.session_state.elisa_df
    batch_col = st.session_state.ELISA_batch_column
    data_col = st.session_state.ELISA_data_column
    condition_col = st.session_state.ELISA_condition_column
    standard_name = st.session_state.ELISA_standard_name
    concentration_col = st.session_state.ELISA_concentration_column
    # Create the spaceholder for the coefficients and get the data
    batch_IDs = elisa_df[batch_col].unique()
    test_coeff = np.empty((len(batch_IDs), 4))
    test_values = np.array(elisa_df[data_col])
    # Per batch calculate the standards coefficients
    for idx, batch in enumerate(batch_IDs):
        standards = np.array(elisa_df.loc[((elisa_df[batch_col]== batch) & (elisa_df[condition_col]==standard_name), data_col)])
        concentrations = np.array(elisa_df.loc[((elisa_df[batch_col]== batch) & (elisa_df[condition_col]==standard_name), concentration_col)])
        concentrations = concentrations[~np.isnan(standards)]
        standards = standards[~np.isnan(standards)]
        # Calculate the standard curve
        p0 = [0, 1, 1, 1]
        temp_coeff = leastsq(residual, p0, args=(standards, concentrations))
        test_coeff[idx,:] = temp_coeff[0]
        # Calculate the concentration values
        values_fltr = (elisa_df[batch_col]==batch) & (elisa_df[condition_col] != standard_name)
        test_values[values_fltr] = fit_concentrations(test_values[values_fltr]*st.session_state.ELISA_dilution_factor, temp_coeff[0])
    # Add the calculated protein concentration to the df
    st.session_state.elisa_df['Protein'] = test_values
    st.session_state.elisa_test_coeff = test_coeff
    

#for key in st.session_state.keys():
#    del st.session_state[key]

# Add a quick start on the sidebar
with st.sidebar:
    logo = Image.open("./img/FGA_Neurospector.png")
    st.image(logo)
    with st.expander("Quick start"):
        st.markdown('''
            This page will giude you through the analysis of ELISA data analysis.  
            The data should be a *.csv :clipboard: files with each row as observation.
                 There has to be one column with an batch indication, one with the ELISA readout, one with the condition, and one with the standard concentrations. Additional column are possible and optional.  
              
            The second step is to indicate if there is any dilution factor :droplet:. If there is no dilution set the dilution factor to "1".  
                    At this stage you can also combined multiple columns containing different categories into one (i.e. one column for genotype and one column for treatment).  
              
            The next step is to calulate the standards:
            - Select the variable with the batch identifier
            - Select the variable with the ELISA readout
            - Select the variable with the condition where the standards are stored, and select the name for the standards
            - Select the variable with the concentrations  
              
            Now the last step is to calculate the standard curves! :sunglasses:  
              
            Finally, to plot the data :bar_chart::
            - Select the variable with the conditions, and the name of the "control" condition (i.e. WT)
            - Select the other conditions that you would like to plot (in the order that you would like them to be plotted)
            - Select if you want to calculate the batch normalized values
            - Select if you want to have the plot as a bar graph or a boxplot
                    
            At this stage you can also download the plot and the data :floppy_disk:!
        ''')

# Main body of the page
if "new_ELISA" not in st.session_state:
    st.session_state.new_ELISA = True

# Data import component: load a csv file
st.subheader('Import data', divider='rainbow')
with st.expander("Data input"):
    st.markdown("""
        Here you can visualize your input data.  
        This could be the outcome of the sCaSpA (*spontaneous Ca2+ Spike Analysis* in MATLAB) or any other analysis tool as long as it is saved as described below.  
        **Important** one of the column in your dataset *must* be your grouping, or condition, variable.
    """)
elisa_input = st.file_uploader("Choose your file", type=["csv"])

# Data preparation component: adjust the data if needed
if elisa_input is not None:
    # Read the data into a pandas dataframe
    if st.session_state.new_ELISA:
        st.session_state.new_ELISA = False
        st.session_state.elisa_df = pd.read_csv(elisa_input)
    var_names = [var for var in st.session_state.elisa_df.keys()]
    name_list = np.append(["##"], var_names)
    # Before doing anything ask if the sample are diluted, and if the conditions are a combination of different columns
    col1, col2, col3 = st.columns([1,3,1])
    with col1:
        st.number_input("ELISA dilution factor", min_value=1, max_value=10, step=0.5, value=2, key="ELISA_dilution_factor")
    with col2:
        conditions_parts = st.multiselect("Select the columns to for a new condition", name_list)
    with col3:
        if(st.button("OK")):
            new_condition = st.session_state.elisa_df[conditions_parts[0]]
            conditions_parts.pop(0)
            for part in conditions_parts:
                new_condition = new_condition + " " + st.session_state.elisa_df[part]
            st.session_state.elisa_df["NewCondition"] = new_condition
    if st.toggle('Show data'):
        st.write(st.session_state.elisa_df)

# Standard preparation component: where is the relevant data
if elisa_input is not None:
    col1, col2, col3, col4 = st.columns(4)    
    with col1:
        idx = 0
        if "ELISA_batch_column" in st.session_state:
            idx = int(np.where(name_list == st.session_state.ELISA_batch_column)[0])
        st.selectbox("Batch variable", name_list, index=idx, key="ELISA_batch_column")
    with col2:
        idx = 0
        if "ELISA_data_column" in st.session_state:
            data_idx = int(np.where(name_list == st.session_state.ELISA_data_column)[0])
        st.selectbox("Data variable", name_list, index=idx, key="ELISA_data_column")
    with col3:
        idx = 0
        if "ELISA_condition_column" in st.session_state:
            condition_idx = int(np.where(name_list == st.session_state.ELISA_condition_column)[0])
        st.selectbox("Condition variable", name_list, index=idx, key="ELISA_condition_column")
        if st.session_state.ELISA_condition_column != "##":
            condition_IDs = st.session_state.elisa_df[st.session_state.ELISA_condition_column].unique()
            standard_ID = st.radio("Select the name of the standard", options=condition_IDs, key="ELISA_standard_name")
    with col4:
        idx = 0
        if "ELISA_concentration_column" in st.session_state:
            concentration_idx = int(np.where(name_list == st.session_state.ELISA_concentration_column)[0])
        concentration_col = st.selectbox("Standard concentrations", name_list, index=idx, key="ELISA_concentration_column")

# Standard plotting components
st.subheader('Standard curve calculation and ploting', divider=True)
with st.expander("Explanation of standard calculation"):
    st.markdown("""
                Now I have to think of an explanation
                """)
if elisa_input is not None:
    elisa_options = [st.session_state.ELISA_batch_column, st.session_state.ELISA_data_column, st.session_state.ELISA_condition_column, st.session_state.ELISA_concentration_column, "##"]
    if (elisa_options.index("##") == 4):
        st.button('Calculate standards', key="ELISA_calculate_standard", on_click=calculate_standards)
        if st.toggle('Show data', key="toggle1"):
            st.write(st.session_state.elisa_df)
        if st.toggle('Plot the standards', key="toggle2"):
            elisa_df = st.session_state.elisa_df
            batch_col = st.session_state.ELISA_batch_column
            data_col = st.session_state.ELISA_data_column
            condition_col = st.session_state.ELISA_condition_column
            standard_name = st.session_state.ELISA_standard_name
            concentration_col = st.session_state.ELISA_concentration_column
            batch_IDs = elisa_df[batch_col].unique()
            # Plot the standard curves
            fig_standards, axs = plt.subplots(1,len(batch_IDs), figsize=(15,5), facecolor='w', edgecolor='k')
            cmap = sns.color_palette(palette='Accent', n_colors=len(batch_IDs))
            for idx, batch in enumerate(batch_IDs):
                standards = np.array(elisa_df.loc[((elisa_df[batch_col]== batch) & (elisa_df[condition_col]==standard_name), data_col)])
                concentrations = np.array(elisa_df.loc[((elisa_df[batch_col]== batch) & (elisa_df[condition_col]==standard_name), concentration_col)])
                concentrations = concentrations[~np.isnan(standards)]
                standards = standards[~np.isnan(standards)]
                test_curve = peval(np.arange(0, 200), st.session_state.elisa_test_coeff[idx,:])
                # Calculate the R-square
                sort_index = np.argsort(concentrations)
                sort_standards = standards[sort_index]
                sort_estimate = peval(concentrations[sort_index], st.session_state.elisa_test_coeff[idx,:])
                r_square = r2_score(sort_standards, sort_estimate)
                # Plot the data
                axs[idx].plot(concentrations, standards, 'ok')
                axs[idx].plot(np.arange(0, 200), test_curve, color=cmap[idx])
                axs[idx].set_xlabel('Concentration (ng/ml)')
                axs[idx].set_ylabel('Absorbance (a.u.)')
                axs[idx].set_title(batch + " - R2: " + "{:.3f}".format(r_square))
            st.pyplot(fig_standards)
            img_standards = io.BytesIO()
            plt.savefig(img_standards, format='pdf')
            st.download_button(label="Download plot", data=img_standards, file_name="ELISA_Standards_Plot.pdf", mime="application/pdf", key="ELISA_save_figure1")

st.subheader('ELISA plot', divider=True)
with st.expander("Explanation for plots"):
    st.markdown("""
                Now I have to think of an explanation
                """)
if elisa_input is not None:
    elisa_options = [st.session_state.ELISA_batch_column, st.session_state.ELISA_data_column, st.session_state.ELISA_condition_column, st.session_state.ELISA_concentration_column, "##"]
    if (elisa_options.index("##") == 4):
        if "elisa_test_coeff" in st.session_state:
            # Get the ELISA data without the standards
            data_df = st.session_state.elisa_df.copy()
            data_df = data_df.loc[(data_df[st.session_state.ELISA_condition_column]!=st.session_state.ELISA_standard_name),:]
            col1, col2, col3 = st.columns([1,3,1])
            with col1:
                st.selectbox("Condition selection", name_list, index=0, key="ELISA_condition_plot")
                if st.session_state.ELISA_condition_plot != "##":
                    condition_IDs = data_df[st.session_state.ELISA_condition_plot].unique()
                    control_ID = st.selectbox("Select control condition", condition_IDs, key="ELISA_control_name")
            with col2:
                if st.session_state.ELISA_condition_plot != "##":
                    conditions_order = st.multiselect("Select the plotting order", condition_IDs, default=control_ID)
            with col3:
                st.write("##")
                b_normalize = st.checkbox("Normalize")
                b_bar = st.checkbox("Bar graph")
            if st.toggle('Plot data', key='toggle3'):
                fig_elisa, ax = plt.subplots(figsize=(15,10))
                batch_col = st.session_state.ELISA_batch_column
                batch_IDs = np.array(data_df[batch_col])
                batches = np.unique(batch_IDs)
                cmap = sns.color_palette(palette='Accent', n_colors=len(batches))
                y_data = 'Protein'
                if b_normalize:
                    condition_col = st.session_state.ELISA_condition_plot
                    condition_IDs = np.array(data_df[condition_col])
                    control_filter = condition_IDs == control_ID
                    temp_values = np.array(data_df['Protein'])
                    for batch in batches:
                        batch_filter = batch_IDs == batch
                        control_value = temp_values[(batch_filter) & (control_filter)]
                        control_value = np.mean(control_value)
                        temp_values[batch_filter] = temp_values[batch_filter] / control_value
                    data_df['NormalizeValues'] = temp_values
                    y_data = 'NormalizeValues'
                if b_bar:
                    sns.barplot(data=data_df, x=st.session_state.ELISA_condition_plot, y=y_data, order=conditions_order, edgecolor='k', facecolor=(.7, .7, .7), errcolor='k')
                else:
                    PROPS = {
                            'boxprops':{'facecolor':(.7, .7, .7), 'edgecolor':'black'},
                        }
                    sns.boxplot(data=data_df, x=st.session_state.ELISA_condition_plot, y=y_data, order=conditions_order, **PROPS)
                sns.swarmplot(data=data_df, x=st.session_state.ELISA_condition_plot, y=y_data, order=conditions_order, hue='BatchID', palette=cmap)
                plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
                if b_normalize:
                    ax.set_ylabel('Concentrations (normalized)')
                else:
                    ax.set_ylabel('Concentration (ng/ml)')
                st.pyplot(fig_elisa)
                img_plot = io.BytesIO()
                plt.savefig(img_plot, format='pdf')
                download_df = data_df.to_csv().encode('utf-8')
                col1, col2 = st.columns(2)
                col1.download_button(label="Download plot", data=img_plot, file_name="ELISA_Plot.pdf", mime="application/pdf", key="ELISA_save_figure2")
                col2.download_button(label="Download data", data=download_df, file_name="ELISA_Results.csv", mime="text/csv")

