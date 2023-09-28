import io
import streamlit as st
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 'truetype'
import seaborn as sns
sns.set_theme(style='ticks', rc={'axes.spines.right':False, 'axes.spines.top':False})
from scipy.stats import ttest_ind
from scipy import linalg
from statsmodels.stats.multitest import multipletests
from PIL import Image

# preprocessing
from sklearn import preprocessing
from sklearn.preprocessing import LabelEncoder, scale

from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
#import umap

# training
from sklearn.utils import resample
from sklearn.model_selection import train_test_split

# metrics
from sklearn.metrics import confusion_matrix, accuracy_score


def process_LDA_data(temp_df, varIdx, columnIdx):
    lda_data = temp_df.iloc[:,varIdx]
    lda_labels = temp_df.loc[:,genotype_col]
    # Returns a scaled version of your sliced data 
    scaler = preprocessing.StandardScaler()
    lda_data = scaler.fit_transform(lda_data)
    lda_data = np.nan_to_num(lda_data, nan=0.0001)
    return lda_data, lda_labels

def plot_LDA_biplot(lda, lda_labels, lda_coef, colors, feat_names, ax=None):
    sum_components = np.sum(abs(lda_coef), axis=1)
    sort_index = np.argsort(sum_components)
    scalex = 1.0 / (lda[:, 0].max() - lda[:, 0].min())
    scaley = 1.0 / (lda[:, 1].max() - lda[:, 1].min())

    if ax is None:
        ax = plt.gca()

    sns.scatterplot(x=lda[:, 0] * scalex, y=lda[:, 1] * scaley, hue=lda_labels, palette=colors, ax=ax)
    for s in sort_index[0:5]:
        ax.arrow(0, 0, lda_coef[s, 0], lda_coef[s, 1], color='r')
        ax.text(lda_coef[s, 0] * 1.5, lda_coef[s, 1] * 1.5, feat_names[s], ha='center', va='center')
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    
    return ax, sort_index

def calculate_fold_change(temp_df, genotype_col, control_group, var_idx):
    '''
    Calculates fold change values for different groups compared to a control group.
    Args:
        temp_df (pandas.DataFrame): Input DataFrame containing the data.
        genotype_col (str): Name of the column in `temp_df` representing the genotype.
        control_group (str): Name of the control group to compare against.
        var_idx (tuple): A tuple specifying the range of indices to consider in the DataFrame.

    Returns:
        pandas.DataFrame: A DataFrame containing the fold change values.
    '''
    # Get the data that we want to use for the dataframe
    temp_df = temp_df.dropna()
    #var_idx = np.arange(var_idx[0], var_idx[1])
    genotypes = temp_df[genotype_col].unique()
    # Calculate the mean values for the control group and save it into a new dataframe
    control_mean = np.mean(temp_df[(temp_df[genotype_col] == control_group)].iloc[:,var_idx], axis=0)
    fold_df = pd.DataFrame()
    fold_df['Features'] = temp_df.columns[var_idx]
    fold_df['Group'] = [control_group] * len(fold_df)
    fold_df['Fold Change'] = control_mean.values
    fold_df['p-value'] = control_mean.values * 0
    # Now loop through the genotypes to calculate the fold-change
    for gen in genotypes:
        if gen != control_group:
            # Step 1: calculate the fold change
            temp_fold_change = np.mean(temp_df[(temp_df[genotype_col] == gen)].iloc[:,var_idx], axis=0) / control_mean
            # Step 2: Perform statistical tests, adjust for multiple testing
            _, p_values = ttest_ind(temp_df[(temp_df[genotype_col] == gen)].iloc[:,var_idx], temp_df[(temp_df[genotype_col] == control_group)].iloc[:,var_idx])
            adjusted_p_values_1 = multipletests(p_values, method='fdr_bh')[1]
            # Step 3: Store the data into the fold_df
            gen_df = pd.DataFrame()
            gen_df['Features'] = temp_df.columns[var_idx]
            gen_df['Group'] = [gen] * len(gen_df)
            gen_df['Fold Change'] = np.log2(temp_fold_change.values)
            gen_df['p-value'] = -np.log10(p_values)
            fold_df = pd.concat([fold_df, gen_df])
    
    st.session_state.groupped_data = fold_df


def convert_df(temp_df):
    return temp_df.to_csv().encode('utf-8')

icon = Image.open("./img/Neurospector.png")

st.set_page_config(
    page_title="FGA ELISA analysis",
    page_icon=icon,
    layout='wide',
)

# Initialize the session
if "new_session" not in st.session_state:
    st.session_state["new_session"] = True
    st.session_state["multiple_group"] = False
    st.session_state.img_df = []
    st.session_state.groupped_data = []
    st.session_state.xMax = 0
    st.session_state.yMax = 0
    st.session_state.biplot = True
    st.session_state.lda_feature_list = False

varIdx = [8, 93]

# Create the basic layout of the page
st.title('Multidimensional data analysis', )
with st.sidebar:
    logo = Image.open("./img/FGA_Neurospector.png")
    st.image(logo)
    with st.expander("Quick start"):
        st.write('''
            This web app will giude you through the analysis of multidimensional data analysis such as for **calcium imaging** data.  
            To use the app you need to load a dataset as a *.csv file, select which variable are important/relevant, and finally select the grouping and filtering variables.  
            The app will automatically calculate the fold change for all the variable that you selected, and plot them as a volcano plot :volcano:.  
            In addition the app with calculate the LDA (*linear discriminant ananlysis*) and plot the data as a scatter plot of the first two components,
            while listing the features in order of absolute weight.
        ''')
    add_input = st.file_uploader("Choose a *.csv file with the data you would like to analyze.", type=["csv"])
    if len(st.session_state.img_df) > 0:
        st.checkbox("Add more grouping factors?", key="multiple_group")
        var_names = [var for var in st.session_state.img_df.keys()]
        options = st.multiselect('Select informative columns', var_names)
        varIdx = st.session_state.img_df.columns.isin(options)
        
    else:
        st.checkbox("Add more grouping factors?", key="multiple_group")
        options = st.multiselect('Select informative columns', [])
        
if add_input is not None:
    st.session_state.img_df = pd.read_csv(add_input)

# Inspect the raw data
st.header('Raw data input', divider=True)
with st.expander("Data input"):
    st.write('''
        Here you can visualize your input data.  
        The data should be a **.csv* files with each raw as observation and each columns as features or descriptors.
        This could be the outcome of the sCaSpA (*spontaneous Ca2+ Spike Analysis* in MATLAB) or any other analysis tool as long as it is saved as described below.  
        **Important** one of the column in your dataset *must* be your grouping, or condition, variable.
    ''')

if len(st.session_state.img_df) > 0:
    img_df = st.session_state.img_df
    if st.toggle("Show raw data", key="Multidimension_RawToggle"):
        st.write(img_df)
    # Get the variables names and ask the user to pick the important ones
    var_names = [var for var in img_df.keys()]
else:
    var_names = []

# Select the gouping and filtering variables
col1, col2, col3 = st.columns(3)
control_group = []
selection_group1 = []
selection_group2 = []
with st.expander("Groups selection"):
    st.write('''
    Here you can select how you would like to group your data.  
    The **first** grouping variable is your conditions, or genotype. This grouping variable define how you would like to divide your data. The radio button will the allow you to select which values is your *control* condition.  
    The **second** and **third** grouping variables are filters, and are available only after enabeling the "Add more grouping factor" in the side bar (:point_left: :white_check_mark:). With those filters you can select which data should be taken into account. Examples of those filters could be:
    - Naive or overexpression
    - Different treatments (inhibitors, diffirent imaging solution)
    - Different culture days
    ''')
    st.write('''
    The radio button for those variable allow you to select which condition you would like to look at.  
    **Important** use the dropdown menu in order of importance, and do not use the *third* drop down before you selected the *second*.
    ''')
with col1:
    genotype_col = st.selectbox("Select the variable with the conditions", np.append(["##"], var_names), index=0, key="sel_genotype")
    if genotype_col != "##":
        genotypes = img_df[genotype_col].unique()
        control_group = st.radio("Select the control condition", options=genotypes, key="radio_genotype")
with col2:
    group1_col = st.selectbox("Select a filter variable", np.append(["##"], var_names), index=0, key="sel_group1", disabled= not st.session_state.multiple_group)
    if group1_col != "##":
        group1 = img_df[group1_col].unique()
        selection_group1 = st.radio("Select the condition that you would like to use", options=group1, key="radio_group1", disabled= not st.session_state.multiple_group)
with col3:
    group2_col = st.selectbox("Select an additional filter varible", np.append(["##"], var_names), index=0, key="sel_group2", disabled= not st.session_state.multiple_group)
    if group2_col != "##":
        group2 = img_df[group2_col].unique()
        selection_group2 = st.radio("Select the condition that you would like to use", options=group2, key="radio_group2", disabled= not st.session_state.multiple_group)

if len(st.session_state.img_df) > 0:
    use_df = img_df.copy()
    if group1_col != "##":
        use_df = img_df.loc[img_df[group1_col]==selection_group1,:]
    if group2_col != "##":
        use_df = use_df[(use_df[group2_col]==selection_group2)]
    st.button("Calculate the fold change", key="fold_change", on_click=calculate_fold_change, args=(use_df, genotype_col, control_group, varIdx), help="You need to select some features from the dropdown menu in the side bar before calculate the fold chage!")

# Show and explain how the fold change is calculated
st.header('Fold change', divider=True)
with st.expander("Fold change explanation"):
    st.markdown("""
                Now that the app knows which features are important, and which grouping factor you want, it will calculate the fold change and the corrected p-value.  
                **Fold change** is a measure describing how much a quantity changes between an original and a subsequent measurement.
                It is defined as the ratio between the two quantities; for quantities A and B the fold change of B with respect to A is: ${\dfrac{B}{A}}$.  
                In other words, a change from 30 to 60 is defined as a fold-change of 2. This is also referred to as a "one fold increase". Similarly, a change from 30 to 15 is referred to as a "0.5-fold decrease".  
                  
                In this app the fold change is calculated with the following steps:
                - Filter the data accordingly to the two user filters
                - Calculate the mean values of the user specified features for the *control* group
                - Calculate the mean values of the user specified features for the other groups.  
                  
                The **p-value** is automatically calculated per combination of feature and group, using the *scipy* library.  
                The values are the corrected for multiple testing using the **fdr_bh** correction from the *statsmodels* library.  
                In the *Benjamini-Hochberg* correction, the level $a$ is calculated so that by plottingh the p-value of the test *k* over *k* and plotting a line with slope ${\dfrac{a}{m}}$, the true discovery are the ones above the line.
                """)

if len(st.session_state.groupped_data) > 0:
    if st.toggle("Show fold change table", key="Multidimension_Foldchange_Data_Toggle"):
        st.session_state.groupped_data
        download_df = convert_df(st.session_state.groupped_data)
        st.download_button(label="Download data", data=download_df, file_name="NetworkAnalysis_FoldChange.csv", mime="text/csv")

# Plot the data indicating the significant threshold as a dotted line
st.subheader('Volcano plot :volcano:')
with st.expander("Volcano plot explanation"):
    st.markdown("""
                A **volcano plot** is a type of scatter-plot that is used to quickly identify changes in large data sets.  
                The volcano plot is constructed by plotting the negative ${\log_{10}(p-value)}$ on the y axis. This results in data points with low p values (highly significant) appearing toward the top of the plot.  
                The x axis is the ${\log_{2}(fold\ change)}$ between the two conditions.
                In this particular case, there are two vertical lines and one :red[horizontal line]:
                - The vertical lines represent the -1 and 1 ${\log_2}$, indicating the values that are half or double of the control respectively;
                - The :red[horizontal line] represent the false discovery rate (**FDR**) threshold.  
                The app will plot the data using the *seaborn* library by default.
                """)
    
if len(st.session_state.groupped_data) > 0:
    if st.toggle("Show volcano plot", key="Multidimension_Volcano_Toggle"):
        volcano_data = st.session_state.groupped_data.copy()
        volcano_data = volcano_data[volcano_data["Group"]!=control_group]
        unique_groups = st.session_state.groupped_data["Group"].unique()
        if  len(unique_groups) == 3:
            cmap=['#252525', '#026842', '#5E9BD1']
        else:
            cmap=['#252525', '#006d2c', '#08519c', '#a50f15', '#54278f',
                    '#636363', '#31a354', '#3182bd', '#de2d26', '#756bb1',
                    '#969696', '#74c476', '#6baed6', '#fb6a4a', '#9e9ac8',
                    '#bdbdbd', '#a1d99b', '#9ecae1', '#fc9272', '#bcbddc']
            cmap = cmap[:len(unique_groups)]
        xlim = [-2, 2]
        ylim = [0, 9]
        fig_volcano = plt.figure()
        sns.scatterplot(data=volcano_data, x='Fold Change', y='p-value', hue='Group', palette=cmap[1:])
        plt.axhline(-np.log10(0.05), color='red', linestyle='dotted')
        plt.axvline(-1, color='black', linestyle='dotted')
        plt.axvline(1, color='black', linestyle='dotted')
        if st.session_state.xMax > 0.0:
            plt.xlim(np.array([-1, 1])*st.session_state.xMax)
        if st.session_state.yMax > 0.0:
            plt.ylim([0, st.session_state.yMax])
        st.pyplot(fig_volcano, use_container_width=True)
        img_volcano = io.BytesIO()
        plt.savefig(img_volcano, format='pdf')
        col1, col2, col3, col4 = st.columns(4)
        st.session_state.xMax = int(col1.number_input('X-axis limit'))
        st.session_state.yMax = int(col2.number_input('Y-axis limit'))
        col3.write("##"); col3.button("Relead")
        col4.write("##"); col4.download_button(label="Download volcano", data=img_volcano, file_name="volcanoPlot.pdf", mime="application/pdf")


st.header('LDA analysis', divider=True)
with st.expander("LDA explanation"):
    st.write("""
            The linear discriminant analysis (**LDA**) is a generalization of *Fisher's linear discriminant*, a method used to find a linear combination of features that characterized, or separate, two or more classes of objects.  
             The LDA is closely related to ANOVA and regression analysis in the attempt to extress one dependent variable as a linear combination of other measurements. It is also related to principal component analysis (**PCA**) in the sense that they are both dimensionality reduction analysis.  
             Here, the LDA is calculatelated using the *sklearn* library:
             - First the data is normalized so that each measurement is represented by a distribution of mean 0 and standard deviation 1.
             - Then all NaN and 0 values are removed from the dataset
             - You then need to select the condition that you would like to plot. Note that the LDA will be calculate on those conditions only
             - Depending if the *bootstrap analysis* is selected the model will be calculated over the entire dataset (no bootstrap) or on 70% of the dataset repeated 1000 times (more accurate, but more computational intensive)
               
            To calculate which features are more important in the model:
             - The scaling factors are retrieved from the model
             - The factors are scaled and sorted on the data  
             The factor can then be plotted on top of the scatter plot of the LDA (default) as a **biplot**, or listed underneath the plot.
             In addition to the scatter plot there is also the option to plot the confusion matrix for the LDA itself, and for the prediction of new conditions on top of the conditions used to calculate the LDA.
             If you opted to perform a bootstrap analysis the number that you see on the confution matrix is the average number of observations in the 1000 repetitions.
             """)
if len(st.session_state.img_df) > 0:    
    if st.toggle("Calculate LDA", key="Multidimension_LDA_Toggle"):
        with st.form("LDA Options"):
            unique_groups = use_df[genotype_col].unique()
            conditions_order = st.multiselect("Select the plotting order", unique_groups, default=control_group)
            col1, col2, col3 = st.columns(3)
            st.session_state.biplot = col1.checkbox("Biplot", key="lda_biplot")
            lda_feature_list = col2.checkbox("Show list of features", key="lda_feature_list")
            do_bootstrap = col3.checkbox("Bootstrap analysis")
            submitted = st.form_submit_button("Update")
            if submitted:
                if  len(conditions_order) == 3:
                    cmap=['#252525', '#026842', '#5E9BD1']
                else:
                    cmap=['#252525', '#006d2c', '#08519c', '#a50f15', '#54278f',
                            '#636363', '#31a354', '#3182bd', '#de2d26', '#756bb1',
                            '#969696', '#74c476', '#6baed6', '#fb6a4a', '#9e9ac8',
                            '#bdbdbd', '#a1d99b', '#9ecae1', '#fc9272', '#bcbddc']
                    cmap = cmap[:len(conditions_order)]
                # Perform an LDA
                if len(unique_groups) > len(conditions_order):
                    lda_df = use_df[use_df[genotype_col].isin(conditions_order)]
                else:
                    lda_df = use_df
                lda_data, lda_labels = process_LDA_data(lda_df, varIdx, genotype_col)
                # set lda model and fit to the data
                linear_discriminant = LDA(solver="svd", store_covariance=True)
                lda = linear_discriminant.fit_transform(lda_data, lda_labels)
                if do_bootstrap:
                    lda_confusion = np.zeros((len(conditions_order), len(conditions_order), 1000))
                    for b in range(1000):
                        data_train, data_test, label_train, label_test = train_test_split(lda_data, lda_labels, test_size=0.3)
                        lda_predict = linear_discriminant.fit(data_train, label_train).predict(data_test)
                        lda_confusion[:,:,b] = confusion_matrix(y_true=label_test, y_pred=lda_predict, labels=conditions_order)
                    lda_confusion = np.mean(lda_confusion, axis=2)
                else:
                    lda_predict = linear_discriminant.fit(lda_data, lda_labels).predict(lda_data)
                    lda_confusion = confusion_matrix(y_true=lda_labels, y_pred=lda_predict, labels=conditions_order)
                tab1, tab2, tab3 = st.tabs(['Scatter', 'Confusion', 'Predict'])
                with tab1:
                    fig_lda = plt.figure()
                    if st.session_state.biplot:
                        # Get the coefficients (scaling factors)
                        lda_coef = linear_discriminant.scalings_
                        feat_names = lda_df.iloc[:,varIdx].columns.tolist()
                        sum_components = np.sum(abs(lda_coef), axis=1)
                        sort_index = np.argsort(sum_components)
                        scalex = 1.0 / (lda[:, 0].max() - lda[:, 0].min())
                        scaley = 1.0 / (lda[:, 1].max() - lda[:, 1].min())
                        sns.scatterplot(x=lda[:, 0] * scalex, y=lda[:, 1] * scaley, hue=lda_labels, palette=cmap)
                        for s in sort_index[0:5]:
                            plt.arrow(0, 0, lda_coef[s, 0], lda_coef[s, 1], color='r')
                            plt.text(lda_coef[s, 0], lda_coef[s, 1], feat_names[s], ha='center', va='center')
                        plt.xlim([-1, 1])
                        plt.ylim([-1, 1])
                    else:
                        if len(conditions_order) > 2:
                            sns.scatterplot(x=lda[:, 0], y=lda[:, 1], hue=lda_labels, palette=cmap)
                            plt.ylabel("LD2")
                        else:
                            sns.histplot(x=lda[:, 0], stat='probability', binwidth=0.1, kde=True, hue=lda_labels, palette=cmap)
                            plt.ylabel("Counts")
                    plt.xlabel("LD1")
                    st.pyplot(fig_lda, use_container_width=True)
                with tab2:
                    fig_confusion = plt.figure()
                    sns.heatmap(data = lda_confusion, cmap='mako', annot=True, xticklabels=conditions_order, yticklabels=conditions_order)
                    plt.xlabel("Predicted")
                    plt.ylabel("True")
                    st.pyplot(fig_confusion)
                with tab3:
                    predict_order = st.multiselect("Select the variables to predict", unique_groups)
                    if len(predict_order) >= 2:
                        predict_df = use_df[use_df[genotype_col].isin(predict_order)]
                        predict_data, predict_labels = process_LDA_data(predict_df, varIdx, genotype_col)
                        predict_test = linear_discriminant.fit(lda_data, lda_labels).predict(predict_data)
                        transform_test = linear_discriminant.fit(lda_data, lda_labels).transform(predict_data)
                        predict_confusion = np.zeros((len(predict_order), len(conditions_order)))
                        for p_i, predict in enumerate(predict_order):
                            temp_fltr = predict_test[predict_labels==predict]
                            for t_i, label in enumerate(conditions_order):
                                predict_confusion[p_i,t_i] = np.sum(temp_fltr == label)
                        fig_confusion = plt.figure()
                        sns.heatmap(data = predict_confusion, cmap='mako', annot=True, xticklabels=conditions_order, yticklabels=predict_order)
                        plt.xlabel("Predicted")
                        plt.ylabel("True")
                        st.pyplot(fig_confusion)
                if lda_feature_list:
                    for s in sort_index:
                        feat_names[s]

                img_lda = io.BytesIO()
                plt.savefig(img_lda, format='pdf')

        if 'img_lda' in locals():
            col1, col2 = st.columns(2)
            col1.download_button(label="Download LDA plot", data=img_lda, file_name="LDA_Plot.pdf", mime="application/pdf")
            download_df = convert_df(pd.DataFrame(lda))
            col2.download_button(label="Download LDA", data=download_df, file_name="Networkanalysis_LDA.csv", mime="text/csv")

st.header('Simple plot', divider=True)
with st.expander('Simple plot explanation'):
    st.markdown("""
                Here you can select one or more features, and the groups that you would like to plot as a boxplot or a bargraph.  
                The plots will follow the order of conditions based on the order that you select. All the boxes, or bars, will be gray, because the colors will be used to plot the individual data points according to their culture batch.
                """)                
if len(st.session_state.img_df) >0:
    if st.toggle('Plot data', key='Multidimensional_Plot_Toggle'):
        with st.form("Plot Options"):
            unique_groups = use_df[genotype_col].unique()
            col1, col2, col3, col4 = st.columns(4)
            feature_plot = col1.multiselect("Select feature", var_names, max_selections=1, key="Multidimension_feature_plot")
            conditions_order = col2.multiselect("Select the plotting order", unique_groups, default=control_group)
            batch_col = col3.selectbox("Select the culutue batch ID", np.append(["##"], var_names), index=0, key="sel_batchID")
            col4.write("##")
            b_normalize = col4.checkbox("Normalize")
            b_bar = col4.checkbox("Bar graph")
            do_plot = st.form_submit_button("Update plot")
            if do_plot:
                fig_plot, ax = plt.subplots(figsize=(15,10))
                batch_IDs = np.array(use_df[batch_col])
                batches = np.unique(batch_IDs)
                cmap = sns.color_palette(palette='Accent', n_colors=len(batches))
                y_data = feature_plot
                if b_normalize:
                    condition_IDs = np.array(use_df[genotype_col])
                    control_filter = condition_IDs == control_group
                    temp_values = np.array(use_df[feature_plot])
                    for batch in batches:
                       batch_filter = batch_IDs == batch
                       control_value = temp_values[(batch_filter) & (control_filter)]
                       control_value = np.mean(control_value)
                       temp_values[batch_filter] = temp_values[batch_filter] / control_value
                    use_df['NormalizeValues'] = temp_values
                    y_data = 'NormalizeValues'
                if b_bar:
                    sns.barplot(data=use_df, x=genotype_col, y=y_data, order=conditions_order, edgecolor='k', facecolor=(.7, .7, .7), errcolor='k')
                else:
                    PROPS = {
                            'boxprops':{'facecolor':(.7, .7, .7), 'edgecolor':'black'},
                        }
                    sns.boxplot(data=use_df, x=genotype_col, y=y_data, order=conditions_order, **PROPS)
                sns.swarmplot(data=use_df, x=genotype_col, y=y_data, order=conditions_order, hue=batch_col, palette=cmap)
                plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
                st.pyplot(fig_plot)
                img_plot = io.BytesIO()
                plt.savefig(img_plot, format='pdf')
        if 'img_plot' in locals():
            st.download_button(label="Download plot", data=img_plot, file_name="Simple_Plot.pdf", mime="application/pdf", key="Multidimension_save_figure2")
                
#if len(st.session_state.groupped_data) > 0:
#    # Try to perform a UMAP analysis
#    reducer = umap.UMAP(n_neighbors=50, min_dist=1)
#    # Use the LDA data
#    embedding = reducer.fit_transform(lda_data)
#    fig_umap = plt.figure()
#    sns.scatterplot(x=embedding[:, 0], y=embedding[:, 1], hue=lda_labels, palette=cmap)
#    st.pyplot(fig_umap, use_container_width=True)
